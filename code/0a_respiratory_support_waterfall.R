# Respiratory Support Waterfall Script using clifpy Python package
# This replaces the manual R implementation with the official Python version

rm(list = ls())

# Load necessary libraries
library(arrow)
library(tidyverse)
library(reticulate)

# Access configuration parameters
source("utils/config.R")
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type
output_path <- config$output_path



# Setup Python env --------------------------------------------------------

print("Checking Python environment setup...")

# Check if virtual environment exists, create if needed
if (!virtualenv_exists("clif-env")) {
  message("Creating Python virtual environment 'clif-env'...")
  virtualenv_create("clif-env")
  message("Installing clifpy package...")
  py_install("clifpy", envname = "clif-env")
  message("clifpy installed successfully!")
} else {
  message("Python virtual environment 'clif-env' already exists")
}

# Activate the virtual environment
use_virtualenv("clif-env")

# Import the specific module from clifpy
print("Loading clifpy.tables.respiratory_support module...")
resp_module <- import("clifpy.tables.respiratory_support")
RespiratorySupport <- resp_module$RespiratorySupport
print("clifpy loaded successfully!")


# Run Waterfall -----------------------------------------------------------

print("Loading respiratory support data via clifpy...")

# Determine file type string (remove leading dot if present)
filetype_clean <- gsub("^\\.", "", file_type)

# Set timezone - adjust this to match your data's timezone
timezone <- "US/Eastern"  # Or hardcode: "UTC", "US/Eastern", etc.

# Load data using clifpy's from_file method
resp_support <- RespiratorySupport$from_file(
  data_directory = tables_path,
  filetype = filetype_clean,
  timezone = timezone
)

print("Running clifpy respiratory waterfall function...")

# Apply waterfall processing
processed <- resp_support$waterfall()

print("Validating processed data...")

# Validate the processed data
processed$validate()

print("Extracting processed DataFrame...")

# Extract the processed DataFrame
py_df <- processed$df

# Define the output path
intermediate_output_path <- file.path(output_path, "intermediate")
dir.create(intermediate_output_path, showWarnings = FALSE, recursive = TRUE)

temp_output_file <- file.path(
  intermediate_output_path, 
  "temp_respiratory_support.parquet"
)

print("Saving processed data from Python...")

# Save directly from Python using pandas
py_df$to_parquet(temp_output_file)

print("Reading data back into R...")

# Read back into R
df_resp_support <- read_parquet(temp_output_file)

print(paste("✓ Processed", nrow(df_resp_support), "records from clifpy waterfall"))

# Clean up temp file
file.remove(temp_output_file)

# Clean up Python objects
rm(resp_support, processed, py_df)


# Additional processing ---------------------------------------------------

print("Processing FiO2 values...")

# Load the device_category to ranges mapping table
category_values <- read_csv("lookup-tables/device_category_to_conversion.csv")

category_values <- category_values %>%
  mutate(device_category = str_trim(device_category), 
         device_category = tolower(device_category))

# Ensure device_category is trimmed of whitespace and merge with mapping
# Convert device_category to lower case
df_resp_support <- df_resp_support %>%
  mutate(device_category = str_trim(device_category), 
         device_category = tolower(device_category)) %>%
  left_join(category_values, by = "device_category")

# Check ranges for fio2_set
df_resp_support_conv <- df_resp_support %>% 
  mutate(fio2_set = case_when(
    !is.na(fio2_set) & !is.na(range_lower) & fio2_set < range_lower ~ range_lower,
    !is.na(fio2_set) & !is.na(range_upper) & fio2_set > range_upper ~ range_upper,
    TRUE ~ fio2_set))

# Impute FiO2 values when fio2_set is NA and lpm_set is not NA
df_resp_support_conv <- df_resp_support_conv %>%
  mutate(fio2_set = case_when(
    is.na(fio2_set) & !is.na(lpm_set) & !is.na(conversion) ~ {
      fio2_imp <- 0.21 + lpm_set * conversion
      pmin(pmax(fio2_imp, range_lower), range_upper)
    },
    TRUE ~ fio2_set))

# Rename fio2_set to fio2_approx
df_resp_support_conv <- df_resp_support_conv %>% 
  rename(fio2_approx = fio2_set)

print("FiO2 summary (original):")
print(summary(df_resp_support$fio2_set))

print("FiO2 summary (after conversion):")
print(summary(df_resp_support_conv$fio2_approx))

# If there are still NA values in fio2_approx, set them to range_lower if available
df_resp_support_conv <- df_resp_support_conv %>%
  mutate(fio2_approx = ifelse(is.na(fio2_approx) & !is.na(range_lower), range_lower, fio2_approx)) %>% 
  # Ensure fio2_approx is 0.21 for Room Air, as lookup table doesn't provide range_lower
  mutate(fio2_approx = ifelse(device_category == "room air", 0.21, fio2_approx)) 

print("FiO2 summary (final):")
print(summary(df_resp_support_conv$fio2_approx))


# Save output -------------------------------------------------------------

# Save the processed data
output_file <- file.path(
  intermediate_output_path, 
  paste0("clif_respiratory_support_processed", file_type))

write_parquet(df_resp_support_conv, output_file)

print(paste("✓ Table exported as parquet to", output_file))
print("✓ Respiratory support waterfall complete!")

# Clean up
rm(df_resp_support, df_resp_support_conv, category_values)