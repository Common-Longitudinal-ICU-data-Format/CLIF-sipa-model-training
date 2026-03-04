# Load necessary libraries
print("Initializing Table1 Script")
library(arrow)
library(tidyverse)
library(stringr)
library(data.table)
library(tictoc)
library(glue)

tic()
# Clear env
rm(list = ls())

# Load data
source("utils/config.R")
output_path <- config$output_path
site_name <- config$site_name

# Create the exportable directory if it doesn't exist
dir.create(file.path(output_path, "exportable"), showWarnings = FALSE, recursive = TRUE)

# Load the SIPA features dataset
data_raw <- read_parquet(file.path(output_path, "final", "sipa_features.parquet"))

# Summarize data to one row per hospitalization_id by averaging numeric clinical variables
data <- data_raw %>%
  group_by(hospitalization_id) %>%
  summarise(
    # Demographics and outcomes - take the first value
    age_at_admission = first(age_at_admission),
    race_category = first(race_category),
    ethnicity_category = first(ethnicity_category),
    sex_category = first(sex_category),
    in_hospital_mortality = first(in_hospital_mortality),

    # SOFA scores - average the two rows
    sofa2_total = mean(sofa2_total, na.rm = TRUE),
    sofa2_brain = mean(sofa2_brain, na.rm = TRUE),
    sofa2_resp = mean(sofa2_resp, na.rm = TRUE),
    sofa2_cv = mean(sofa2_cv, na.rm = TRUE),
    sofa2_liver = mean(sofa2_liver, na.rm = TRUE),
    sofa2_kidney = mean(sofa2_kidney, na.rm = TRUE),
    sofa2_hemo = mean(sofa2_hemo, na.rm = TRUE),

    # Clinical variables - average the two rows
    gcs = mean(gcs, na.rm = TRUE),
    map = mean(map, na.rm = TRUE),
    creatinine = mean(creatinine, na.rm = TRUE),
    bilirubin = mean(bilirubin, na.rm = TRUE),
    potassium = mean(potassium, na.rm = TRUE),
    ph = mean(ph, na.rm = TRUE),
    bicarbonate = mean(bicarbonate, na.rm = TRUE),
    pao2 = mean(pao2, na.rm = TRUE),
    spo2 = mean(spo2, na.rm = TRUE),
    fio2 = mean(fio2, na.rm = TRUE),
    pf = mean(pf, na.rm = TRUE),
    sf = mean(sf, na.rm = TRUE)
  ) %>%
  ungroup()


# Helper for median (IQR), 10th/90th deciles, and NA count
summary_stats_na <- function(x) {
  non_na_count <- sum(!is.na(x))
  na_count <- sum(is.na(x))
  if (non_na_count == 0) {
    return(paste0("NA: ", na_count))
  }
  paste0(formatC(median(x, na.rm = TRUE), digits = 2, format = "f"),
         " (", 
         formatC(quantile(x, 0.25, na.rm = TRUE), digits = 2, format = "f"), ", ",
         formatC(quantile(x, 0.75, na.rm = TRUE), digits = 2, format = "f"), ")",
         "; ",
         formatC(quantile(x, 0.1, na.rm = TRUE), digits = 2, format = "f"), "-",
         formatC(quantile(x, 0.9, na.rm = TRUE), digits = 2, format = "f"),
         "; ", na_count)
}


# Define new variable display names
var_display_names <- list(
  "age_at_admission" = "Age (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_total" = "SOFA-2 Total Score (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_brain" = "SOFA-2 Brain Score (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_resp" = "SOFA-2 Respiratory Score (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_cv" = "SOFA-2 Cardiovascular Score (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_liver" = "SOFA-2 Liver Score (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_kidney" = "SOFA-2 Kidney Score (Median, IQR; 10th-90th Decile; NA)",
  "sofa2_hemo" = "SOFA-2 Hematology Score (Median, IQR; 10th-90th Decile; NA)",
  "gcs" = "Glasgow Coma Scale (Median, IQR; 10th-90th Decile; NA)",
  "map" = "Mean Arterial Pressure (Median, IQR; 10th-90th Decile; NA)",
  "creatinine" = "Creatinine (Median, IQR; 10th-90th Decile; NA)",
  "bilirubin" = "Bilirubin (Median, IQR; 10th-90th Decile; NA)",
  "potassium" = "Potassium (Median, IQR; 10th-90th Decile; NA)",
  "ph" = "pH (Median, IQR; 10th-90th Decile; NA)",
  "bicarbonate" = "Bicarbonate (Median, IQR; 10th-90th Decile; NA)",
  "pao2" = "PaO2 (Median, IQR; 10th-90th Decile; NA)",
  "spo2" = "SpO2 (Median, IQR; 10th-90th Decile; NA)",
  "fio2" = "FiO2 (Median, IQR; 10th-90th Decile; NA)",
  "pf" = "PaO2/FiO2 Ratio (Median, IQR; 10th-90th Decile; NA)",
  "sf" = "SpO2/FiO2 Ratio (Median, IQR; 10th-90th Decile; NA)")

table1 <- list()

# Number of ICU Encounters (N)
table1[["ICU Encounters (N)"]] <- nrow(data)

# In-Hospital Mortality Rate
mortality_n <- sum(data$in_hospital_mortality, na.rm = TRUE)
mortality_pct <- 100 * mean(data$in_hospital_mortality, na.rm = TRUE)
table1[["In-Hospital Mortality Rate (N, %)"]] <- paste0(mortality_n, " (", formatC(mortality_pct, digits = 1, format = "f"), "%)")

# Sex (N, % female)
n_female <- sum(data$sex_category == "Female", na.rm = TRUE)
pct_female <- 100 * mean(data$sex_category == "Female", na.rm = TRUE)
table1[["Sex (N, % Female)"]] <- paste0(n_female, " (", formatC(pct_female, digits = 1, format = "f"), "%)")

# Blank Race header row
table1[["Race (N, %)"]] <- ""
race_levels <- c("White", "Asian", "Native Hawaiian or Other Pacific Islander", 
                 "Black or African American", "Unknown", "Other", "American Indian or Alaska Native")
race_tab <- data %>%
  filter(race_category %in% race_levels) %>%
  group_by(race_category) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Total = sum(N),
         Percent = 100 * N / Total,
         Race = paste0(N, " (", formatC(Percent, digits = 1, format = "f"), "%)")) %>%
  select(race_category, Race)

for (level in race_levels) {
  val <- race_tab$Race[race_tab$race_category == level]
  table1[[paste0("  ", level)]] <- if (length(val) > 0) val else "0 (0.0%)"
}

# Blank Ethnicity header row
table1[["Ethnicity (N, %)"]] <- ""
ethnicity_levels <- c("Hispanic", "Non-Hispanic", "Unknown")
eth_tab <- data %>%
  filter(ethnicity_category %in% ethnicity_levels) %>%
  group_by(ethnicity_category) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Total = sum(N),
         Percent = 100 * N / Total,
         Eth = paste0(N, " (", formatC(Percent, digits = 1, format = "f"), "%)")) %>%
  select(ethnicity_category, Eth)

for (level in ethnicity_levels) {
  val <- eth_tab$Eth[eth_tab$ethnicity_category == level]
  table1[[paste0("  ", level)]] <- if (length(val) > 0) val else "0 (0.0%)"
}

# Numeric variables
numeric_vars <- names(var_display_names)
for (v in numeric_vars) {
  display_name <- var_display_names[[v]]
  table1[[display_name]] <- summary_stats_na(data[[v]])
}

# Build the table in the desired order
desired_order <- c(
  "ICU Encounters (N)",
  "In-Hospital Mortality Rate (N, %)",
  "Sex (N, % Female)",
  "Race (N, %)",
  paste0("  ", race_levels),
  "Ethnicity (N, %)",
  paste0("  ", ethnicity_levels),
  unname(unlist(var_display_names))
)

# Defensive extraction: always return blank if not found
table1_df <- tibble::tibble(
  Variable = desired_order,
  !!site_name := sapply(desired_order, function(x) table1[[x]] %||% "")
)

# Export as a separate CSV file
write.csv(table1_df, file.path(output_path, "exportable", glue("table1_{site_name}.csv")), row.names = FALSE)
print("Table 1 exported as CSV to output/exportable")
toc()