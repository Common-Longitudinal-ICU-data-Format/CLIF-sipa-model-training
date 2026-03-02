# Feature set creation for model development
# Requires: clif_sofa2_scores.parquet created by sofa2_calculation.py
# Requires: sipa_clif_cohort.parquet created by 01_cohort_identification.R

# Load necessary libraries
print("Initialized Feature Set Processing Script")

library(arrow)
library(tidyverse)
library(stringr)
library(data.table)
library(tictoc)
library(duckdb)

# Clear env
rm(list = ls())

# Path
source("utils/config.R")
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type
output_path <- config$output_path


# Load data
cohort <- read_parquet(file.path(output_path, "intermediate", "sipa_clif_cohort.parquet"))
sofa2_results <- read_parquet(file.path(output_path, "intermediate", "clif_sofa2_scores.parquet"))

tic()
# Only certain variables are important from sofa2_results
sofa2_selected <- sofa2_results %>%
  select(hospitalization_id, start_dttm, end_dttm, sofa2_total, sofa2_brain, sofa2_resp, sofa2_cv, 
         sofa2_liver, sofa2_kidney, sofa2_hemo, gcs_min, has_sedation, 
         has_delirium_drug, pf_ratio, sf_ratio, device_category, 
         pao2_at_worst, spo2_at_worst, fio2_at_worst, has_ecmo, 
         map_min, norepi_epi_maxsum, dopa_max, has_other_non_dopa, 
         has_other_vaso, has_mechanical_cv_support, bilirubin_total, 
         creatinine, potassium, ph, bicarbonate)

# Rename specific variables for clarity and consistency
sofa2_selected <- sofa2_selected %>% 
  rename(gcs = gcs_min,
         map = map_min,
         pf = pf_ratio,
         sf = sf_ratio,
         pao2 = pao2_at_worst,
         spo2 = spo2_at_worst,
         fio2 = fio2_at_worst,
         bilirubin = bilirubin_total)

# Drop duplicates from cohort based on hospitalization_id
# Select appropriate patient demographic variables and vasoactive medication variables
cohort <- cohort %>%
  distinct(hospitalization_id, .keep_all = TRUE) %>% 
  select(hospitalization_id, age_at_admission, race_category, 
         ethnicity_category, sex_category, in_hospital_mortality, milrinone, 
         norepinephrine, phenylephrine, vasopressin, dopamine,
         epinephrine, dobutamine, angiotensin, metaraminol, fio2_filled)

# Merge cohort with sofa2_selected on hospitalization_id
final_data <- cohort %>%
  left_join(sofa2_selected, by = "hospitalization_id")

# Create the directory if it doesn't exist
dir.create(file.path(output_path, "final"), showWarnings = FALSE, recursive = TRUE)

# Export final data
write_parquet(final_data, file.path(output_path, "final", "sipa_features.parquet"))
toc()
print(paste("Feature set exported to", file.path(output_path, "final")))
