# Load libraries
library(tidyverse)
library(knitr)
library(caret)
library(pROC)
library(ggplot2)
library(gridExtra)
library(lightgbm)
library(mgcv)
library(glmnet)
library(arrow)
library(wesanderson)

rm(list = ls())

# Load data
source("utils/config.R")
output_path <- config$output_path
models_path <- file.path(output_path, "exportable", "models")
visuals_path <- file.path(models_path, "model_visualizations")
dir.create(visuals_path, showWarnings = FALSE, recursive = TRUE)

# --- Data Loading and Reshaping (Mirrors 04_model_training.R) ---
# Load the long-format data
data_long <- read_parquet(file.path(output_path, "final", "sipa_features.parquet"))

# Rank rows within each hospitalization to identify pre/post periods
data_ranked <- data_long %>%
  group_by(hospitalization_id) %>%
  arrange(start_dttm) %>%
  mutate(period_rank = row_number()) %>%
  ungroup()

# Define which columns are static vs. which need to be suffixed
static_cols <- c("hospitalization_id", "age_at_admission", "race_category", "ethnicity_category", "sex_category", "in_hospital_mortality")
feature_cols <- setdiff(names(data_long), c(static_cols, "start_dttm", "end_dttm", "period_rank"))

# Create the 'pre' dataframe
pre_data <- data_ranked %>%
  filter(period_rank == 1) %>%
  select(all_of(static_cols), all_of(feature_cols))
names(pre_data) <- c(static_cols, paste0(feature_cols, "_pre"))

# Create the 'post' dataframe
post_data <- data_ranked %>%
  filter(period_rank == 2) %>%
  select(hospitalization_id, all_of(feature_cols))
names(post_data) <- c("hospitalization_id", paste0(feature_cols, "_post"))

# Join pre and post dataframes to create the final wide format
data <- pre_data %>%
  full_join(post_data, by = "hospitalization_id")


# --- Load Models and Results ---
# Load the single summary RDS file which contains all results
model_summary <- readRDS(file.path(models_path, "model_summary.rds"))

# Load all individual model objects
glm_sofa_score <- readRDS(file.path(models_path, "glm_sofa_score.rds"))
glm_sofa_only_pre <- readRDS(file.path(models_path, "glm_sofa_only_pre.rds"))
glm_sofa_only_all <- readRDS(file.path(models_path, "glm_sofa_only_all.rds"))
glm_sofa_age_pre <- readRDS(file.path(models_path, "glm_sofa_age_pre.rds"))
glm_sofa_age_all <- readRDS(file.path(models_path, "glm_sofa_age_all.rds"))
gam_sofa_only_pre <- readRDS(file.path(models_path, "gam_sofa_only_pre.rds"))
gam_sofa_only_all <- readRDS(file.path(models_path, "gam_sofa_only_all.rds"))
gam_sofa_age_pre <- readRDS(file.path(models_path, "gam_sofa_age_pre.rds"))
gam_sofa_age_all <- readRDS(file.path(models_path, "gam_sofa_age_all.rds"))
glmnet_sofa_only_pre <- readRDS(file.path(models_path, "glmnet_sofa_only_pre.rds"))
glmnet_sofa_only_all <- readRDS(file.path(models_path, "glmnet_sofa_only_all.rds"))
glmnet_sofa_age_pre <- readRDS(file.path(models_path, "glmnet_sofa_age_pre.rds"))
glmnet_sofa_age_all <- readRDS(file.path(models_path, "glmnet_sofa_age_all.rds"))
lightgbm_sofa_only_pre <- lgb.load(file.path(models_path, "lightgbm_sofa_only_pre.txt"))
lightgbm_sofa_only_all <- lgb.load(file.path(models_path, "lightgbm_sofa_only_all.txt"))
lightgbm_sofa_age_pre <- lgb.load(file.path(models_path, "lightgbm_sofa_age_pre.txt"))
lightgbm_sofa_age_all <- lgb.load(file.path(models_path, "lightgbm_sofa_age_all.txt"))

# --- AUC Table Creation ---
# Create and save a data frame for the AUCs using the loaded summary object
auc_table <- data.frame(
  feature_set = c(
    "SOFA Score Only",
    "SOFA Variables before life support",
    "SOFA Variables before and after life support",
    "SOFA Variables + Age before life support",
    "SOFA Variables + Age before and after life support"
  ),
  GLM = c(
    max(model_summary$glm_sofa_score$ROC),
    max(model_summary$glm_sofa_only_pre$ROC),
    max(model_summary$glm_sofa_only_all$ROC),
    max(model_summary$glm_sofa_age_pre$ROC),
    max(model_summary$glm_sofa_age_all$ROC)
  ),
  GAM = c(
    NA, # No GAM model for SOFA score only
    model_summary$gam_sofa_only_pre_auc,
    model_summary$gam_sofa_only_all_auc,
    model_summary$gam_sofa_age_pre_auc,
    model_summary$gam_sofa_age_all_auc
  ),
  Elastic_Net = c(
    NA, # No Elastic Net model for SOFA score only
    model_summary$glmnet_sofa_only_pre$ROC,
    model_summary$glmnet_sofa_only_all$ROC,
    model_summary$glmnet_sofa_age_pre$ROC,
    model_summary$glmnet_sofa_age_all$ROC
  ),
  LightGBM = c(
    NA, # No LightGBM model for SOFA score only
    model_summary$lgbm_sofa_only_pre_auc,
    model_summary$lgbm_sofa_only_all_auc,
    model_summary$lgbm_sofa_age_pre_auc,
    model_summary$lgbm_sofa_age_all_auc
  )
)

print(auc_table)
write.csv(auc_table, file.path(visuals_path, "auc_table.csv"), row.names = FALSE)


# --- Feature Set Recreation (Mirrors 04_model_training.R) ---
sofa_vars_list <- c("milrinone", "norepinephrine", "phenylephrine", "vasopressin", 
                    "dopamine", "epinephrine", "dobutamine", "angiotensin", 
                    "metaraminol", "gcs", "sf", "map", "bilirubin", "potassium", 
                    "bicarbonate", "creatinine", "ph", "pf")

hosp_sofa_score <- data %>%
  select(sofa2_total_pre, sofa2_total_post) %>%
  rowwise() %>%
  mutate(worst_sofa_score = max(sofa2_total_pre, sofa2_total_post, na.rm = TRUE)) %>%
  ungroup() %>%
  select(worst_sofa_score)

sofa_only_pre <- data %>%
  select(all_of(paste0(sofa_vars_list, "_pre"))) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sofa_only_all <- data %>%
  select(all_of(paste0(sofa_vars_list, "_pre")), all_of(paste0(sofa_vars_list, "_post"))) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sofa_age_pre <- data %>%
  select(all_of(paste0(sofa_vars_list, "_pre")), age_at_admission) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sofa_age_all <- data %>%
   select(all_of(paste0(sofa_vars_list, "_pre")), all_of(paste0(sofa_vars_list, "_post")), age_at_admission) %>%
   mutate(across(everything(), ~replace_na(.x, 0)))

output <- factor(data$in_hospital_mortality, levels = c(0, 1), labels = c("Alive", "Dead"))


# --- Visualizations ---

# Initialize lists to store results
confusion_matrices <- list()
calibration_plots <- list()

models <- list(
  "GLM_SOFA_Score" = glm_sofa_score,
  "GLM_SOFA_Vars_Pre" = glm_sofa_only_pre,
  "GLM_SOFA_Vars_All" = glm_sofa_only_all,
  "GLM_SOFA_Age_Pre" = glm_sofa_age_pre,
  "GLM_SOFA_Age_All" = glm_sofa_age_all,
  "GAM_SOFA_Vars_Pre" = gam_sofa_only_pre,
  "GAM_SOFA_Vars_All" = gam_sofa_only_all,
  "GAM_SOFA_Age_Pre" = gam_sofa_age_pre,
  "GAM_SOFA_Age_All" = gam_sofa_age_all,
  "Elastic_Net_SOFA_Vars_Pre" = glmnet_sofa_only_pre,
  "Elastic_Net_SOFA_Vars_All" = glmnet_sofa_only_all,
  "Elastic_Net_SOFA_Age_Pre" = glmnet_sofa_age_pre,
  "Elastic_Net_SOFA_Age_All" = glmnet_sofa_age_all,
  "LightGBM_SOFA_Vars_Pre" = lightgbm_sofa_only_pre,
  "LightGBM_SOFA_Vars_All" = lightgbm_sofa_only_all,
  "LightGBM_SOFA_Age_Pre" = lightgbm_sofa_age_pre,
  "LightGBM_SOFA_Age_All" = lightgbm_sofa_age_all
)

feature_sets <- list(
  "GLM_SOFA_Score" = hosp_sofa_score,
  "GLM_SOFA_Vars_Pre" = sofa_only_pre,
  "GLM_SOFA_Vars_All" = sofa_only_all,
  "GLM_SOFA_Age_Pre" = sofa_age_pre,
  "GLM_SOFA_Age_All" = sofa_age_all,
  "GAM_SOFA_Vars_Pre" = sofa_only_pre,
  "GAM_SOFA_Vars_All" = sofa_only_all,
  "GAM_SOFA_Age_Pre" = sofa_age_pre,
  "GAM_SOFA_Age_All" = sofa_age_all,
  "Elastic_Net_SOFA_Vars_Pre" = sofa_only_pre,
  "Elastic_Net_SOFA_Vars_All" = sofa_only_all,
  "Elastic_Net_SOFA_Age_Pre" = sofa_age_pre,
  "Elastic_Net_SOFA_Age_All" = sofa_age_all,
  "LightGBM_SOFA_Vars_Pre" = as.matrix(sofa_only_pre),
  "LightGBM_SOFA_Vars_All" = as.matrix(sofa_only_all),
  "LightGBM_SOFA_Age_Pre" = as.matrix(sofa_age_pre),
  "LightGBM_SOFA_Age_All" = as.matrix(sofa_age_all)
)

data$race_group <- ifelse(data$race_category == "Black or African American", "Black", "Non-Black")

for (model_name in names(models)) {
  model <- models[[model_name]]
  features <- feature_sets[[model_name]]

  # Confusion Matrix
  if (inherits(model, "lgb.Booster")) {
    pred_probs_cm <- predict(model, features)
    pred_char <- ifelse(pred_probs_cm > 0.5, "Dead", "Alive")
  } else if (inherits(model, "gam")) {
    pred_probs_cm <- predict(model, newdata = as.data.frame(features), type = "response")
    pred_char <- ifelse(pred_probs_cm > 0.5, "Dead", "Alive")
  } else { # This now only handles caret models
    pred_char <- as.character(predict(model, newdata = as.data.frame(features)))
  }
  predictions <- factor(pred_char, levels = levels(output))
  cm <- confusionMatrix(predictions, output, positive = "Dead")

  # Extract metrics
  metrics_df <- as.data.frame(t(c(
    cm$overall,
    cm$byClass
  )))
  confusion_matrices[[model_name]] <- metrics_df

  # Calibration Plot
  if (inherits(model, "lgb.Booster")) {
    pred_probs_cal <- predict(model, features)
  } else if (inherits(model, "gam")) {
    pred_probs_cal <- predict(model, newdata = as.data.frame(features), type = "response")
  } else {
    pred_probs_cal <- predict(model, newdata = as.data.frame(features), type = "prob")$Dead
  }

  cal_data <- data.frame(
    prob = pred_probs_cal,
    y = output,
    race = data$race_group
  )

  cal_obj_black <- calibration(y ~ prob, data = subset(cal_data, race == "Black"), class = "Dead")
  cal_obj_non_black <- calibration(y ~ prob, data = subset(cal_data, race == "Non-Black"), class = "Dead")

  plot_data <- rbind(
    data.frame(cal_obj_black$data, race = "Black Patients"),
    data.frame(cal_obj_non_black$data, race = "Non-Black Patients")
  )

  p <- ggplot(plot_data, aes(x = midpoint, y = Percent, color = race)) +
    geom_line() +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      title = paste("Calibration Plot for", model_name),
      x = "Predicted Probability",
      y = "Observed Frequency"
    ) +
    theme_classic() +
    scale_color_manual(values = wes_palette("GrandBudapest2", n = 2))

  calibration_plots[[model_name]] <- p
}

# Save the results
saveRDS(confusion_matrices, file.path(visuals_path, "confusion_matrices.rds"))
saveRDS(calibration_plots, file.path(visuals_path, "calibration_plots.rds"))

print("Finished. Model visuals and tables saved to output/exportable/models/model_visualizations")
