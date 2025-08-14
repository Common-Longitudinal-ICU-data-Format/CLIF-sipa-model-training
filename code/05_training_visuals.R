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
models_path <- file.path(output_path, "models")
visuals_path <- file.path(models_path, "model_visualizations")
dir.create(visuals_path, showWarnings = FALSE, recursive = TRUE)

load(file.path(models_path, "trained_models_summary.RData"))
data <- read_parquet(paste0(output_path, "/sipa_features.parquet"))

# Load all models from their individual files
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


# Create and save a data frame for the AUCs
auc_table <- data.frame(
  feature_set = c(
    "SOFA Score Only",
    "SOFA Variables before life support",
    "SOFA Variables before and after life support",
    "SOFA Variables + Age before life support",
    "SOFA Variables + Age before and after life support"
  ),
  GLM = c(
    max(glm_sofa_score$results$ROC),
    max(glm_sofa_only_pre$results$ROC),
    max(glm_sofa_only_all$results$ROC),
    max(glm_sofa_age_pre$results$ROC),
    max(glm_sofa_age_all$results$ROC)
  ),
  GAM = c(
    NA, # No GAM model for SOFA score only
    gam_sofa_only_pre_auc,
    gam_sofa_only_all_auc,
    gam_sofa_age_pre_auc,
    gam_sofa_age_all_auc
  ),
  Elastic_Net = c(
    NA, # No Elastic Net model for SOFA score only
    max(glmnet_sofa_only_pre$results$ROC),
    max(glmnet_sofa_only_all$results$ROC),
    max(glmnet_sofa_age_pre$results$ROC),
    max(glmnet_sofa_age_all$results$ROC)
  ),
  LightGBM = c(
    NA, # No LightGBM model for SOFA score only
    lightgbm_sofa_only_pre_results$best_score,
    lightgbm_sofa_only_all_results$best_score,
    lightgbm_sofa_age_pre_results$best_score,
    lightgbm_sofa_age_all_results$best_score
  )
)

print(auc_table)
write.csv(auc_table, file.path(visuals_path, "auc_table.csv"), row.names = FALSE)

# Initialize lists to store results
confusion_matrices <- list()
calibration_plots <- list()

# Recreate feature sets
hosp_sofa_score <- data %>%
  select(sofa_score_pre, sofa_score_post) %>%
  rowwise() %>%
  mutate(worst_sofa_score = max(sofa_score_pre, sofa_score_post)) %>%
  select(worst_sofa_score)

sofa_only_pre <- data %>%
  select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sofa_only_all <- data %>%
  select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre,p_f_post, s_f_post, platelets_post, bilirubin_post, map_post, gcs_post, creatinine_post, phenylephrine_post, norepinephrine_post, vasopressin_post, dopamine_post, dobutamine_post, milrinone_post, epinephrine_post, angiotensin_post) %>%
    mutate(across(everything(), ~replace_na(.x, 0)))

sofa_age_pre <- data %>%
  select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre, age_at_admission) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sofa_age_all <- data %>%
   select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre,p_f_post, s_f_post, platelets_post, bilirubin_post, map_post, gcs_post, creatinine_post, phenylephrine_post, norepinephrine_post, vasopressin_post, dopamine_post, dobutamine_post, milrinone_post, epinephrine_post, angiotensin_post, age_at_admission) %>%
    mutate(across(everything(), ~replace_na(.x, 0)))

output <- factor(data$in_hospital_mortality, levels = c(0, 1), labels = c("Alive", "Dead"))

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

data$race_group <- ifelse(data$race == "Black or African American", "Black", "Non-Black")

for (model_name in names(models)) {
  model <- models[[model_name]]
  features <- feature_sets[[model_name]]

  # Confusion Matrix
  if (inherits(model, "lgb.Booster")) {
    pred_probs_cm <- predict(model, features)
    pred_char <- ifelse(pred_probs_cm > 0.5, "Dead", "Alive")
  } else if (inherits(model, "gam")) {
    pred_probs_cm <- predict(model, newdata = features, type = "response")
    pred_char <- ifelse(pred_probs_cm > 0.5, "Dead", "Alive")
  } else { # This now only handles caret models
    pred_char <- as.character(predict(model, newdata = features))
  }
  predictions <- factor(pred_char, levels = levels(output))
  cm <- confusionMatrix(predictions, output, positive = "Dead")

  # Extract metrics and store in a data frame
  metrics_df <- data.frame(
    TP = cm$table[2, 2],
    TN = cm$table[1, 1],
    FP = cm$table[1, 2],
    FN = cm$table[2, 1],
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]),
    PPV = unname(cm$byClass["Pos Pred Value"])
  )
  rownames(metrics_df) <- NULL
  
  confusion_matrices[[model_name]] <- metrics_df

  # Calibration Plot
  if (inherits(model, "lgb.Booster")) {
    pred_probs_cal <- predict(model, features)
  } else if (inherits(model, "gam")) {
    pred_probs_cal <- predict(model, newdata = features, type = "response")
  } else {
    pred_probs_cal <- predict(model, newdata = features, type = "prob")$Dead
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

print("Finished. Model visuals and tables saved to output/models/model_visualizations")

# --- Allocation Efficiency Analysis ---

# 1. Find the best model based on AUC
candidate_models <- list(
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

candidate_aucs <- c(
  max(glm_sofa_only_pre$results$ROC),
  max(glm_sofa_only_all$results$ROC),
  max(glm_sofa_age_pre$results$ROC),
  max(glm_sofa_age_all$results$ROC),
  gam_sofa_only_pre_auc,
  gam_sofa_only_all_auc,
  gam_sofa_age_pre_auc,
  gam_sofa_age_all_auc,
  max(glmnet_sofa_only_pre$results$ROC),
  max(glmnet_sofa_only_all$results$ROC),
  max(glmnet_sofa_age_pre$results$ROC),
  max(glmnet_sofa_age_all$results$ROC),
  lightgbm_sofa_only_pre_results$best_score,
  lightgbm_sofa_only_all_results$best_score,
  lightgbm_sofa_age_pre_results$best_score,
  lightgbm_sofa_age_all_results$best_score
)

best_model_idx <- which.max(candidate_aucs)
best_model_name <- names(candidate_models)[best_model_idx]
best_model_obj <- candidate_models[[best_model_idx]]
best_model_features <- feature_sets[[best_model_name]]

# 2. Get models and predictions
pred_glm_sofa <- predict(glm_sofa_score, newdata = hosp_sofa_score, type = "prob")$Dead

if (inherits(best_model_obj, "lgb.Booster")) {
  pred_best_model <- predict(best_model_obj, best_model_features)
} else if (inherits(best_model_obj, "gam")) {
  pred_best_model <- predict(best_model_obj, newdata = best_model_features, type = "response")
} else {
  pred_best_model <- predict(best_model_obj, newdata = best_model_features, type = "prob")$Dead
}

# 3. Define and calculate efficiency
calculate_efficiency <- function(predicted_mortality, true_mortality, s_grid = seq(0, 1, by = 0.01)) {
  N <- length(true_mortality)
  survivors <- 1 - true_mortality
  res <- purrr::map_dfr(s_grid, function(s) {
    V <- floor(s * N)
    if (V == 0) {
      eff <- NA_real_
    } else {
      selected_idx <- order(predicted_mortality)[1:V]
      survivors_allocated <- sum(survivors[selected_idx])
      eff <- survivors_allocated / V
    }
    tibble(s = s, efficiency = eff)
  })
  res
}

true_mortality_vector <- data$in_hospital_mortality
eff_random <- calculate_efficiency(runif(length(true_mortality_vector)), true_mortality_vector)
eff_glm_sofa <- calculate_efficiency(pred_glm_sofa, true_mortality_vector)
eff_best_model <- calculate_efficiency(pred_best_model, true_mortality_vector)

# 4. Combine for plotting
eff_all <- bind_rows(
  eff_random %>% mutate(model = "Random"),
  eff_glm_sofa %>% mutate(model = "GLM SOFA Score"),
  eff_best_model %>% mutate(model = best_model_name)
)

# 5. Plot and save
allocation_plot <- ggplot(eff_all, aes(x = s, y = efficiency, color = model)) +
  geom_line(size = 1.2, na.rm = TRUE) +
  scale_x_continuous(name = "Life support supply, s", breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(name = "Allocation Efficiency", limits = c(0, 1)) +
  ggtitle("Comparison of Allocation Efficiency Across Supply Levels") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  scale_color_manual(values = wes_palette("GrandBudapest2", n = 3))

saveRDS(allocation_plot, file.path(visuals_path, "allocation_efficiency_plot.rds"))

print("Finished Allocation Efficiency Analysis.")
