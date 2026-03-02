# Load data
library(arrow)
library(ggplot2)
library(caret)
library(pROC)
library(tictoc)
library(mgcv)
library(lightgbm)
library(tidyverse)
library(glmnet)

rm(list = ls())

source("utils/config.R")
output_path <- config$output_path

set.seed(42) # the meaning of life


# Data Loading and Reshaping ----------------------------------------------
tic("Loading and reshaping data")

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
# Automatically determine feature columns, excluding keys and temporary rank column
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

toc()


# Input and Output Vector Creation ----------------------------------------------------
tic("Creating input and output vectors.")

# Define the output vector
output <- data$in_hospital_mortality

# Model 1 Input: Worst SOFA Score
# Using the new sofa2_total_pre/post columns
hosp_sofa_score <- data %>%
  select(sofa2_total_pre, sofa2_total_post) %>%
  rowwise() %>%
  mutate(worst_sofa_score = max(sofa2_total_pre, sofa2_total_post, na.rm = TRUE)) %>%
  ungroup() %>%
  select(worst_sofa_score)

# List of all SOFA variables to include
sofa_vars_list <- c("milrinone", "norepinephrine", "phenylephrine", "vasopressin", 
                    "dopamine", "epinephrine", "dobutamine", "angiotensin", 
                    "metaraminol", "gcs", "sf", "map", "bilirubin", "potassium", 
                    "bicarbonate", "creatinine", "ph", "pf")

# Model 2 Input: SOFA Variables before life support
sofa_only_pre <- data %>%
  select(all_of(paste0(sofa_vars_list, "_pre"))) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# Model 3 Input: SOFA variables before and after life support
sofa_only_all <- data %>%
  select(all_of(paste0(sofa_vars_list, "_pre")),
         all_of(paste0(sofa_vars_list, "_post"))) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# Model 4 Input: SOFA variables and age before life support
sofa_age_pre <- data %>%
  select(all_of(paste0(sofa_vars_list, "_pre")), age_at_admission) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# Model 5 Input: SOFA variables and age before and after life support
sofa_age_all <- data %>%
   select(all_of(paste0(sofa_vars_list, "_pre")),
          all_of(paste0(sofa_vars_list, "_post")),
          age_at_admission) %>%
   mutate(across(everything(), ~replace_na(.x, 0)))

toc()


# Elastic Net Feature Analysis -------------------------------------------------------------
tic("Determine feature importance using elastic net.")
x_sofa_only_pre <- as.matrix(sofa_only_pre)
x_sofa_only_all <- as.matrix(sofa_only_all)
x_sofa_age_pre <- as.matrix(sofa_age_pre)
x_sofa_age_all <- as.matrix(sofa_age_all)
y_glmnet <- data$in_hospital_mortality

cvfit_sofa_only_pre <- cv.glmnet(x_sofa_only_pre, y_glmnet, family = "binomial", type.measure = "auc")
sofa_only_pre_coef <- coef(cvfit_sofa_only_pre, s = "lambda.min")

cvfit_sofa_only_all <- cv.glmnet(x_sofa_only_all, y_glmnet, family = "binomial", type.measure = "auc")
sofa_only_all_coef <- coef(cvfit_sofa_only_all, s = "lambda.min")

cvfit_sofa_age_pre <- cv.glmnet(x_sofa_age_pre, y_glmnet, family = "binomial", type.measure = "auc")
sofa_age_pre_coef <- coef(cvfit_sofa_age_pre, s = "lambda.min")

cvfit_sofa_age_all <- cv.glmnet(x_sofa_age_all, y_glmnet, family = "binomial", type.measure = "auc")
sofa_age_all_coef <- coef(cvfit_sofa_age_all, s = "lambda.min")

print(sofa_only_pre_coef)
print(sofa_only_all_coef)
toc()


# Standard Logistic Regression --------------------------------------------

tic("Training Logistic Regression Models")
print("-----LOGISTIC REGRESSION-----")

# Set up 5-fold cross-validation
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Model 1: SOFA Score
model_df_1 <- data.frame(
  hosp_sofa_score,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead"))
)
glm_sofa_score <- train(
  output ~ worst_sofa_score,
  data = model_df_1,
  method = "glm", family = binomial(link = "logit"), metric = "ROC", trControl = control
)
print(glm_sofa_score)

# Model 2: SOFA Variables Before
model_df_2 <- data.frame(
  sofa_only_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead"))
)
glm_sofa_only_pre <- train(
  output ~ .,
  data = model_df_2,
  method = "glm", family = binomial(link = "logit"), metric = "ROC", trControl = control
)
print(glm_sofa_only_pre)

# Model 3: SOFA Variables Before and After
model_df_3 <- data.frame(
  sofa_only_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead"))
)
glm_sofa_only_all <- train(
  output ~ .,
  data = model_df_3,
  method = "glm", family = binomial(link = "logit"), metric = "ROC", trControl = control
)
print(glm_sofa_only_all)

# Model 4: SOFA Variables + Age Before
model_df_4 <- data.frame(
  sofa_age_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead"))
)
glm_sofa_age_pre <- train(
  output ~ .,
  data = model_df_4,
  method = "glm", family = binomial(link = "logit"), metric = "ROC", trControl = control
)
print(glm_sofa_age_pre)

# Model 5: SOFA Variables + Age Before and After
model_df_5 <- data.frame(
  sofa_age_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead"))
)
glm_sofa_age_all <- train(
  output ~ .,
  data = model_df_5,
  method = "glm", family = binomial(link = "logit"), metric = "ROC", trControl = control
)
print(glm_sofa_age_all)
toc()


# Generalized Additive Models ---------------------------------------------
tic("Training Generalized Additive Models")
print("-----GENERALIZED ADDITIVE MODELS-----")

# Helper to construct GAM formula from variable list
create_gam_formula <- function(data, vars, add_age = FALSE) {
  terms <- sapply(vars, function(var) {
    # Use a linear term if variable has fewer than 5 unique values, otherwise use a spline
    if (length(unique(data[[var]])) < 5) {
      return(var)
    } else {
      return(paste0("s(", var, ", k=3)"))
    }
  })
  
  full_formula_str <- paste(terms, collapse = " + ")
  
  if (add_age && "age_at_admission" %in% names(data)) {
    full_formula_str <- paste(full_formula_str, "+ age_at_admission")
  }
  
  return(as.formula(paste("output ~", full_formula_str)))
}

# Model 1: SOFA Variables Before
model_df <- data.frame(sofa_only_pre, output = factor(output, levels = c(0, 1)))
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)
for(i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])
  test_idx <- folds[[i]]
  
  # Create formula based on the actual training data for this fold
  gam_formula <- create_gam_formula(model_df[train_idx, ], vars = paste0(sofa_vars_list, "_pre"))
  
  gam_model <- gam(gam_formula, data = model_df[train_idx, ], family = binomial())
  models[[i]] <- gam_model
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  auc_values[i] <- as.numeric(pROC::roc(model_df$output[test_idx], probs)$auc)
}
best_idx <- which.max(auc_values)
gam_sofa_only_pre_auc <- auc_values[best_idx]
gam_sofa_only_pre <- models[[best_idx]]
print(summary(gam_sofa_only_pre))
cat("Best AUC (SOFA Pre):", gam_sofa_only_pre_auc, "\n")

# Model 2: SOFA Variables Before and After
model_df <- data.frame(sofa_only_all, output = factor(output, levels = c(0, 1)))
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)
for(i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])
  test_idx <- folds[[i]]
  
  gam_formula <- create_gam_formula(model_df[train_idx, ], vars = c(paste0(sofa_vars_list, "_pre"), paste0(sofa_vars_list, "_post")))
  
  gam_model <- gam(gam_formula, data = model_df[train_idx, ], family = binomial())
  models[[i]] <- gam_model
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  auc_values[i] <- as.numeric(pROC::roc(model_df$output[test_idx], probs)$auc)
}
best_idx <- which.max(auc_values)
gam_sofa_only_all_auc <- auc_values[best_idx]
gam_sofa_only_all <- models[[best_idx]]
print(summary(gam_sofa_only_all))
cat("Best AUC (SOFA All):", gam_sofa_only_all_auc, "\n")

# Model 3: SOFA Variables + Age Before
model_df <- data.frame(sofa_age_pre, output = factor(output, levels = c(0, 1)))
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)
for(i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])
  test_idx <- folds[[i]]

  gam_formula <- create_gam_formula(model_df[train_idx, ], vars = paste0(sofa_vars_list, "_pre"), add_age = TRUE)

  gam_model <- gam(gam_formula, data = model_df[train_idx, ], family = binomial())
  models[[i]] <- gam_model
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  auc_values[i] <- as.numeric(pROC::roc(model_df$output[test_idx], probs)$auc)
}
best_idx <- which.max(auc_values)
gam_sofa_age_pre_auc <- auc_values[best_idx]
gam_sofa_age_pre <- models[[best_idx]]
print(summary(gam_sofa_age_pre))
cat("Best AUC (SOFA + Age Pre):", gam_sofa_age_pre_auc, "\n")

# Model 4: SOFA Variables + Age Before and After
model_df <- data.frame(sofa_age_all, output = factor(output, levels = c(0, 1)))
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)
for(i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])
  test_idx <- folds[[i]]
  
  gam_formula <- create_gam_formula(model_df[train_idx, ], vars = c(paste0(sofa_vars_list, "_pre"), paste0(sofa_vars_list, "_post")), add_age = TRUE)
  
  gam_model <- gam(gam_formula, data = model_df[train_idx, ], family = binomial())
  models[[i]] <- gam_model
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  auc_values[i] <- as.numeric(pROC::roc(model_df$output[test_idx], probs)$auc)
}
best_idx <- which.max(auc_values)
gam_sofa_age_all_auc <- auc_values[best_idx]
gam_sofa_age_all <- models[[best_idx]]
print(summary(gam_sofa_age_all))
cat("Best AUC (SOFA + Age All):", gam_sofa_age_all_auc, "\n")
toc()


# Elastic Net Model Training ----------------------------------------------
tic("Training Elastic Net Models")
print("-----ELASTIC NET-----")
# Re-use trainControl from GLM section
glmnet_sofa_only_pre <- train(output ~ ., data = model_df_2, method = "glmnet", trControl = control, metric = "ROC", tuneLength = 10)
# Function to get best results, needs to be defined if not already in scope
get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  return(best_result)
}
print(get_best_result(glmnet_sofa_only_pre))

glmnet_sofa_only_all <- train(output ~ ., data = model_df_3, method = "glmnet", trControl = control, metric = "ROC", tuneLength = 10)
print(get_best_result(glmnet_sofa_only_all))

glmnet_sofa_age_pre <- train(output ~ ., data = model_df_4, method = "glmnet", trControl = control, metric = "ROC", tuneLength = 10)
print(get_best_result(glmnet_sofa_age_pre))

glmnet_sofa_age_all <- train(output ~ ., data = model_df_5, method = "glmnet", trControl = control, metric = "ROC", tuneLength = 10)
print(get_best_result(glmnet_sofa_age_all))
toc()


# LightGBM Training -------------------------------------------------------
tic("Training LightGBM Models")
print("-----LIGHTGBM-----")
y <- data$in_hospital_mortality
params <- list(objective = "binary", metric = "auc")

# Model 1
train_df_lgb1 <- lgb.Dataset(data = as.matrix(sofa_only_pre), label = y)
lightgbm_sofa_only_pre_results <- lgb.cv(params = params, data = train_df_lgb1, nfold = 5, nrounds = 100, stratified = TRUE, early_stopping_rounds = 10, verbose = 0)
lightgbm_sofa_only_pre <- lgb.train(params = params, data = train_df_lgb1, nrounds = lightgbm_sofa_only_pre_results$best_iter)

# Model 2
train_df_lgb2 <- lgb.Dataset(data = as.matrix(sofa_only_all), label = y)
lightgbm_sofa_only_all_results <- lgb.cv(params = params, data = train_df_lgb2, nfold = 5, nrounds = 100, stratified = TRUE, early_stopping_rounds = 10, verbose = 0)
lightgbm_sofa_only_all <- lgb.train(params = params, data = train_df_lgb2, nrounds = lightgbm_sofa_only_all_results$best_iter)

# Model 3
train_df_lgb3 <- lgb.Dataset(data = as.matrix(sofa_age_pre), label = y)
lightgbm_sofa_age_pre_results <- lgb.cv(params = params, data = train_df_lgb3, nfold = 5, nrounds = 100, stratified = TRUE, early_stopping_rounds = 10, verbose = 0)
lightgbm_sofa_age_pre <- lgb.train(params = params, data = train_df_lgb3, nrounds = lightgbm_sofa_age_pre_results$best_iter)

# Model 4
train_df_lgb4 <- lgb.Dataset(data = as.matrix(sofa_age_all), label = y)
lightgbm_sofa_age_all_results <- lgb.cv(params = params, data = train_df_lgb4, nfold = 5, nrounds = 100, stratified = TRUE, early_stopping_rounds = 10, verbose = 0)
lightgbm_sofa_age_all <- lgb.train(params = params, data = train_df_lgb4, nrounds = lightgbm_sofa_age_all_results$best_iter)

toc()


# Model Saving ------------------------------------------------------------
tic("Saving models and comparing results.")
models_path <- file.path(output_path, "exportable", "models")
dir.create(models_path, showWarnings = FALSE, recursive = TRUE)

# Remove training data from caret models before saving
glm_sofa_score$trainingData <- NULL
glm_sofa_only_pre$trainingData <- NULL
glm_sofa_only_all$trainingData <- NULL
glm_sofa_age_pre$trainingData <- NULL
glm_sofa_age_all$trainingData <- NULL
glmnet_sofa_only_pre$trainingData <- NULL
glmnet_sofa_only_all$trainingData <- NULL
glmnet_sofa_age_pre$trainingData <- NULL
glmnet_sofa_age_all$trainingData <- NULL

# Save models
saveRDS(glm_sofa_score, file.path(models_path, "glm_sofa_score.rds"))
saveRDS(glm_sofa_only_pre, file.path(models_path, "glm_sofa_only_pre.rds"))
saveRDS(glm_sofa_only_all, file.path(models_path, "glm_sofa_only_all.rds"))
saveRDS(glm_sofa_age_pre, file.path(models_path, "glm_sofa_age_pre.rds"))
saveRDS(glm_sofa_age_all, file.path(models_path, "glm_sofa_age_all.rds"))
saveRDS(gam_sofa_only_pre, file.path(models_path, "gam_sofa_only_pre.rds"))
saveRDS(gam_sofa_only_all, file.path(models_path, "gam_sofa_only_all.rds"))
saveRDS(gam_sofa_age_pre, file.path(models_path, "gam_sofa_age_pre.rds"))
saveRDS(gam_sofa_age_all, file.path(models_path, "gam_sofa_age_all.rds"))
saveRDS(glmnet_sofa_only_pre, file.path(models_path, "glmnet_sofa_only_pre.rds"))
saveRDS(glmnet_sofa_only_all, file.path(models_path, "glmnet_sofa_only_all.rds"))
saveRDS(glmnet_sofa_age_pre, file.path(models_path, "glmnet_sofa_age_pre.rds"))
saveRDS(glmnet_sofa_age_all, file.path(models_path, "glmnet_sofa_age_all.rds"))
lgb.save(lightgbm_sofa_only_pre, file.path(models_path, "lightgbm_sofa_only_pre.txt"))
lgb.save(lightgbm_sofa_only_all, file.path(models_path, "lightgbm_sofa_only_all.txt"))
lgb.save(lightgbm_sofa_age_pre, file.path(models_path, "lightgbm_sofa_age_pre.txt"))
lgb.save(lightgbm_sofa_age_all, file.path(models_path, "lightgbm_sofa_age_all.txt"))

# Save summary of results
model_summary <- list(
  glm_sofa_score = glm_sofa_score$results,
  glm_sofa_only_pre = glm_sofa_only_pre$results,
  glm_sofa_only_all = glm_sofa_only_all$results,
  glm_sofa_age_pre = glm_sofa_age_pre$results,
  glm_sofa_age_all = glm_sofa_age_all$results,
  gam_sofa_only_pre_auc = gam_sofa_only_pre_auc,
  gam_sofa_only_all_auc = gam_sofa_only_all_auc,
  gam_sofa_age_pre_auc = gam_sofa_age_pre_auc,
  gam_sofa_age_all_auc = gam_sofa_age_all_auc,
  glmnet_sofa_only_pre = get_best_result(glmnet_sofa_only_pre),
  glmnet_sofa_only_all = get_best_result(glmnet_sofa_only_all),
  glmnet_sofa_age_pre = get_best_result(glmnet_sofa_age_pre),
  glmnet_sofa_age_all = get_best_result(glmnet_sofa_age_all),
  lgbm_sofa_only_pre_auc = lightgbm_sofa_only_pre_results$best_score,
  lgbm_sofa_only_all_auc = lightgbm_sofa_only_all_results$best_score,
  lgbm_sofa_age_pre_auc = lightgbm_sofa_age_pre_results$best_score,
  lgbm_sofa_age_all_auc = lightgbm_sofa_age_all_results$best_score
)
saveRDS(model_summary, file.path(models_path, "model_summary.rds"))

toc()
print("Finished.")
