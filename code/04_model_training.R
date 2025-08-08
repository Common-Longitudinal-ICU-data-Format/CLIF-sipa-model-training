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

# Clear environment

# Load data
data <- read_parquet(paste0(output_path, "/sipa_features.parquet"))

tic("Creating input and output vectors.")
# Define output vector and input matrix. Input matrices:

# 1.  SOFA Score Only
# 2.  SOFA Variables before life support
# 3.  SOFA Variables before and after life support
# 4.  SOFA Variables + Age before life support
# 5.  SOFA Variables + Age before and after life support

# Create output and input matrices
output <- data$in_hospital_mortality

# Worst SOFA Score for the hospitalization
hosp_sofa_score <- data %>%
  select(sofa_score_pre, sofa_score_post) %>%
  rowwise() %>%
  mutate(worst_sofa_score = max(sofa_score_pre, sofa_score_post)) %>%
  select(worst_sofa_score)

# SOFA variables before life support
sofa_only_pre <- data %>%
  select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

# SOFA variables before and after life support
sofa_only_all <- data %>%
  select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre,p_f_post, s_f_post, platelets_post, bilirubin_post, map_post, gcs_post, creatinine_post, phenylephrine_post, norepinephrine_post, vasopressin_post, dopamine_post, dobutamine_post, milrinone_post, epinephrine_post, angiotensin_post) %>%
    mutate(across(everything(), ~replace_na(.x, 0)))

# SOFA variables and age before life support
sofa_age_pre <- data %>%
  select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre, age_at_admission) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

sofa_age_all <- data %>%
   select(p_f_pre, s_f_pre, platelets_pre, bilirubin_pre, map_pre, gcs_pre, creatinine_pre, phenylephrine_pre, norepinephrine_pre, vasopressin_pre, dopamine_pre, dobutamine_pre, milrinone_pre, epinephrine_pre, angiotensin_pre,p_f_post, s_f_post, platelets_post, bilirubin_post, map_post, gcs_post, creatinine_post, phenylephrine_post, norepinephrine_post, vasopressin_post, dopamine_post, dobutamine_post, milrinone_post, epinephrine_post, angiotensin_post, age_at_admission) %>%
    mutate(across(everything(), ~replace_na(.x, 0)))

toc()

# Analysis of Variables

# Run an elastic net model to see which features are selected across all the variables proposed above.

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

# Logistic Regression

tic()

print("-----LOGISTIC REGRESSION-----")

# For each of the four feature sets described above, run a logistic regression model with 5-fold cross validation.

## Model 1: SOFA Score

# Set up model dataframe
model_df <- data.frame(
  hosp_sofa_score = hosp_sofa_score,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# Set up 5-fold cross-validation with AUC as the metric
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final")

# Train
glm_sofa_score <- train(
  output ~ worst_sofa_score,
  data = model_df,
  method = "glm",
  family = binomial(link = "logit"),
  metric = "ROC",
  trControl = control
)

# View results
print(glm_sofa_score)

## Model 2: SOFA Variables Before Life Support Initiation 

# Set up model dataframe
model_df <- data.frame(
  sofa_only_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# Set up 5-fold cross-validation with AUC as the metric
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final")

# Train
glm_sofa_only_pre <- train(
  output ~ .,
  data = model_df,
  method = "glm",
  family = binomial(link = "logit"),
  metric = "ROC",
  trControl = control)

print(glm_sofa_only_pre)

## Model 3: SOFA Variables Before and After Life Support Initiation 

set.seed(42) # The meaning of life

# Set up model dataframe
model_df <- data.frame(
  sofa_only_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# Set up 5-fold cross-validation with AUC as the metric
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final")

# Train
glm_sofa_only_all <- train(
  output ~ .,
  data = model_df,
  method = "glm",
  family = binomial(link = "logit"),
  metric = "ROC",
  trControl = control)

print(glm_sofa_only_all)

## Model 4: SOFA Variables + Age Before Life Support Initiation 

# Set up model dataframe
model_df <- data.frame(
  sofa_age_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# Set up 5-fold cross-validation with AUC as the metric
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final")

# Train
glm_sofa_age_pre <- train(
  output ~ .,
  data = model_df,
  method = "glm",
  family = binomial(link = "logit"),
  metric = "ROC",
  trControl = control)

print(glm_sofa_age_pre)

## Model 5: SOFA Variables + Age Before and After Life Support Initiation
# Set up model dataframe
model_df <- data.frame(
  sofa_age_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))
# Set up 5-fold cross-validation with AUC as the metric
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final")
# Train
glm_sofa_age_all <- train(
  output ~ .,
  data = model_df,
  method = "glm",
  family = binomial(link = "logit"),
  metric = "ROC",
  trControl = control)
print(glm_sofa_age_all)

toc()

# Generalized Additive Models

print("-----GENERALIZED ADDITIVE MODELS-----")

tic()
# Constructs a set of GAMs for each of the four feature sets. Uses splines. 

# 1.  SOFA Variables before life support
# 2.  SOFA Variables before and after life support
# 3.  SOFA Variables + Age before life support
# 4.  SOFA Variables + Age before and after life support

## Model 1: SOFA Variables Before Life Support Initiation 

# Set up model dataframe
model_df <- data.frame(
  sofa_only_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# 5-fold cross validation
folds <- createFolds(model_df$output, k = 5, list = TRUE)

auc_values <- numeric(5)
models <- vector("list", 5)

# Drops angiotensin
for(i in seq_along(folds)) {
  train_idx <- setdiff(seq_len(nrow(model_df)), folds[[i]])
  test_idx <- folds[[i]]
  
  gam_model <- gam(
    output ~ s(p_f_pre, k=3) + s(s_f_pre, k=3) + s(platelets_pre, k=3) + s(bilirubin_pre, k=3) +
             s(map_pre, k=3) + s(gcs_pre, k=3) + s(creatinine_pre, k=3) + s(phenylephrine_pre, k=3) +
             s(norepinephrine_pre, k=3) + s(vasopressin_pre, k=3) + s(dopamine_pre, k=3) +
             s(dobutamine_pre, k=3) + s(milrinone_pre, k=3) + s(epinephrine_pre, k=3),
    data = model_df[train_idx, ],
    family = binomial()
  )
  
  # Store the model
  models[[i]] <- gam_model
  
  # Predict probabilities on test fold
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  # Calculate AUC
  auc <- tryCatch({
    pROC::roc(model_df$output[test_idx], probs)$auc
  }, error = function(e) NA)
  auc_values[i] <- auc
}

# Find the best AUC and corresponding model
best_idx <- which.max(auc_values)
gam_sofa_only_pre_auc <- auc_values[best_idx]
gam_sofa_only_pre <- models[[best_idx]]

print(gam_sofa_only_pre)
cat("Best AUC:", gam_sofa_only_pre_auc, "\n")

## Model 2: SOFA Variables Before and After Life Support Initiation 

# Set up model dataframe
model_df <- data.frame(
  sofa_only_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))
# 5-fold cross validation
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)

# Drops angiotensin
for(i in seq_along(folds)) {
  train_idx <- setdiff(seq_len(nrow(model_df)), folds[[i]])
  test_idx <- folds[[i]]
  
  gam_model <- gam(
    output ~ s(p_f_pre, k=3) + s(s_f_pre, k=3) + s(platelets_pre, k=3) + s(bilirubin_pre, k=3) +
             s(map_pre, k=3) + s(gcs_pre, k=3) + s(creatinine_pre, k=3) + s(phenylephrine_pre, k=3) +
             s(norepinephrine_pre, k=3) + s(vasopressin_pre, k=3) + s(dopamine_pre, k=3) +
             s(dobutamine_pre, k=3) + s(milrinone_pre, k=3) + s(epinephrine_pre, k=3) +
             s(p_f_post, k=3) + s(s_f_post, k=3) + s(platelets_post, k=3) + s(bilirubin_post, k=3) +
             s(map_post, k=3) + s(gcs_post, k=3) + s(creatinine_post, k=3) +
             s(phenylephrine_post, k=3) + s(norepinephrine_post, k=3) +
             s(vasopressin_post, k=3) + s(dopamine_post, k=3) +
             s(dobutamine_post, k=3) + s(milrinone_post, k=3) +
             s(epinephrine_post, k=3),
    data = model_df[train_idx, ],
    family = binomial()
  )
  
  # Store the model
  models[[i]] <- gam_model
  
  # Predict probabilities on test fold
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  # Calculate AUC
  auc <- tryCatch({
    pROC::roc(model_df$output[test_idx], probs)$auc
  }, error = function(e) NA)
  auc_values[i] <- auc
}
# Find the best AUC and corresponding model
best_idx <- which.max(auc_values)
gam_sofa_only_all_auc <- auc_values[best_idx]
gam_sofa_only_all <- models[[best_idx]]

# Output
print(gam_sofa_only_all)
cat("Best AUC:", gam_sofa_only_all_auc, "\n")

## Model 3: SOFA Variables + Age Before Life Support Initiation

# Set up model dataframe
model_df <- data.frame(
  sofa_age_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# 5-fold cross validation
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)

# Drops angiotensin
for(i in seq_along(folds)) {
  train_idx <- setdiff(seq_len(nrow(model_df)), folds[[i]])
  test_idx <- folds[[i]]
  
  gam_model <- gam(
    output ~ s(p_f_pre, k=3) + s(s_f_pre, k=3) + s(platelets_pre, k=3) + s(bilirubin_pre, k=3) +
             s(map_pre, k=3) + s(gcs_pre, k=3) + s(creatinine_pre, k=3) + s(phenylephrine_pre, k=3) +
             s(norepinephrine_pre, k=3) + s(vasopressin_pre, k=3) + s(dopamine_pre, k=3) +
             s(dobutamine_pre, k=3) + s(milrinone_pre, k=3) + s(epinephrine_pre, k=3) +
             age_at_admission,
    data = model_df[train_idx, ],
    family = binomial()
  )
  
  # Store the model
  models[[i]] <- gam_model
  
  # Predict probabilities on test fold
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  # Calculate AUC
  auc <- tryCatch({
    pROC::roc(model_df$output[test_idx], probs)$auc
  }, error = function(e) NA)
  auc_values[i] <- auc
}
# Find the best AUC and corresponding model
best_idx <- which.max(auc_values)
gam_sofa_age_pre_auc <- auc_values[best_idx]
gam_sofa_age_pre <- models[[best_idx]]

# Output
print(gam_sofa_age_pre)
cat("Best AUC:", gam_sofa_age_pre_auc, "\n")

## Model 4: SOFA Variables + Age Before and After Life Support Initiation

# Set up model dataframe
model_df <- data.frame(
  sofa_age_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))
# 5-fold cross validation
folds <- createFolds(model_df$output, k = 5, list = TRUE)
auc_values <- numeric(5)
models <- vector("list", 5)

# Drops angiotensin from pre and post
for(i in seq_along(folds)) {
  train_idx <- setdiff(seq_len(nrow(model_df)), folds[[i]])
  test_idx <- folds[[i]]
  
  gam_model <- gam(
    output ~ s(p_f_pre, k=3) + s(s_f_pre, k=3) + s(platelets_pre, k=3) + s(bilirubin_pre, k=3) +
             s(map_pre, k=3) + s(gcs_pre, k=3) + s(creatinine_pre, k=3) + s(phenylephrine_pre, k=3) +
             s(norepinephrine_pre, k=3) + s(vasopressin_pre, k=3) + s(dopamine_pre, k=3) +
             s(dobutamine_pre, k=3) + s(milrinone_pre, k=3) + s(epinephrine_pre, k=3) +
             age_at_admission +
             s(p_f_post, k=3) + s(s_f_post, k=3) + s(platelets_post, k=3) + s(bilirubin_post, k=3) +
             s(map_post, k=3) + s(gcs_post, k=3) + s(creatinine_post, k=3) +
             s(phenylephrine_post, k=3) + s(norepinephrine_post, k=3) +
             s(vasopressin_post, k=3) + s(dopamine_post, k=3) +
             s(dobutamine_post, k=3) + s(milrinone_post, k=3) +
             s(epinephrine_post, k=3),
    data = model_df[train_idx, ],
    family = binomial()
  )
  
  # Store the model
  models[[i]] <- gam_model
  
  # Predict probabilities on test fold
  probs <- predict(gam_model, newdata = model_df[test_idx, ], type = "response")
  # Calculate AUC
  auc <- tryCatch({
    pROC::roc(model_df$output[test_idx], probs)$auc
  }, error = function(e) NA)
  auc_values[i] <- auc
}
# Find the best AUC and corresponding model
best_idx <- which.max(auc_values)
gam_sofa_age_all_auc <- auc_values[best_idx]
gam_sofa_age_all <- models[[best_idx]]

# Output
print(gam_sofa_age_all)
cat("Best AUC:", gam_sofa_age_all_auc, "\n")

toc()

# Elastic Net

print("-----ELASTIC NET-----")

## Model 1: SOFA Variables Before Life Support Initiation 

tic()

# Function that returns the best result for elastic net
get_best_result = function(caret_fit) {
  # Find the best row in caret_fit$results
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  
  # Get best alpha and lambda
  best_alpha = caret_fit$bestTune$alpha
  best_lambda = caret_fit$bestTune$lambda
  
  # Get coefficients from the final glmnet model at best lambda
  coefs = as.matrix(coef(caret_fit$finalModel, s = best_lambda))
  
  # Return as a named list
  list(
    AUC = best_result$ROC,
    lambda = best_lambda,
    alpha = best_alpha,
    coefficients = coefs
  )
}

# Define model matrix for elastic net training
model_df <- data.frame(
  sofa_only_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead"))
)

# 5-fold cross validation for training
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
)

glmnet_sofa_only_pre <- train(
  output ~.,
  data = model_df,
  method = "glmnet",
  trControl = control,
  metric = "ROC",
  tuneLength = 10
)

get_best_result(glmnet_sofa_only_pre)

## Model 2: SOFA Variables Before and After Life Support Initiation

# Define model matrix for elastic net training
model_df <- data.frame(
  sofa_only_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# 5-fold cross validation for training
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

glmnet_sofa_only_all <- train(
  output ~.,
  data = model_df,
  method = "glmnet",
  trControl = control,
  metric = "ROC",
  tuneLength = 10)

get_best_result(glmnet_sofa_only_all)

## Model 3: SOFA Variables + Age Before Life Support Initiation

# Define model matrix for elastic net training
model_df <- data.frame(
  sofa_age_pre,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# 5-fold cross validation for training
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

# Train
glmnet_sofa_age_pre <- train(
  output ~.,
  data = model_df,
  method = "glmnet",
  trControl = control,
  metric = "ROC",
  tuneLength = 10)

get_best_result(glmnet_sofa_age_pre)

## Model 4: SOFA Variables + Age Before and After Life Support Initiation

# Define model matrix for elastic net training
model_df <- data.frame(
  sofa_age_all,
  output = factor(output, levels = c(0, 1), labels = c("Alive", "Dead")))

# 5-fold cross validation for training
control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary)

glmnet_sofa_age_all <- train(
  output ~.,
  data = model_df,
  method = "glmnet",
  trControl = control,
  metric = "ROC",
  tuneLength = 10)

get_best_result(glmnet_sofa_age_all)

toc()

# LightGBM

tic()

print("-----LIGHTGBM-----")

## Model 1: SOFA Variables Before Life Support Initiation

# Convert output to a 1/0 vector
y <- data$in_hospital_mortality

# Train LightGBM model with SOFA variables before life support initiation
train_df <- lgb.Dataset(data = as.matrix(sofa_only_pre), label = y)

params <- list(
  objective = "binary",
  metric = "auc"
)

lightgbm_sofa_only_pre_results <- lgb.cv(
  params = params,
  data = train_df,
  nfold = 5,
  nrounds = 100,
  stratified = TRUE,
  eval = "auc",
  early_stopping_rounds = 10,
  verbose = 0)

best_iter <- lightgbm_sofa_only_pre_results$best_iter
lightgbm_sofa_only_pre_results$best_score

# Train
lightgbm_sofa_only_pre <- lgb.train(
  params = params,
  data = train_df,
  nrounds = best_iter)

## Model 2: SOFA Variables Before and After Life Support Initiation

# Convert output to a 1/0 vector
y <- data$in_hospital_mortality

# Train LightGBM model with SOFA variables before  and after life support initiation
train_df <- lgb.Dataset(data = as.matrix(sofa_only_all), label = y)

params <- list(
  objective = "binary",
  metric = "auc")

lightgbm_sofa_only_all_results <- lgb.cv(
  params = params,
  data = train_df,
  nfold = 5,
  nrounds = 100,
  stratified = TRUE,
  eval = "auc",
  early_stopping_rounds = 10,
  verbose = 0)

best_iter <- lightgbm_sofa_only_all_results$best_iter
lightgbm_sofa_only_all_results$best_score

# Train
lightgbm_sofa_only_all <- lgb.train(
  params = params,
  data = train_df,
  nrounds = best_iter)

## Model 3: SOFA Variables + Age Before Life Support Initiation

# Convert output to a 1/0 vector
y <- data$in_hospital_mortality
# Train LightGBM model with SOFA variables and age before life support initiation
train_df <- lgb.Dataset(data = as.matrix(sofa_age_pre), label = y)

params <- list(
  objective = "binary",
  metric = "auc")

lightgbm_sofa_age_pre_results <- lgb.cv(
  params = params,
  data = train_df,
  nfold = 5,
  nrounds = 100,
  stratified = TRUE,
  eval = "auc",
  early_stopping_rounds = 10,
  verbose = 0)

best_iter <- lightgbm_sofa_age_pre_results$best_iter
lightgbm_sofa_age_pre_results$best_score
# Train
lightgbm_sofa_age_pre <- lgb.train(
  params = params,
  data = train_df,
  nrounds = best_iter)

## Model 4: SOFA Variables + Age Before and After Life Support Initiation

# Convert output to a 1/0 vector
y <- data$in_hospital_mortality

# Train LightGBM model with SOFA variables and age before life support initiation
train_df <- lgb.Dataset(data = as.matrix(sofa_age_all), label = y)

params <- list(
  objective = "binary",
  metric = "auc")

lightgbm_sofa_age_all_results <- lgb.cv(
  params = params,
  data = train_df,
  nfold = 5,
  nrounds = 100,
  stratified = TRUE,
  eval = "auc",
  early_stopping_rounds = 10,
  verbose = 0)

best_iter <- lightgbm_sofa_age_all_results$best_iter
lightgbm_sofa_age_all_results$best_score

# Train
lightgbm_sofa_age_all <- lgb.train(
  params = params,
  data = train_df,
  nrounds = best_iter)

toc()

# Save models

tic("Saving models.")
# Create a directory to save the models
models_path <- file.path(output_path, "models")
dir.create(models_path, showWarnings = FALSE)

# Remove training data from caret models before saving
glm_sofa_score$trainingData <- NULL
glm_sofa_only_pre$trainingData <- NULL
glm_sofa_only_all$trainingData <- NULL
glm_sofa_age_pre$trainingData <- NULL
glm_sofa_age_all$trainingData <- NULL

glm_sofa_score$finalModel$data <- NULL
glm_sofa_only_pre$finalModel$data <- NULL
glm_sofa_only_all$finalModel$data <- NULL
glm_sofa_age_pre$finalModel$data <- NULL
glm_sofa_age_all$finalModel$data <- NULL

glm_sofa_score$finalModel$model <- NULL
glm_sofa_only_pre$finalModel$model <- NULL
glm_sofa_only_all$finalModel$model <- NULL
glm_sofa_age_pre$finalModel$model <- NULL
glm_sofa_age_all$finalModel$model <- NULL

glmnet_sofa_only_pre$trainingData <- NULL
glmnet_sofa_only_all$trainingData <- NULL
glmnet_sofa_age_pre$trainingData <- NULL
glmnet_sofa_age_all$trainingData <- NULL

# Remove training data from GAM models before saving
gam_sofa_only_pre$data <- NULL
gam_sofa_only_all$data <- NULL
gam_sofa_age_pre$data <- NULL
gam_sofa_age_all$data <- NULL

# Save GLM models
saveRDS(glm_sofa_score, file.path(models_path, "glm_sofa_score.rds"))
saveRDS(glm_sofa_only_pre, file.path(models_path, "glm_sofa_only_pre.rds"))
saveRDS(glm_sofa_only_all, file.path(models_path, "glm_sofa_only_all.rds"))
saveRDS(glm_sofa_age_pre, file.path(models_path, "glm_sofa_age_pre.rds"))
saveRDS(glm_sofa_age_all, file.path(models_path, "glm_sofa_age_all.rds"))

# Save GAM models
saveRDS(gam_sofa_only_pre, file.path(models_path, "gam_sofa_only_pre.rds"))
saveRDS(gam_sofa_only_all, file.path(models_path, "gam_sofa_only_all.rds"))
saveRDS(gam_sofa_age_pre, file.path(models_path, "gam_sofa_age_pre.rds"))
saveRDS(gam_sofa_age_all, file.path(models_path, "gam_sofa_age_all.rds"))

# Save GLMNET models
saveRDS(glmnet_sofa_only_pre, file.path(models_path, "glmnet_sofa_only_pre.rds"))
saveRDS(glmnet_sofa_only_all, file.path(models_path, "glmnet_sofa_only_all.rds"))
saveRDS(glmnet_sofa_age_pre, file.path(models_path, "glmnet_sofa_age_pre.rds"))
saveRDS(glmnet_sofa_age_all, file.path(models_path, "glmnet_sofa_age_all.rds"))

# Save LightGBM models
lgb.save(lightgbm_sofa_only_pre, file.path(models_path, "lightgbm_sofa_only_pre.txt"))
lgb.save(lightgbm_sofa_only_all, file.path(models_path, "lightgbm_sofa_only_all.txt"))
lgb.save(lightgbm_sofa_age_pre, file.path(models_path, "lightgbm_sofa_age_pre.txt"))
lgb.save(lightgbm_sofa_age_all, file.path(models_path, "lightgbm_sofa_age_all.txt"))

# Also save the AUCs and model results in a single file for reference
save(
  # GLM models' results are in the objects
  glm_sofa_score, glm_sofa_only_pre, glm_sofa_only_all, glm_sofa_age_pre, glm_sofa_age_all,
  # GAM models' AUCs
  gam_sofa_only_pre_auc, gam_sofa_only_all_auc, gam_sofa_age_pre_auc, gam_sofa_age_all_auc,
  # GLMNET models' results are in the objects
  glmnet_sofa_only_pre, glmnet_sofa_only_all, glmnet_sofa_age_pre, glmnet_sofa_age_all,
  # LightGBM results
  lightgbm_sofa_only_pre_results, lightgbm_sofa_only_all_results, lightgbm_sofa_age_pre_results, lightgbm_sofa_age_all_results,
  file = file.path(models_path, "trained_models_summary.RData"))

# Compare models and select the best one

# Create a list of all models and their AUCs
models_and_aucs <- list(
  list(model = glm_sofa_score, auc = max(glm_sofa_score$results$ROC), name = "GLM SOFA Score"),
  list(model = glm_sofa_only_pre, auc = max(glm_sofa_only_pre$results$ROC), name = "GLM SOFA Variables Before"),
  list(model = glm_sofa_only_all, auc = max(glm_sofa_only_all$results$ROC), name = "GLM SOFA Variables Before and After"),
  list(model = glm_sofa_age_pre, auc = max(glm_sofa_age_pre$results$ROC), name = "GLM SOFA Variables + Age Before"),
  list(model = glm_sofa_age_all, auc = max(glm_sofa_age_all$results$ROC), name = "GLM SOFA Variables + Age Before and After"),
  list(model = gam_sofa_only_pre, auc = gam_sofa_only_pre_auc, name = "GAM SOFA Variables Before"),
  list(model = gam_sofa_only_all, auc = gam_sofa_only_all_auc, name = "GAM SOFA Variables Before and After"),
  list(model = gam_sofa_age_pre, auc = gam_sofa_age_pre_auc, name = "GAM SOFA Variables + Age Before"),
  list(model = gam_sofa_age_all, auc = gam_sofa_age_all_auc, name = "GAM SOFA Variables + Age Before and After"),
  list(model = glmnet_sofa_only_pre, auc = max(glmnet_sofa_only_pre$results$ROC), name = "GLMNET SOFA Variables Before"),
  list(model = glmnet_sofa_only_all, auc = max(glmnet_sofa_only_all$results$ROC), name = "GLMNET SOFA Variables Before and After"),
  list(model = glmnet_sofa_age_pre, auc = max(glmnet_sofa_age_pre$results$ROC), name = "GLMNET SOFA Variables + Age Before"),
  list(model = glmnet_sofa_age_all, auc = max(glmnet_sofa_age_all$results$ROC), name = "GLMNET SOFA Variables + Age Before and After"),
  list(model = lightgbm_sofa_only_pre, auc = lightgbm_sofa_only_pre_results$best_score, name = "LightGBM SOFA Variables Before"),
  list(model = lightgbm_sofa_only_all, auc = lightgbm_sofa_only_all_results$best_score, name = "LightGBM SOFA Variables Before and After"),
  list(model = lightgbm_sofa_age_pre, auc = lightgbm_sofa_age_pre_results$best_score, name = "LightGBM SOFA Variables + Age Before"),
  list(model = lightgbm_sofa_age_all, auc = lightgbm_sofa_age_all_results$best_score, name = "LightGBM SOFA Variables + Age Before and After")
)

toc()

# Find the best model
best_model_info <- models_and_aucs[[which.max(sapply(models_and_aucs, function(x) x$auc))]]

# Print the best model's information
cat("Best Model:", best_model_info$name, "\n")
cat("Best AUC:", best_model_info$auc, "\n")

# Save the best model
if (grepl("LightGBM", best_model_info$name)) {
  lgb.save(best_model_info$model, file.path(output_path, "best_model.txt"))
} else {
  saveRDS(best_model_info$model, file.path(output_path, "best_model.rds"))
}

# Print hyperparameters if they exist
if (!is.null(best_model_info$model$bestTune)) {
  cat("Hyperparameters:\n")
  print(best_model_info$model$bestTune)
}

print("Finished. Check the output folder for best_model.")
