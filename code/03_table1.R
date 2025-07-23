# Load necessary libraries
library(arrow)
library(tidyverse)
library(stringr)
library(kableExtra)
library(stringr)
library(htmltools)

# Clear env
rm(list = ls())

# Load data
source("utils/config.R")
output_path <- config$output_path
site_name <- config$site_name

data <- read_parquet(paste0(output_path, "/sipa_features.parquet"))

# Correct column names by removing trailing underscore
setnames(data, names(data), sub("__", "_", names(data)))

# Average _pre and _post columns
pre_cols <- names(data)[grepl("_pre$", names(data))]
post_cols <- names(data)[grepl("_post$", names(data))]

base_vars <- sub("_pre$", "", pre_cols)

for (var in base_vars) {
  pre_col <- paste0(var, "_pre")
  post_col <- paste0(var, "_post")
  if (pre_col %in% names(data) && post_col %in% names(data)) {
    data <- data %>%
      mutate(!!var := rowMeans(select(., all_of(c(pre_col, post_col))), na.rm = TRUE))}}

# This script outputs a Table 1 as an HTML.

# Helper for median (IQR)
med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  paste0(formatC(median(x), digits=1, format="f"),
         " (", 
         formatC(quantile(x, 0.25), digits=1, format="f"), ", ",
         formatC(quantile(x, 0.75), digits=1, format="f"), ")")}

# Summarize to hospitalization_id (ICU encounter) level
data_summary <- data %>%
  group_by(hospitalization_id) %>%
  summarise(
    patient_id = first(patient_id),
    in_hospital_mortality = first(in_hospital_mortality),
    sex_category = first(sex_category),
    race_category = first(race_category),
    ethnicity_category = first(ethnicity_category),
    age_at_admission = first(age_at_admission),
    p_f = median(p_f, na.rm = TRUE),
    s_f = median(s_f, na.rm = TRUE),
    platelets = median(platelets, na.rm = TRUE),
    bilirubin = median(bilirubin, na.rm = TRUE),
    map = median(map, na.rm = TRUE),
    dobutamine = median(dobutamine, na.rm = TRUE),
    dopamine = median(dopamine, na.rm = TRUE),
    norepinephrine = median(norepinephrine, na.rm = TRUE),
    phenylephrine = median(phenylephrine, na.rm = TRUE),
    epinephrine = median(epinephrine, na.rm = TRUE),
    gcs = median(gcs, na.rm = TRUE),
    creatinine = median(creatinine, na.rm = TRUE),
    sofa_score = median(sofa_score, na.rm = TRUE),
  ) %>%
  ungroup() %>%
  mutate(Survival = ifelse(in_hospital_mortality == 0, "Survivor", "Non-Survivor"))

# Patient-level summary for demographics and survival
patient_summary <- data_summary %>%
  group_by(patient_id) %>%
  summarise(
    race_category = first(race_category),
    ethnicity_category = first(ethnicity_category),
    sex_category = first(sex_category),
    Survival = ifelse(any(in_hospital_mortality == 1), "Non-Survivor", "Survivor")
  )

# Define new variable display names
var_display_names <- list(
  "p_f" = "PaO2/FiO2 (median, IQR)",
  "s_f" = "SpO2/FiO2 (median, IQR)",
  "platelets" = "Platelet Count (median, IQR)",
  "bilirubin" = "Bilirubin (median, IQR)",
  "map" = "Mean Arterial Pressure (median, IQR)",
  "dobutamine" = "Dobutamine (median, IQR)",
  "dopamine" = "Dopamine (median, IQR)",
  "norepinephrine" = "Norepinephrine (median, IQR)",
  "phenylephrine" = "Phenylephrine (median, IQR)",
  "epinephrine" = "Epinephrine (median, IQR)",
  "gcs" = "Glasgow Coma Scale (median, IQR)",
  "creatinine" = "Creatinine (median, IQR)",
  "sofa_score" = "SOFA Score (median, IQR)")

table1 <- list()

# Number of ICU Encounters (N)
icu_encounters <- data_summary %>%
  count(Survival)
table1[["ICU Encounters (N)"]] <- setNames(icu_encounters$n, icu_encounters$Survival)

# Number of Patients (N)
n_patients <- patient_summary %>%
  count(Survival)
table1[["Patients (N)"]] <- setNames(n_patients$n, n_patients$Survival)

# Sex (N, % female) at patient level
sex_pat <- patient_summary %>%
  group_by(Survival) %>%
  summarise(
    N_female = sum(sex_category == "Female", na.rm = TRUE),
    Percent_female = 100 * mean(sex_category == "Female", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Sex = paste0(N_female, " (", formatC(Percent_female, digits=1, format="f"), "%)"))
table1[["Sex (N, % Female)"]] <- setNames(sex_pat$Sex, sex_pat$Survival)

# Age (median, IQR) -- at encounter level, placed after Sex
age <- data_summary %>%
  group_by(Survival) %>%
  summarise(stat = med_iqr(age_at_admission), .groups = "drop")
table1[["Age (median, IQR)"]] <- setNames(
  c(ifelse("Survivor" %in% age$Survival, age$stat[age$Survival == "Survivor"], ""),
    ifelse("Non-Survivor" %in% age$Survival, age$stat[age$Survival == "Non-Survivor"], "")),
  c("Survivor", "Non-Survivor")
)

# Blank Race header row
table1[["Race (N, %)"]] <- c("", "")

race_levels <- c("White", "Asian", "Native Hawaiian or Other Pacific Islander", 
                 "Black or African American", "Unknown", "Other", "American Indian or Alaska Native")
race_tab <- patient_summary %>%
  filter(race_category %in% race_levels) %>%
  group_by(Survival, race_category) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Survival) %>%
  mutate(Total = sum(N),
         Percent = 100 * N / Total,
         Race = paste0(N, " (", formatC(Percent, digits=1, format="f"), "%)")) %>%
  ungroup() %>%
  select(Survival, race_category, Race) %>%
  pivot_wider(names_from = Survival, values_from = Race) %>%
  rename(Variable = race_category)

for (i in seq_along(race_levels)) {
  row <- race_tab[race_tab$Variable == race_levels[i],]
  # Just the formatted value, not the variable name
  table1[[paste0("  ", race_levels[i])]] <- setNames(
    c(ifelse("Survivor" %in% names(row), row[["Survivor"]], ""),
      ifelse("Non-Survivor" %in% names(row), row[["Non-Survivor"]], "")),
    c("Survivor", "Non-Survivor")
  )
}

# Blank Ethnicity header row
table1[["Ethnicity (N, %)"]] <- c("", "")
ethnicity_levels <- c("Hispanic", "Non-Hispanic", "Unknown")
eth_tab <- patient_summary %>%
  filter(ethnicity_category %in% ethnicity_levels) %>%
  group_by(Survival, ethnicity_category) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Survival) %>%
  mutate(Total = sum(N),
         Percent = 100 * N / Total,
         Eth = paste0(N, " (", formatC(Percent, digits=1, format="f"), "%)")) %>%
  ungroup() %>%
  select(Survival, ethnicity_category, Eth) %>%
  pivot_wider(names_from = Survival, values_from = Eth) %>%
  rename(Variable = ethnicity_category)

for (i in seq_along(ethnicity_levels)) {
  row <- eth_tab[eth_tab$Variable == ethnicity_levels[i],]
  table1[[paste0("  ", ethnicity_levels[i])]] <- setNames(
    c(ifelse("Survivor" %in% names(row), row[["Survivor"]], ""),
      ifelse("Non-Survivor" %in% names(row), row[["Non-Survivor"]], "")),
    c("Survivor", "Non-Survivor")
  )
}

# Numeric variables (median, IQR) with new names (still at encounter level)
numeric_vars <- c("p_f", "s_f", "platelets", "bilirubin", "map",
                  "dobutamine", "dopamine", "norepinephrine", 
                  "phenylephrine", "epinephrine", "gcs", "creatinine", "sofa_score")
for (v in numeric_vars) {
  res <- data_summary %>%
    group_by(Survival) %>%
    summarise(stat = med_iqr(.data[[v]]), .groups = "drop")
  display_name <- var_display_names[[v]]
  table1[[display_name]] <- setNames(
    c(ifelse("Survivor" %in% res$Survival, res$stat[res$Survival == "Survivor"], ""),
      ifelse("Non-Survivor" %in% res$Survival, res$stat[res$Survival == "Non-Survivor"], "")),
    c("Survivor", "Non-Survivor"))}

# Build the table in the desired order
desired_order <- c(
  "ICU Encounters (N)",
  "Patients (N)",
  "Sex (N, % Female)",
  "Age (median, IQR)",
  "Race (N, %)",
  paste0("  ", race_levels),
  "Ethnicity (N, %)",
  paste0("  ", ethnicity_levels),
  unname(unlist(var_display_names)))

# Defensive extraction: always return blank if not found or not named
table1_df <- tibble::tibble(
  Variable = desired_order,
  Survivor = sapply(desired_order, function(x) {
    if (!is.null(table1[[x]]) && "Survivor" %in% names(table1[[x]])) table1[[x]][["Survivor"]] else
      if (!is.null(table1[[x]]) && length(table1[[x]]) == 2 && all(table1[[x]] == "")) "" else ""}),
  `Non-Survivor` = sapply(desired_order, function(x) {
    if (!is.null(table1[[x]]) && "Non-Survivor" %in% names(table1[[x]])) table1[[x]][["Non-Survivor"]] else
      if (!is.null(table1[[x]]) && length(table1[[x]]) == 2 && all(table1[[x]] == "")) "" else ""}))

# Export as a separate HTML file
table1_html <- table1_df %>%
  mutate(Variable = ifelse(stringr::str_detect(Variable, "^  "), 
                           paste0("&nbsp;&nbsp;", stringr::str_trim(Variable)), 
                           Variable)) %>%
  kable("html", escape = FALSE, col.names = c("Variable", "Survivors", "Non-Survivors")) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
    full_width = FALSE) %>%
  scroll_box(width = "100%")

# Save to HTML file
save_html(table1_html, file = "/Users/cdiaz/Desktop/SRP/SRP SOFA/table1.html")