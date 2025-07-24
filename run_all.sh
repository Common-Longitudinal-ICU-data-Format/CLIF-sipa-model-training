#!/bin/bash

# This script runs the R scripts in the correct order.
# Make sure to run this script from the root directory of the project.
Rscript code/0a_respiratory_support_waterfall.R
Rscript code/01_cohort_identification.R
Rscript code/02_feature_set_processing.R
Rscript code/03_table1.R
Rscript code/04_model_training.R
