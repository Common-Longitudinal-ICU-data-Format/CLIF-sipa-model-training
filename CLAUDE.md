# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the CLIF SIPA project - a research codebase for developing an improved in-hospital mortality predictor based on SOFA (Sequential Organ Failure Assessment) score. The project analyzes critical care data from the CLIF consortium to build a more robust statistical model for mortality prediction.

## Essential Development Commands

### Environment Setup
- `Rscript code/00_renv_restore.R` - Set up R environment with required packages
- Configure `config/config.json` by copying from `config_template.json` and updating with site-specific settings

### Running the Analysis Pipeline
Execute scripts in this exact order:
```bash
./run_all.sh
```

Or run individually:
```bash
# Run from project root directory
Rscript code/0a_respiratory_support_waterfall.R  # Respiratory device categorization
Rscript code/01_cohort_identification.R          # Create study cohort
Rscript code/02_feature_set_processing.R         # Process features for modeling
Rscript code/03_table1.R                         # Generate descriptive statistics
Rscript code/04_model_training.R                 # Train mortality prediction models
```

## Project Architecture

### Data Flow
1. **Raw CLIF Tables**: Six required tables (patient, hospitalization, vitals, labs, medication_admin_continuous, respiratory_support)
2. **Cohort Creation**: Adult ICU patients on life support ≥6 hours (2018-2023)
3. **Feature Engineering**: SOFA-based variables with temporal aggregation
4. **Model Training**: Classification models for in-hospital mortality

### Key Components
- `utils/config.R`: Configuration loader using `config.json`
- `utils/outlier_handler.R`: Data quality and outlier detection
- `lookup-tables/`: Device mapping tables for respiratory support categorization
- `outlier-thresholds/`: Predefined thresholds for data validation

### Database Strategy
Uses DuckDB for in-memory analytics with support for multiple file formats:
- Parquet (preferred)
- CSV
- FST

### R Environment
- Uses `renv` for reproducible package management
- Key packages: tidyverse, arrow, duckdb, gtsummary, data.table
- R version: 4.4.2

## Configuration Requirements

Before running any analysis:
1. Copy `config/config_template.json` to `config/config.json`
2. Update required fields:
   - `site_name`: Institution identifier
   - `tables_path`: Path to CLIF data tables
   - `file_type`: Data format (.parquet, .csv, or .fst)
   - `output_path`: Results output directory

## Data Requirements

Must have CLIF-formatted tables with specific fields:
- **patient**: demographics and identifiers
- **hospitalization**: admission/discharge times, age, mortality
- **vitals**: MAP, SpO2, weight measurements
- **labs**: bilirubin, creatinine, platelets, PaO2
- **medication_admin_continuous**: vasoactive medications
- **respiratory_support**: ventilation and oxygen therapy details