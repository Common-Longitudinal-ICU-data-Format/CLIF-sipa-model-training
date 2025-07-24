# *Development of an improved in-hospital mortality predictor based on SOFA*

## CLIF VERSION

2.1

## Objective

*The purpose of this project is to develop a more robust, statistically sound model for in-hospital mortality using the variables outlined by the Sequential Organ Failure Assessment (SOFA) score.*

## Required CLIF tables and fields

Please refer to the online [CLIF data dictionary](https://clif-consortium.github.io/website/data-dictionary.html), [ETL tools](https://github.com/clif-consortium/CLIF/tree/main/etl-to-clif-resources), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields.

1\. **patient**: `patient_id`, `race_category`, `ethnicity_category`, `sex_category`

2\. **hospitalization**: `patient_id`, `hospitalization_id`, `admission_dttm`, `discharge_dttm`, `age_at_admission`, `death_dttm`

3\. **vitals**: `hospitalization_id`, `recorded_dttm`, `vital_value` - `vital_category == 'map'|'spo2'|'weight_kg'`

4\. **labs**: `hospitalization_id`, `lab_result_dttm`, `lab_category`, `lab_value` - `lab_category == 'bilirubin_total'|'creatinine'|'platelet_count'|'pao2_arterial'`

5\. **medication_admin_continuous**: `hospitalization_id`, `admin_dttm`, `med_name`, `med_category`, `med_dose`, `med_dose_unit` - `med_group == 'vasoactives'`

6\. **respiratory_support**: `hospitalization_id`, `recorded_dttm`, `device_category`, `mode_category`, `tracheostomy`, `fio2_set`, `lpm_set`, `resp_rate_set`, `peep_set`, `resp_rate_obs`

## Cohort identification

The study population included all adults (age \>= 18 years) that were admitted to the intensive care unit and identified as having been on life support for at least six hours. Life support was defined as receiving vasoactive medications, invasive or non-invasive mechanical ventilation, or high-flow/facemask oxygen therapy for hypoxic respiratory failure. In addition, the years selected are **beginning of 2018 through the end of 2023**.

## Expected Results

A well-calibrated in-hospital mortality classification model that performs better than SOFA for critical resource allocation.

## Detailed Instructions for running the project

## 1. Update `config/config.json`

Follow instructions in the [config/README.md](config/README.md) file for detailed configuration steps.

## 2. Set up the project environment

Run `00_renv_restore.R` to set up the project environment

## 3. Run code

Please read items 1-5 carefully. Run code in the following order:

```{bash}
#!/bin/bash

# This script runs the R scripts in the correct order.
Rscript code/0a_respiratory_support_waterfall.R
Rscript code/01_cohort_identification.R
Rscript code/02_feature_set_processing.R
Rscript code/03_table1.R
Rscript code/04_model_training.R
```

1.  `0a_respiratory_support_waterfall.R`. This script runs Nick Ingraham's respiratory waterfall algorithm which will horizontally fill in various device categories. Requires lookup-tables `device_conversion_table_updated` and `device_name_mapper`.

2.  `01_cohort_identification.R`. This script creates the cohort dataframe. This script also outputs data needed to create a STROBE diagram. Make sure to specify the correct dates to select from `clif_hospitalization`.

3.  `02_feature_set_processing.R`. This script creates the feature set needed to train/test the model.

4.  `03_table1.R`. This script outputs data needed to create a Table 1. It also outputs a STROBE diagram.

5.  `04_model_training.R` This script trains the models and outputs the results.
