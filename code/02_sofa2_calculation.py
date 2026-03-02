"""
SOFA-2 Calculation Script
Calculates SOFA-2 scores for a cohort using clifpy
"""

import pandas as pd
import json
from pathlib import Path
from clifpy.utils.sofa2 import calculate_sofa2
import yaml

# ============================================================================
# LOAD CONFIGURATION
# ============================================================================

print("\nLoading configuration from config/config.json...")

# Get the script's directory and go up one level to project root
script_dir = Path(__file__).parent
project_root = script_dir.parent
config_file = project_root / "config" / "config.json"

with open(config_file, 'r') as f:
    config = json.load(f)

SITE_NAME = config['site_name']
TABLES_PATH = config['tables_path']
OUTPUT_PATH = config['output_path']
FILE_TYPE = config['file_type'].lstrip('.')

# Change this based on your data
TIMEZONE = "US/Eastern" 

print(f"✓ Configuration loaded:")
print(f"  Site: {SITE_NAME}")
print(f"  Tables path: {TABLES_PATH}")
print(f"  Output path: {OUTPUT_PATH}")
print(f"  File type: {FILE_TYPE}")
print(f"  Timezone: {TIMEZONE}")

# ============================================================================
# CREATE CLIF CONFIG
# ============================================================================

print("\nSetting up CLIF configuration...")

# CLIF config file
clif_config = {
    'data_directory': TABLES_PATH,
    'filetype': FILE_TYPE,
    'timezone': TIMEZONE
}

# Save config to YAML file (clifpy expects YAML)
clif_config_path = Path("config/clif_config.yaml")
clif_config_path.parent.mkdir(exist_ok=True)

with open(clif_config_path, 'w') as f:
    yaml.dump(clif_config, f)

print(f"✓ CLIF config created at: {clif_config_path}")

# ============================================================================
# LOAD COHORT DATA
# ============================================================================

print("\nLoading cohort data...")

# Load the sipa_clif_cohort_sofa file
cohort_file = Path(OUTPUT_PATH) / "intermediate" / f"sipa_clif_cohort_sofa.{FILE_TYPE}"

if not cohort_file.exists():
    raise FileNotFoundError(f"Cohort file not found: {cohort_file}")

if FILE_TYPE == "parquet":
    cohort_df = pd.read_parquet(cohort_file)
elif FILE_TYPE == "csv":
    cohort_df = pd.read_csv(cohort_file)
else:
    raise ValueError(f"Unsupported file type: {FILE_TYPE}")

print(f"✓ Loaded {len(cohort_df)} records from {cohort_file.name}")
print(f"  Columns: {list(cohort_df.columns)}")

print(f"\nFirst few rows:")
print(cohort_df.head())

# ============================================================================
# CALCULATE SOFA-2
# ============================================================================

print("\n" + "="*70)
print("Running SOFA-2 calculation...")
print("This may take a while depending on cohort size...")
print("="*70 + "\n")

try:
    sofa_results = calculate_sofa2(
        cohort_df=cohort_df,
        clif_config_path=str(clif_config_path),
        return_rel=False  # Return pandas DataFrame
    )
    
    print("\n✓ SOFA-2 calculation complete!")
    print(f"✓ Results contain {len(sofa_results)} records")
    
except Exception as e:
    print(f"\n✗ Error during SOFA-2 calculation:")
    print(f"  {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
    raise

# ============================================================================
# VIEW RESULTS
# ============================================================================

print("\n" + "="*70)
print("SOFA-2 RESULTS SUMMARY")
print("="*70)

# Get SOFA-2 score columns
sofa_cols = [col for col in sofa_results.columns if col.startswith('sofa2_')]

print("\nSOFA-2 Score Columns:")
for col in sofa_cols:
    print(f"  - {col}")

print("\nSOFA-2 Score Statistics:")
print(sofa_results[sofa_cols].describe())

print("\nSample of results:")
display_cols = ['hospitalization_id', 'start_dttm', 'end_dttm'] + sofa_cols
available_display_cols = [col for col in display_cols if col in sofa_results.columns]
print(sofa_results[available_display_cols].head(10))

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("\n" + "="*70)
print("Saving results...")
print("="*70 + "\n")

# Create output directory
output_dir = Path(OUTPUT_PATH) / "intermediate"
output_dir.mkdir(parents=True, exist_ok=True)

# Save the results
output_file = output_dir / f"clif_sofa2_scores.{FILE_TYPE}"

if FILE_TYPE == "parquet":
    sofa_results.to_parquet(output_file, index=False)
elif FILE_TYPE == "csv":
    sofa_results.to_csv(output_file, index=False)

print(f"✓ SOFA-2 scores saved to: {output_file}")
print(f"✓ File size: {output_file.stat().st_size / (1024*1024):.2f} MB")

# Also save a summary CSV for quick review
summary_file = output_dir / "sofa2_summary.csv"
summary_df = sofa_results[available_display_cols]
summary_df.to_csv(summary_file, index=False)
print(f"✓ Summary saved to: {summary_file}")

print("\n" + "="*70)
print("SOFA-2 CALCULATION COMPLETE!")
print("="*70)
print(f"\nOutput files:")
print(f"  Full results: {output_file}")
print(f"  Summary: {summary_file}")
