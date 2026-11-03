# Reproducibility package

This package contains analysis scripts, output tables, the locally frozen protocol and notes. It does not include participant-level data or install software.

Download DEMO_L.xpt, DIQ_L.xpt, GHB_L.xpt, BMX_L.xpt, BPXO_L.xpt and BPQ_L.xpt from CDC into `03_data/raw/`, using URLs like `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.xpt` with the filename substituted. Verify the six downloads against `03_data/file_hashes.json` by running `python 04_analysis/verify_input_hashes.py` from the package directory. Python 3.11 or newer is required for this check. The manifest paths are relative to the package directory; any missing file or mismatch returns a nonzero exit status.

Prerequisites: R 4.5.3 with `survey` and `foreign`. Run `Rscript 04_analysis/run.R <absolute-paper-directory> <R-library-directory>`, then `Rscript 04_analysis/quality_review.R <absolute-paper-directory> <R-library-directory>` and `Rscript 04_analysis/tables_review.R <absolute-paper-directory> <R-library-directory>`. Both review scripts source and rerun `run.R`, overwriting generated tables; use a fresh directory if existing outputs must be preserved. Data files must be in `<absolute-paper-directory>/03_data/raw/`. The bootstrap uses 2,000 replicates.

The plotting scripts were written for the original workspace and require path edits elsewhere. The Python checks use paths relative to this package and require pandas and numpy. Run `python 04_analysis/independent_check.py`, `python 04_analysis/check_survey_variance.py`, `python 04_analysis/missingness_check.py` and `python 04_analysis/auc_linearization_check.py` after downloading the source data; the missingness and AUC checks overwrite their generated CSVs. The separately written Python checks verify selected point estimates and Taylor-linearised standard errors for five proportions. The checks were conducted within the same project and do not validate the AUC bootstrap or replace independent human review.

Bootstrap results use analyst-generated PSU replicates, not CDC-provided replicate weights. An independent statistical review remains necessary. The protocol was frozen locally and was not publicly registered. A single elevated HbA1c is not a clinical diagnosis.

The AUC linearization calculation is exploratory and post hoc. It is a separately coded first-order cross-check within the same project, not an independent human review and not a validation of the bootstrap confidence interval.

Run `python 04_analysis/selection_comparison.py` to reproduce the post hoc unweighted age/sex selection comparison; it overwrites its generated CSV and does not correct selection bias.

Run `Rscript 04_analysis/outcome_characteristics.R <absolute-paper-directory> <R-library-directory>` for post hoc outcome-stratified characteristics. It sources run.R and overwrites generated main outputs; use a fresh directory. Sparse cells and NCHS reliability need independent review.

Run `Rscript 04_analysis/auc_linearization_survey_check.R <absolute-paper-directory> <R-library-directory>` after the Python AUC linearization check. It sources run.R and overwrites some main tables. This is the same influence formula checked using R survey, not an independent human review.

Recreate the journal-oriented figures with `python 04_analysis/plot_acta_figures.py --output-dir <directory>` after generating flow.csv and screening_metrics.csv. Outputs are EPS vectors and preview PNGs; captions are supplied separately in the review package.
