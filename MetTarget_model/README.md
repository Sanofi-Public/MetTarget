MetTarget Score Prediction for IBD
============================================

This folder provides example input data, pre-trained MetTarget models, and scripts to
reproduce MetTarget score predictions for inflammatory bowel disease (IBD),
including a single-gene SPHK1 example evaluated in the paper.

The prediction script is self-contained and requires only base R and openxlsx.


Contents
--------

| File | Description |
|------|-------------|
| input.IBD.data.xlsx | Feature data for all example genes (4 worksheets; see below) |
| plain.MetTarget.pretrained.UC.RDS | Plain MetTarget model for ulcerative colitis (UC) |
| plain.MetTarget.pretrained.CD.RDS | Plain MetTarget model for Crohn's disease (CD) |
| predict.MetTarget.score.plain.R | Self-contained R script for plain-model prediction |
| run.sh | Bash wrapper that runs all four example predictions |
| example.UC.MetTarget.output.txt | Example UC output (reference for reproducibility) |
| example.CD.MetTarget.output.txt | Example CD output (reference for reproducibility) |
| example.UC.SPHK1.MetTarget.output.txt | Example UC.SPHK1 output (reference for reproducibility) |
| example.CD.SPHK1.MetTarget.output.txt | Example CD.SPHK1 output (reference for reproducibility) |

Input data (input.IBD.data.xlsx)
--------------------------------

The workbook has four worksheets:

| Sheet | Format | Description |
|-------|--------|-------------|
| UC | wide (genes × features) | Main UC input for batch prediction |
| CD | wide (genes × features) | Main CD input for batch prediction |
| UC.SPHK1 | long (features × SPHK1) | Example UC feature profile for SPHK1 |
| CD.SPHK1 | long (features × SPHK1) | Example CD feature profile for SPHK1 |

The UC and CD sheets contain feature values for genes drawn from public databases:

- UC drug targets in clinical trial phase 2 or above (n = 48)  
- CD drug targets in clinical trial phase 2 or above (n = 34)  
- Non-IBD immune targets of FDA-approved therapies (n = 151)  

All model input features are included (134 for UC, 76 for CD).

The UC.SPHK1 and CD.SPHK1 sheets provide the complete feature profile for SPHK1,
a gene evaluated in the paper, in long format (one row per feature). predict.MetTarget.score.plain.R
accepts both wide and long input formats and auto-transposes long sheets before prediction.


Models
------

Two plain MetTarget models are provided as R serialized objects (*.RDS):

- plain.MetTarget.pretrained.UC.RDS — ensemble of 5 radial SVM models stored as base-R
  lists; final score is the median of the five model probabilities
- plain.MetTarget.pretrained.CD.RDS — same structure for Crohn's disease

Each model stores preprocessing parameters, feature names, support vectors, SVM
coefficients, kernel sigma, and Platt calibration values. The models can be loaded
with base R only.

Use plain.MetTarget.pretrained.UC.RDS for UC and UC.SPHK1 sheets; use
plain.MetTarget.pretrained.CD.RDS for CD and CD.SPHK1 sheets.

Requirements
------------

- R ≥ 4.0.0
- R package: openxlsx

Tested environment (R 4.2.0, CentOS Linux 7):

  openxlsx_4.2.5.2

Usage
-----

To run a single sheet manually:

  Rscript predict.MetTarget.score.plain.R <excel_file> <sheet_name> <model.RDS>

Arguments:  
  excel_file   — path to input.IBD.data.xlsx  
  sheet_name   — worksheet name (UC, CD, UC.SPHK1, or CD.SPHK1)  
  model.RDS    — plain.MetTarget.pretrained.UC.RDS (for UC, UC.SPHK1) or  
                 plain.MetTarget.pretrained.CD.RDS (for CD, CD.SPHK1)  

Or, from this directory, run all four example predictions:

  bash run.sh

This executes:

  - Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx UC plain.MetTarget.pretrained.UC.RDS  
  - Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx CD plain.MetTarget.pretrained.CD.RDS  
  - Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx UC.SPHK1 plain.MetTarget.pretrained.UC.RDS  
  - Rscript predict.MetTarget.score.plain.R input.IBD.data.xlsx CD.SPHK1 plain.MetTarget.pretrained.CD.RDS  

Logs are written to task.<sheet_name>.log.


Output
------

Each run produces a tab-delimited file:

  <sheet_name>.MetTarget.output.txt

For example, the UC.SPHK1 run writes UC.SPHK1.MetTarget.output.txt.

Columns:
  gene                    — gene symbol  
  predicted.MetTarget.score   — median predicted probability across the 5 models (0–1)

Batch runs (UC, CD) return one row per gene. SPHK1 runs (UC.SPHK1, CD.SPHK1) return
a single row for SPHK1.

Missing feature values are imputed as 0 before prediction.


Reproducibility
---------------

The provided example.*.MetTarget.output.txt files are the expected results from running
the scripts on input.IBD.data.xlsx:

  - example.UC.MetTarget.output.txt        ↔ UC.MetTarget.output.txt  
  - example.CD.MetTarget.output.txt        ↔ CD.MetTarget.output.txt  
  - example.UC.SPHK1.MetTarget.output.txt  ↔ UC.SPHK1.MetTarget.output.txt  
  - example.CD.SPHK1.MetTarget.output.txt  ↔ CD.SPHK1.MetTarget.output.txt  

After execution, compare your output files to these reference files. For example:

  md5sum UC.SPHK1.MetTarget.output.txt example.UC.SPHK1.MetTarget.output.txt  

Batch outputs (UC, CD) match the reference scores to machine precision; md5 checksums
may differ slightly due to floating-point formatting.

Expected run time for the full demo (4 runs): under one minute.
