# MetTarget

MetTarget is a machine learning framework for prioritizing immunometabolic therapeutic targets in inflammatory bowel disease (IBD). It integrates public genetics, bulk RNA-seq, and single-cell RNA-seq data to score and rank candidate targets in Crohn's disease (CD) and ulcerative colitis (UC).

## What's in this repository

This repository includes trained MetTarget models (R objects) for CD and UC, prediction code, and code and data to reproduce manuscript figures. All features for model training and testing are provided (76 for CD, 134 for UC). Target sets covered include clinical phase 2+ IBD drug targets (34 CD, 48 UC), FDA-approved non-IBD immune targets (n = 151), and SPHK1, a new target identified by MetTarget.

## Requirements

- R (>= 4.x recommended)
- R packages listed in `DESCRIPTION`

## Usage

See `MetTarget_model/README.md`.

## License

See `LICENSE.txt`.
