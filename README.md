# GRMTree for Longitudinal Measurement Invariance in PROMs

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/1373758225.svg)](https://doi.org/10.5281/zenodo.22803364)

<!-- badges: end -->

Reproducibility bundle for the manuscript:

**"Tree-based item response theory model for assessing longitudinal measurement invariance in patient-reported outcome measures"**
Sajobi, T. T., Arimoro, O. I., Bohm, E., Lix, L. M.
*Quality of Life Research* (under revision).

Contact: olayinka.arimoro@ucalgary.ca

---

## What is here

This repository contains R code to reproduce the simulation study and the empirical Oxford Knee Score analysis reported in the manuscript.

```
├── README.md                  # This file
├── LICENSE                    # MIT
├── 01_helpers.R               # Shared helper functions used by simulation scripts
├── 02_grmtree_type1.R         # Type I error simulation for GRMTree (6 and 12 items)
├── 03_grmtree_power.R         # Power simulation for GRMTree (6 and 12 items)
├── 04_lrt_type1.R             # Type I error simulation for LRT (6 and 12 items)
├── 05_lrt_power.R             # Power simulation for LRT (6 and 12 items)
└── 06_oks_analysis.R          # Empirical OKS analysis (GRMTree + LRT)
```

## Requirements

- R (>= 4.5.0)
- Packages: `grmtree`, `mirt`, `tidyverse`, `readxl` (for OKS analysis only)

The `grmtree` package is available on CRAN:

```r
install.packages(c("grmtree", "mirt", "tidyverse", "readxl", "ggplot2", "EnvStats"))
```

## Running the simulation

Each simulation script is self-contained and can be run independently:

```r
source("01_helpers.R")       # load once per session
source("02_grmtree_type1.R") # or any of 03-05
```

Or from the command line:

```bash
Rscript 02_grmtree_type1.R
```

**Runtime.** Each script fits many models sequentially (500 replications × many conditions × two scale lengths). Full runs take several hours on a standard laptop. The scripts save their output as `.rds` files at the end so results can be inspected without re-running.

## Empirical OKS analysis

Script `06_oks_analysis.R` assumes access to the cleaned OKS dataset (`cleaned_oks_data.xlsx`), which contains patient responses to the 12 OKS items at pre- and 1-year post-total-knee-arthroplasty, together with patient demographics. This dataset cannot be shared publicly due to data-use restrictions of the Winnipeg Regional Health Authority Joint Replacement Registry. Investigators seeking access should contact the registry directly.

The script demonstrates:
- Baseline GRM fit and fit indices (RMSEA, CFI, SRMR)
- GRMTree fit with time as the sole partitioning covariate
- Two-factor longitudinal GRM with LRT (global and item-level with Bonferroni correction)

## Simulation design

- **Items:** 6 and 12 (mirroring short and moderate-length PROMs)
- **Categories:** 5 (ordinal)
- **Sample sizes:** 200, 500, 750, 1500
- **Latent trait correlation across occasions (ρ):** 0.3, 0.6, 0.9
- **Non-invariance items:** 3 items (items 1, 3, 6 for the 6-item scale; items 1, 6, 12 for the 12-item scale)
- **Non-invariance type:** uniform (thresholds shift), non-uniform (discrimination shifts), both
- **Non-invariance magnitude (Δ):** 0.5 (medium), 1.0 (large)
- **Latent mean shift (μ):** 0 or 0.5
- **Replications:** 500 per condition
- **Seed:** `set.seed(2026)`

## Citation

If you use this code, please cite the manuscript above and the underlying `grmtree` package:

Arimoro, O. I., Sajobi, T. T., & Lix, L. M. (2026). grmtree: Recursive Partitioning for Graded Response Models (Version 0.3.0) [R package]. Comprehensive R Archive Network (CRAN). 
