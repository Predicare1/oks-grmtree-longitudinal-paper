# GRMTree for Longitudinal Measurement Invariance in PROMs

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22803365.svg)](https://doi.org/10.5281/zenodo.22803365)

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
├── sim_supp_null/            # Null simulation: change in latent mean or variance
└── sim_supp_null_B_lrt/      # LRT rerun for the variance conditions (variance freely estimated)
```

## Requirements

- R (>= 4.5.0)
- Packages: `grmtree`, `mirt`, `tidyverse`; for the OKS analysis also `readxl` and `EnvStats`

```r
install.packages(c("grmtree", "mirt", "tidyverse", "readxl", "EnvStats"))
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

## Supplementary null simulation (Supplementary Table A7)

Examines Type I error when item parameters are identical at both occasions but the latent trait distribution changes at follow-up: mean shift of 0.25 or 0.5 (variance 1), or variance of 0.5 or 1.5 (mean 0). 6 and 12 items; N = 200, 500, 750, 1,500; ρ = 0.6; 500 replications per condition.

- `sim_supp_null/`: GRMTree for all conditions, and the LRT for the mean-shift conditions.
- `sim_supp_null_B_lrt/`: the LRT for the variance conditions, with the follow-up variance freely estimated.

These were run on a computing cluster (SLURM), one job per condition, with seed `2026 + task ID`. Any single condition can also be run locally from inside its folder:

```bash
Rscript scripts/sim_supp_null_run.R 9            # task IDs 1 to 32
Rscript scripts/sim_supp_null_B_lrt_run.R 17     # task IDs 17 to 32
```

The reported results are in `supplementary_null_results.csv` and `B_rerun_lrt_freevar.csv`.

## Empirical OKS analysis

Script `06_oks_analysis.R` assumes access to the cleaned OKS dataset (`cleaned_oks_data.xlsx`), which contains patient responses to the 12 OKS items at pre- and 1-year post-total-knee-arthroplasty, together with patient demographics. This dataset cannot be shared publicly due to data-use restrictions of the Winnipeg Regional Health Authority Joint Replacement Registry. Investigators seeking access should contact the registry directly.

The script reproduces:

- Descriptive statistics (Table 3, Supplementary Table A5, Supplementary Figure A1)
- Baseline GRM fit (RMSEA, CFI, SRMSR)
- GRMTree with time as the only partitioning covariate (Figure 1, Table 4)
- Two-factor longitudinal GRM with likelihood ratio tests, with the post-surgery latent mean and variance freely estimated: global test and item-level tests of overall, uniform, and non-uniform non-invariance with Bonferroni correction (Supplementary Table A6)

## Citation

If you use this code, please cite the manuscript above and the underlying `grmtree` package:

Arimoro, O. I., Sajobi, T. T., & Lix, L. M. (2026). grmtree: Recursive Partitioning for Graded Response Models (Version 0.3.0) [R package]. Comprehensive R Archive Network (CRAN). 
