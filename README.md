# These Four First

Replication code for:

> Gibson, B. (2025). *These Four First: A Robust Endogenous Baseline for Relational Event Sequences*. SocArXiv. https://osf.io/preprints/socarxiv/mu7j3_v1

## Overview

This repository contains the analysis code and pre-computed results for a study identifying the four endogenous terms that consistently enter greedy forward AIC selection within the first four steps across 217 empirical interaction sequences spanning ten social domains. The robust baseline is: **SndInitiation**, **IOviolation**, **Repetition**, and **RecRecencyRank**.

## Files

| File | Description |
|------|-------------|
| `unified_analysis.R` | Main analysis script: greedy forward AIC selection, take-one-out analyses, all figures and tables |
| `fit_final_models.R` | Fits the per-network "final" greedy model (last step with negative ΔAICs), extracts coefficients and SEs for the appendix table |
| `gen_table.R` | Generates the appendix coefficient table from saved `final_model_coefs.Rdata` |
| `stability_results.Rdata` | Per-network greedy paths, TOO-4, and TOO-12 results (output of `unified_analysis.R`) |
| `greedy_results.Rdata` | Greedy paths in legacy format for compatibility with figure scripts |
| `final_model_coefs.Rdata` | Per-network final-model coefficients and SEs (output of `fit_final_models.R`) |
| `paper/` | LaTeX source, bibliography, and figures for the manuscript |

## Requirements

- R (≥ 4.0)
- Packages: `llrem`, `survival`, `parallel`, `ggplot2`, `dplyr`, `tidyr`, `scales`, `gridExtra`

The `llrem` package implements the separable relational event model estimator used throughout.

## Usage

### Reproduce the main analysis

Edit the configuration block near the top of `unified_analysis.R`:

```r
PROJ_DIR    <- "/path/to/FirstFour"   # project root
DATA_DIR    <- "/path/to/DataP"       # root of network data files
N_NONEVENTS <- 120L
N_CORES     <- parallel::detectCores() - 1L
```

Then run:

```r
source("unified_analysis.R")
```

This writes all figures and LaTeX tables to `paper/` and saves `stability_results.Rdata` and `greedy_results.Rdata`.

### Reproduce the appendix coefficient table

`fit_final_models.R` also requires the raw network data files. Edit `PROJ_DIR` at the top and run:

```r
source("fit_final_models.R")
```

Or, if you only want to regenerate the table from the saved results:

```r
source("gen_table.R")
```

## Data

The raw network data files are not distributed with this repository. The 217 sequences span ten social domains: gang violence, international conflict, email, private messages, glance networks, phone calls (DoD), phone calls (organized crime), SMS, Twitter/EMO, and classroom interactions. See the Data Sources appendix in the paper for full provenance and citations.

## Citation

```
@misc{gibson2025firstfour,
  author    = {Gibson, B.},
  title     = {These Four First: A Robust Endogenous Baseline for Relational Event Sequences},
  year      = {2025},
  publisher = {SocArXiv},
  url       = {https://osf.io/preprints/socarxiv/mu7j3_v1}
}
```
