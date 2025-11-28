# Metabolic state shapes cortisol reactivity to acute stress: A systematic review and meta-analysis of metabolic and hormonal modulators

A Bayesian meta-analysis examining how glucose, progesterone, and estradiol influence cortisol reactivity to acute stress.

## Authors
- Madeleine Kördel
- Maria Meier
- Anne Kühnel
- Nils B. Kroemer

## Overview

This repository contains the complete code and data for a Bayesian meta-analysis investigating the effects of metabolic and hormonal modulators (glucose, progesterone, and estradiol) on cortisol stress reactivity. 
The analysis employs Bayesian methods to aggregate effect sizes, assess heterogeneity, and evaluate evidence strength through Bayes factors.

## Repository Structure

### Data Files
- `df.RDa` - Prepared data frame (R binary format)
- `df.csv` - Data frame in CSV format
- `master_table.xlsx` - Original extracted data from systematic review

### Scripts
- `meta_stress_dataprep.R` - Data preparation script that processes the original dataset
- `MainEffectMAstress.R` - Main analysis script containing all Bayesian meta-analysis procedures

## Analysis Features

### Primary Analyses
- **Bayesian meta-analysis** for three conditions: glucose, progesterone, and estradiol
- **Effect size aggregation** using the BHHR method with correlation = 0.5
- **Half-Cauchy priors** on heterogeneity parameter (τ) with scale = 0.5
- **Normal priors** on overall effect size (μ) with mean = 0, SD = 1.5

### Robustness Checks
- **Bayes factor sensitivity analysis** across different prior specifications
- **Outlier detection** using interquartile range (IQR) method
- **Sensitivity analyses** excluding identified outliers
- **Publication bias assessment** through funnel plots

### Visualizations
The script generates several publication-ready figures:
- Forest plots for each meta-analysis (Figures 3, 5)
- Funnel plot for publication bias assessment (Figure 4)
- Bayes factor robustness plot (Figure 6)
- Posterior distribution plots (Figure 7)

## Key Findings

The analysis reveals a robust modulatory role of metabolic state, specifically glucose availability, on cortisol stress reactivity, while evidence for sex hormone effects remains inconclusive.

## Requirements

### R Packages
The script automatically installs and loads required packages via `pacman`:
- `bayesmeta` - Bayesian meta-analysis
- `metafor` - Meta-analysis methods
- `ggplot2` - Data visualization
- `dplyr`, `tidyr` - Data manipulation
- `ggdist` - Distribution visualization
- Additional packages: `compute.es`, `cowplot`, `data.table`, `esc`, `forestplot`, `knitr`, `MAd`, `reactable`, `readr`, `readxl`, `rmarkdown`, `R.rsp`, `stringr`, `scales`, `writexl`

## Usage

1. Ensure all data files are in the working directory
2. Run `meta_stress_dataprep.R` to prepare the dataset (if starting from raw data)
3. Execute `MainEffectMAstress.R` for the complete analysis pipeline

The script will automatically:
- Load and process the data
- Perform Bayesian meta-analyses for all three conditions
- Generate all figures and statistical outputs
- Conduct robustness checks and sensitivity analyses

## Statistical Methods

- **Bayesian framework** with informative priors
- **Random-effects meta-analysis** accounting for between-study heterogeneity
- **Bayes factors** for evidence evaluation
- **Credible intervals** for effect size estimation
- **Outlier-robust analyses** for sensitivity testing

## Citation

If you use this code or data, please cite the original publication. (https://doi.org/10.1016/j.ynstr.2025.100764)


## Contact

For questions or issues, please contact the corresponding authors or open an issue in this repository.
