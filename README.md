# Improving Comparability in Network Meta-Analysis: Dose, Heterogeneity, and Mixed Treatments in Depression Research

**Katherine Elise Herder**
Doctor of Philosophy, Department of Epidemiology and Biostatistics, 2026
Mel and Enid Zuckerman College of Public Health, The University of Arizona

---

## Overview

This repository contains all code, data, and supporting files necessary to reproduce the analyses presented in the dissertation. The work evaluates how dose is incorporated into network meta-analyses (NMA) of depression treatments, tests model-based dose-response methods (MBNMA) in both nonpharmacologic and pharmacologic settings, and provides practical guidance for applied researchers.

## Repository Structure

Each chapter (2–5) has its own folder containing up to four subfolders (when applicable): data, code, figures, results.

### Chapter 2: Current Practices and Evidence Gaps

An empirical scoping review of 22 depression NMAs (2020–2025) examining how dose is handled in practice.

**Data:** 
- .csv file containing all extracted data from the 22 articles.
  
**Code:**
- `chapter2_figures.R` — Produces Figures 2.2 and 2.3 (dose specification heatmap; primary analysis approach by intervention category)

**Dependencies:** `dplyr`, `ggplot2`

---

### Chapter 3: Evaluating Dose-Response Methods for Nonpharmacologic Interventions

Reanalysis of the Tian et al. (2024) exercise-for-depression network (45 RCTs, 4 exercise modalities). Compares six dose-handling approaches: lumping, splitting, and four MBNMA functional forms (linear, Emax, quadratic, natural cubic spline). Includes sensitivity analysis using total minutes/week as an alternative dose metric.

**Data:** Tian_Dataset.xlsx - Supplemental dataset from Tian et al. (2024), loaded as `Tian`.

**Code (run in order):**
1. `01_load_data.R - Load in .xlsx
2. `02_clean_data.R` — Adds the sensitivity dose metric to the Tian dataset
3. `03_describe_data.R` — Produces all descriptive statistics reported in the Network Characteristics section
4. `04_nonpharm.R` — Full analysis: model fitting, effect comparison (splitting vs Emax), minimum effective dose estimation, consistency assessment, sensitivity analysis, and cross-method comparison

**Dependencies:** `MBNMAdose` (≥ 0.5.0), `gemtc`, `tidyverse`, `coda`, `magick`

---

### Chapter 4: Evaluating Dose-Response Pooling Assumptions in Pharmacologic Networks

Simulation study and empirical analyses of antidepressant dose-response networks. Tests whether the default shared-ED50 specification in MBNMAdose is adequate across single-class (SSRI) and multi-class (GRISELDA) networks. Includes class-effect extensions (Models C1 and C2) and prior sensitivity analysis.

**Data:**
- `ssri` — Built-in dataset from the `MBNMAdose` package (60 RCTs, 5 SSRIs)
- `Griselda_Fixed.xlsx` — Derived from the full GRISELDA network (https://data.mendeley.com/datasets/83rthbp8ys/2) -> all trials with fixed doses were manually labeled.
  
**Code (run in order):**
1. `01_load_data.R` - Load Griselda_Fixed.xlsx and ssri
2. `02_clean_data.R` — Cleans the full GRISELDA dataset: restricts to fixed-dose trials, computes response/remission rates, classifies outcome types, assigns drug classes, computes fluoxetine-equivalent doses, and produces `Griselda_clean`
3. `03_describe_data.R` — Produces all descriptive statistics for the SSRI and GRISELDA network characteristics sections, plus the cross-network overlap analysis
4. `04_ssri.R` — SSRI empirical analysis: fits Models A, B, Ref, and Lumped; model comparison (DIC, LPML); dose-response curves; agent rankings (unified d_hat approach); prior sensitivity analysis
5. `05_simulation.R` — Simulation study: data-generating mechanism calibrated to SSRI Model B posteriors; 100 replications; evaluates bias, MSE, coverage, ranking accuracy (Spearman ρ, SUCRA, per-agent rank match rate), and model selection accuracy (DIC, LPML)
6. `06_griselda.R` — GRISELDA empirical analysis: fits all six models (A, B, C1, C2, Ref, Lumped); class-effect extensions; expanded model comparison; prior sensitivity analysis

**Dependencies:** `MBNMAdose` (≥ 0.5.0), `R2jags`, `dplyr`, `tidyr`, `ggplot2`, `coda`, `stringr`, `patchwork`, `kableExtra`

**Note:** The simulation script (step 5) takes approximately 25 hours to complete at 100 replications. Checkpoints are saved every 5 replications and the script supports resumption from the most recent checkpoint if interrupted.

---

### Chapter 5: Exploratory Simulations and Practical Guidance for Researchers

Exploratory simulations varying three dimensions of the SSRI network to identify boundary conditions for reliable treatment rankings: heterogeneity magnitude (τ), evidence density, and network connectivity.

**Prerequisites:** Chapter 4's simulation script (sections 0–6) must be run first, as Chapter 5 relies on the following objects from that session: `network_template`, `true_params`, `true_tau`, `true_mu_mean`, `true_mu_sd`, `generate_ssri_dataset()`, `fit_all_models()`, `extract_metrics()`, `compute_lpml()`, and MCMC settings.

**Data:** Simulated from the SSRI network structure (no external data required beyond what is loaded in Chapter 4).

**Code:**
- `chapter5.R` — Runs all simulation scenarios (heterogeneity sweep, evidence density, network connectivity) and compiles results across dimensions. Includes rank match rate and SUCRA compilation. Produces network diagrams for evidence density and connectivity scenarios (Figures 5.1 and 5.2)

**Dependencies:** Same as Chapter 4.

**Note:** At 100 replications across 11 scenarios (3 τ values + 3 density + 4 connectivity + 1 dose-enriched), the full simulation takes approximately 4–5 days. Checkpoints are saved every 5 replications per scenario. The script supports automatic resumption from the most recent checkpoint.

---

## Software Environment

All analyses were conducted in R (version 4.5.1) with JAGS (version 4.3.1) as the MCMC sampler. Key packages:

| Package | Version | Role |
|---------|---------|------|
| `MBNMAdose` | 0.5.0 | Dose-response MBNMA and standard NMA |
| `R2jags` | — | R interface to JAGS |
| `gemtc` | — | Network plots (Chapter 3 only) |
| `coda` | — | MCMC diagnostics |
| `dplyr` / `tidyr` | — | Data manipulation |
| `ggplot2` | — | Visualization |

## Reproducibility

- Random seeds are set at the beginning of each script and passed to JAGS where possible.
- Minor variations in results may occur across computing environments due to differences in numerical precision or JAGS version.
- Complete session information is saved to `results/session_info.txt` in each chapter folder upon script completion.

## Contact

Katherine Elise Herder
keherder@arizona.edu
