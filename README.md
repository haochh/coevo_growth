# Growth rate structures the dynamics of bacterial trait evolution

**Data and R code accompanying the study _“Growth rate structures the dynamics of bacterial trait evolution”_.**

## Overview

Bacterial species vary widely in their growth rates, but the evolutionary relationships between growth rate and other bacterial traits remain incompletely understood. This study combines experimentally measured minimal doubling times for **283 bacterial species** with **11 genomic, metabolic, ecological, and environmental traits** to investigate their evolutionary tempo, adaptive regimes, and patterns of coevolution.

We use phylogenetic comparative analyses and Ornstein–Uhlenbeck (OU) models to examine whether growth rate and associated traits evolve independently, show directional evolutionary coupling, or exhibit reciprocal coupling. The analyses distinguish evolutionary dynamics at instantaneous and longer macroevolutionary timescales.

## Repository contents

```text
coevo_growth/
├── Data/                             # Data inputs (see files in this directory)
├── univariate_tempo_mode_experiGR.R  # Phylogenetic signal and univariate evolutionary models
├── ExperiGR_bayou.R                  # Bayesian single- versus multiple-optimum OU analyses
├── mvOU_Bivariate_ExperiGR.R         # Bivariate evolutionary models and coevolutionary coupling
└── README.md                         # Project documentation
```

The directory and script names above reflect the public repository. Consult the scripts for the exact input filenames, output paths, and run parameters. Not all downstream analyses or figure-generation workflows are necessarily represented by these three scripts.

## Data

The main dataset comprises experimentally measured growth rates and 11 associated traits:

| Trait group | Traits |
| --- | --- |
| Growth | Maximum growth rate, represented as the negative logarithm of minimal doubling time |
| Genomic | Genome size, GC percentage, ribosomal RNA operon (`rrn`) copy number |
| Metabolic | Metabolic capacity, metabolic versatility, metabolic completeness |
| Environmental | Optimum temperature, optimum salinity, optimum pH, oxygen tolerance |
| Ecological | Niche breadth |

Minimal doubling times were compiled from the Madin et al. bacterial trait database. Representative genomes and the reference phylogeny were obtained using **GTDB release 220**. Other sources and derivation methods include rrnDB, KEGG/KofamScan, GenomeSPOT, and ProkAtlas, as detailed in the accompanying manuscript.

For comparative modelling, traits were standardized to a mean of zero and a standard deviation of one. The GTDB phylogeny was pruned to the study species; phylogenetic tree height was rescaled to one for multivariate OU fitting.

**Data provenance:** Please consult the manuscript Methods and original data providers when reusing source-derived measurements. Verify the exact file-level contents of `Data/` before treating any individual file as raw, processed, or analysis-ready data.

## Analytical workflow

### 1. Univariate trait evolution

**Script:** [`univariate_tempo_mode_experiGR.R`](univariate_tempo_mode_experiGR.R)

Quantifies phylogenetic signal and compares Brownian motion (BM), Ornstein–Uhlenbeck (OU), and early-burst (EB) models. Model support is assessed using corrected Akaike information criterion (AICc); evolutionary half-lives are evaluated for fitted OU models.

### 2. Adaptive regime inference

**Script:** [`ExperiGR_bayou.R`](ExperiGR_bayou.R)

Uses the R package `bayou` to compare single-optimum (OU1) and multiple-optimum (OUM) evolutionary models. The manuscript describes reversible-jump MCMC inference of adaptive regime shifts and Bayes-factor comparisons using marginal likelihood estimates.

### 3. Bivariate coevolution

**Script:** [`mvOU_Bivariate_ExperiGR.R`](mvOU_Bivariate_ExperiGR.R)

Fits phylogenetic bivariate models of growth rate and individual associated traits. The manuscript compares BM and OU models, alternative adaptive-regime definitions, different adaptation-matrix structures (independent, unidirectional, or bidirectional coupling), and alternative assumptions about stochastic covariance and measurement error. Model selection uses AICc.

The manuscript further interprets fitted adaptation matrices through their eigenvalues and eigenvectors to compare the tempo of joint evolutionary modes.

## Software requirements

The analyses are written in **R**. Major packages reported in the manuscript include:

| R package | Reported version | Role |
| --- | --- | --- |
| `ape` | 5.7.1 | Phylogenetic tree processing |
| `phytools` | 2.5.2 | Phylogenetic signal |
| `geiger` | 2.0.11 | Univariate BM, OU, and EB models |
| `bayou` | 2.3.2 | Bayesian adaptive-regime inference |
| `mvSLOUCH` | 2.7.7 | Multivariate evolutionary model fitting |

These are manuscript-reported versions, not a verified complete dependency lockfile. Additional packages may be loaded by individual scripts. For exact reproducibility, record your R version, all package versions, and operating-system details (for example, with `sessionInfo()` or an `renv.lock` file).

## Getting started

Clone the repository or download it as a ZIP archive:

```bash
git clone https://github.com/haochh/coevo_growth.git
cd coevo_growth
```

Open the repository directory as your working directory in R or RStudio. Before running an analysis, inspect the corresponding script for its package requirements, input file paths, output directories, random seeds, and computational settings.

The conceptual analysis order is:

1. `univariate_tempo_mode_experiGR.R`
2. `ExperiGR_bayou.R`
3. `mvOU_Bivariate_ExperiGR.R`

**Reproducibility note:** This is an analytical overview, not a verified one-command pipeline. Some analyses may be computationally intensive, and successful execution may require adjusting paths or supplying intermediate files. In particular, Bayesian MCMC and repeated multivariate optimizations may take substantial time. Do not assume these scripts can be run sequentially without checking their inputs and dependencies.

## Interpretation

The analyses are designed to evaluate the direction and tempo of *evolutionary coupling* inferred from comparative phylogenetic models. Directional coefficients should not be interpreted as direct experimental evidence of mechanistic causation.

## Citation

If you use these data or scripts, please cite the accompanying manuscript:

> _Growth rate structures the dynamics of bacterial trait evolution_. [Full author list, publication year, journal, and DOI to be added when available.]

Please also cite the original databases and software relevant to the data and methods reused.

## Data and code availability

The data and analysis code maintained for this study are available at:

**https://github.com/haochh/coevo_growth**

For a permanent, citable snapshot, archive a tagged release in Zenodo and add its DOI here once issued.

## License and contact

A reuse license has not been specified in this README. Please add an appropriate `LICENSE` file after confirming rights to redistribute third-party data and code.

For questions, corrections, or reproducibility issues, please [open a GitHub issue](https://github.com/haochh/coevo_growth/issues).
