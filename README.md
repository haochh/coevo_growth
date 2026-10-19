# Growth rate structures the dynamics of bacterial trait evolution

This repository contains the data, R scripts, and intermediate model results associated with the study **“Growth rate structures the dynamics of bacterial trait evolution.”** The analyses examine the coevolution between experimentally measured growth rates and 11 genomic, metabolic, ecological, and environmental traits across 283 bacterial species using multivariate Ornstein–Uhlenbeck (mvOU) models. The scripts focus on three stages of analysis: **univariate evolution**, **Bayesian regime inference**, and **bivariate coevolution**. The manuscript provides the complete trait definitions, data sources, and detailed methods.


## Repository contents

```text
coevo_growth/
├── README.md
├── Data/
│   ├── GR_Madin_imputated_RRN_11_factor_scaled.RData
│   ├── tree_madin.tree
│   ├── data_bivariate_mvOU_experiGR.RData
│   ├── tree_mvSLOUCH.RData
│   └── Regime_list.RData
├── Results/
│   ├──bayou/
│      ├──Bayou_OU1_ExperiGR_11factors.RData
│      ├──Bayou_OUM_ExperiGR_11factors.RData
│   ├──mvOU/
│      ├──mvOU_BAYOU_SyyDiagonal_niter100.RData
│      ├──mvOU_BAYOU_SyyDiagonal_niter100_ME.RData
│      ├──mvOU_BAYOU_SyyUpperTri_niter100.RData                                             
├── univariate_tempo_mode_experiGR.R
├── ExperiGR_bayou.R
├── mvOU_Bivariate_ExperiGR.R
└── Model_comparison_mvOU_Biva_ExperiGR.R
```

The **`Data/`** directory contains the species-level trait data, phylogenetic trees, and evolutionary-regime assignments used for model fitting. The R scripts cover four parts of the analysis: univariate evolutionary models, Bayesian regime inference, bivariate mvOU fitting, and comparison of fitted mvOU models. The **`Results/`** directory groups intermediate `.RData` files from the Bayesian and mvOU analyses. The principal saved model fits are described under **Intermediate results**.

## Data

The primary dataset combines experimentally measured minimum doubling times from the Madin et al. bacterial trait database with 11 associated traits. Representative genomes and the phylogenetic reference tree were derived from **GTDB release 220**. Additional trait sources and processing methods are described in the manuscript's Methods section.

| File | Contents | Used in |
| --- | --- | --- |
| [`GR_Madin_imputated_RRN_11_factor_scaled.RData`](Data/GR_Madin_imputated_RRN_11_factor_scaled.RData) | Species-level table containing experimentally derived growth rates and the 11 standardized associated traits. Species are identified by `gtdb_accession`. | Univariate analyses and Bayesian OU models |
| [`tree_madin.tree`](Data/tree_madin.tree) | Newick-format phylogeny of the species included in the experimental-growth-rate dataset, pruned from the GTDB reference tree. | Univariate analyses and Bayesian OU models |
| [`data_bivariate_mvOU_experiGR.RData`](Data/data_bivariate_mvOU_experiGR.RData) | Species-by-trait matrix prepared for mvOU analysis, with `Growth_rate_experiment` and the 11 associated traits ordered to match the phylogeny. | Bivariate mvOU fitting |
| [`tree_mvSLOUCH.RData`](Data/tree_mvSLOUCH.RData) | Phylogenetic tree prepared for `mvSLOUCH`, with tree height scaled to one. | Bivariate mvOU fitting |
| [`Regime_list.RData`](Data/Regime_list.RData) | Branch-level regime assignments inferred from Bayesian OU analyses, including growth-rate-defined, trait-defined, and combined growth-rate–trait regimes. | Bivariate mvOU fitting |

The analysed variables are:

| Trait group | Variables in the R scripts |
| --- | --- |
| Growth rate | `Growth_rate_experiment` |
| Genomic | `Genome_size`, `GC_percentage`, `rRNA_operon_copy_number` |
| Metabolic | `Metabolic_capacity`, `Metabolic_versatility`, `Metabolic_completeness` |
| Environmental | `Temperature`, `Salinity`, `pH`, `Oxygen` |
| Ecological | `Niche_breadth` |

Growth rate is defined as the **negative logarithm of minimum doubling time**; higher values correspond to faster growth. Traits were standardized to mean zero and standard deviation one for comparative analyses. `Temperature`, `Salinity`, and `pH` refer to estimated growth optima; `Oxygen` represents oxygen tolerance, and `Niche_breadth` represents the number of habitats occupied by a species.

## R scripts

### [`univariate_tempo_mode_experiGR.R`](univariate_tempo_mode_experiGR.R)

Analyses phylogenetic signal and fits univariate evolutionary models for growth rate and each of the 11 associated traits. It estimates Blomberg's *K* and Pagel's *λ*, compares Brownian motion (BM), Ornstein–Uhlenbeck (OU), and early-burst (EB) models using AICc, and extracts evolutionary parameters and computes OU half-lives to compare the tempos of evolution. The code uses `geiger::fitContinuous` with 100 optimization iterations for each model.

**Inputs:** `GR_Madin_imputated_RRN_11_factor_scaled.RData`, `tree_madin.tree`  
**Manuscript outputs:** Table 1; Figure S1.

### [`ExperiGR_bayou.R`](ExperiGR_bayou.R)

Fits single-optimum (OU1) and multiple-optimum (OUM) OU models using `bayou`. The OUM models use reversible-jump Markov chain Monte Carlo (rjMCMC) to infer shifts in adaptive optima. The script compares OU1 and OUM through stepping-stone marginal likelihoods and Bayes factors, summarizes the posterior number of regimes, identifies regime shifts, and constructs the growth-rate, trait-specific, and combined regime assignments used in subsequent mvOU analyses.

**Inputs:** `GR_Madin_imputated_RRN_11_factor_scaled.RData`, `tree_madin.tree`; saved Bayesian fit files for model comparison and regime summaries.  
**Manuscript outputs:** Figure 3; Figures S2–S3; Tables S2–S3.

### [`mvOU_Bivariate_ExperiGR.R`](mvOU_Bivariate_ExperiGR.R)

Fits bivariate evolutionary models for growth rate paired with each associated trait using `mvSLOUCH`. Model configurations vary in their selective-regime definition, adaptation-matrix structure, and stochastic covariance structure. The adaptation matrix allows independent (`Diagonal`), unidirectional (`UpperTri` or `LowerTri`), or bidirectional (`Any`) evolutionary coupling. Models are evaluated with no regime shift, growth-rate-defined regimes, trait-defined regimes, or combined regimes.

For the stochastic covariance matrix (`Syy`), the script fits both `Diagonal` and `UpperTri` structures. The diagonal-`Syy` analyses additionally consider measurement-error settings of 0, 0.01, 0.05, and 0.10. Each model configuration is fitted over 100 optimization runs, using a two-stage procedure to obtain starting values and then estimate the constrained model.

**Inputs:** `data_bivariate_mvOU_experiGR.RData`, `tree_mvSLOUCH.RData`, `Regime_list.RData`.  
**Outputs:** Fitted mvOU models used for model selection, parameter extraction, and Figures 4–5 and S4–S10.

### [`Model_comparison_mvOU_Biva_ExperiGR.R`](Model_comparison_mvOU_Biva_ExperiGR.R)

Processes the fitted mvOU models and identifies the best-supported evolutionary configurations. For each growth-rate–trait pair, adaptation-matrix structure, and regime definition, it extracts likelihood-based statistics (AICc, AIC, BIC, SIC, and log-likelihood) and retains the fit with the **lowest AICc across 100 optimization runs**. It compares alternative measurement-error settings and stochastic covariance structures, then identifies the overall best-supported model for each trait pair.

The selected models are used to extract adaptation matrices and characterize directional coevolution over instantaneous evolutionary timescales. Eigen-decomposition of these matrices yields eigenvalues, eigenvectors, and evolutionary half-lives describing the tempo of coevolution at longer timescales.

**Inputs:** `mvOU_BAYOU_SyyDiagonal_niter100.RData`, `mvOU_BAYOU_SyyDiagonal_niter100_ME.RData`, `mvOU_BAYOU_SyyUpperTri_niter100.RData`.  
**Manuscript outputs:** Figure 4 (AICc comparison of regimes and coupling modes); Figure 5a (direction and strength of coupling); Figure 5b (evolutionary axes and half-lives); Figure S5 (measurement-error comparisons); Figure S8 (`Syy` comparisons); Figure S10 (eigen-decomposition summaries).

## Intermediate results

The following saved `.RData` files contain computationally intensive model-fitting results used by the analysis scripts. They allow model comparisons and plots to be generated from existing fits without repeating the complete Bayesian MCMC or mvOU optimization procedures.

| File | Description | Associated script |
| --- | --- | --- |
| `Bayou_OUM_ExperiGR_11factors.RData` | Bayesian multiple-optimum OU fits, including MCMC samples and stepping-stone marginal-likelihood estimates for growth rate and associated traits. | `ExperiGR_bayou.R` |
| `Bayou_OU1_ExperiGR_11factors.RData` | Bayesian single-optimum OU fits used as the comparison with OUM models. | `ExperiGR_bayou.R` |
| `mvOU_BAYOU_SyyDiagonal_niter100.RData` | Bivariate mvOU model fits with diagonal stochastic covariance and no measurement error, covering alternative regime definitions and adaptation-matrix structures across 100 optimization runs. | `mvOU_Bivariate_ExperiGR.R`; `Model_comparison_mvOU_Biva_ExperiGR.R` |
| `mvOU_BAYOU_SyyDiagonal_niter100_ME.RData` | Bivariate mvOU fits with diagonal stochastic covariance and measurement-error levels of 0.01, 0.05, and 0.10. | `mvOU_Bivariate_ExperiGR.R`; `Model_comparison_mvOU_Biva_ExperiGR.R` |
| `mvOU_BAYOU_SyyUpperTri_niter100.RData` | Bivariate mvOU fits allowing correlated stochastic perturbations through an upper-triangular covariance structure, without measurement error. | `mvOU_Bivariate_ExperiGR.R`; `Model_comparison_mvOU_Biva_ExperiGR.R` |

## Software

All analyses were performed in **R**. Core packages include `ape`, `phytools`, `geiger`, `bayou`, `mvSLOUCH`, `PCMBaseCpp`, `tidyverse`, `doParallel`, `ggpubr`, `ggh4x`, `patchwork`, and `gt`. 
