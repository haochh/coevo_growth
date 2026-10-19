# Growth rate structures the dynamics of bacterial trait evolution

**Data, R scripts, and computational workflow accompanying the manuscript _“Growth rate structures the dynamics of bacterial trait evolution”_.**

## Overview

This repository contains materials for phylogenetic comparative analyses of bacterial growth rate and 11 associated genomic, metabolic, ecological, and environmental traits across **283 bacterial species**. The analyses assess phylogenetic signal, fit univariate evolutionary models, identify adaptive-regime shifts with Bayesian Ornstein–Uhlenbeck (OU) models, and fit bivariate multivariate-OU (mvOU) models to compare alternative evolutionary coupling structures.

The scripts focus on three analysis stages: **univariate evolution**, **Bayesian regime inference**, and **bivariate coevolution**. The associated manuscript provides the complete trait definitions, data provenance, and statistical methods.

## Repository structure

```text
coevo_growth/
├── README.md
├── Data/
│   ├── GR_Madin_imputated_RRN_11_factor_scaled.RData
│   ├── tree_madin.tree
│   ├── data_bivariate_mvOU_experiGR.RData
│   ├── tree_mvSLOUCH.RData
│   └── Regime_list.RData
├── univariate_tempo_mode_experiGR.R
├── ExperiGR_bayou.R
└── mvOU_Bivariate_ExperiGR.R
```

The five files shown under `Data/` are **confirmed analysis inputs**; the directory may contain other materials not enumerated here. The public repository currently does not show a top-level `Results/` directory. The scripts also refer to intermediate result objects/files that may have been generated separately and are described below.

## Data files

| File | Contents and purpose | Used by |
| --- | --- | --- |
| [`Data/GR_Madin_imputated_RRN_11_factor_scaled.RData`](Data/GR_Madin_imputated_RRN_11_factor_scaled.RData) | Standardized species-by-trait data. The scripts expect an object called `GR_Madin_imputated_RRN_11_factor_scaled` with a `gtdb_accession` identifier, experimentally derived growth rate (`Growth_rate_experiment`), and 11 trait columns. | Univariate and `bayou` analyses |
| [`Data/tree_madin.tree`](Data/tree_madin.tree) | Newick phylogeny for the study species, with tips corresponding to GTDB accessions; rescaled to tree height 1 within the univariate and `bayou` scripts. | Univariate and `bayou` analyses |
| [`Data/data_bivariate_mvOU_experiGR.RData`](Data/data_bivariate_mvOU_experiGR.RData) | Prearranged standardized data for bivariate mvOU models. The script expects an object `data_bivariate_mvOU_experiGR` whose row names match phylogenetic tip labels. | Bivariate mvOU analyses |
| [`Data/tree_mvSLOUCH.RData`](Data/tree_mvSLOUCH.RData) | Tree object (`tree_mvSLOUCH`) used in mvSLOUCH modelling; its height was already rescaled to 1. | Bivariate mvOU analyses |
| [`Data/Regime_list.RData`](Data/Regime_list.RData) | Regime annotations (`Regime_list`) indexed by focal trait. Each entry is expected to contain `regime_trait1` (growth-rate-defined regime), `regime_trait2` (focal-trait-defined regime), and `regime_trait1_trait2` (combined regime) for phylogenetic edges. | Bivariate mvOU analyses |

The `Regime_list` object can be constructed from `bayou` posterior summaries using the final section of `ExperiGR_bayou.R`. The standardized data and tree used in bivariate analyses are distributed as separate preprocessed objects.

### Trait variables used in the scripts

| Variable name in R | Trait |
| --- | --- |
| `Growth_rate_experiment` | Growth rate, defined in the manuscript as the negative logarithm of minimum doubling time |
| `Genome_size`, `GC_percentage`, `rRNA_operon_copy_number` | Genomic traits |
| `Metabolic_capacity`, `Metabolic_versatility`, `Metabolic_completeness` | Metabolic traits |
| `Temperature`, `Salinity`, `pH`, `Oxygen` | Environmental traits |
| `Niche_breadth` | Ecological trait |

The manuscript describes experimental doubling times from Madin et al.; representative genomes and the phylogeny from GTDB release 220; and supporting data derived from rrnDB, KEGG/KofamScan, GenomeSPOT, and ProkAtlas. All trait columns are standardized for the comparative analyses. See the manuscript Methods for original measurements, transformations, definitions, and source citations.

## R scripts and their outputs

### 1. [`univariate_tempo_mode_experiGR.R`](univariate_tempo_mode_experiGR.R) — univariate evolutionary tempo and mode

**Inputs:** `GR_Madin_imputated_RRN_11_factor_scaled.RData`, `tree_madin.tree`.

**Analysis:**
- Aligns trait observations to the phylogenetic tip order and rescales tree height to 1.
- Fits Brownian motion (`BM`), Ornstein–Uhlenbeck (`OU`), and early-burst (`EB`) models to growth rate and each of 11 traits using `geiger::fitContinuous()` with `niter = 100`.
- Compares models using AICc.
- Estimates Blomberg's K and Pagel's λ with `phytools::phylosig()`.
- Extracts fitted parameters and calculates OU evolutionary half-lives (`log(2) / alpha`), including half-life as a percentage of tree height.

**Named results created in the R session:**

| R object | Meaning | Manuscript output |
| --- | --- | --- |
| `Fit_Continuous_ExperiGR_11factors` | Full BM/OU/EB fitted-model objects, grouped by trait and model | Underlying fitted models |
| `AICc_Fit_Continuous_ExperiGR_11factors` | AICc for each trait × evolutionary model | Figure S1 |
| `PhySignal_ExperiGR_11factors` | Full phylogenetic-signal test objects | Table 1 |
| `PhySignal_ExperiGR_11factors_result` | Summary of K, λ, and p-values | Table 1 |
| `Param_Fit_Continuous` | Best-fitting models, parameter estimates, OU half-lives | Table 1 |
| `Univariate_FitContinuous_final_output`, `gt_table` | Merged results and publication-formatted summary | Table 1 |

The script draws the Figure S1 comparison plot and constructs a Table 1-style `gt` table; it does **not** explicitly export these displays or save its fitted objects to files.

### 2. [`ExperiGR_bayou.R`](ExperiGR_bayou.R) — adaptive-regime inference

**Inputs:** `GR_Madin_imputated_RRN_11_factor_scaled.RData`, `tree_madin.tree`. The result-summary portions additionally use previously computed `Bayou_OUM_ExperiGR_11factors.RData` and `Bayou_OU1_ExperiGR_11factors.RData` files.

**Analysis:**
- Fits a multi-optimum OU (`OUM`) model for each of the 12 variables with `bayou` reversible-jump MCMC.
- Fits corresponding single-optimum (`OU1`) models.
- Uses 100,000 MCMC iterations, 30% burn-in, and 10,000 stepping-stone iterations for marginal-likelihood estimation (as specified in the script).
- Compares OUM and OU1 log marginal likelihoods and computes log Bayes factors.
- Summarizes posterior numbers and locations of regime shifts; draws regime maps.
- Constructs two-regime annotations for growth rate and each focal trait, including a combined-regime column for downstream mvOU fitting.

**Named results created or loaded:**

| R object / referenced file | Meaning | Manuscript output |
| --- | --- | --- |
| `Bayou_OUM_ExperiGR_11factors` / `Bayou_OUM_ExperiGR_11factors.RData` | Per-trait OUM MCMC chains and stepping-stone summaries (`MCMC`, `SS`) | Figures 3a, S2, S3; Table S2 |
| `Bayou_OU1_ExperiGR_11factors` / `Bayou_OU1_ExperiGR_11factors.RData` | Per-trait OU1 MCMC chains and stepping-stone summaries | Figure 3a; Table S2 |
| `bayou_OUM_regime` | Posterior summaries of regime counts | Figure 3a; Table S2 |
| `Bayes_factor_OU1_OUM`, `Bayou_results` | Marginal likelihoods, log Bayes factors, and merged regime summaries | Figure 3a |
| `Regime_list` / `Data/Regime_list.RData` | Growth-rate, trait-specific, and combined phylogenetic regime assignments | Input to bivariate mvOU models |

`Regime_list.RData` is available as an input file, but the script does not contain an explicit `save()` command for its newly calculated regime list or for the MCMC collections. The two `Bayou_*.RData` names above are **referenced intermediate files**, not confirmed files in `Data/`.

### 3. [`mvOU_Bivariate_ExperiGR.R`](mvOU_Bivariate_ExperiGR.R) — bivariate coevolution model fits

**Inputs:** `data_bivariate_mvOU_experiGR.RData`, `tree_mvSLOUCH.RData`, `Regime_list.RData`.

For each pair of `Growth_rate_experiment` and one of the 11 associated traits, the script fits `mvSLOUCH::ouchModel()` under alternative settings:

- **Adaptive-regime definitions:** no shift (`null`/OU1), growth-rate regime (`regime_trait1`), associated-trait regime (`regime_trait2`), and combined regime (`regime_trait1_trait2`).
- **Adaptation matrix A:** `Diagonal` (uncoupled), `UpperTri`, `LowerTri` (unidirectional coupling), and `Any` (fully coupled).
- **Stochastic covariance structure Syy:** `Diagonal` (independent perturbations) and `UpperTri` (correlated perturbations).
- **Measurement error:** none for both Syy types; additionally, `0.01` (`low`), `0.05` (`mid`), and `0.10` (`high`) for `Syy = Diagonal`.
- **Optimization:** 100 repeated starts; initial `DecomposablePositive` fits followed by constrained fits seeded from the initial parameter estimates. Repeated fitting is parallelized.

**Named result collections:**

| R object | Model set | Main use |
| --- | --- | --- |
| `mvOU_BAYOU_SyyDiagonal_niter100` | `Syy = Diagonal`, no measurement error | Baseline bivariate OU and regime/model comparisons |
| `mvOU_BAYOU_SyyDiagonal_niter100_ME` | `Syy = Diagonal`, three nonzero measurement-error levels | Measurement-error sensitivity analysis |
| `mvOU_BAYOU_SyyUpperTri_niter100` | `Syy = UpperTri`, no measurement error | Correlated-versus-independent stochastic perturbation comparison |

Each collection stores results by focal trait; the error-sensitive collection has an additional `low`/`mid`/`high` level. Within each optimization replicate, `start` contains unconstrained initial model fits, and `final` contains fitted models grouped by adaptation-matrix type and regime definition. Individual final fits also record their seeds. These detailed objects are the **intermediate results** underlying comparisons of AICc, regime definitions, coupling modes, and measurement-error robustness (Figures 4, S4–S8 and Tables S4–S5 in the manuscript).

**Important:** This script performs the primary model fitting, but does not contain a complete downstream script for selecting best fits across all repetitions, tabulating all comparisons, plotting Figure 4, or calculating the eigenvector/eigenvalue summaries in Figure 5. It also does not explicitly save its three large result collections to `.RData` or `.rds` files. Those post-processing outputs should not be presented as automatically reproducible from this script alone.

## Intermediate results: how to read them

The largest objects in this project are fitted model collections rather than additional raw measurements. The table below distinguishes the **data inputs confirmed in the repository** from **model objects named by the scripts**.

| Type | Location or object | How to use it |
| --- | --- | --- |
| Analysis-ready inputs | `Data/*.RData`, `Data/tree_madin.tree` (listed above) | Load before model fitting; align accession identifiers and phylogeny |
| Univariate fitted models | `Fit_Continuous_ExperiGR_11factors` | Extract model likelihoods, AICc, and evolutionary parameters without refitting |
| Bayesian chains and marginal likelihoods | `Bayou_OUM_ExperiGR_11factors`, `Bayou_OU1_ExperiGR_11factors` | Recreate Bayes-factor summaries and regime maps when saved model objects are available |
| Regime annotations | `Data/Regime_list.RData` | Supply `regime_trait1`, `regime_trait2`, and combined regimes to bivariate fitting |
| Bivariate mvOU fit ensembles | `mvOU_BAYOU_SyyDiagonal_niter100`, `mvOU_BAYOU_SyyDiagonal_niter100_ME`, `mvOU_BAYOU_SyyUpperTri_niter100` | Retrieve replicate fits, AICc, adaptation matrices, covariance matrices, and starting seeds |

For any `.RData` file, identify the objects it contains before analysis:

```r
file <- "Data/Regime_list.RData"
objects_loaded <- load(file)
print(objects_loaded)
str(Regime_list, max.level = 2)
```

The presence and placement of any additional large intermediate result files should be documented separately if those objects are uploaded or archived outside this repository.

## Reproducing the analyses

### Requirements

R and the following packages are used across the scripts: `tidyverse`, `ape`, `geiger`, `phytools`, `gt`, `bayou`, `patchwork`, `mvSLOUCH`, `PCMBaseCpp`, `doParallel`, and `foreach`. The manuscript reports versions **ape 5.7.1**, **phytools 2.5.2**, **geiger 2.0.11**, **bayou 2.3.2**, and **mvSLOUCH 2.7.7**. The exact R version and complete dependency snapshot are not included in the repository.

### Obtain the files

```bash
git clone https://github.com/haochh/coevo_growth.git
cd coevo_growth
```

The three R scripts are currently in the repository root, whereas the verified input files are inside `Data/`. Their `load()` and `read.tree()` calls use **bare filenames** (for example, `load("Regime_list.RData")`) rather than `Data/...` paths. Therefore, either run the script blocks with `Data/` as the working directory or update these calls to point to `Data/` before executing them. Do not run the scripts unchanged from the repository root and expect input paths to resolve.

### Suggested order of analysis

1. **Univariate models:** Run the univariate script's preparation, BM/OU/EB fitting, AICc, phylogenetic-signal, and table/figure sections. This reproduces the analyses for Table 1 and Figure S1, subject to the code adjustments below.
2. **Bayesian adaptive regimes:** Run the OUM and OU1 fitting sections in `ExperiGR_bayou.R`, retain the MCMC result collections, then run its Bayes-factor and regime-annotation sections. The precomputed `Data/Regime_list.RData` allows readers to proceed to bivariate analyses without rerunning the MCMC.
3. **Bivariate mvOU models:** Load the three processed inputs and run the `Syy = Diagonal` (without/with measurement error) and `Syy = UpperTri` fitting blocks in `mvOU_Bivariate_ExperiGR.R`. Save the three fit collections before closing R. For each competing model, compare the retained best-fitting replicate using AICc and then extract the fitted A and Syy matrices as needed for the manuscript's summaries.

### Practical reproducibility notes

Several script-level details require attention before an independent researcher can run all three scripts end-to-end:

- `univariate_tempo_mode_experiGR.R` calls `phylosig()` without attaching `phytools`, and its parameter-extraction section uses publication-formatted trait labels to index a list created with original variable names. The latter lookup needs to use consistent names.
- `ExperiGR_bayou.R` references `Bk` in its stepping-stone calls without defining it in the script; it also contains an absolute `/Users/.../` path to an existing OUM result file. The MCMC collections and `Regime_list` need explicit saving if refitted.
- `mvOU_Bivariate_ExperiGR.R` uses `regime_df` in an example fit before assigning it, sets `ncores <- 100`, and hardcodes `.libPaths("/home/haochh/Rlib")` on parallel workers. Adjust these to the available computational environment. Its results are held in memory unless saved explicitly.
- The repository does not currently include a fully automated pipeline for all final figures, tables, and downstream eigen-analyses. Models with 100 optimization replicates and Bayesian MCMC can be computationally intensive.

When rerunning and archiving analyses, retain package versions (`sessionInfo()`), random seeds, fitted-model objects, and intermediate summaries. Bayesian chains also require convergence diagnostics before interpretation or model comparison.
