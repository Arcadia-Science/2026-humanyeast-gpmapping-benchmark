# yeastandhuman-gpmapping

[A benchmark for the next generation of genotype-phenotype mapping](https://doi.org/10.57844/arcadia-27pw-kx6m)

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

Genomic prediction (GP) mapping across simulated yeast and human traits. This
repository contains the code used to benchmark several penalized linear
regression methods — ridge, lasso, LARS, and elastic net (via scikit-learn) and
a PyTorch ridge with SGD — as well as standard GWAS via PLINK, across a range
of simulated genetic architectures.

Note that the UK Biobank data is private and requires an application for
access. See [https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/](https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/).

## Folder structure

| Folder | Description |
| --- | --- |
| `scripts/` | Python and R scripts and shell wrappers for each analysis step |
| `envs/` | Conda environment YAML files |
| `input_data/` | Please download input data from Zenodo (instructions below) |
| `figure_intermediates/yeast_simulated_data_1510_ukbb_simulated_traits_1105/` | Please download aggregated results files for plotting from Zenodo (instructions below) |

## Required files

The following files must be present to run the pipeline end-to-end.

### Pipeline configuration

| File | Description |
| --- | --- |
| `Snakefile` | Pipeline definition; reads all settings from `config.yaml` |
| `config.yaml` | Seeds, prefixes, file paths, and parameter grids for the full run |

### Scripts

| File | Used by |
| --- | --- |
| `scripts/generate_phenotypes.r` | `simulate_phenotypes`, `split_phenotypes` |
| `scripts/regress.py` | `regress` (ridge, lasso, elasticnet, lars) |
| `scripts/plink.sh` | `plink_gwas` |
| `scripts/run_tuning.sh` | `tune_pytorch` |
| `scripts/fit_linear_sgd_optuna.py` | called by `run_tuning.sh` |
| `scripts/run_final_fit.sh` | `final_fit_pytorch` |
| `scripts/fit_linear_sgd_cli.py` | called by `run_final_fit.sh` |
| `scripts/aggregate_regression_outputs.py` | `aggregate_sklearn`, `aggregate_lars`, `aggregate_pytorch` |
| `scripts/build_gwas_matrices.py` | `build_gwas_matrices`, `gwas_predict` |
| `scripts/pub_figures_snakemake.r` | `pub_figures` |

### Conda environments

| File | Used by |
| --- | --- |
| `envs/scikit.yml` | `regress`, `aggregate_sklearn`, `aggregate_lars`, `aggregate_pytorch`, `build_gwas_matrices`, `gwas_predict` |
| `envs/pytorch.yml` | `tune_pytorch`, `final_fit_pytorch` |
| `envs/plink.yml` | `plink_gwas`, `clump_gwas` |
| `envs/r_env.yml` | `simulate_phenotypes`, `split_phenotypes`, `pub_figures` |

### Input data

The following files should be placed in `input_data/`. See [Download data](#1-download-data) below for instructions.

| File | Description |
| --- | --- |
| `genotypes_binarized.feather` | Binarized genotype matrix |
| `genotypes_binarized.transposed.tsv` | Transposed TSV genotype matrix |
| `genotypes.vcf` | VCF genotype matrix (plink input) |
| `SNP_list_pos_corrected.txt` | SNP metadata |
| `plink_bfile_prefix.bed` | Plink binary genotype file |
| `plink_bfile_prefix.bim` | Plink SNP information file |
| `plink_bfile_prefix.fam` | Plink sample information file |

To reproduce publication figures only (skipping the full pipeline), place the figure intermediate files in `figure_intermediates/yeast_simulated_data_1510_ukbb_simulated_traits_1105/`. See [Download data](#1-download-data) for the full list.

## Installation

This project uses conda environments. Install
[miniconda](https://docs.conda.io/projects/miniconda/en/latest/).

Then install [mamba](https://mamba.readthedocs.io/en/latest/), which you can do using conda:
```bash
conda install -n base -c conda-forge mamba
```

Then install snakemake into your base environment:
```bash
mamba install -n base -c conda-forge -c bioconda snakemake
```

## Analysis Pipeline

### 1. Download data

Please download data from [Zenodo](https://doi.org/10.5281/zenodo.19860006). The full-size yeast data is included, as well as a subset which is easier to run locally. To run the analysis using a subset of the yeast data, please place the following files into a folder titled `input_data/`:

| File | Description |
| --- | --- |
| `genotypes_binarized.feather` | Yeast genotype matrix in feather format |
| `genotypes_binarized.transposed.tsv` | Yeast genotype matrix in TSV format |
| `genotypes.vcf` | Yeast genotype matrix in VCF format |
| `genotypes_binarized_subset_5000_snpsubset_1000_centered.feather` | Subsetted centered yeast genotype matrix for running pipeline locally |
| `genotypes_binarized_subset_5000_snpsubset_1000_uncentered.feather` | Subsetted uncentered yeast genotype matrix for running pipeline locally |
| `genotypes_binarized_subset_5000_snpsubset_1000_ids.feather` | Subsetted list of yeast strain IDs |
| `yeast_subset_5000_snpsubset_1000_allele_frequencies.Rda` | Allele frequencies for subsetted yeast genotypes for running pipeline locally |
| `SNP_list_pos_corrected.txt` | List of yeast SNP IDs |
| `plink_bfile_prefix.fam` | Yeast genotype matrix in plink format (.fam file) |
| `plink_bfile_prefix.bed` | Yeast genotype matrix in plink format (.bed file) |
| `plink_bfile_prefix.bim` | Yeast genotype matrix in plink format (.bim file) |

To recreate the plots from our publication, place the following files into a folder titled `figure_intermediates/yeast_simulated_data_1510_ukbb_simulated_traits_1105/`:

| File | Description |
| --- | --- |
| `combined_all_betas_human_1105.feather` | Estimated SNP effects for all methods and traits in human |
| `combined_all_betas_yeast_1510.feather` | Estimated SNP effects for all methods and traits in yeast |
| `human_littlelonger_with_fullinfo_1105.feather` | Simulated true and estimated effect sizes for all methods and traits in human |
| `yeast_littlelonger_with_fullinfo_1510.feather` | Simulated true and estimated effect sizes for all methods and traits in yeast |
| `human_true_avg_distance_1105.feather` | Average distances between simulated true and top estimates for all methods and traits in human |
| `yeast_true_avg_distance_1510.feather` | Average distances between simulated true and top estimates for all methods and traits in yeast |
| `yeast_cumulative_1510.feather` | Cumulative true and false positives for ROC plot in yeast |
| `yeast_roc_1510.feather` | Summarized true and false positives for ROC plot in yeast |
| `yeast_roc_approx_250_1510.feather` | Summarized true and false positives for approximate ROC plot in yeast |
| `yeast_snp_pairs_1510.feather` | Distances between SNP pairs in yeast |
| `all_methods_effects_correlations_yeast_simulated_data_1510_ukbb_simulated_traits_1105.feather` | Variant effect r^2 for all methods and traits for human and yeast in feather format |
| `all_methods_max_prediction_correlations_yeast_simulated_data_1510_ukbb_simulated_traits_1105.feather` | Maximum phenotype prediction r^2 for all methods and traits for human and yeast in feather format |
| `all_methods_parameters_yeast_simulated_data_1510_ukbb_simulated_traits_1105.feather` | Training parameters for relevant methods and traits for human and yeast in feather format |
| `all_methods_prediction_correlations_yeast_simulated_data_1510_ukbb_simulated_traits_1105.feather` | Phenotype prediction r^2 for all methods and traits for human and yeast in feather format |

### 2. Create conda environments

| Environment file | Used for |
| --- | --- |
| `envs/scikit.yml` | sklearn regression (`regress.py`) |
| `envs/pytorch.yml` | PyTorch ridge regression and Optuna tuning |
| `envs/plink.yml` | GWAS with PLINK2 |
| `envs/r_env.yml` | phenotype generation and figure making |

Create the environments by running the following in your shell:
```bash
mamba env create -n scikit --file envs/scikit.yml
mamba env create -n pytorch --file envs/pytorch.yml
mamba env create -n plink --file envs/plink.yml
mamba env create -n r_env --file envs/r_env.yml
```

### 3. Edit `config.yaml` (optional)

If you want to change any of the settings, modify `config.yaml`. The Snakefile reads all values from `config.yaml` — you should not need to edit scripts directly for a standard run. Keep all the parameters the same to run the analysis using a subset of the yeast data. Running the scikit-learn and pytorch methods using the full dataset requires increased computational power.

For the PyTorch ridge step, tuning hyperparameters (trial count, alpha/LR range, timeout, etc.) can be overridden via environment variables without editing the script. See `scripts/run_tuning.sh` for the full list of available variables and their defaults.

### 4. Run the full pipeline

This command will run the analysis pipeline on a subset of the yeast data:
```bash
snakemake --cores 8 --use-conda
```

The pipeline runs these steps (see the Snakefile header for the full DAG):

| Step | Rules | Runs in parallel? |
| --- | --- | --- |
| Simulate traits | `simulate_phenotypes` → `split_phenotypes` | — |
| Regression | `regress {ridge,lasso,elasticnet,lars}` | Yes (per method) |
| GWAS | `plink_gwas {test,train}` | Yes (per split) |
| PyTorch tuning | `tune_pytorch` | Yes (with regression/GWAS) |
| PyTorch final fit | `final_fit_pytorch` | After tuning |
| Aggregate | `aggregate_sklearn`, `aggregate_lars`, `aggregate_pytorch`, `build_gwas_matrices` | Yes |
| Polygenic scores | `gwas_predict {test,train}` | Yes (per split) |
| Plot | `pub_figures` | — |

To run just the figures step (reproduce published plots only):
```bash
snakemake --cores 8 logs/pub_figures_1510_1105.done --use-conda
```

## Notes

- In the regression step, the "test" split of genotypes and phenotypes is used for training (instead of the "train" split) due to computational limitations. This is not a bug, but the test/train split can also be switched if you have a powerful machine and want to try it.
- plink2 is not available on conda for some computers. If you cannot install it via conda, install it [here](https://www.cog-genomics.org/plink/2.0/)
- Data/code for plotting figures 9 (ROC plots for human) and 11 (variant effect prediction by minor allele frequency) from our [publication](https://doi.org/10.57844/arcadia-27pw-kx6m0) are not included for human data security reasons.
