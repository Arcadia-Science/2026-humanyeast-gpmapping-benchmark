# yeastandhuman-gpmapping

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

Genomic prediction (GP) mapping across simulated yeast and human traits. This
repository contains the code used to benchmark several penalized linear
regression methods — ridge, lasso, LARS, and elastic net (via scikit-learn) and
a PyTorch ridge with SGD — as well as standard GWAS via PLINK, across a range
of simulated genetic architectures.

Note that the UK Biobank data is private and requires an application for
access. See [https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/](https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/).


## Folder structure

```
scripts/          Python and R scripts and shell wrappers for each analysis step
envs/             Conda environment YAML files
input_data/       Please download from Zenodo (instructions below)

```

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

## Reproducibility

Input data are centred genotype matrices and normalised phenotype matrices
stored in [Apache Feather](https://arrow.apache.org/docs/python/feather.html)
format. The raw yeast genotype array has already been downloaded from S3 and binarized by
`scripts/save_input_data.py`. Follow the steps to download files from Zenodo.

### 1. Download data

Please download from [Zenodo](https://doi.org/10.5281/zenodo.19860006).

### 2. Create conda environments

| Environment file | Used for |
|---|---|
| `envs/scikit.yml` | sklearn regression (`regress.py`) |
| `envs/pytorch.yml` | PyTorch ridge regression and Optuna tuning |
| `envs/plink.yml` | GWAS with PLINK2 |
| `envs/r_env.yml` | phenotype generation and figure making |

```bash
mamba env create -n scikit --file envs/scikit.yml
mamba env create -n pytorch --file envs/pytorch.yml
mamba env create -n plink --file envs/plink.yml
mamba env create -n r_env --file envs/r_env.yml
```

### 3. Edit `config.yaml` (optional)

Set the `seed`, `prefix`, phenotype simulation parameters, and any method
settings. The Snakefile reads all values from `config.yaml` — you should not
need to edit scripts directly for a standard run.

For the PyTorch ridge step, the phenotype list and tuning hyperparameters live
in `scripts/run_tuning.sh` and `scripts/run_final_fit.sh` and must be edited
there directly.

### 4. Run the full pipeline

```bash
snakemake --cores 8 --use-conda
```

The pipeline runs these steps (see the Snakefile header for the full DAG):

| Step | Rules | Runs in parallel? |
|---|---|---|
| Simulate traits | `simulate_phenotypes` → `split_phenotypes` | — |
| Regression | `regress {ridge,lasso,elasticnet,lars}` | Yes (per method) |
| GWAS | `plink_gwas {test,train}` | Yes (per split) |
| PyTorch tuning | `tune_pytorch` | Yes (with regression/GWAS) |
| PyTorch final fit | `final_fit_pytorch` | After tuning |
| Aggregate | `aggregate_sklearn`, `aggregate_lars`, `aggregate_pytorch`, `build_gwas_matrices` | Yes |
| Polygenic scores | `gwas_predict {test,train}` | Yes (per split) |

Note that in the regression step, the "test" split of genotypes and phenotypes is used for training (instead of the "train" split) due to computational limitations. This is not a bug, but the test/train split can also be switched if you have a powerful machine and want to try it.
