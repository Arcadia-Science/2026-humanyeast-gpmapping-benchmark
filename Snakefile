# Snakefile — yeastandhuman-gpmapping genomic prediction pipeline
#
# Runs the full analysis from genotype data through phenotype simulation,
# regression, GWAS, and polygenic score prediction.
#
# ── Quick start ───────────────────────────────────────────────────────────────
#
#   Edit config.yaml to set your seed, prefix, and parameter grids, then:
#
#   # Run the full pipeline
#   snakemake --cores N
#
#   # Dry run: print rules and filenames without executing
#   snakemake -n
#
#   # Per-rule conda environments (requires envs/*.yml)
#   snakemake --cores N --use-conda
#
#
# ── Pipeline order ────────────────────────────────────────────────────────────
#
#   1. simulate_phenotypes          — simulate quantitative traits (phase 1)
#      split_phenotypes             — write train/test feather files (phase 2)
#   2. (parallel)
#      regress          {method}    — sklearn CV regression (ridge/lasso/elasticnet/lars)
#      plink_gwas       {split}     — GWAS via plink2
#      tune_pytorch                 — Optuna hyperparameter search (pytorch ridge)
#   3. final_fit_pytorch            — fit pytorch ridge with tuned hyperparameters
#   4. aggregate_sklearn {method}   — collect per-trait sklearn outputs into matrices
#      aggregate_lars               — same for lars (separate: different filename tag)
#      aggregate_pytorch            — collect pytorch ridge outputs into matrices
#      build_gwas_matrices {split}  — assemble plink GWAS results into beta/p-value matrices
#      clump_gwas        {split}    — LD-clump GWAS results (plink2 --clump, per split)
#   5. gwas_predict     {split}     — compute polygenic scores from filtered betas
#   6. plot             {split}     — creates plots seen in our publication
#
# ─────────────────────────────────────────────────────────────────────────────

import os

configfile: "config.yaml"

# ── Constants derived from config ─────────────────────────────────────────────

SEED   = config["seed"]
PREFIX = config["prefix"]

# Formatted p-value string matching build_gwas_matrices.py's fmt_threshold():
#   e.g. 1.2e-6  →  "1.20e-06"
P_STR = f"{config['gwas_p_threshold']:.2e}".replace("+", "")

# Key directories
TT_DIR     = f"test_train_seed_{SEED}"         # train/test feather files
INT_DIR    = f"intermediates_seed_{SEED}"       # per-trait sklearn CSV outputs
TUNING_DIR = config["pytorch_tuning_dir"]       # Optuna JSON results
FINAL_DIR  = f"final_fit_results_seed_{SEED}"  # pytorch ridge final-fit outputs

# Genotype basename for --geno arg (e.g. "yeast_genotypes_binarized.feather")
GENO_BASENAME = os.path.basename(config["geno_file"])
# Genotype prefix without extension (used by split_phenotypes / phase 2)
GENO_PREFIX = os.path.splitext(GENO_BASENAME)[0]

# Suffix appended to genotype files by generate_phenotypes.r when subsetting
_geno_suffix = ""
if config.get("subset_geno"):
    _geno_suffix += f"_subset_{config['subset_geno']}"
if config.get("subset_snps"):
    _geno_suffix += f"_snpsubset_{config['subset_snps']}"
# Genotype prefix including any subsetting suffix — used as the --geno arg in
# split_phenotypes so the R script can find the files simulate_phenotypes wrote
GENO_PREFIX_SUBSETTED = GENO_PREFIX + _geno_suffix
# Full-sample centered genotype file produced as a side effect of simulate_phenotypes
GENO_CENTERED_FILE = os.path.join(
    config["geno_dir"],
    GENO_PREFIX_SUBSETTED + "_centered.feather"
)
# Allele frequency file produced by simulate_phenotypes; basename only (R loads via geno_dir)
ALLELE_FREQ_BASENAME = PREFIX + _geno_suffix + "_allele_frequencies.Rda"
ALLELE_FREQ_FILE     = os.path.join(config["geno_dir"], ALLELE_FREQ_BASENAME)

# Phenotype simulation parameters as space-separated strings for R list args
VAVG_STR   = ",".join(str(v) for v in config["vavg_ratios"])
QTL_STR    = ",".join(str(q) for q in config["qtl_numbers"])
BSENSE_STR = ",".join(str(b) for b in config["broad_sense"])

# Optional --pheno-filter flag for regress.py; empty string if not configured
PHENO_FILTER_ARG = (
    "--pheno-filter " + " ".join(config["pheno_filter"])
    if config.get("pheno_filter")
    else ""
)

SKLEARN_METHODS = config["sklearn_methods"]
NON_LARS        = [m for m in SKLEARN_METHODS if m != "lars"]
GWAS_SPLITS     = config["gwas_splits"]
LARS_MAXITER    = config["lars_maxiter"]

CLUMP_R2   = config["clump_r2"]
CLUMP_P1   = config["clump_p1"]
CLUMP_BINS = " ".join(str(b) for b in config["clump_bins"])

YEAST_PLOT_SEED   = config["yeast_figure_seed"]
HUMAN_PLOT_SEED   = config.get("human_figure_seed", "")

YEAST_FIG_PREFIX  = config["yeast_figure_prefix"]
HUMAN_FIG_PREFIX  = config.get("human_figure_prefix", "human")

# Pre-computed figure intermediates directory (downloaded from Zenodo).
# Path mirrors the naming used by pub_figures_snakemake.r:
#   figure_intermediates/{yeast_prefix}_{yeast_seed}[_{human_prefix}_{human_seed}]
_fig_int_suffix = f"{YEAST_FIG_PREFIX}_{YEAST_PLOT_SEED}"
if config.get("human_figure_seed"):
    _fig_int_suffix += f"_{HUMAN_FIG_PREFIX}_{HUMAN_PLOT_SEED}"
FIG_INT_DIR = os.path.join("figure_intermediates", _fig_int_suffix)

# ── rule all — top-level targets ──────────────────────────────────────────────
# Snakemake works backwards from these targets to determine which rules to run.

rule all:
    input:
        # ── sklearn regression aggregates ──────────────────────────────────────
        # Non-lars methods: one coefficient matrix per method
        expand(
            "aggregate_outputs_{method}_test_seed_{seed}/"
            "coeff_matrix_{method}_test_seed_{seed}.feather",
            method=NON_LARS,
            seed=SEED,
        ),
        # Lars: filenames include maxiter tag (e.g. lars_maxiter1000)
        (
            f"aggregate_outputs_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}/"
            f"coeff_matrix_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}.feather"
        ),
        # ── PyTorch ridge aggregate ────────────────────────────────────────────
        (
            f"aggregate_outputs_pytorch_ridge_train_seed_{SEED}/"
            f"coeff_matrix_pytorch_ridge_train_seed_{SEED}.feather"
        ),
        # ── GWAS outputs ───────────────────────────────────────────────────────
        # Beta and p-value matrices from plink GWAS results
        expand(
            "aggregate_outputs_{split}_seed_{seed}/"
            "beta_matrix_{split}_seed_{seed}.feather",
            split=GWAS_SPLITS,
            seed=SEED,
        ),
        # Polygenic scores computed from p-value-filtered betas
        expand(
            "aggregate_outputs_{split}_seed_{seed}/"
            "polygenic_scores_{split}_seed_{seed}_p{p_str}.feather",
            split=GWAS_SPLITS,
            seed=SEED,
            p_str=[P_STR],
        ),
        # ── LD clumping ────────────────────────────────────────────────────────
        expand(
            "logs/clump_gwas_{split}_{seed}.done",
            split=GWAS_SPLITS,
            seed=SEED,
        ),
        # ── Publication figures ────────────────────────────────────────────────
        f"logs/pub_figures_{YEAST_PLOT_SEED}_{HUMAN_PLOT_SEED}.done",




# ── Step 1a: Simulate quantitative traits ────────────────────────────────────
# Reads the binarized genotype matrix and simulates traits for all combinations
# of vavg_ratio × qtl_number × broad_sense. Writes Rda files and an allele
# frequency file. A sentinel (.done) file marks completion because R produces
# many output files with no single canonical path.

rule simulate_phenotypes:
    """Simulate quantitative traits (generate_phenotypes.r, phase 1)."""
    input:
        geno     = config["geno_file"],
        snp_file = config["snp_file"],
    output:
        done        = touch(f"logs/simulate_phenotypes_{SEED}.done"),
        centered    = GENO_CENTERED_FILE,
        allele_freqs = ALLELE_FREQ_FILE,
    log:
        f"logs/simulate_phenotypes_{SEED}.log",
    conda:
        "envs/r_env.yml",
    params:
        geno_basename    = GENO_BASENAME,
        vavg_str         = VAVG_STR,
        qtl_str          = QTL_STR,
        bsense_str       = BSENSE_STR,
        subset_geno_arg  = (
            f"--subset_geno {config['subset_geno']}" if config.get("subset_geno") else ""
        ),
        subset_snps_arg  = (
            f"--subset_snps {config['subset_snps']}" if config.get("subset_snps") else ""
        ),
    shell:
        """
        Rscript scripts/generate_phenotypes.r \\
            --snp_file {input.snp_file} \\
            --geno_dir {config[geno_dir]} \\
            --geno {params.geno_basename} \\
            {params.subset_geno_arg} \\
            {params.subset_snps_arg} \\
            --geno_seed """ + SEED + """ \\
            --ploidy {config[ploidy]} \\
            --file_save_prefix """ + PREFIX + """ \\
            --single_snp_trait FALSE \\
            --remake_traits TRUE \\
            --calculate_phenos TRUE \\
            --save_test_train FALSE \\
            --vavg_ratios {params.vavg_str} \\
            --QTL_numbers {params.qtl_str} \\
            --broad_sense {params.bsense_str} \\
            --seed """ + SEED + """ \\
        > {log} 2>&1
        """

# ── Step 1b: Save train/test splits ──────────────────────────────────────────
# Centers genotypes, normalizes phenotypes, and writes train/test feather files
# to test_train_seed_{SEED}/. These files are the primary inputs for all
# downstream regression and GWAS steps.

rule split_phenotypes:
    """Write centered train/test genotypes and normalized phenotypes (phase 2)."""
    input:
        sim_done     = rules.simulate_phenotypes.output.done,
        centered     = rules.simulate_phenotypes.output.centered,
        allele_freqs = rules.simulate_phenotypes.output.allele_freqs,
        snp_file     = config["snp_file"],
    output:
        train_geno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_genotypes_centered.feather",
        test_geno   = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_genotypes_centered.feather",
        train_pheno = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_phenotypes_normalized.feather",
        test_pheno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_phenotypes_normalized.feather",
        all_pheno   = f"{TT_DIR}/{PREFIX}_seed_{SEED}_all_phenotypes_normalized.csv",
        train_ids   = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_ids.txt",
        test_ids    = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_ids.txt",
        snp_ids     = f"{TT_DIR}/{PREFIX}_seed_{SEED}_snp_ids.txt",
    log:
        f"logs/split_phenotypes_{SEED}.log",
    conda:
        "envs/r_env.yml",
    params:
        geno_prefix = GENO_PREFIX_SUBSETTED,
    shell:
        """
        Rscript scripts/generate_phenotypes.r \\
            --snp_file {input.snp_file} \\
            --geno_dir {config[geno_dir]} \\
            --geno {params.geno_prefix} \\
            --ploidy {config[ploidy]} \\
            --file_save_prefix """ + PREFIX + """ \\
            --single_snp_trait FALSE \\
            --remake_traits FALSE \\
            --calculate_phenos FALSE \\
            --save_test_train TRUE \\
            --train_prop {config[train_prop]} \\
            --seed """ + SEED + """ \\
            --normalize TRUE \\
            --allele_freqs """ + ALLELE_FREQ_BASENAME + """ \\
        >> {log} 2>&1
        Rscript -e "
            library(arrow)
            df      <- arrow::read_feather('{output.train_geno}')
            snp_ids <- setdiff(names(df), c('IID', 'FID'))
            writeLines(snp_ids, '{output.snp_ids}')
            cat('SNP IDs written:', length(snp_ids), '\\n')
        " >> {log} 2>&1
        """


# ── Step 2: sklearn regression (all methods run in parallel) ──────────────────
# Runs regress.py for each configured sklearn method. All methods use the same
# rule via the {method} wildcard. LARS output filenames automatically include the
# maxiter tag when --method lars is passed (handled inside regress.py). A
# sentinel marks completion because regress.py writes one CSV per trait.

rule regress:
    """Cross-validated sklearn regression for one method."""
    input:
        train_geno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_genotypes_centered.feather",
        test_geno   = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_genotypes_centered.feather",
        train_pheno = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_phenotypes_normalized.feather",
        test_pheno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_phenotypes_normalized.feather",
    output:
        touch(f"logs/regress_{{method}}_{SEED}.done"),
    wildcard_constraints:
        method="|".join(SKLEARN_METHODS),
    log:
        f"logs/regress_{{method}}_{SEED}.log",
    conda:
        "envs/scikit.yml",
    params:
        pheno_filter = PHENO_FILTER_ARG,
        # Pass --max-iter only for lars (other methods default to 10000 inside regress.py)
        max_iter_arg = lambda wildcards: (
            f"--max-iter {LARS_MAXITER}" if wildcards.method == "lars" else ""
        ),
    shell:
        """
        python scripts/regress.py \\
            --method {wildcards.method} \\
            --seed """ + SEED + """ \\
            --prefix """ + PREFIX + """ \\
            --n-folds {config[n_folds]} \\
            {params.max_iter_arg} \\
            {params.pheno_filter} \\
        > {log} 2>&1
        """


# ── Step 2 (pre): Install plink2 ─────────────────────────────────────────────
# Runs ensure_plink2.sh once and writes the resolved binary to bin/plink2 so
# both plink_gwas and clump_gwas can declare it as an explicit DAG dependency
# instead of relying on an implicit PATH side-effect.

rule setup_plink2:
    """Install plink2 to bin/plink2 (conda or direct download via ensure_plink2.sh)."""
    output:
        "bin/plink2",
    log:
        "logs/setup_plink2.log",
    conda:
        "envs/plink.yml",
    shell:
        """
        bash scripts/ensure_plink2.sh > {log} 2>&1
        if [[ ! -x bin/plink2 ]]; then
            mkdir -p bin
            ln -sf "$(command -v plink2)" bin/plink2
        fi
        """


# ── Step 2: GWAS via plink2 (per split, runs in parallel with regress) ────────
# Calls plink.sh gwas which (1) converts the VCF to plink binary format and
# (2) runs a linear GWAS keeping only the individuals in the given split.
# Produces one *.glm.linear file per trait in plink_outputs_{split}_seed_{SEED}/.

rule plink_gwas:
    """Run linear GWAS for one data split using plink2."""
    input:
        plink2_bin = "bin/plink2",
        vcf       = config["plink_vcf"],
        train_ids = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_ids.txt",
        test_ids  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_ids.txt",
        all_pheno = f"{TT_DIR}/{PREFIX}_seed_{SEED}_all_phenotypes_normalized.csv",
    output:
        touch(f"logs/plink_gwas_{{split}}_{SEED}.done"),
    wildcard_constraints:
        split="|".join(GWAS_SPLITS),
    log:
        f"logs/plink_gwas_{{split}}_{SEED}.log",
    conda:
        "envs/plink.yml"
    params:
        bfile       = config["plink_bfile"],
        tt_dir      = TT_DIR,
    shell:
        """
        bash scripts/plink.sh gwas \\
            --vcf {input.vcf} \\
            --bfile {params.bfile} \\
            --prefix """ + PREFIX + """ \\
            --seed """ + SEED + """ \\
            --split {wildcards.split} \\
            --test-train-dir {params.tt_dir} \\
            --chr-set {config[chr_set]} \\
        > {log} 2>&1
        """


# ── Step 4: LD clumping of plink GWAS results (per split) ────────────────────
# Loops over all *.glm.linear files in plink_outputs_{split}_seed_{SEED}/ and
# runs plink2 --clump on each, producing per-trait *.clumps files in the same
# directory. A sentinel marks completion because output count equals trait count.

rule clump_gwas:
    """LD-clump plink GWAS results for one data split."""
    input:
        gwas_done  = f"logs/plink_gwas_{{split}}_{SEED}.done",
        plink2_bin = "bin/plink2",
    output:
        touch(f"logs/clump_gwas_{{split}}_{SEED}.done"),
    wildcard_constraints:
        split="|".join(GWAS_SPLITS),
    log:
        f"logs/clump_gwas_{{split}}_{SEED}.log",
    conda:
        "envs/plink.yml"
    params:
        bfile      = config["plink_bfile"],
        gwas_dir   = lambda wildcards: f"plink_outputs_{wildcards.split}_seed_{SEED}",
        clump_r2   = CLUMP_R2,
        clump_p1   = CLUMP_P1,
        clump_bins = CLUMP_BINS,
    shell:
        """
        export PATH="${{PWD}}/bin:${{PATH}}"
        (
        ls {params.gwas_dir}/*.glm.linear \\
            | awk '{{split($1,a,".glm.linear"); print a[1]}}' \\
            | while read fle ; do
                plink2 --bfile {params.bfile} \\
                       --clump "${{fle}}.glm.linear" \\
                       --clump-unphased \\
                       --clump-r2 {params.clump_r2} \\
                       --clump-bins {params.clump_bins} \\
                       --clump-log10 output-only \\
                       --clump-p1 {params.clump_p1} \\
                       --chr-set {config[chr_set]} \\
                       --out "${{fle}}_r2_{params.clump_r2}"
              done
        ) > {log} 2>&1
        """


# ── Step 2: Optuna hyperparameter tuning for PyTorch ridge ───────────────────
# Calls run_tuning.sh, which derives the phenotype list at runtime from the
# train phenotype feather file and calls fit_linear_sgd_optuna.py for each.
# Results are JSON files in TUNING_DIR.

rule tune_pytorch:
    """Optuna hyperparameter search for PyTorch ridge regression."""
    input:
        train_geno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_genotypes_centered.feather",
        test_geno   = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_genotypes_centered.feather",
        train_pheno = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_phenotypes_normalized.feather",
        test_pheno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_phenotypes_normalized.feather",
    output:
        touch(f"logs/tune_pytorch_{SEED}.done"),
    log:
        f"logs/tune_pytorch_{SEED}.log",
    conda:
        "envs/pytorch.yml"
    shell:
        """
        N_TRIALS={config[optuna_n_trials]} bash scripts/run_tuning.sh \\
            --seed """ + SEED + """ \\
            --prefix """ + PREFIX + """ \\
            --test-train-dir """ + TT_DIR + """ \\
            --output-dir """ + TUNING_DIR + """ \\
            --val-seed """ + SEED + """ \\
        > {log} 2>&1
        """


# ── Step 3: Final PyTorch ridge fit ──────────────────────────────────────────
# Reads the best alpha and learning rate from each Optuna JSON result and runs
# a final fit for each phenotype. Writes per-trait weights and predictions CSVs
# to final_fit_results_seed_{SEED}/.

rule final_fit_pytorch:
    """Final PyTorch ridge fit using Optuna-tuned hyperparameters."""
    input:
        tuning_done = f"logs/tune_pytorch_{SEED}.done",
        train_geno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_genotypes_centered.feather",
        test_geno   = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_genotypes_centered.feather",
        train_pheno = f"{TT_DIR}/{PREFIX}_seed_{SEED}_train_phenotypes_normalized.feather",
        test_pheno  = f"{TT_DIR}/{PREFIX}_seed_{SEED}_test_phenotypes_normalized.feather",
    output:
        touch(f"logs/final_fit_pytorch_{SEED}.done"),
    log:
        f"logs/final_fit_pytorch_{SEED}.log",
    conda:
        "envs/pytorch.yml"
    shell:
        """
        bash scripts/run_final_fit.sh \\
            --seed """ + SEED + """ \\
            --prefix """ + PREFIX + """ \\
            --test-train-dir """ + TT_DIR + """ \\
            --tuning-dir """ + TUNING_DIR + """ \\
            --output-dir """ + FINAL_DIR + """ \\
            --torch-seed """ + SEED + """ \\
        > {log} 2>&1
        """


# ── Step 4: Aggregate sklearn outputs (non-lars) ──────────────────────────────
# Collects per-trait coefficient and prediction CSVs from intermediates_seed_{SEED}/
# and stacks them into SNPs × traits and individuals × traits Feather matrices.

rule aggregate_sklearn:
    """Aggregate per-trait sklearn outputs into coefficient and prediction matrices."""
    input:
        f"logs/regress_{{method}}_{SEED}.done",
    output:
        coeffs = (
            "aggregate_outputs_{method}_test_seed_" + SEED + "/"
            "coeff_matrix_{method}_test_seed_" + SEED + ".feather"
        ),
        preds  = (
            "aggregate_outputs_{method}_test_seed_" + SEED + "/"
            "prediction_matrix_{method}_test_seed_" + SEED + ".feather"
        ),
        params = (
            "aggregate_outputs_{method}_test_seed_" + SEED + "/"
            "params_{method}_test_seed_" + SEED + ".csv"
        ),
    wildcard_constraints:
        # Lars is excluded — handled separately because its filenames have a maxiter tag
        method="|".join(NON_LARS),
    log:
        f"logs/aggregate_sklearn_{{method}}_{SEED}.log",
    conda:
        "envs/scikit.yml"
    shell:
        """
        python scripts/aggregate_regression_outputs.py \\
            --method {wildcards.method} \\
            --seed """ + SEED + """ \\
            --prefix """ + PREFIX + """ \\
        > {log} 2>&1
        """


# ── Step 4: Aggregate lars outputs ────────────────────────────────────────────
# Lars output filenames include a maxiter tag (e.g. lars_maxiter1000), so this
# rule is separate from aggregate_sklearn to pass --maxiter explicitly.

rule aggregate_lars:
    """Aggregate per-trait LARS outputs (filenames include maxiter tag)."""
    input:
        f"logs/regress_lars_{SEED}.done",
    output:
        coeffs = (
            f"aggregate_outputs_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}/"
            f"coeff_matrix_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}.feather"
        ),
        preds  = (
            f"aggregate_outputs_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}/"
            f"prediction_matrix_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}.feather"
        ),
        params = (
            f"aggregate_outputs_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}/"
            f"params_lars_maxiter{LARS_MAXITER}_test_seed_{SEED}.csv"
        ),
    log:
        f"logs/aggregate_lars_{SEED}.log",
    conda:
        "envs/scikit.yml"
    shell:
        """
        python scripts/aggregate_regression_outputs.py \\
            --method lars \\
            --seed """ + SEED + """ \\
            --prefix """ + PREFIX + """ \\
            --maxiter """ + str(LARS_MAXITER) + """ \\
        > {log} 2>&1
        """


# ── Step 4: Aggregate PyTorch ridge outputs ───────────────────────────────────
# Collects ridge_weights_*.csv and ridge_predictions_*.csv from FINAL_DIR and
# assembles coefficient and prediction matrices. pytorch_ridge uses the "train"
# data split label (due to train/test labelling convention in this pipeline).

rule aggregate_pytorch:
    """Aggregate PyTorch ridge weights and predictions into matrices."""
    input:
        f"logs/final_fit_pytorch_{SEED}.done",
    output:
        coeffs = (
            f"aggregate_outputs_pytorch_ridge_train_seed_{SEED}/"
            f"coeff_matrix_pytorch_ridge_train_seed_{SEED}.feather"
        ),
        preds  = (
            f"aggregate_outputs_pytorch_ridge_train_seed_{SEED}/"
            f"prediction_matrix_pytorch_ridge_train_seed_{SEED}.feather"
        ),
        params = (
            f"aggregate_outputs_pytorch_ridge_train_seed_{SEED}/"
            f"params_pytorch_ridge_train_seed_{SEED}.csv"
        ),
    log:
        f"logs/aggregate_pytorch_{SEED}.log",
    conda:
        "envs/scikit.yml"
    shell:
        """
        python scripts/aggregate_regression_outputs.py \\
            --method pytorch_ridge \\
            --seed """ + SEED + """ \\
            --prefix """ + PREFIX + """ \\
            --pytorch-dir """ + FINAL_DIR + """ \\
            --pytorch-params-dir """ + TUNING_DIR + """ \\
        > {log} 2>&1
        """


# ── Step 4: Build GWAS beta and p-value matrices ──────────────────────────────
# Reads all *.glm.linear files from plink_outputs_{split}_seed_{SEED}/ and
# stacks them into SNPs × traits beta and p-value matrices. Also writes a
# p-value-filtered beta matrix used for polygenic score prediction.

rule build_gwas_matrices:
    """Stack plink GWAS results into beta and p-value matrices."""
    input:
        f"logs/plink_gwas_{{split}}_{SEED}.done",
    output:
        beta          = f"aggregate_outputs_{{split}}_seed_{SEED}/beta_matrix_{{split}}_seed_{SEED}.feather",
        pvalue        = f"aggregate_outputs_{{split}}_seed_{SEED}/pvalue_matrix_{{split}}_seed_{SEED}.feather",
        beta_filtered = (
            f"aggregate_outputs_{{split}}_seed_{SEED}/"
            f"beta_matrix_filtered_{{split}}_seed_{SEED}_p{P_STR}.feather"
        ),
    wildcard_constraints:
        split="|".join(GWAS_SPLITS),
    log:
        f"logs/build_gwas_matrices_{{split}}_{SEED}.log",
    conda:
        "envs/scikit.yml"
    shell:
        """
        python scripts/build_gwas_matrices.py build \\
            --data-split {wildcards.split} \\
            --seed """ + SEED + """ \\
            --p-threshold {config[gwas_p_threshold]} \\
        > {log} 2>&1
        """


# ── Step 5: Compute polygenic scores ─────────────────────────────────────────
# Loads the p-value-filtered beta matrix and multiplies by the full centered
# genotype matrix to produce polygenic scores (individuals × traits).

rule gwas_predict:
    """Compute polygenic scores from p-value-filtered GWAS betas."""
    input:
        beta_filtered = (
            f"aggregate_outputs_{{split}}_seed_{SEED}/"
            f"beta_matrix_filtered_{{split}}_seed_{SEED}_p{P_STR}.feather"
        ),
        geno_centered = GENO_CENTERED_FILE,
    output:
        scores = (
            f"aggregate_outputs_{{split}}_seed_{SEED}/"
            f"polygenic_scores_{{split}}_seed_{SEED}_p{P_STR}.feather"
        ),
    wildcard_constraints:
        split="|".join(GWAS_SPLITS),
    log:
        f"logs/gwas_predict_{{split}}_{SEED}.log",
    conda:
        "envs/scikit.yml"
    shell:
        """
        python scripts/build_gwas_matrices.py predict \\
            --data-split {wildcards.split} \\
            --seed """ + SEED + """ \\
            --p-threshold {config[gwas_p_threshold]} \\
            --geno-matrix {input.geno_centered} \\
        > {log} 2>&1
        """


# ── Step 6: Publication figures ───────────────────────────────────────────────
# Reads pre-computed figure intermediates from FIG_INT_DIR (downloaded from
# Zenodo — see README step 1) and generates publication-quality SVG figures.
# Because inputs come from FIG_INT_DIR rather than the pipeline's aggregate
# outputs, this rule can be run standalone after downloading the intermediates
# without executing the full pipeline. A sentinel marks completion because the
# script produces many output SVG files.

rule pub_figures:
    """Generate publication figures from pre-computed figure intermediates."""
    input:
        # Yeast figure intermediates — must be downloaded from Zenodo (README step 1)
        f"{FIG_INT_DIR}/combined_all_betas_yeast_{YEAST_PLOT_SEED}.feather",
        f"{FIG_INT_DIR}/yeast_littlelonger_with_fullinfo_{YEAST_PLOT_SEED}.feather",
        f"{FIG_INT_DIR}/yeast_cumulative_{YEAST_PLOT_SEED}.feather",
        f"{FIG_INT_DIR}/yeast_roc_{YEAST_PLOT_SEED}.feather",
        # Human figure intermediates — private UK Biobank data, never produced by
        # this pipeline. Must be downloaded from Zenodo before running pub_figures.
        # Skipped when human_figure_seed is not set in config.
        (
            [
                f"{FIG_INT_DIR}/combined_all_betas_human_{HUMAN_PLOT_SEED}.feather",
                f"{FIG_INT_DIR}/human_littlelonger_with_fullinfo_{HUMAN_PLOT_SEED}.feather",
                f"{FIG_INT_DIR}/human_true_avg_distance_{HUMAN_PLOT_SEED}.feather",
                f"{FIG_INT_DIR}/all_methods_effects_correlations_{YEAST_FIG_PREFIX}_{YEAST_PLOT_SEED}_{HUMAN_FIG_PREFIX}_{HUMAN_PLOT_SEED}.feather",
                f"{FIG_INT_DIR}/all_methods_max_prediction_correlations_{YEAST_FIG_PREFIX}_{YEAST_PLOT_SEED}_{HUMAN_FIG_PREFIX}_{HUMAN_PLOT_SEED}.feather",
                f"{FIG_INT_DIR}/all_methods_parameters_{YEAST_FIG_PREFIX}_{YEAST_PLOT_SEED}_{HUMAN_FIG_PREFIX}_{HUMAN_PLOT_SEED}.feather",
                f"{FIG_INT_DIR}/all_methods_prediction_correlations_{YEAST_FIG_PREFIX}_{YEAST_PLOT_SEED}_{HUMAN_FIG_PREFIX}_{HUMAN_PLOT_SEED}.feather",
            ]
            if config.get("human_figure_seed")
            else []
        ),
    output:
        touch(f"logs/pub_figures_{YEAST_PLOT_SEED}_{HUMAN_PLOT_SEED}.done"),
    log:
        f"logs/pub_figures_{YEAST_PLOT_SEED}_{HUMAN_PLOT_SEED}.log",
    conda:
        "envs/r_env.yml",
    params:
        output_dir     = config.get("figures_output_dir", "plots"),
        yeast_seed     = config["yeast_figure_seed"],
        yeast_prefix   = config["yeast_figure_prefix"],
        yeast_p_thr    = P_STR,
        human_seed_arg = (
            f"--human-seed {config['human_figure_seed']}"
            if config.get("human_figure_seed")
            else ""
        ),
        human_prefix   = config.get("human_figure_prefix", "human"),
        human_bim_arg  = (
            f"--human-bim-file {config['human_bim_file']}"
            if config.get("human_bim_file")
            else ""
        ),
    shell:
        """
        Rscript scripts/pub_figures_snakemake.r \\
            --base-dir . \\
            --output-dir {params.output_dir} \\
            --yeast-seed {params.yeast_seed} \\
            --yeast-prefix {params.yeast_prefix} \\
            --yeast-plink-threshold {params.yeast_p_thr} \\
            {params.human_seed_arg} \\
            --human-prefix {params.human_prefix} \\
            {params.human_bim_arg} \\
            --lars-maxiter {config[lars_maxiter]} \\
        > {log} 2>&1
        """
