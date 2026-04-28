#!/bin/bash
# plink.sh — PLINK2 tasks for genomic prediction
#
# Three subcommands, each wrapping one or more plink2 calls:
#
#   gwas   Convert a VCF to binary format and run linear GWAS
#   ld     Compute LD: pruning, LD with lead SNPs, and unfiltered pairwise LD
#   clump  Clump GWAS results for every trait in a results directory
#
# Run from the repository root with the plink conda environment active:
#   conda activate plink
#
# ─────────────────────────────────────────────────────────────────────────────
# Usage
# ─────────────────────────────────────────────────────────────────────────────
#
#   gwas
#     bash scripts/plink.sh gwas \
#       --vcf input_data/yeast_genotypes_binarized.vcf \
#       --bfile input_data/yeast_genotypes_binarized \
#       --prefix yeast_simulated_data \
#       --seed 1510 \
#       --split test \
#       --test-train-dir test_train_seed_1510 \
#       --chr-set -16
#
#   ld
#     bash scripts/plink.sh ld \
#       --bfile input_data/yeast_genotypes_binarized \
#       --out-dir yeast_plink_ld
#
#   clump
#     bash scripts/plink.sh clump \
#       --bfile input_data/yeast_genotypes_binarized \
#       --seed 1510 \
#       --split train \
#       --clump-r2 0.5
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Ensure plink2 is installed ────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "${SCRIPT_DIR}/ensure_plink2.sh"
export PATH="${PWD}/bin:${PATH}"   # pick up binary if it was downloaded locally

# ── Helpers ───────────────────────────────────────────────────────────────────

usage() {
    sed -n '/^# Usage/,/^# ───/p' "$0" | grep -v '^# ───' | sed 's/^# \?//'
    exit 1
}

die() { echo "[ERROR] $*" >&2; exit 1; }

require_arg() {
    # require_arg <value> <flag-name>
    [[ -n "$1" ]] || die "$2 is required for this subcommand."
}

# ── Parse subcommand ──────────────────────────────────────────────────────────

[[ $# -ge 1 ]] || usage
subcommand="$1"
shift

# ─────────────────────────────────────────────────────────────────────────────
# gwas — VCF → binary + linear GWAS
# ─────────────────────────────────────────────────────────────────────────────

cmd_gwas() {
    local vcf="" bfile="" prefix="yeast_simulated_data" seed="1510"
    local split="test" test_train_dir="" chr_set="-16" results_dir=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --vcf)            vcf="$2";            shift 2 ;;
            --bfile)          bfile="$2";          shift 2 ;;
            --prefix)         prefix="$2";         shift 2 ;;
            --seed)           seed="$2";           shift 2 ;;
            --split)          split="$2";          shift 2 ;;
            --test-train-dir) test_train_dir="$2"; shift 2 ;;
            --chr-set)        chr_set="$2";        shift 2 ;;
            --results-dir)    results_dir="$2";    shift 2 ;;
            -h|--help)        usage ;;
            *) die "Unknown argument: $1" ;;
        esac
    done

    require_arg "$vcf"    "--vcf"
    require_arg "$bfile"  "--bfile"

    local tt_dir="${test_train_dir:-test_train_seed_${seed}}"
    results_dir="${results_dir:-plink_outputs_${split}_seed_${seed}}"

    local inds_file="./${tt_dir}/${prefix}_seed_${seed}_${split}_ids.txt"
    local pheno_file="./${tt_dir}/${prefix}_seed_${seed}_all_phenotypes_normalized.csv"

    echo "─────────────────────────────────────────────────────────────────────"
    echo "GWAS"
    echo "  VCF           : $vcf"
    echo "  Binary plink  : $bfile"
    echo "  Individuals   : $inds_file"
    echo "  Phenotypes    : $pheno_file"
    echo "  Results dir   : $results_dir"
    echo "  Seed          : $seed   Split: $split   Chr-set: $chr_set"
    echo "─────────────────────────────────────────────────────────────────────"

    [[ -f "$inds_file" ]]  || die "Individuals file not found: $inds_file"
    [[ -f "$pheno_file" ]] || die "Phenotypes file not found: $pheno_file"

    # Convert VCF to binary format (skip if already done)
    if [[ ! -f "${bfile}.bed" ]]; then
        # plink2 on macOS exits non-zero because it cannot delete its own
        # intermediate temp files (a known macOS-specific cleanup bug).
        # The .bed/.bim/.fam outputs are written successfully beforehand.
        # Filter the spurious deletion errors from the log; verify outputs exist.
        plink2 --vcf "$vcf" --make-bed --out "$bfile" 2>&1 \
            | grep -v "^Error: Failed to delete.*-temporary" \
            || true
        [[ -f "${bfile}.bed" ]] || die "plink2 --make-bed failed to create ${bfile}.bed"
    else
        echo "Binary plink files already exist at ${bfile}.* — skipping VCF conversion."
    fi

    # plink2 writes -9 in the phenotype column of the .fam file by default;
    # replace with 0 so downstream tools see "no phenotype" rather than a
    # missing-value code. Use $$ (PID) in the temp filename to avoid a race
    # condition when train and test splits run in parallel — both use the same
    # bfile, so they must not share the same .fam.tmp path.
    local fam_tmp="${bfile}.fam.$$.tmp"
    awk '{$NF = 0; print}' "${bfile}.fam" > "${fam_tmp}" \
        && mv "${fam_tmp}" "${bfile}.fam"

    # Linear GWAS (no covariate filtering)
    mkdir -p "$results_dir"

    # plink2 needs "#IID" (with leading #) to recognise it as a single-column
    # ID file and match by IID only (FID defaults to 0, matching the 0 FID that
    # plink2 assigns when importing from VCF).
    # Note: BSD sed on macOS does not interpret \t as a tab, so we prepend '#'
    # to line 1 without any tab matching.
    local pheno_plink="${results_dir}/pheno_for_plink.txt"
    sed '1s/^/#/' "${pheno_file}" > "${pheno_plink}"

    plink2 --bfile "$bfile" \
        --no-input-missing-phenotype \
        --seed "$seed" \
        --keep "$inds_file" \
        --glm allow-no-covars \
        --pheno "${pheno_plink}" \
        --chr-set "$chr_set" \
        --out "${results_dir}/plink_nofiltering"

    echo "Done. Results in: $results_dir"
}

# ─────────────────────────────────────────────────────────────────────────────
# ld — LD pruning, LD with lead SNPs, unfiltered pairwise LD
# ─────────────────────────────────────────────────────────────────────────────

cmd_ld() {
    local bfile="" out_dir="yeast_plink_ld" window_kb="250" unfiltered_kb="300"

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --bfile)          bfile="$2";          shift 2 ;;
            --out-dir)        out_dir="$2";        shift 2 ;;
            --window-kb)      window_kb="$2";      shift 2 ;;
            --unfiltered-kb)  unfiltered_kb="$2";  shift 2 ;;
            -h|--help)        usage ;;
            *) die "Unknown argument: $1" ;;
        esac
    done

    require_arg "$bfile" "--bfile"
    mkdir -p "$out_dir"

    local stem="${out_dir}/$(basename "$out_dir")"

    echo "─────────────────────────────────────────────────────────────────────"
    echo "LD"
    echo "  Binary plink  : $bfile"
    echo "  Output dir    : $out_dir"
    echo "  Prune window  : ${window_kb}kb   Unfiltered window: ${unfiltered_kb}kb"
    echo "─────────────────────────────────────────────────────────────────────"

    # 1. Pruning — identify lead SNPs
    echo "Step 1/3: LD pruning (${window_kb}kb, r2 < 0.8)..."
    plink2 --bfile "$bfile" \
        --indep-pairwise "${window_kb}kb" 0.8 \
        --out "${stem}_leadsnps"

    # 2. LD between all SNPs and lead SNPs
    echo "Step 2/3: LD with lead SNPs (${window_kb}kb window)..."
    plink2 --bfile "$bfile" \
        --r2-unphased \
        --ld-snp-list "${stem}_leadsnps.prune.in" \
        --ld-window-kb "$window_kb" \
        --ld-window-r2 0 \
        --out "${stem}_ld"

    # 3. Unfiltered pairwise LD
    echo "Step 3/3: Unfiltered pairwise LD (${unfiltered_kb}kb window)..."
    plink2 --bfile "$bfile" \
        --r2-unphased \
        --ld-window-kb "$unfiltered_kb" \
        --ld-window-r2 0 \
        --out "${stem}_ld_unfiltered"

    echo "Done. LD files in: $out_dir"
}

# ─────────────────────────────────────────────────────────────────────────────
# clump — Clump GWAS results for every trait
# ─────────────────────────────────────────────────────────────────────────────

cmd_clump() {
    local bfile="" seed="1510" split="train" clump_r2="0.5" results_dir=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --bfile)       bfile="$2";       shift 2 ;;
            --seed)        seed="$2";        shift 2 ;;
            --split)       split="$2";       shift 2 ;;
            --clump-r2)    clump_r2="$2";    shift 2 ;;
            --results-dir) results_dir="$2"; shift 2 ;;
            -h|--help)     usage ;;
            *) die "Unknown argument: $1" ;;
        esac
    done

    require_arg "$bfile" "--bfile"

    results_dir="${results_dir:-plink_outputs_${split}_seed_${seed}}"

    [[ -d "$results_dir" ]] || die "Results directory not found: $results_dir"

    local glm_files=("${results_dir}"/*.glm.linear)
    [[ -e "${glm_files[0]}" ]] || die "No .glm.linear files found in: $results_dir"

    echo "─────────────────────────────────────────────────────────────────────"
    echo "Clump"
    echo "  Binary plink  : $bfile"
    echo "  Results dir   : $results_dir"
    echo "  Seed          : $seed   Split: $split   Clump r2: $clump_r2"
    echo "  Traits        : ${#glm_files[@]}"
    echo "─────────────────────────────────────────────────────────────────────"

    local n=0
    for glm_file in "${glm_files[@]}"; do
        # Strip .glm.linear to get the per-trait output prefix
        local out_prefix="${glm_file%.glm.linear}_r2_${clump_r2}"
        n=$((n + 1))
        echo "[$n/${#glm_files[@]}] $(basename "$glm_file")"
        plink2 --bfile "$bfile" \
            --clump "$glm_file" \
            --clump-unphased \
            --clump-r2 "$clump_r2" \
            --clump-bins 1.2e-6 1e-6 1e-5 1e-4 1e-3 1e-2 .05 \
            --clump-log10 output-only \
            --clump-p1 1.2e-6 \
            --out "$out_prefix"
    done

    echo "Done. Clumped $n trait(s)."
}

# ─────────────────────────────────────────────────────────────────────────────
# Dispatch
# ─────────────────────────────────────────────────────────────────────────────

case "$subcommand" in
    gwas)  cmd_gwas  "$@" ;;
    ld)    cmd_ld    "$@" ;;
    clump) cmd_clump "$@" ;;
    -h|--help) usage ;;
    *) die "Unknown subcommand: '$subcommand'. Choose from: gwas, ld, clump" ;;
esac
