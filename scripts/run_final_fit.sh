#!/bin/bash
# run_final_fit.sh — Final PyTorch ridge fit using Optuna-tuned hyperparameters.
#
# Reads every phenotype column from the train phenotype Feather file and calls
# fit_linear_sgd_cli.py once per phenotype, picking up the best alpha and
# learning rate that were saved by run_tuning.sh / fit_linear_sgd_optuna.py.
#
# Run from the repository root with the pytorch conda environment active:
#   conda activate pytorch
#
# ─────────────────────────────────────────────────────────────────────────────
# Usage
# ─────────────────────────────────────────────────────────────────────────────
#
#   bash scripts/run_final_fit.sh \
#       --seed SEED \
#       --prefix PREFIX \
#       [--test-train-dir DIR] \
#       [--tuning-dir DIR] \
#       [--output-dir DIR] \
#       [--torch-seed INT]
#
# ─────────────────────────────────────────────────────────────────────────────
# Flags
# ─────────────────────────────────────────────────────────────────────────────
#
#   --seed           STR   Seed identifier embedded in input file names (required)
#   --prefix         STR   Filename prefix for Feather data files (required)
#   --test-train-dir DIR   Directory containing train/test Feather files
#                          (default: test_train_seed_<seed>)
#   --tuning-dir     DIR   Directory containing Optuna tuning results produced
#                          by run_tuning.sh (default: pytorch_tuning_results)
#   --output-dir     DIR   Output directory for per-phenotype fit results
#                          (default: final_fit_results_seed_<seed>)
#   --torch-seed     INT   Random seed for PyTorch weight initialisation and
#                          DataLoader shuffling (default: 42)
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────

SEED=""
PREFIX=""
TEST_TRAIN_DIR=""
TUNING_DIR=""
OUTPUT_DIR=""
TORCH_SEED=42

while [[ $# -gt 0 ]]; do
    case "$1" in
        --seed)           SEED="$2";           shift 2 ;;
        --prefix)         PREFIX="$2";         shift 2 ;;
        --test-train-dir) TEST_TRAIN_DIR="$2"; shift 2 ;;
        --tuning-dir)     TUNING_DIR="$2";     shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2";     shift 2 ;;
        --torch-seed)     TORCH_SEED="$2";     shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$SEED" ]]   || { echo "ERROR: --seed is required"   >&2; exit 1; }
[[ -n "$PREFIX" ]] || { echo "ERROR: --prefix is required" >&2; exit 1; }

TEST_TRAIN_DIR="${TEST_TRAIN_DIR:-test_train_seed_${SEED}}"
TUNING_DIR="${TUNING_DIR:-pytorch_tuning_results}"
OUTPUT_DIR="${OUTPUT_DIR:-final_fit_results_seed_${SEED}}"

# ── File paths ────────────────────────────────────────────────────────────────

TRAIN_GENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_train_genotypes_centered.feather"
TEST_GENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_test_genotypes_centered.feather"
TRAIN_PHENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_train_phenotypes_normalized.feather"
TEST_PHENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_test_phenotypes_normalized.feather"

# ── Phenotype list (derived from the feather file at runtime) ─────────────────

[[ -f "$TRAIN_PHENO" ]] || { echo "ERROR: phenotype file not found: $TRAIN_PHENO" >&2; exit 1; }

# mapfile is bash 4+ only; macOS ships bash 3.2, so use a while-read loop.
# Use 'python' (not 'python3') — conda envs may not create the python3 symlink.
PHENOTYPES=()
while IFS= read -r line; do
    [[ -n "$line" ]] && PHENOTYPES+=("$line")
done < <(python -c "
import pandas as pd
df = pd.read_feather('${TRAIN_PHENO}')
id_cols = {'IID', 'FID', '#IID', '#FID'}
for c in df.columns:
    if c not in id_cols:
        print(c)
")

[[ ${#PHENOTYPES[@]} -gt 0 ]] || { echo "ERROR: no phenotype columns found in $TRAIN_PHENO" >&2; exit 1; }
echo "Found ${#PHENOTYPES[@]} phenotype(s) in $TRAIN_PHENO"

# ── Run ───────────────────────────────────────────────────────────────────────

n=${#PHENOTYPES[@]}
for i in "${!PHENOTYPES[@]}"; do
    pheno="${PHENOTYPES[$i]}"
    echo "[$(( i + 1 ))/${n}] ${pheno}"

    python scripts/fit_linear_sgd_cli.py \
        --train-geno     "${TRAIN_GENO}" \
        --test-geno      "${TEST_GENO}" \
        --train-pheno    "${TRAIN_PHENO}" \
        --test-pheno     "${TEST_PHENO}" \
        --phenotype-name "${pheno}" \
        --tuning-dir     "${TUNING_DIR}" \
        --output-dir     "${OUTPUT_DIR}" \
        --which-seed     "${SEED}" \
        --torch-seed     "${TORCH_SEED}"

    echo "---"
done

echo "Done. Final fit complete for ${n} phenotype(s)."
