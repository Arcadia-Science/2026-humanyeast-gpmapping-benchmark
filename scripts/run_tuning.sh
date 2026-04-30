#!/bin/bash
# run_tuning.sh — Optuna hyperparameter tuning for PyTorch ridge regression.
#
# Reads every phenotype column from the train phenotype Feather file and calls
# fit_linear_sgd_optuna.py once per phenotype to search for the best alpha and
# learning rate via Bayesian optimisation.  Each trial is wrapped with `timeout`
# so a slow phenotype cannot block the whole run indefinitely.
#
# Tuning hyperparameters (alpha range, learning-rate range, trial count, etc.)
# are controlled by environment variables so they can be overridden without
# editing the script.
#
# Run from the repository root with the pytorch conda environment active:
#   conda activate pytorch
#
# ─────────────────────────────────────────────────────────────────────────────
# Usage
# ─────────────────────────────────────────────────────────────────────────────
#
#   bash scripts/run_tuning.sh \
#       --seed SEED \
#       --prefix PREFIX \
#       [--test-train-dir DIR] \
#       [--output-dir DIR] \
#       [--val-seed INT]
#
# ─────────────────────────────────────────────────────────────────────────────
# Flags
# ─────────────────────────────────────────────────────────────────────────────
#
#   --seed           STR   Seed identifier embedded in input file names (required)
#   --prefix         STR   Filename prefix for Feather data files (required)
#   --test-train-dir DIR   Directory containing train/test Feather files
#                          (default: test_train_seed_<seed>)
#   --output-dir     DIR   Output directory for Optuna tuning results
#                          (default: pytorch_tuning_results)
#   --val-seed       INT   Random seed for val split, TPE sampler, and torch
#                          (default: 42)
#
# ─────────────────────────────────────────────────────────────────────────────
# Environment variables (all optional — defaults shown)
# ─────────────────────────────────────────────────────────────────────────────
#
#   N_TRIALS     INT    Optuna trials per phenotype (default: 30)
#   N_JOBS       INT    Parallel Optuna workers per phenotype (default: 1)
#   MAX_EPOCHS   INT    Maximum SGD epochs per trial (default: 12)
#   ALPHA_MIN    FLOAT  Lower bound of ridge alpha search range (default: 0.001)
#   ALPHA_MAX    FLOAT  Upper bound of ridge alpha search range (default: 3)
#   LR_MIN       FLOAT  Lower bound of learning-rate search range (default: 1e-6)
#   LR_MAX       FLOAT  Upper bound of learning-rate search range (default: 1e-3)
#   TIMEOUT      INT    Seconds before Optuna stops accepting new trials;
#                       the process is hard-killed TIMEOUT+60 s after start
#                       (default: 7500)
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Arguments ─────────────────────────────────────────────────────────────────

SEED=""
PREFIX=""
TEST_TRAIN_DIR=""
OUTPUT_DIR=""
VAL_SEED=42

while [[ $# -gt 0 ]]; do
    case "$1" in
        --seed)           SEED="$2";           shift 2 ;;
        --prefix)         PREFIX="$2";         shift 2 ;;
        --test-train-dir) TEST_TRAIN_DIR="$2"; shift 2 ;;
        --output-dir)     OUTPUT_DIR="$2";     shift 2 ;;
        --val-seed)       VAL_SEED="$2";       shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$SEED" ]]   || { echo "ERROR: --seed is required"   >&2; exit 1; }
[[ -n "$PREFIX" ]] || { echo "ERROR: --prefix is required" >&2; exit 1; }

TEST_TRAIN_DIR="${TEST_TRAIN_DIR:-test_train_seed_${SEED}}"
OUTPUT_DIR="${OUTPUT_DIR:-pytorch_tuning_results}"

# ── Tuning hyperparameters ────────────────────────────────────────────────────

N_TRIALS=${N_TRIALS:-30}
N_JOBS=${N_JOBS:-1}
MAX_EPOCHS=${MAX_EPOCHS:-12}
ALPHA_MIN=${ALPHA_MIN:-0.001}
ALPHA_MAX=${ALPHA_MAX:-3}
LR_MIN=${LR_MIN:-1e-6}
LR_MAX=${LR_MAX:-1e-3}
# TIMEOUT controls both Optuna (graceful stop) and the shell hard-kill wrapper.
# Optuna stops accepting new trials after TIMEOUT seconds; the shell kills the
# process TIMEOUT + 60 seconds after it starts.
TIMEOUT=${TIMEOUT:-7500}

# ── File paths ────────────────────────────────────────────────────────────────

# NOTE: train/test file labels are intentionally swapped here — see regress.py
# for explanation of this convention.
TRAIN_GENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_test_genotypes_centered.feather"
TRAIN_PHENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_test_phenotypes_normalized.feather"
TEST_GENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_train_genotypes_centered.feather"
TEST_PHENO="${TEST_TRAIN_DIR}/${PREFIX}_seed_${SEED}_train_phenotypes_normalized.feather"

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

    timeout --kill-after=60s "${TIMEOUT}s" \
        python scripts/fit_linear_sgd_optuna.py \
            --train-geno      "${TRAIN_GENO}" \
            --test-geno       "${TEST_GENO}" \
            --train-pheno     "${TRAIN_PHENO}" \
            --test-pheno      "${TEST_PHENO}" \
            --phenotype-name  "${pheno}" \
            --reg-mode        ridge \
            --output-dir      "${OUTPUT_DIR}" \
            --alpha-min       "${ALPHA_MIN}" \
            --alpha-max       "${ALPHA_MAX}" \
            --lr-min          "${LR_MIN}" \
            --lr-max          "${LR_MAX}" \
            --n-trials        "${N_TRIALS}" \
            --n-jobs          "${N_JOBS}" \
            --max-epochs      "${MAX_EPOCHS}" \
            --timeout         "${TIMEOUT}" \
            --val-seed        "${VAL_SEED}" \
            --pruning \
            --verbose

    echo "---"
done

echo "Done. Tuning complete for ${n} phenotype(s)."
