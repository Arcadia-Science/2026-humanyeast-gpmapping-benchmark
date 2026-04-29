#!/usr/bin/env python3
"""
aggregate_regression_outputs.py

Aggregates coefficient and prediction matrices from regression methods
(ridge, lasso, elasticnet, lars, pytorch_ridge) across all traits for a
given seed, and saves them as feather files.

─── SUPPORTED METHODS ──────────────────────────────────────────────────────

  sklearn methods  : ridge, lasso, elasticnet, lars
    Input dir      : intermediates_seed_{seed}/
    File prefixes  :
      ridge        : ridge_coeffs_*, ridge_predictions_*, ridge_alphas_*
      lasso        : lasso_coeffs_*, lasso_predictions_*, lasso_alphas_*
      elasticnet  : elasticnet_coeffs_*, elasticnet_predictions_*, elasticnet_alphas_*
      lars         : lars_maxiter{N}_coeffs_*, lars_maxiter{N}_predictions_*, lars_maxiter{N}_alphas_*
    Data split     : test
    Coeffs files   : 1 col (coefficient) or 2 cols (SNP, coefficient) with/without header
    Predictions    : 1 col (prediction) or 2 cols (IID, prediction) with/without header

  pytorch_ridge    :
    Weights/preds  : final_fit_results_seed_{seed}/
      ridge_weights_seed_{seed}_{trait}.csv  — cols: snp_id, weight (with header)
      ridge_predictions_seed_{seed}_{trait}.csv — cols: IID, true_phenotype, predicted_phenotype
    Params dir     : yeast_tuning_results_test/
      optuna_best_pheno_{idx}_ridge.json  — keys: phenotype_name, best_alpha, best_learning_rate
    Data split     : train

─── OUTPUT FILES ───────────────────────────────────────────────────────────

All outputs are saved to: aggregate_outputs_{method}_seed_{seed}/

  coeff_matrix_{method}_{split}_seed_{seed}[_maxiter{N}].feather
      SNPs x traits matrix of coefficients/weights

  prediction_matrix_{method}_{split}_seed_{seed}[_maxiter{N}].feather
      individuals x traits matrix of predictions

  params_{method}_{split}_seed_{seed}[_maxiter{N}].csv
      One row per trait with columns: trait, alpha, [l1_ratio, learning_rate]
      (NA where not applicable)

─── USAGE ──────────────────────────────────────────────────────────────────

  python aggregate_regression_outputs.py --method ridge --seed 1105
  python aggregate_regression_outputs.py --method lars --seed 1105 --maxiter 1000
  python aggregate_regression_outputs.py \\
      --method pytorch_ridge --seed 1510 \\
      --pytorch-dir final_fit_results_seed_1510 \\
      --pytorch-params-dir yeast_tuning_results_test

Options:
  --method          STR    One of: ridge, lasso, elasticnet, lars, pytorch_ridge (required)
  --seed            INT    Random seed (required)
  --base-dir        DIR    Parent dir containing intermediates_seed_* (default: .)
  --output-base-dir DIR    Parent dir for aggregate_outputs_* (default: same as --base-dir)
  --prefix          STR    Filename prefix for ID/SNP files (default: yeast_simulated_data)

  --maxiter         INT    Max iterations suffix for lars (e.g. 1000 -> maxiter1000)
  --pytorch-dir     DIR    Directory containing pytorch_ridge weights/predictions files
                           (required for pytorch_ridge)
  --pytorch-params-dir DIR Directory containing optuna JSON param files
                           (required for pytorch_ridge)
  --threads         INT    Parallel threads for loading files (default: 4)
"""

import argparse
import glob
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import pyarrow.feather as feather


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SKLEARN_METHODS  = {'ridge', 'lasso', 'elasticnet', 'lars'}
PYTORCH_METHODS  = {'pytorch_ridge'}
ALL_METHODS      = SKLEARN_METHODS | PYTORCH_METHODS
SKLEARN_SPLIT    = 'test'
PYTORCH_SPLIT    = 'train'


# ---------------------------------------------------------------------------
# Naming helpers
# ---------------------------------------------------------------------------

def make_label(method: str, seed: int, maxiter: int = None) -> str:
    """Return the core label used in filenames.

    Standard : 'ridge_test_seed_1105'
    Lars     : 'lars_maxiter1000_test_seed_1105'
    """
    split = PYTORCH_SPLIT if method in PYTORCH_METHODS else SKLEARN_SPLIT
    if maxiter is not None:
        return f"{method}_maxiter{maxiter}_{split}_seed_{seed}"
    return f"{method}_{split}_seed_{seed}"


def output_dir(base: str, method: str, seed: int, maxiter: int = None) -> str:
    split = PYTORCH_SPLIT if method in PYTORCH_METHODS else SKLEARN_SPLIT
    label = f"{method}_maxiter{maxiter}" if maxiter is not None else method
    return os.path.join(base, f"aggregate_outputs_{label}_{split}_seed_{seed}")


def maxiter_tag(maxiter: int = None) -> str:
    return f"maxiter{maxiter}_" if maxiter is not None else ""


# ---------------------------------------------------------------------------
# Generic file-loading helpers
# ---------------------------------------------------------------------------

def load_feather_or_csv(path: str, id_col: str = None) -> pd.DataFrame:
    """Load feather/csv/tsv, optionally setting id_col as index."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.feather':
        df = feather.read_feather(path)
    elif ext == '.csv':
        df = pd.read_csv(path, low_memory=False)
    elif ext in ('.tsv', '.txt'):
        df = pd.read_csv(path, sep='\t', low_memory=False)
    else:
        try:
            df = feather.read_feather(path)
        except Exception:
            df = pd.read_csv(path, low_memory=False)
    if id_col:
        if id_col not in df.columns:
            sys.exit(f"[ERROR] Column '{id_col}' not found in {path}")
        df = df.set_index(id_col)
    return df


def read_snp_order(seed: int, base_dir: str, file_prefix: str) -> list:
    """Read SNP order from test_train_seed_{seed}/{file_prefix}_seed_{seed}_snp_ids.txt"""
    path = os.path.join(base_dir, f"test_train_seed_{seed}",
                        f"{file_prefix}_seed_{seed}_snp_ids.txt")
    if not os.path.isfile(path):
        sys.exit(f"[ERROR] SNP IDs file not found: {path}")
    with open(path) as f:
        snps = [line.strip() for line in f if line.strip()]
    print(f"  SNP order loaded: {len(snps)} SNPs from {path}")
    return snps


def read_iid_order(seed: int, split: str, base_dir: str, file_prefix: str) -> list:
    """Read individual order from test_train_seed_{seed}/.

    split='test' -> {file_prefix}_seed_{seed}_train_ids.txt (held-out from train)
    split='train' -> {file_prefix}_seed_{seed}_test_ids.txt (held-out from test)
    """
    id_file = 'train_ids.txt' if split == SKLEARN_SPLIT else 'test_ids.txt'
    path = os.path.join(base_dir, f"test_train_seed_{seed}",
                        f"{file_prefix}_seed_{seed}_{id_file}")
    if not os.path.isfile(path):
        sys.exit(f"[ERROR] IID file not found: {path}")
    df = pd.read_csv(path, sep=' ')
    # Header is either "#FID" + "IID" (two cols) or "#IID" (one col)
    df.columns = [c.lstrip('#') for c in df.columns]

    if 'IID' in df.columns:
        iids = list(df['IID'].astype(str))
    else:
        # Single column named something like "IID" after stripping #
        iids = list(df.iloc[:, 0].astype(str))
    print(f"  IID order loaded: {len(iids)} individuals from {path}")
    return iids


# ---------------------------------------------------------------------------
# Coeffs / predictions file parsers
# ---------------------------------------------------------------------------

def parse_coeffs_file(path: str, expected_snps: list) -> pd.Series:
    """
    Parse a coefficients CSV file.
    Handles:
      - 2 columns with header (SNP, coefficient)
      - 1 column with no header (coefficient only, order matches expected_snps)
    Returns a Series indexed by SNP.
    """
    # Peek at first line to detect header
    with open(path) as f:
        first_line = f.readline().strip()

    n_cols = len(first_line.split(','))

    if n_cols >= 2:
        # Has header + SNP column
        df = pd.read_csv(path)
        # Find the coefficient column (last numeric column)
        snp_col   = df.columns[0]
        coeff_col = df.columns[1]
        series = pd.to_numeric(df[coeff_col], errors='coerce')
        series.index = df[snp_col].values
    else:
        # No header, single column — use expected_snps for index
        df = pd.read_csv(path, header=None, names=['coeff'])
        if len(df) != len(expected_snps):
            print(f"  [WARN] {os.path.basename(path)}: "
                  f"got {len(df)} coefficients, expected {len(expected_snps)} SNPs — "
                  f"returning NaN column.", file=sys.stderr)
            return pd.Series(np.nan, index=expected_snps)
        series = pd.to_numeric(df['coeff'], errors='coerce')
        series.index = expected_snps

    return series


def extract_iid(raw_id: str) -> str:
    """Extract just the IID from a raw ID string.

    If the value is "FID IID" (space-separated pair), return only the IID.
    Otherwise return the value as-is.
    """
    parts = str(raw_id).strip().split()
    return parts[-1] if len(parts) == 2 else str(raw_id).strip()


def parse_predictions_file(path: str, expected_iids: list) -> pd.Series:
    """
    Parse a predictions CSV file.
    Handles:
      - 2 columns with header (IID, prediction)
      - 1 column with no header (prediction only, order matches expected_iids)
    IID values of the form "FID IID" are reduced to just the IID part.
    Returns a Series indexed by IID.
    """
    with open(path) as f:
        first_line = f.readline().strip()

    n_cols = len(first_line.split(','))

    if n_cols >= 2:
        df = pd.read_csv(path)
        iid_col  = df.columns[0]
        pred_col = df.columns[1]
        series = pd.to_numeric(df[pred_col], errors='coerce')
        series.index = df[iid_col].apply(extract_iid).values
    else:
        df = pd.read_csv(path, header=None, names=['pred'])
        if len(df) != len(expected_iids):
            print(f"  [WARN] {os.path.basename(path)}: "
                  f"got {len(df)} predictions, expected {len(expected_iids)} individuals — "
                  f"returning NaN column.", file=sys.stderr)
            return pd.Series(np.nan, index=expected_iids)
        series = pd.to_numeric(df['pred'], errors='coerce')
        series.index = [str(i) for i in expected_iids]

    return series


def parse_alphas_file(path: str, method: str) -> dict:
    """
    Parse an alphas text file.
    ridge/lasso/lars : 2 lines  -> trait_name, alpha
    elasticnet      : 3 lines  -> trait_name, alpha, l1_ratio
    Returns dict with keys: trait, alpha, [l1_ratio]
    """
    with open(path) as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    result = {'alpha': np.nan, 'l1_ratio': np.nan, 'learning_rate': np.nan}
    if len(lines) >= 1:
        result['trait'] = lines[0]
    if len(lines) >= 2:
        try:
            result['alpha'] = float(lines[1])
        except ValueError:
            pass
    if method == 'elasticnet' and len(lines) >= 3:
        try:
            result['l1_ratio'] = float(lines[2])
        except ValueError:
            pass
    return result


# ---------------------------------------------------------------------------
# Trait name extraction
# ---------------------------------------------------------------------------

def extract_trait_from_sklearn_filename(filename: str, method: str,
                                         output_type: str, maxiter: int = None) -> str:
    """
    Extract trait name from sklearn output filename.
    e.g. ridge_coeffs_trait_overlap_vavg_0.8_numQTL_5000_numChr_1_H2_0.75.csv
         lars_predictions_maxiter1000_trait_overlap_...csv
    """
    base = os.path.basename(filename)
    # Remove extension
    base = re.sub(r'\.(csv|txt)$', '', base)
    # Build prefix to strip
    if maxiter is not None:
        prefix = f"{method}_maxiter{maxiter}_{output_type}_"
    else:
        prefix = f"{method}_{output_type}_"
    if base.startswith(prefix):
        return base[len(prefix):]
    return base


# ---------------------------------------------------------------------------
# sklearn aggregation
# ---------------------------------------------------------------------------

def aggregate_sklearn(args):
    """Collect per-trait sklearn CSV outputs and write coefficient, prediction, and params feathers."""
    method  = args.method
    seed    = args.seed
    maxiter = args.maxiter
    mtag    = maxiter_tag(maxiter)

    intermediates_dir = os.path.join(args.base_dir, f"intermediates_seed_{seed}")
    if not os.path.isdir(intermediates_dir):
        sys.exit(f"[ERROR] Intermediates directory not found: {intermediates_dir}")

    out_dir = output_dir(args.output_base_dir, method, seed, maxiter)
    os.makedirs(out_dir, exist_ok=True)
    label = make_label(method, seed, maxiter)

    print(f"Method         : {method}")
    print(f"Seed           : {seed}")
    print(f"Intermediates  : {intermediates_dir}")
    print(f"Output dir     : {out_dir}")

    # ------------------------------------------------------------------
    # Load reference SNP and IID orders
    # ------------------------------------------------------------------
    snp_order = read_snp_order(seed, args.base_dir, args.prefix)
    iid_order = read_iid_order(seed, SKLEARN_SPLIT, args.base_dir, args.prefix)

    # ------------------------------------------------------------------
    # Discover trait files
    # ------------------------------------------------------------------
    file_prefix    = f"{method}_{mtag}"   # "lars_maxiter1000_" or "ridge_"
    coeffs_pattern = os.path.join(intermediates_dir, f"{file_prefix}coeffs_*.csv")
    preds_pattern  = os.path.join(intermediates_dir, f"{file_prefix}predictions_*.csv")
    alpha_ext      = 'txt'
    alpha_pattern  = os.path.join(intermediates_dir, f"{file_prefix}alphas_*.{alpha_ext}")

    coeffs_files = sorted(glob.glob(coeffs_pattern))
    preds_files  = sorted(glob.glob(preds_pattern))
    alpha_files  = sorted(glob.glob(alpha_pattern))

    print(f"\nFound {len(coeffs_files)} coeffs files")
    print(f"Found {len(preds_files)} predictions files")
    print(f"Found {len(alpha_files)} alphas files")

    if not coeffs_files and not preds_files:
        sys.exit(f"[ERROR] No files found matching patterns in {intermediates_dir}")

    # Build trait list from coeffs files (fall back to preds files)
    source_files = coeffs_files if coeffs_files else preds_files
    trait_names = sorted(set(
        extract_trait_from_sklearn_filename(f, method, 'coeffs' if coeffs_files else 'predictions', maxiter)
        for f in source_files
    ))
    print(f"Traits found   : {len(trait_names)}")

    # Map trait -> file paths
    def file_map(files, output_type):
        return {
            extract_trait_from_sklearn_filename(f, method, output_type, maxiter): f
            for f in files
        }

    coeffs_map = file_map(coeffs_files, 'coeffs')
    preds_map  = file_map(preds_files, 'predictions')
    alpha_map  = file_map(alpha_files, 'alphas')

    # ------------------------------------------------------------------
    # Load IID reference from first predictions file that has IID column
    # (for traits where predictions have no IID column)
    # ------------------------------------------------------------------
    iid_ref = None
    for trait, path in preds_map.items():
        with open(path) as f:
            first_line = f.readline().strip()
        if len(first_line.split(',')) >= 2:
            df_ref = pd.read_csv(path)
            iid_ref = [extract_iid(x) for x in df_ref.iloc[:, 0].astype(str)]
            if len(iid_ref) != len(iid_order):
                print(f"  [WARN] IID reference from {os.path.basename(path)}: "
                      f"{len(iid_ref)} individuals but iid order file has {len(iid_order)}. "
                      f"Using IID order file.", file=sys.stderr)
                iid_ref = [str(i) for i in iid_order]
            else:
                print(f"  IID reference loaded from {os.path.basename(path)}: {len(iid_ref)} individuals")
            break
    if iid_ref is None:
        print("  No predictions file with IID column found; using IID order file.")
        iid_ref = [str(i) for i in iid_order]

    # ------------------------------------------------------------------
    # Assemble coefficient matrix
    # ------------------------------------------------------------------
    print("\nAssembling coefficient matrix ...")
    coeff_arr = np.full((len(snp_order), len(trait_names)), np.nan, dtype=np.float64)

    def load_coeffs(trait):
        if trait not in coeffs_map:
            print(f"  [WARN] No coeffs file for trait: {trait}", file=sys.stderr)
            return trait, None
        series = parse_coeffs_file(coeffs_map[trait], snp_order)
        return trait, series

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(load_coeffs, t): t for t in trait_names}
        for i, future in enumerate(as_completed(futures), 1):
            trait, series = future.result()
            j = trait_names.index(trait)
            if series is not None:
                # Align by SNP order
                aligned = series.reindex(snp_order)
                coeff_arr[:, j] = aligned.values
            if i % 50 == 0 or i == len(trait_names):
                print(f"  Loaded {i}/{len(trait_names)} coeffs files...")

    coeff_df = pd.DataFrame(coeff_arr, columns=trait_names)
    coeff_df.insert(0, 'SNP', snp_order)
    coeff_out = os.path.join(out_dir, f"coeff_matrix_{label}.feather")
    print(f"Writing coeff_matrix_{label}.feather ...")
    feather.write_feather(coeff_df, coeff_out)

    # ------------------------------------------------------------------
    # Assemble prediction matrix
    # ------------------------------------------------------------------
    print("\nAssembling prediction matrix ...")
    pred_arr = np.full((len(iid_ref), len(trait_names)), np.nan, dtype=np.float64)

    def load_preds(trait):
        if trait not in preds_map:
            print(f"  [WARN] No predictions file for trait: {trait}", file=sys.stderr)
            return trait, None
        series = parse_predictions_file(preds_map[trait], iid_ref)
        return trait, series

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(load_preds, t): t for t in trait_names}
        for i, future in enumerate(as_completed(futures), 1):
            trait, series = future.result()
            j = trait_names.index(trait)
            if series is not None:
                aligned = series.reindex([str(x) for x in iid_ref])
                pred_arr[:, j] = aligned.values
            if i % 50 == 0 or i == len(trait_names):
                print(f"  Loaded {i}/{len(trait_names)} predictions files...")

    pred_df = pd.DataFrame(pred_arr, columns=trait_names)
    pred_df.insert(0, 'IID', iid_ref)
    pred_out = os.path.join(out_dir, f"prediction_matrix_{label}.feather")
    print(f"Writing prediction_matrix_{label}.feather ...")
    feather.write_feather(pred_df, pred_out)

    # ------------------------------------------------------------------
    # Assemble params table
    # ------------------------------------------------------------------
    print("\nAssembling params table ...")
    split = SKLEARN_SPLIT
    param_rows = []
    for trait in trait_names:
        row = {
            'method': method,
            'maxiter': maxiter if maxiter is not None else np.nan,
            'split': split,
            'seed': seed,
            'trait': trait,
            'alpha': np.nan,
            'l1_ratio': np.nan,
            'learning_rate': np.nan,
        }
        if trait in alpha_map:
            parsed = parse_alphas_file(alpha_map[trait], method)
            row['alpha'] = parsed.get('alpha', np.nan)
            if method == 'elasticnet':
                row['l1_ratio'] = parsed.get('l1_ratio', np.nan)
        param_rows.append(row)

    params_df = pd.DataFrame(param_rows, columns=['method', 'maxiter', 'split', 'seed',
                                                    'trait', 'alpha', 'l1_ratio', 'learning_rate'])
    params_out = os.path.join(out_dir, f"params_{label}.csv")
    params_df.to_csv(params_out, index=False)
    print(f"Writing params_{label}.csv ...")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_coeff_filled = int(np.sum(~np.isnan(coeff_arr)))
    n_pred_filled  = int(np.sum(~np.isnan(pred_arr)))
    print("\n=== Done ===")
    print(f"  Traits              : {len(trait_names)}")
    print(f"  SNPs                : {len(snp_order)}")
    print(f"  Individuals         : {len(iid_ref)}")
    print(f"  Coeff non-NA        : {n_coeff_filled:,} / {len(snp_order)*len(trait_names):,}")
    print(f"  Pred non-NA         : {n_pred_filled:,} / {len(iid_ref)*len(trait_names):,}")
    print(f"  coeff_matrix        -> {coeff_out}")
    print(f"  prediction_matrix   -> {pred_out}")
    print(f"  params              -> {params_out}")


# ---------------------------------------------------------------------------
# PyTorch-ridge aggregation
# ---------------------------------------------------------------------------

def aggregate_pytorch_ridge(args):
    """Collect per-trait PyTorch ridge CSVs and Optuna JSONs and write coefficient, prediction, and params feathers."""
    seed = args.seed

    if not args.pytorch_dir:
        sys.exit("[ERROR] --pytorch-dir is required for pytorch_ridge.")
    if not args.pytorch_params_dir:
        sys.exit("[ERROR] --pytorch-params-dir is required for pytorch_ridge.")

    pytorch_dir   = args.pytorch_dir
    params_dir    = args.pytorch_params_dir
    out_dir       = output_dir(args.output_base_dir, 'pytorch_ridge', seed)
    os.makedirs(out_dir, exist_ok=True)
    label         = make_label('pytorch_ridge', seed)

    print(f"Method         : pytorch_ridge")
    print(f"Seed           : {seed}")
    print(f"PyTorch dir    : {pytorch_dir}")
    print(f"Params dir     : {params_dir}")
    print(f"Output dir     : {out_dir}")

    # ------------------------------------------------------------------
    # Load param JSON files to build trait list and param table
    # ------------------------------------------------------------------
    json_files = sorted(glob.glob(os.path.join(params_dir, "optuna_best_pheno_*_ridge.json")))
    print(f"\nFound {len(json_files)} param JSON files")
    if not json_files:
        sys.exit(f"[ERROR] No JSON files found in {params_dir}")

    param_rows  = []
    trait_names = []
    trait_to_json = {}

    for jf in json_files:
        with open(jf) as f:
            data = json.load(f)
        trait = data.get('phenotype_name', os.path.basename(jf))
        best_alpha = data.get('best_alpha', np.nan)
        best_lr    = data.get('best_learning_rate', np.nan)
        param_rows.append({
            'method':        'pytorch_ridge',
            'maxiter':       np.nan,
            'split':         PYTORCH_SPLIT,
            'seed':          seed,
            'trait':         trait,
            'alpha':         best_alpha,
            'l1_ratio':      np.nan,
            'learning_rate': best_lr,
        })
        trait_names.append(trait)
        trait_to_json[trait] = jf

    trait_names = sorted(trait_names)
    print(f"Traits found   : {len(trait_names)}")

    # ------------------------------------------------------------------
    # Discover weights and predictions files; map to traits
    # ------------------------------------------------------------------
    weights_files = sorted(glob.glob(
        os.path.join(pytorch_dir, f"ridge_weights_seed_{seed}_*.csv")
    ))
    preds_files = sorted(glob.glob(
        os.path.join(pytorch_dir, f"ridge_predictions_seed_{seed}_*.csv")
    ))
    print(f"Found {len(weights_files)} weights files")
    print(f"Found {len(preds_files)} predictions files")

    def pytorch_trait_from_filename(path: str, seed: int, ftype: str) -> str:
        base = os.path.basename(path)
        base = re.sub(r'\.csv$', '', base)
        prefix = f"ridge_{ftype}_seed_{seed}_"
        if base.startswith(prefix):
            return base[len(prefix):]
        return base

    weights_map = {
        pytorch_trait_from_filename(f, seed, 'weights'): f
        for f in weights_files
    }
    preds_map = {
        pytorch_trait_from_filename(f, seed, 'predictions'): f
        for f in preds_files
    }

    # ------------------------------------------------------------------
    # Load SNP order from first weights file
    # ------------------------------------------------------------------
    if not weights_files:
        sys.exit("[ERROR] No weights files found — cannot determine SNP order.")
    first_weights = pd.read_csv(weights_files[0])
    snp_order = list(first_weights.iloc[:, 0].astype(str))
    print(f"  SNP order from weights file: {len(snp_order)} SNPs")

    # ------------------------------------------------------------------
    # Load IID order from first predictions file
    # ------------------------------------------------------------------
    if not preds_files:
        sys.exit("[ERROR] No predictions files found — cannot determine IID order.")
    first_preds = pd.read_csv(preds_files[0])
    iid_order = list(first_preds.iloc[:, 0].astype(str))
    print(f"  IID order from predictions file: {len(iid_order)} individuals")

    # ------------------------------------------------------------------
    # Assemble weights matrix
    # ------------------------------------------------------------------
    print("\nAssembling weights matrix ...")
    coeff_arr = np.full((len(snp_order), len(trait_names)), np.nan, dtype=np.float64)

    def load_weights(trait):
        if trait not in weights_map:
            print(f"  [WARN] No weights file for trait: {trait}", file=sys.stderr)
            return trait, None
        df = pd.read_csv(weights_map[trait])
        snp_col    = df.columns[0]
        weight_col = df.columns[1]
        # Verify length
        if len(df) != len(snp_order):
            print(f"  [WARN] {trait}: weights file has {len(df)} SNPs, "
                  f"expected {len(snp_order)} — aligning by SNP ID.", file=sys.stderr)
        series = pd.to_numeric(df[weight_col], errors='coerce')
        series.index = df[snp_col].astype(str).values
        return trait, series

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(load_weights, t): t for t in trait_names}
        for i, future in enumerate(as_completed(futures), 1):
            trait, series = future.result()
            j = trait_names.index(trait)
            if series is not None:
                aligned = series.reindex(snp_order)
                coeff_arr[:, j] = aligned.values
            if i % 50 == 0 or i == len(trait_names):
                print(f"  Loaded {i}/{len(trait_names)} weights files...")

    coeff_df = pd.DataFrame(coeff_arr, columns=trait_names)
    coeff_df.insert(0, 'SNP', snp_order)
    coeff_out = os.path.join(out_dir, f"coeff_matrix_{label}.feather")
    print(f"Writing coeff_matrix_{label}.feather ...")
    feather.write_feather(coeff_df, coeff_out)

    # ------------------------------------------------------------------
    # Assemble predictions matrix
    # ------------------------------------------------------------------
    print("\nAssembling predictions matrix ...")
    pred_arr = np.full((len(iid_order), len(trait_names)), np.nan, dtype=np.float64)

    def load_preds(trait):
        if trait not in preds_map:
            print(f"  [WARN] No predictions file for trait: {trait}", file=sys.stderr)
            return trait, None
        df = pd.read_csv(preds_map[trait])
        iid_col  = df.columns[0]
        pred_col = df.columns[2]  # predicted_phenotype is 3rd column
        # Verify individual count
        if len(df) != len(iid_order):
            print(f"  [WARN] {trait}: predictions file has {len(df)} individuals, "
                  f"expected {len(iid_order)} — aligning by IID.", file=sys.stderr)
        series = pd.to_numeric(df[pred_col], errors='coerce')
        series.index = df[iid_col].astype(str).values
        return trait, series

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(load_preds, t): t for t in trait_names}
        for i, future in enumerate(as_completed(futures), 1):
            trait, series = future.result()
            j = trait_names.index(trait)
            if series is not None:
                aligned = series.reindex(iid_order)
                pred_arr[:, j] = aligned.values
            if i % 50 == 0 or i == len(trait_names):
                print(f"  Loaded {i}/{len(trait_names)} predictions files...")

    pred_df = pd.DataFrame(pred_arr, columns=trait_names)
    pred_df.insert(0, 'IID', iid_order)
    pred_out = os.path.join(out_dir, f"prediction_matrix_{label}.feather")
    print(f"Writing prediction_matrix_{label}.feather ...")
    feather.write_feather(pred_df, pred_out)

    # ------------------------------------------------------------------
    # Write params table
    # ------------------------------------------------------------------
    params_df  = pd.DataFrame(param_rows, columns=['method', 'maxiter', 'split', 'seed',
                                                    'trait', 'alpha', 'l1_ratio', 'learning_rate'])
    params_out = os.path.join(out_dir, f"params_{label}.csv")
    params_df.to_csv(params_out, index=False)
    print(f"Writing params_{label}.csv ...")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_coeff_filled = int(np.sum(~np.isnan(coeff_arr)))
    n_pred_filled  = int(np.sum(~np.isnan(pred_arr)))
    print("\n=== Done ===")
    print(f"  Traits              : {len(trait_names)}")
    print(f"  SNPs                : {len(snp_order)}")
    print(f"  Individuals         : {len(iid_order)}")
    print(f"  Coeff non-NA        : {n_coeff_filled:,} / {len(snp_order)*len(trait_names):,}")
    print(f"  Pred non-NA         : {n_pred_filled:,} / {len(iid_order)*len(trait_names):,}")
    print(f"  coeff_matrix        -> {coeff_out}")
    print(f"  prediction_matrix   -> {pred_out}")
    print(f"  params              -> {params_out}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate regression outputs into coefficient and prediction matrices.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument('--method', required=True, choices=sorted(ALL_METHODS),
                        help='Regression method to aggregate')
    parser.add_argument('--seed', required=True, type=int,
                        help='Random seed')
    parser.add_argument('--base-dir', default='.',
                        help='Parent directory containing intermediates_seed_* (default: .)')
    parser.add_argument('--output-base-dir', default=None,
                        help='Parent directory for aggregate_outputs_* '
                             '(default: same as --base-dir)')
    parser.add_argument('--prefix', default='yeast_simulated_data',
                        help='Filename prefix for train/test ID and SNP files, e.g. '
                             '"yeast_simulated_data" or "ukbb_simulated_traits" '
                             '(default: yeast_simulated_data)')

    parser.add_argument('--maxiter', type=int, default=None,
                        help='Max iterations suffix for lars (e.g. --maxiter 1000)')
    parser.add_argument('--pytorch-dir', default=None,
                        help='Directory containing pytorch_ridge weights/predictions files')
    parser.add_argument('--pytorch-params-dir', default=None,
                        help='Directory containing optuna JSON param files')
    parser.add_argument('--threads', type=int, default=4,
                        help='Parallel threads for file loading (default: 4)')

    args = parser.parse_args()

    if args.output_base_dir is None:
        args.output_base_dir = args.base_dir

    # Validate lars maxiter
    if args.method == 'lars' and args.maxiter is None:
        parser.error("--maxiter is required for lars.")

    # Validate pytorch args
    if args.method == 'pytorch_ridge':
        if not args.pytorch_dir:
            parser.error("--pytorch-dir is required for pytorch_ridge.")
        if not args.pytorch_params_dir:
            parser.error("--pytorch-params-dir is required for pytorch_ridge.")

    if args.method in SKLEARN_METHODS:
        aggregate_sklearn(args)
    elif args.method in PYTORCH_METHODS:
        aggregate_pytorch_ridge(args)


if __name__ == '__main__':
    main()
