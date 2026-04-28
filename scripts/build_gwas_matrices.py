#!/usr/bin/env python3
"""
build_gwas_matrices.py

Two modes of operation. Both require --data-split and --seed, which
determine all input and output directory/file names:

  input  directory : plink_outputs_{data_split}_seed_{which_seed}/
  output directory : aggregate_outputs_{data_split}_seed_{which_seed}/

─── MODE 1: build ──────────────────────────────────────────────────────────
Reads all plink2 GWAS result files and builds two parallel matrices
(SNPs x traits):
  1. beta_matrix_{data_split}_seed_{which_seed}.feather
  2. pvalue_matrix_{data_split}_seed_{which_seed}.feather

Optionally also writes a p-value-filtered beta matrix:
  3. beta_matrix_filtered_p{p_threshold}_{data_split}_seed_{which_seed}.feather

Sign correction rule:
  - If ALT == A1  -> keep BETA as-is
  - If REF == A1  -> multiply BETA by -1

Usage:
  python build_gwas_matrices.py build --data-split {test|train} --seed N [OPTIONS]

Options:
  --data-split  STR      "test" or "train" (required)
  --seed  INT      Seed number (required)
  --base-dir    DIR      Parent directory containing plink_outputs_* folders
                         (default: current directory)
  --pattern     GLOB     Glob pattern to match files within the input dir
                         (default: "*.glm.linear")
  --p-threshold FLOAT    Optional p-value threshold for a filtered BETA matrix
  --threads     INT      Number of CPU threads for parallel loading (default: 4)

─── MODE 2: predict ────────────────────────────────────────────────────────
Loads pre-built beta and p-value matrices, filters betas by a p-value
threshold, aligns SNPs with a genotype matrix, then computes:

  polygenic_scores = genotype_matrix  @  beta_matrix_filtered
                     (individuals x SNPs)   (SNPs x traits)
  -> result: individuals x traits

Output file:
  polygenic_scores_p{p_threshold}_{data_split}_seed_{which_seed}.feather

Usage:
  python build_gwas_matrices.py predict --data-split {test|train} --seed N [OPTIONS]
  python scripts/build_gwas_matrices.py predict --data-split train --seed 1510 \
    --p-threshold 1.2e-6 --geno-matrix input_data/yeast_genotypes_binarized_centered.feather

Options:
  --data-split    STR    "test" or "train" (required)
  --seed    INT    Seed number (required)
  --base-dir      DIR    Parent directory containing aggregate_outputs_* folders
                         (default: current directory)
  --geno-matrix   FILE   Path to genotype matrix (.feather, .tsv, or .csv)
                         Rows = individuals, columns = SNPs (required)
  --p-threshold   FLOAT  P-value threshold for filtering betas (required)
  --id-col        STR    Individual ID column in genotype matrix (default: first col)
  --missing-geno  FLOAT  Fill value for SNPs absent from genotype matrix (default: 0)
"""

import argparse
import glob
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import pyarrow.feather as feather


# ---------------------------------------------------------------------------
# Naming helpers
# ---------------------------------------------------------------------------

def make_suffix(data_split: str, which_seed: int, num_covariates: int = None) -> str:
    """Return the shared filename/directory suffix.

    Standard : 'train_seed_42'
    Biobank  : 'cov6_train_seed_42'
    """
    base = f"{data_split}_seed_{which_seed}"
    if num_covariates is not None:
        return f"cov{num_covariates}_{base}"
    return base


def input_dir_name(base_dir: str, suffix: str, num_covariates: int = None) -> str:
    base = os.path.join(base_dir, f"plink_outputs_{suffix}")
    if num_covariates is not None:
        return os.path.join(base, f"cov{num_covariates}")
    return base


def output_dir_name(base_dir: str, suffix: str) -> str:
    return os.path.join(base_dir, f"aggregate_outputs_{suffix}")


def fmt_threshold(p: float) -> str:
    """Format a p-value threshold into a filename-safe string, e.g. 5e-08."""
    return f"{p:.2e}".replace("+", "")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def extract_trait_name(filepath: str, prefix: str = 'plink_nofiltering') -> str:
    """Derive a short trait label from the filename.

    Expected filename format: {prefix}.{trait_name}.glm.linear

    Standard run : prefix = 'plink_nofiltering'
    Biobank run  : prefix = 'cov{num_covariates}nofiltering'
    """
    base = os.path.basename(filepath)
    base = re.sub(rf'^{re.escape(prefix)}\.', '', base)
    base = re.sub(r'\.glm\.linear$', '', base)
    return base


def parse_p_to_log10(p_series: pd.Series) -> pd.Series:
    """Convert a string Series of p-values to -log10(P) without float64 underflow.

    Parses scientific notation manually so that a value like 1.1925e-325
    becomes -log10(1.1925) - (-325) ~ 324.924, rather than underflowing to 0.0.
    Non-numeric / NA values are returned as NaN.
    """
    def _convert(s):
        if pd.isna(s):
            return np.nan
        s = str(s).strip().lower()
        if s in ('na', 'nan', '.', ''):
            return np.nan
        try:
            if 'e' in s:
                coeff_str, exp_str = s.split('e')
                coeff = float(coeff_str)
                exp   = int(exp_str)
                if coeff <= 0:
                    return np.nan
                return -np.log10(coeff) - exp
            else:
                val = float(s)
                if val <= 0:
                    return np.nan
                return -np.log10(val)
        except (ValueError, TypeError):
            return np.nan

    return p_series.map(_convert)


def load_file(filepath: str, prefix: str = 'plink_nofiltering'):
    """
    Load a single .glm.linear file.
    Returns (trait_name, DataFrame indexed by SNP ID with BETA_* and P_* columns).
    """
    trait = extract_trait_name(filepath, prefix=prefix)
    try:
        df = pd.read_csv(filepath, sep='\t', comment=None, low_memory=False, dtype={'P': str})
        df.columns = [c.lstrip('#') for c in df.columns]

        required = {'CHROM', 'POS', 'ID', 'REF', 'ALT', 'A1', 'BETA', 'P'}
        missing = required - set(df.columns)
        if missing:
            print(f"[WARN] {trait}: missing columns {missing} — skipping", file=sys.stderr)
            return trait, None

        # Include ERRCODE so we can identify CONST_OMITTED_ALLELE rows
        errcode_present = 'ERRCODE' in df.columns
        cols = ['ID', 'REF', 'ALT', 'A1', 'TEST', 'BETA', 'P']
        if errcode_present:
            cols.append('ERRCODE')
        df = df[cols].copy()
        df = df[df['TEST'] == 'ADD'].copy()
        df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')
        # Parse P as strings -> -log10(P) to avoid float64 underflow on tiny p-values
        df['P'] = parse_p_to_log10(df['P'])

        # Flag CONST_OMITTED_ALLELE rows: BETA=0, -log10(P)=0 (i.e. P=1)
        if errcode_present:
            const_mask = df['ERRCODE'].str.contains('CONST_OMITTED_ALLELE', na=False)
            n_const = const_mask.sum()
            if n_const > 0:
                print(f"[INFO] {trait}: {n_const} SNP(s) with ERRCODE=CONST_OMITTED_ALLELE "
                      f"— setting BETA=0, P=1.", file=sys.stderr)
                df.loc[const_mask, 'BETA'] = 0.0
                df.loc[const_mask, 'P']    = 0.0  # -log10(1) = 0

        # Drop any remaining rows where BETA or P are still NA
        n_before = len(df)
        dropped_mask = df['BETA'].isna() | df['P'].isna()
        if dropped_mask.any():
            dropped_df = df[dropped_mask][['ID', 'BETA', 'P', 'ERRCODE'] if errcode_present
                                          else ['ID', 'BETA', 'P']].head(5)
            n_dropped = dropped_mask.sum()
            print(f"[WARN] {trait}: {n_dropped} SNP(s) dropped due to unparseable "
                  f"BETA or P values. First up to 5:", file=sys.stderr)
            print(dropped_df.to_string(index=False), file=sys.stderr)
        df.dropna(subset=['BETA', 'P'], inplace=True)

        # Sign correction: flip BETA when REF is the effect allele
        ref_is_a1 = df['REF'].str.upper() == df['A1'].str.upper()
        df['BETA'] = np.where(ref_is_a1, -df['BETA'], df['BETA'])

        df.drop_duplicates(subset='ID', keep='first', inplace=True)
        result = df[['ID', 'BETA', 'P']].rename(
            columns={'BETA': f'BETA_{trait}', 'P': f'P_{trait}'}
        ).set_index('ID')

        return trait, result

    except Exception as exc:
        print(f"[ERROR] {trait}: {exc}", file=sys.stderr)
        return trait, None


def load_matrix(path: str, id_col: str = None) -> pd.DataFrame:
    """Load a matrix from .feather, .tsv, or .csv. Returns DataFrame with ID as index."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.feather':
        df = feather.read_feather(path)
    elif ext in ('.tsv', '.txt'):
        df = pd.read_csv(path, sep='\t', low_memory=False)
    elif ext == '.csv':
        df = pd.read_csv(path, low_memory=False)
    else:
        try:
            df = feather.read_feather(path)
        except Exception:
            df = pd.read_csv(path, sep='\t', low_memory=False)

    col = id_col if id_col else df.columns[0]
    if col not in df.columns:
        sys.exit(f"[ERROR] ID column '{col}' not found in {path}")
    return df.set_index(col)


# ---------------------------------------------------------------------------
# Build mode
# ---------------------------------------------------------------------------

def run_build(args):
    suffix     = make_suffix(args.data_split, args.seed, args.num_covariates if args.biobank else None)
    input_dir  = input_dir_name(args.base_dir, make_suffix(args.data_split, args.seed), args.num_covariates if args.biobank else None)
    output_dir = output_dir_name(args.base_dir, suffix)

    if not os.path.isdir(input_dir):
        sys.exit(f"[ERROR] Input directory not found: {input_dir}")
    os.makedirs(output_dir, exist_ok=True)

    print(f"data_split : {args.data_split}")
    print(f"which_seed : {args.seed}")
    print(f"Input  dir : {input_dir}")
    print(f"Output dir : {output_dir}")

    # Discover files
    pattern = os.path.join(input_dir, args.pattern)
    files   = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"[ERROR] No files matched: {pattern}")
    print(f"\nFound {len(files)} files. Loading with {args.threads} thread(s)...")

    # Determine filename prefix based on biobank mode
    if args.biobank:
        file_prefix = f"cov{args.num_covariates}_nofiltering"
    else:
        file_prefix = "plink_nofiltering"
    print(f"File prefix: {file_prefix}")

    # Load files in parallel
    trait_dfs = {}
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {pool.submit(load_file, f, file_prefix): f for f in files}
        for i, future in enumerate(as_completed(futures), 1):
            trait, df = future.result()
            if df is not None:
                trait_dfs[trait] = df
            if i % 50 == 0 or i == len(files):
                print(f"  Loaded {i}/{len(files)} files...")

    if not trait_dfs:
        sys.exit("[ERROR] No files loaded successfully.")
    print(f"\nSuccessfully loaded {len(trait_dfs)} traits.")

    # Build unified SNP index
    print("Building SNP index...")
    all_snps = sorted(
        set().union(*[set(df.index) for df in trait_dfs.values()]),
        key=lambda s: int(re.sub(r'\D', '', s)) if re.search(r'\d', s) else s
    )
    print(f"Total unique SNPs: {len(all_snps)}")
    snp_list_out = os.path.join(output_dir, f"snp_list_{suffix}.txt")
    with open(snp_list_out, 'w') as fh:
        fh.write(os.linesep.join(all_snps) + os.linesep)
    print(f"SNP list written to: {snp_list_out}")

    # Assemble matrices
    print("Assembling matrices (this may take a moment)...")
    trait_names = sorted(trait_dfs.keys())
    n_snp   = len(all_snps)
    n_trait = len(trait_names)
    snp_idx = {s: i for i, s in enumerate(all_snps)}

    beta_arr = np.full((n_snp, n_trait), np.nan, dtype=np.float64)
    pval_arr = np.full((n_snp, n_trait), np.nan, dtype=np.float64)

    for j, trait in enumerate(trait_names):
        df       = trait_dfs[trait]
        beta_col = f'BETA_{trait}'
        p_col    = f'P_{trait}'
        # Align df to the master SNP index order, then assign in one vectorised step
        aligned  = df[[beta_col, p_col]].reindex(all_snps)
        beta_arr[:, j] = aligned[beta_col].values
        pval_arr[:, j] = aligned[p_col].values

    # Write outputs
    def write_feather_with_snp(arr, path):
        out_df = pd.DataFrame(arr, columns=trait_names)
        out_df.insert(0, 'SNP', all_snps)
        feather.write_feather(out_df, path)

    beta_name = f"beta_matrix_{suffix}.feather"
    pval_name = f"pvalue_matrix_{suffix}.feather"
    beta_out  = os.path.join(output_dir, beta_name)
    pval_out  = os.path.join(output_dir, pval_name)

    print(f"Writing {beta_name} ...")
    write_feather_with_snp(beta_arr, beta_out)
    print(f"Writing {pval_name} ...")
    # pval_arr already contains -log10(P) as computed in load_file
    write_feather_with_snp(pval_arr, pval_out)

    # Optional filtered BETA matrix
    filt_out = None
    if args.p_threshold is not None:
        p_str     = fmt_threshold(args.p_threshold)
        filt_name = f"beta_matrix_filtered_{suffix}_p{p_str}.feather"
        filt_out  = os.path.join(output_dir, filt_name)
        print(f"Writing {filt_name} (P <= {args.p_threshold}) ...")
        beta_filt = beta_arr.copy()
        beta_filt[pval_arr < -np.log10(args.p_threshold)] = np.nan
        write_feather_with_snp(beta_filt, filt_out)

    # Summary
    total_entries = n_snp * n_trait
    filled_beta   = int(np.sum(~np.isnan(beta_arr)))
    filled_pct    = 100 * filled_beta / total_entries if total_entries else 0

    print("\n=== Done ===")
    print(f"  Matrix dimensions : {n_snp} SNPs x {n_trait} traits")
    print(f"  Non-NA entries    : {filled_beta:,} / {total_entries:,} ({filled_pct:.1f}%)")
    print(f"  {beta_name} -> {beta_out}")
    print(f"  {pval_name} -> {pval_out}")
    if filt_out:
        sig_entries = int(np.sum(~np.isnan(beta_filt)))
        print(f"  {filt_name} ({sig_entries:,} entries) -> {filt_out}")


# ---------------------------------------------------------------------------
# Predict mode
# ---------------------------------------------------------------------------

def run_predict(args):
    suffix     = make_suffix(args.data_split, args.seed, args.num_covariates if args.biobank else None)
    output_dir = output_dir_name(args.base_dir, suffix)
    os.makedirs(output_dir, exist_ok=True)

    p_str = fmt_threshold(args.p_threshold)

    print(f"data_split : {args.data_split}")
    print(f"which_seed : {args.seed}")
    print(f"Output dir : {output_dir}")

    # Resolve beta and p-value matrix paths (look in aggregate_outputs dir by default)
    beta_path = args.beta_matrix or os.path.join(
        output_dir, f"beta_matrix_{suffix}.feather"
    )
    pval_path = args.pvalue_matrix or os.path.join(
        output_dir, f"pvalue_matrix_{suffix}.feather"
    )

    for path, label in [(beta_path, 'beta'), (pval_path, 'p-value')]:
        if not os.path.isfile(path):
            sys.exit(f"[ERROR] Could not find {label} matrix: {path}\n"
                     f"        Run 'build' first, or supply the path explicitly.")

    print(f"\nLoading beta matrix    : {beta_path}")
    beta_df = load_matrix(beta_path, id_col='SNP')

    print(f"Loading p-value matrix : {pval_path}")
    log10_pval_df = load_matrix(pval_path, id_col='SNP')
    # Matrix stores -log10(P); convert threshold to -log10 scale for comparison
    pval_df = log10_pval_df  # comparison done in log10 space below

    if not beta_df.index.equals(pval_df.index):
        sys.exit("[ERROR] Beta and p-value matrices have different SNP indices.")
    if not beta_df.columns.equals(pval_df.columns):
        sys.exit("[ERROR] Beta and p-value matrices have different trait columns.")

    # Filter betas by p-value threshold
    # Set non-significant and NaN p-value entries to 0 so they contribute nothing
    print(f"\nFiltering betas at P <= {args.p_threshold} ...")
    beta_filt = beta_df.fillna(0.0).copy()
    log10_threshold = -np.log10(args.p_threshold)
    fail_threshold = (pval_df < log10_threshold) | pval_df.isna()
    beta_filt[fail_threshold] = 0.0

    n_sig = (beta_filt != 0.0).sum()
    print(f"  Significant SNP-trait pairs — "
          f"min={n_sig.min()}, median={int(n_sig.median())}, "
          f"max={n_sig.max()}, total={n_sig.sum():,}")

    # Load genotype matrix (individuals x SNPs)
    print(f"\nLoading genotype matrix: {args.geno_matrix}")
    geno_df = load_matrix(args.geno_matrix, id_col=args.id_col)
    n_indiv, n_geno_snps = geno_df.shape
    print(f"  Genotype matrix: {n_indiv} individuals x {n_geno_snps} SNPs")

    # Strip trailing allele suffix from genotype SNP IDs (e.g. "rs123_A" -> "rs123")
    geno_df.columns = geno_df.columns.str.replace(r'_[ACGT]+$', '', regex=True)

    # Align SNP columns of genotype matrix to SNP rows of beta matrix
    beta_snps = beta_filt.index
    geno_snps = geno_df.columns

    snps_in_both   = beta_snps.intersection(geno_snps)
    snps_only_beta = beta_snps.difference(geno_snps)
    snps_only_geno = geno_snps.difference(beta_snps)

    print(f"\nSNP alignment:")
    print(f"  SNPs in beta matrix              : {len(beta_snps):,}")
    print(f"  SNPs in genotype matrix          : {len(geno_snps):,}")
    print(f"  SNPs in common                   : {len(snps_in_both):,}")
    if len(snps_only_beta) > 0:
        print(f"  SNPs only in beta (filled {args.missing_geno}): {len(snps_only_beta):,}")
    if len(snps_only_geno) > 0:
        print(f"  SNPs only in genotype (dropped)  : {len(snps_only_geno):,}")
    print(f"  First 5 SNP IDs in beta matrix   : {list(beta_snps[:5])}")
    print(f"  First 5 SNP IDs in geno matrix   : {list(geno_snps[:5])}")
    if len(snps_in_both) == 0:
        sys.exit("[ERROR] No SNPs in common between beta and genotype matrices. "
                 "Check that SNP ID formats match (see samples above).")

    # Reindex genotype columns to exactly match beta SNP order
    geno_aligned = geno_df.reindex(columns=beta_snps, fill_value=args.missing_geno)

    # All non-significant betas already set to 0; fillna catches any remaining NaNs
    beta_for_mult = beta_filt.fillna(0).values   # (n_snps, n_traits)
    geno_values   = geno_aligned.values           # (n_indiv, n_snps)

    # Matrix multiplication
    print("\nComputing polygenic scores ...")
    scores = geno_values @ beta_for_mult          # (n_indiv, n_traits)

    # Save result
    scores_name = f"polygenic_scores_{suffix}_p{p_str}.feather"
    scores_out  = os.path.join(output_dir, scores_name)
    print(f"Writing {scores_name} ...")
    scores_df = pd.DataFrame(scores,
                             index=geno_df.index,
                             columns=beta_filt.columns)
    scores_df.index.name = 'individual'
    scores_df = scores_df.reset_index()
    feather.write_feather(scores_df, scores_out)

    print("\n=== Done ===")
    print(f"  Polygenic score matrix : {n_indiv} individuals x {len(beta_filt.columns)} traits")
    print(f"  {scores_name} -> {scores_out}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def add_shared_args(p):
    """Add common arguments to a subparser."""
    p.add_argument('--data-split', required=True, choices=['test', 'train'],
                   help='"test" or "train"')
    p.add_argument('--seed', required=True, type=int,
                   help='Numeric seed identifier')
    p.add_argument('--base-dir',   default='.',
                   help='Parent directory containing plink_outputs_* / '
                        'aggregate_outputs_* folders (default: .)')
    p.add_argument('--biobank',    action='store_true', default=False,
                   help='Enable biobank mode. Requires --num-covariates.')
    p.add_argument('--num-covariates', type=int, default=None,
                   help='Number of covariates used in biobank GWAS files '
                        '(required when --biobank is set). Sets the file prefix '
                        'to "cov{num_covariates}nofiltering".')


def main():
    parser = argparse.ArgumentParser(
        description="Build GWAS beta/p-value matrices or compute polygenic scores.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # ── build ──────────────────────────────────────────────────────────
    build_p = subparsers.add_parser(
        'build',
        help='Build SNP x trait BETA and P-value matrices from plink2 GWAS files.'
    )
    add_shared_args(build_p)
    build_p.add_argument('--pattern',     default='*.glm.linear',
                         help='Glob pattern to match files (default: *.glm.linear)')
    build_p.add_argument('--p-threshold', type=float, default=None,
                         help='Optional p-value threshold for a filtered BETA matrix '
                              '(e.g. 5e-8)')
    build_p.add_argument('--threads',     type=int,   default=4,
                         help='Parallel threads for file loading (default: 4)')

    # ── predict ────────────────────────────────────────────────────────
    pred_p = subparsers.add_parser(
        'predict',
        help='Compute polygenic scores: genotype_matrix @ filtered_beta_matrix.'
    )
    add_shared_args(pred_p)
    pred_p.add_argument('--geno-matrix',   required=True,
                        help='Genotype matrix file (.feather/.tsv/.csv). '
                             'Rows=individuals, columns=SNPs.')
    pred_p.add_argument('--p-threshold',   type=float, required=True,
                        help='P-value threshold for filtering betas (e.g. 5e-8)')
    pred_p.add_argument('--beta-matrix',   default=None,
                        help='Path to beta matrix feather (default: auto-resolved '
                             'from aggregate_outputs dir)')
    pred_p.add_argument('--pvalue-matrix', default=None,
                        help='Path to p-value matrix feather (default: auto-resolved '
                             'from aggregate_outputs dir)')
    pred_p.add_argument('--id-col',        default=None,
                        help='Individual ID column in genotype matrix (default: first col)')
    pred_p.add_argument('--missing-geno',  type=float, default=0.0,
                        help='Fill value for SNPs missing from genotype matrix (default: 0)')

    args = parser.parse_args()

    # Validate biobank arguments
    if args.biobank and args.num_covariates is None:
        parser.error("--num-covariates is required when --biobank is set.")
    if not args.biobank and args.num_covariates is not None:
        parser.error("--num-covariates was supplied but --biobank was not set.")

    if args.command == 'build':
        run_build(args)
    elif args.command == 'predict':
        run_predict(args)


if __name__ == '__main__':
    main()
