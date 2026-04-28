#!/usr/bin/env python3
"""
regress.py — Cross-validated sklearn regression for genomic prediction.

Fits a penalized linear regression model (ridge, lasso, lars, or elasticnet)
to centered genotype data and writes per-phenotype coefficients, best
hyperparameters, and predictions to an intermediates directory.

Expected data layout (Feather format, produced by save_input_data.py):
  <test-train-dir>/
    {prefix}_seed_{seed}_train_genotypes_centered.feather
    {prefix}_seed_{seed}_test_genotypes_centered.feather
    {prefix}_seed_{seed}_test_phenotypes_normalized.feather

Run from the repository root.  Example usage:

  # Ridge on all phenotypes, seed 1510
  python scripts/regress.py --method ridge --seed 1510

  # Lasso on phenotypes matching two substrings
  python scripts/regress.py --method lasso --seed 1510 \\
      --pheno-filter numQTL_100_ _numChr_16_

  # ElasticNet with a custom alpha grid
  python scripts/regress.py --method elasticnet --seed 1510 \\
      --alphas logspace:-5:2:15 --l1-ratios 0.5,0.9,0.99

  # LARS on an explicit phenotype list
  python scripts/regress.py --method lars --seed 1510 \\
      --phenos trait_overlap_vavg_1_numQTL_100_numChr_16_H2_0.25 \\
               trait_overlap_vavg_0.5_numQTL_100_numChr_16_H2_0.5

  # UK Biobank mode: upload each output file to DNAnexus after saving
  python scripts/regress.py --method elasticnet --seed 1105 \\
      --prefix ukbb_simulated_traits --biobank \\
      --dx-project exploration --dx-folder /intermediates_seed_1105/ \\
      --phenos <phenotype_name>
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.feather as feather
from sklearn.linear_model import (
    ElasticNet,
    ElasticNetCV,
    Lasso,
    LassoCV,
    LarsCV,
    Ridge,
    RidgeCV,
)


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────


def load_feather(path: str, id_col: str = None) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Load a Feather file, drop the ID column, return (data, col_names, ids)."""
    df = feather.read_feather(path)
    if id_col is None:
        id_col = df.columns[0]
    ids = df[id_col].values
    col_names = [c for c in df.columns if c != id_col]
    data = df[col_names].to_numpy()
    print(
        f"  → Loaded {os.path.basename(path)}: "
        f"{len(ids)} individuals, {data.shape[1]} features/phenotypes"
    )
    return data, col_names, ids


# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    parser.add_argument(
        "--method",
        required=True,
        choices=["ridge", "lasso", "lars", "elasticnet"],
        help="Regression method.",
    )

    # Data paths / identifiers
    parser.add_argument(
        "--seed",
        default="1510",
        help="Seed identifier embedded in input file names (default: 1510).",
    )
    parser.add_argument(
        "--prefix",
        default="yeast_simulated_data",
        help="Filename prefix for Feather data files (default: yeast_simulated_data).",
    )
    parser.add_argument(
        "--test-train-dir",
        default=None,
        help="Directory containing train/test Feather files. "
        "Default: test_train_seed_<seed>/",
    )
    parser.add_argument(
        "--intermediates-dir",
        default=None,
        help="Output directory for coefficients, alphas, and predictions. "
        "Default: intermediates_seed_<seed>/",
    )

    # CV / solver settings
    parser.add_argument(
        "--n-folds",
        type=int,
        default=5,
        help="Cross-validation folds (default: 5).",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=None,
        help="Maximum solver iterations. Defaults: lars=1000, lasso/elasticnet=10000. "
        "Ignored for ridge.",
    )
    parser.add_argument(
        "--alphas",
        default=None,
        help="Alpha regularisation grid. Accepts a comma-separated list of floats "
        "or logspace notation 'logspace:start:stop:n'. "
        "Defaults: ridge=logspace:-3:7:16, lasso/elasticnet=logspace:-7:2:10. "
        "Ignored for lars (path determined automatically).",
    )
    parser.add_argument(
        "--l1-ratios",
        default="0.8,0.9,0.95,0.99",
        help="Comma-separated L1-ratio grid for elasticnet (default: 0.8,0.9,0.95,0.99). "
        "1.0 = pure lasso, 0.0 = pure ridge.",
    )

    # Phenotype selection
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--pheno-filter",
        nargs="+",
        metavar="SUBSTR",
        help="One or more substrings; only phenotypes whose names contain ALL "
        "substrings are included. Example: --pheno-filter numQTL_100_ _numChr_16_",
    )
    group.add_argument(
        "--phenos",
        nargs="+",
        metavar="PHENOTYPE",
        help="Explicit phenotype names to run. Mutually exclusive with --pheno-filter.",
    )

    # DNAnexus / UK Biobank
    parser.add_argument(
        "--biobank",
        action="store_true",
        help="Enable UK Biobank mode: upload each output file to DNAnexus after saving. "
        "Requires dxpy and the dx CLI to be installed and authenticated.",
    )
    parser.add_argument(
        "--dx-project",
        default=None,
        help="DNAnexus project name to upload files to (required when --biobank is set).",
    )
    parser.add_argument(
        "--dx-folder",
        default=None,
        help="Destination folder path within the DNAnexus project "
        "(default: /intermediates_seed_<seed>/).",
    )

    return parser.parse_args()


def _parse_alphas(alphas_str: str | None, method: str) -> np.ndarray:
    if alphas_str is None:
        return np.logspace(-3, 7, 16) if method == "ridge" else np.logspace(-7, 2, 10)
    if alphas_str.startswith("logspace:"):
        parts = alphas_str.split(":")
        _, start, stop, n = parts
        return np.logspace(float(start), float(stop), int(n))
    return np.array([float(x) for x in alphas_str.split(",")])


# ─────────────────────────────────────────────────────────────────────────────
# DNAnexus helpers (only used when --biobank is set)
# ─────────────────────────────────────────────────────────────────────────────


def setup_dnanexus(dx_project: str) -> None:
    """Clear the workspace environment variable and select the DNAnexus project."""
    import dxpy  # noqa: F401 — imported here so the rest of the script works without dxpy

    print("Configuring DNAnexus environment...")
    os.environ.pop("DX_WORKSPACE_ID", None)
    result = subprocess.run(["dx", "select", dx_project], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  dx select failed: {result.stderr.strip()}")
        sys.exit(1)
    print(f"  Selected project: {dx_project}")


def upload_to_dnanexus(local_path: str, dx_project: str, dx_folder: str) -> None:
    """Upload a local file to the configured DNAnexus project and folder."""
    import dxpy

    filename = os.path.basename(local_path)
    print(f"  Uploading {filename} to {dx_project}:{dx_folder} ...")
    project = dxpy.find_one_project(name=dx_project, return_handler=True)
    dxpy.upload_local_file(
        filename=local_path,
        project=project.get_id(),
        folder=dx_folder,
        wait_on_close=True,
    )
    print(f"  Upload complete → {dx_project}:{dx_folder}/{filename}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    args = parse_args()

    method = args.method
    seed = args.seed
    prefix = args.prefix
    test_train_dir = args.test_train_dir or f"test_train_seed_{seed}/"
    intermediates_dir = args.intermediates_dir or f"intermediates_seed_{seed}/"
    n_folds = args.n_folds
    max_iter = args.max_iter or (1000 if method == "lars" else 10000)
    alphas = _parse_alphas(args.alphas, method)
    l1_ratios = [float(x) for x in args.l1_ratios.split(",")]

    # ── DNAnexus setup ─────────────────────────────────────────────────────────
    biobank = args.biobank
    if biobank:
        if not args.dx_project:
            print("Error: --dx-project is required when --biobank is set.")
            sys.exit(1)
        dx_project = args.dx_project
        dx_folder = args.dx_folder or f"/intermediates_seed_{seed}/"
        setup_dnanexus(dx_project)

    # ── Load data ──────────────────────────────────────────────────────────────
    # NOTE: "train" and "test" file labels are intentionally swapped here to
    # preserve the convention used in the original per-method scripts: the file
    # called *train_genotypes* is used as the held-out test set, and vice versa.
    print("\nLoading data...")
    training_geno, snp_ids, test_iids = load_feather(
        f"{test_train_dir}{prefix}_seed_{seed}_train_genotypes_centered.feather"
    )

    testing_geno, _, _ = load_feather(
        f"{test_train_dir}{prefix}_seed_{seed}_test_genotypes_centered.feather"
    )
    testing_pheno, pheno_names_list, _ = load_feather(
        f"{test_train_dir}{prefix}_seed_{seed}_test_phenotypes_normalized.feather"
    )

    print(f"\n  → {len(pheno_names_list)} phenotypes available, {len(snp_ids)} SNPs")

    # ── Select phenotypes ──────────────────────────────────────────────────────
    if args.phenos:
        for p in args.phenos:
            if p not in pheno_names_list:
                print(f"Error: '{p}' not found in phenotype columns.")
                sys.exit(1)
        some_phenos = args.phenos
        some_phenos_indices = [pheno_names_list.index(p) for p in some_phenos]
    elif args.pheno_filter:
        some_phenos_indices = [
            i
            for i, name in enumerate(pheno_names_list)
            if all(f in name for f in args.pheno_filter)
        ]
        some_phenos = [pheno_names_list[i] for i in some_phenos_indices]
    else:
        some_phenos_indices = list(range(len(pheno_names_list)))
        some_phenos = list(pheno_names_list)

    print(f"  → {len(some_phenos)} phenotype(s) selected")

    # ── Prepare matrices ───────────────────────────────────────────────────────
    pheno_train = testing_pheno[:, some_phenos_indices]
    geno_train = testing_geno
    geno_test = training_geno

    p_train = np.mean(geno_train, axis=0)
    geno_train_c = geno_train - p_train
    geno_test_c = geno_test - p_train

    print(f"\nGenotype shapes — train: {geno_train.shape}, test: {geno_test.shape}")
    print("Genotype matrices centred using training-set allele frequencies.")

    # ── Fit models ─────────────────────────────────────────────────────────────
    n_outcomes = pheno_train.shape[1]
    best_alphas = np.zeros(n_outcomes)
    best_l1_ratios = np.zeros(n_outcomes) if method == "elasticnet" else None

    print(f"\n{'─' * 70}")
    print(f"{method.capitalize()} regression — {n_outcomes} phenotype(s), {n_folds}-fold CV")
    if method != "lars":
        print(f"Alpha grid: {alphas[0]:.2e} → {alphas[-1]:.2e} ({len(alphas)} values)")
    if method == "elasticnet":
        print(f"L1-ratio grid: {l1_ratios}")
    if method in ("lasso", "lars", "elasticnet"):
        print(f"Max iterations: {max_iter}")
    print(f"{'─' * 70}\n")

    os.makedirs(intermediates_dir, exist_ok=True)

    for i, pheno_name in enumerate(some_phenos):
        print(f"[{i + 1}/{n_outcomes}] {pheno_name}")

        # LARS output filenames include max_iter to match the original naming convention
        tag = f"lars_maxiter{max_iter}" if method == "lars" else method
        coef_file = f"{intermediates_dir}{tag}_coeffs_{pheno_name}.csv"
        alpha_file = f"{intermediates_dir}{tag}_alphas_{pheno_name}.txt"
        pred_file = f"{intermediates_dir}{tag}_predictions_{pheno_name}.csv"

        if os.path.exists(coef_file):
            print("  Skipping — coefficients already exist.\n")
            continue

        y_train = pheno_train[:, i]

        # ── Method-specific CV + final fit ─────────────────────────────────────
        if method == "lasso":
            cv_model = LassoCV(alphas=alphas, cv=n_folds, n_jobs=-1, verbose=0, max_iter=max_iter)
            cv_model.fit(geno_train_c, y_train)
            best_alpha = cv_model.alpha_
            model = Lasso(alpha=best_alpha, max_iter=max_iter)
            model.fit(geno_train_c, y_train)

        elif method == "lars":
            # LarsCV traverses the full regularisation path and selects the best
            # alpha via CV — the fitted object is already the final model.
            cv_model = LarsCV(cv=n_folds, max_iter=max_iter, n_jobs=-1)
            cv_model.fit(geno_train_c, y_train)
            best_alpha = cv_model.alpha_
            model = cv_model

        elif method == "ridge":
            cv_model = RidgeCV(alphas=alphas, cv=n_folds)
            cv_model.fit(geno_train_c, y_train)
            best_alpha = cv_model.alpha_
            model = Ridge(alpha=best_alpha)
            model.fit(geno_train_c, y_train)

        elif method == "elasticnet":
            cv_model = ElasticNetCV(
                alphas=alphas,
                l1_ratio=l1_ratios,
                cv=n_folds,
                n_jobs=-1,
                verbose=0,
                max_iter=max_iter,
            )
            cv_model.fit(geno_train_c, y_train)
            best_alpha = cv_model.alpha_
            best_l1_ratio = cv_model.l1_ratio_
            best_l1_ratios[i] = best_l1_ratio
            model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1_ratio, max_iter=max_iter)
            model.fit(geno_train_c, y_train)

        best_alphas[i] = best_alpha
        print(f"  Best alpha: {best_alpha:.6e}")
        if method == "elasticnet":
            print(f"  Best l1_ratio: {best_l1_ratio:.4f}")

        # ── Save hyperparameters ───────────────────────────────────────────────
        with open(alpha_file, "w") as f:
            f.write(f"{pheno_name}\n{best_alpha}\n")
            if method == "elasticnet":
                f.write(f"{best_l1_ratio}\n")
        print(f"  Hyperparameters → {alpha_file}")
        if biobank:
            upload_to_dnanexus(alpha_file, dx_project, dx_folder)

        # ── Save coefficients ──────────────────────────────────────────────────
        n_nonzero = int(np.sum(model.coef_ != 0))
        print(f"  Non-zero coefficients: {n_nonzero} / {len(model.coef_)}")
        pd.DataFrame({"SNP": snp_ids, "coefficient": model.coef_}).to_csv(coef_file, index=False)
        print(f"  Coefficients     → {coef_file}")
        if biobank:
            upload_to_dnanexus(coef_file, dx_project, dx_folder)

        # ── Save predictions ───────────────────────────────────────────────────
        predictions = model.predict(geno_test_c)
        pd.DataFrame({"IID": test_iids, "prediction": predictions}).to_csv(
            pred_file, index=False
        )
        print(f"  Predictions      → {pred_file}\n")
        if biobank:
            upload_to_dnanexus(pred_file, dx_project, dx_folder)

    # ── Summary ────────────────────────────────────────────────────────────────
    fitted = np.sum(best_alphas != 0)
    print("─" * 70)
    print(f"Done. {fitted}/{n_outcomes} phenotype(s) fitted (rest skipped).")
    if fitted > 0:
        used = best_alphas[best_alphas != 0]
        print(f"Best alphas — min: {used.min():.3e}, max: {used.max():.3e}, mean: {used.mean():.3e}")
    if method == "elasticnet" and fitted > 0:
        used_l1 = best_l1_ratios[best_l1_ratios != 0]
        print(
            f"Best l1_ratios — min: {used_l1.min():.3f}, "
            f"max: {used_l1.max():.3f}, mean: {used_l1.mean():.3f}"
        )


if __name__ == "__main__":
    main()
