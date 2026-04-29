#!/usr/bin/env python3
"""
fit_linear_sgd_optuna.py

Optuna hyperparameter tuning for PyTorch ridge (or lasso) regression,
run once per phenotype. Searches over alpha and learning_rate using
TPE sampling; results are written as JSON and read by fit_linear_sgd_cli.py
for the final fit.
Called once per phenotype by run_tuning.sh during the tune_pytorch step.

─── INPUTS ─────────────────────────────────────────────────────────────────

  Genotype / phenotype feather files (train + test), produced by split_phenotypes:
    test_train_seed_{seed}/{prefix}_seed_{seed}_{train|test}_{genotypes|phenotypes}_*.feather

─── OUTPUTS ────────────────────────────────────────────────────────────────

  All written to --output-dir (default: optuna_results/):

    optuna_trials_pheno_{index}_{reg_mode}.csv  — all trial results
    optuna_best_pheno_{index}_{reg_mode}.json   — best alpha and learning_rate
                                                   (read by fit_linear_sgd_cli.py)
    optuna_studies.db                           — Optuna SQLite database

─── USAGE ──────────────────────────────────────────────────────────────────

  python scripts/fit_linear_sgd_optuna.py \\
      --train-geno  test_train_seed_1510/yeast_simulated_data_seed_1510_train_genotypes_centered.feather \\
      --test-geno   test_train_seed_1510/yeast_simulated_data_seed_1510_test_genotypes_centered.feather \\
      --train-pheno test_train_seed_1510/yeast_simulated_data_seed_1510_train_phenotypes_normalized.feather \\
      --test-pheno  test_train_seed_1510/yeast_simulated_data_seed_1510_test_phenotypes_normalized.feather \\
      --phenotype-name trait_name \\
      --n-trials 25 \\
      --output-dir pytorch_tuning_results

Options:
  --train-geno      FILE   Training genotype feather (IID + SNP columns) (required)
  --test-geno       FILE   Test genotype feather (IID + SNP columns) (required)
  --train-pheno     FILE   Training phenotype feather (IID + trait columns) (required)
  --test-pheno      FILE   Test phenotype feather (IID + trait columns) (required)
  --phenotype-name  STR    Trait column to tune; must exist in both pheno files (required)
  --reg-mode        STR    ridge (L2) or lasso (smoothed L1) (default: ridge)
  --huber-beta      FLOAT  Huber transition point for lasso mode (default: 1e-4)
  --n-trials        INT    Number of Optuna trials (default: 25)
  --n-jobs          INT    Synchronous parallel Optuna jobs (default: 5)
  --timeout         INT    Time limit in seconds (default: None)
  --study-name      STR    Optuna study name (default: ridge_pheno_{index})
  --alpha-min       FLOAT  Minimum alpha search bound (default: 1e-4)
  --alpha-max       FLOAT  Maximum alpha search bound (default: 1)
  --lr-min          FLOAT  Minimum learning rate search bound (default: 1e-6)
  --lr-max          FLOAT  Maximum learning rate search bound (default: 1e-2)
  --max-epochs      INT    Maximum epochs per trial (default: 15)
  --batch-size      INT    Mini-batch size (default: 128)
  --patience        INT    Early stopping patience in epochs (default: 3)
  --min-delta       FLOAT  Minimum validation MSE improvement to reset patience (default: 0.0001)
  --val-fraction    FLOAT  Fraction of training data held out for validation (default: 0.15)
  --val-seed        INT    Random seed for the validation split (default: 42)
  --output-dir      DIR    Where to write all output files (default: optuna_results)
  --db-path         FILE   Optuna SQLite database path (default: output_dir/optuna_studies.db)
  --device          STR    auto | cpu | cuda (default: auto)
  --verbose                Print detailed trial information
  --pruning                Enable Optuna MedianPruner to stop poor trials early
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
from scipy.stats import pearsonr
import pandas as pd
from pathlib import Path
import argparse
import optuna
from optuna.pruners import MedianPruner
from optuna.samplers import TPESampler
from sklearn.metrics import r2_score, mean_squared_error
from datetime import datetime
import json

#########################################################################################
# ARGUMENT PARSING
#########################################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Optuna hyperparameter tuning for ridge regression'
    )

    # Data arguments — four separate feather files
    parser.add_argument('--train-geno', type=str, required=True,
                       help='Path to training genotype feather file (IID + SNP columns)')
    parser.add_argument('--test-geno', type=str, required=True,
                       help='Path to test genotype feather file (IID + SNP columns)')
    parser.add_argument('--train-pheno', type=str, required=True,
                       help='Path to training phenotype feather file (IID + trait columns)')
    parser.add_argument('--test-pheno', type=str, required=True,
                       help='Path to test phenotype feather file (IID + trait columns)')
    parser.add_argument('--phenotype-name', type=str, required=True,
                       help='Name of the phenotype column to tune (must match a column in the pheno files)')

    # Model arguments
    parser.add_argument('--reg-mode', type=str, default='ridge',
                   choices=['ridge', 'lasso'],
                   help='Regularization mode: ridge (L2) or lasso (smoothed L1) (default: ridge)')
    parser.add_argument('--huber-beta', type=float, default=1e-4,
                   help='Huber transition point for lasso mode (default: 1e-4)')

    # Optuna arguments
    parser.add_argument('--n-trials', type=int, default=25,
                       help='Number of Optuna trials (default: 25)')
    parser.add_argument('--n-jobs', type=int, default=5,
                       help='Number of synchronous Optuna jobs (default: 5)')
    parser.add_argument('--timeout', type=int, default=None,
                       help='Time limit in seconds (default: None)')
    parser.add_argument('--study-name', type=str, default=None,
                       help='Optuna study name (default: ridge_pheno_{index})')

    # Search space
    parser.add_argument('--alpha-min', type=float, default=1e-4,
                       help='Minimum alpha value (default: 1e-4)')
    parser.add_argument('--alpha-max', type=float, default=1,
                       help='Maximum alpha value (default: 1)')
    parser.add_argument('--lr-min', type=float, default=1e-6,
                       help='Minimum learning rate (default: 1e-6)')
    parser.add_argument('--lr-max', type=float, default=1e-2,
                       help='Maximum learning rate (default: 1e-2)')

    # Training arguments
    parser.add_argument('--max-epochs', type=int, default=15,
                       help='Maximum epochs per trial (default: 15)')
    parser.add_argument('--batch-size', type=int, default=128,
                       help='Batch size (default: 128)')
    parser.add_argument('--patience', type=int, default=3,
                       help='Early stopping patience (default: 3)')
    parser.add_argument('--min-delta', type=float, default=0.0001,
                       help='Minimum improvement for early stopping (default: 0.0001)')

    # Validation split
    parser.add_argument('--val-fraction', type=float, default=0.15,
                       help='Fraction of training data for validation (default: 0.15)')
    parser.add_argument('--val-seed', type=int, default=42,
                       help='Random seed for validation split (default: 42)')

    # Output arguments
    parser.add_argument('--output-dir', type=str, default='optuna_results',
                       help='Directory to save results (default: optuna_results)')
    parser.add_argument('--db-path', type=str, default=None,
                       help='Path to Optuna database (default: output_dir/optuna.db)')

    # Behavior flags
    parser.add_argument('--device', type=str, default='auto',
                       choices=['auto', 'cpu', 'cuda'],
                       help='Device to use (default: auto)')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed trial information')
    parser.add_argument('--pruning', action='store_true',
                       help='Enable Optuna pruning (stop bad trials early)')

    return parser.parse_args()

#########################################################################################
# SETUP
#########################################################################################

def setup_device(device_arg):
    """Setup computation device"""
    if device_arg == 'auto':
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(device_arg)
    return device

#########################################################################################
# DATA LOADING
#########################################################################################

def load_feather_data(geno_path, pheno_path):
    """
    Load and inner-join genotype and phenotype feather files on IID.

    Returns:
        geno_cols: list of SNP column names
        pheno_cols: list of trait column names
        merged_df: DataFrame with IID, all SNP columns, all trait columns
    """
    geno_df = pd.read_feather(geno_path)
    pheno_df = pd.read_feather(pheno_path)

    # Identify SNP and trait columns (everything except IID)
    geno_cols = [c for c in geno_df.columns if c != 'IID']
    pheno_cols = [c for c in pheno_df.columns if c != 'IID']

    # Inner join on IID to align samples
    merged_df = pd.merge(geno_df, pheno_df, on='IID', how='inner')

    if len(merged_df) == 0:
        raise ValueError(
            f"No matching IIDs found between genotype file ({geno_path}) "
            f"and phenotype file ({pheno_path}). Check that IID values match."
        )

    return geno_cols, pheno_cols, merged_df

#########################################################################################
# DATASET CLASS
#########################################################################################

class FeatherPhenotypeDataset(torch.utils.data.Dataset):
    """
    Dataset backed by in-memory DataFrames (already merged geno + pheno).
    Accepts a subset of row indices for train/val splitting.
    """

    def __init__(self, merged_df, geno_cols, phenotype_col, indices):
        """
        Args:
            merged_df:      Full merged DataFrame (geno + pheno, aligned by IID)
            geno_cols:      List of SNP column names to use as features
            phenotype_col:  Name of the single trait column to predict
            indices:        Array of integer row indices for this split
        """
        # Subset once up front — avoids repeated indexing inside __getitem__
        subset = merged_df.iloc[indices].reset_index(drop=True)

        self.geno = torch.tensor(
            subset[geno_cols].values, dtype=torch.float32
        )
        self.pheno = torch.tensor(
            subset[[phenotype_col]].values, dtype=torch.float32
        )

    def __len__(self):
        return len(self.pheno)

    def __getitem__(self, idx):
        return self.pheno[idx], self.geno[idx]

#########################################################################################
# MODEL
#########################################################################################

class RidgeRegression(nn.Module):
    """Linear model for ridge regression"""
    def __init__(self, n_loci, n_phen=1):
        super(RidgeRegression, self).__init__()
        self.linear = nn.Linear(n_loci, n_phen)

    def forward(self, x):
        return self.linear(x)

def l2_penalty(model):
    """Squared L2 norm of all non-bias parameters (ridge regularisation term)."""
    penalty = 0
    for name, param in model.named_parameters():
        if 'bias' not in name:
            penalty += torch.sum(param ** 2)
    return penalty

def l1_penalty(model, beta=1e-4):
    """Smoothed L1 norm of all non-bias parameters via Huber loss (lasso regularisation term)."""
    penalty = 0
    for name, param in model.named_parameters():
        if 'bias' not in name:
            penalty += F.huber_loss(
                param,
                torch.zeros_like(param),
                delta=beta,
                reduction='sum'
            )
    return penalty

#########################################################################################
# TRAINING WITH PRUNING
#########################################################################################

def train_with_pruning(model, train_loader, val_loader,
                       alpha, learning_rate,
                       max_epochs, min_delta, patience,
                       device, trial=None,
                       reg_mode='ridge', huber_beta=1e-4):
    """
    Train model with optional Optuna pruning.

    Args:
        trial: Optuna trial object (if None, no pruning)
    """

    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=3
    )

    best_val_mse = float('inf')
    best_epoch = 0
    best_model_state = None
    patience_counter = 0

    for epoch in range(max_epochs):
        # Training
        model.train()
        train_loss = 0

        for phens, gens in train_loader:
            phens = phens.to(device)
            gens = gens.to(device)

            output = model(gens)
            mse_loss = F.mse_loss(output, phens)

            reg_term = l2_penalty(model) if reg_mode == 'ridge' else l1_penalty(model, beta=huber_beta)
            total_loss = mse_loss + alpha * reg_term

            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            train_loss += mse_loss.item()

        # Validation
        model.eval()
        val_mse = 0

        with torch.no_grad():
            for phens, gens in val_loader:
                phens = phens.to(device)
                gens = gens.to(device)
                output = model(gens)
                mse_loss = F.mse_loss(output, phens)
                val_mse += mse_loss.item()

        val_mse = val_mse / len(val_loader)

        # Optuna pruning (if trial provided)
        if trial is not None:
            trial.report(val_mse, epoch)
            if trial.should_prune():
                raise optuna.TrialPruned()

        # Learning rate scheduling
        scheduler.step(val_mse)

        # Early stopping check
        if val_mse < (best_val_mse - min_delta):
            best_val_mse = val_mse
            best_epoch = epoch
            patience_counter = 0
            best_model_state = {k: v.cpu().detach().clone() for k, v in model.state_dict().items()}
        else:
            patience_counter += 1

        if patience_counter >= patience:
            break

    # Restore best model
    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    return model, best_val_mse, epoch + 1

#########################################################################################
# EVALUATION
#########################################################################################

def evaluate_model(model, loader, device):
    """Evaluate model and return predictions and metrics"""
    model.eval()

    all_true = []
    all_pred = []

    with torch.no_grad():
        for phens, gens in loader:
            phens = phens.to(device)
            gens = gens.to(device)
            predictions = model(gens)

            all_true.append(phens.cpu().numpy())
            all_pred.append(predictions.cpu().numpy())

    y_true = np.concatenate(all_true)
    y_pred = np.concatenate(all_pred)

    # Calculate metrics
    corr, p_val = pearsonr(y_true.flatten(), y_pred.flatten())
    r2 = r2_score(y_true, y_pred)
    mse = mean_squared_error(y_true, y_pred)

    return {
        'mse': float(mse),
        'correlation': float(corr),
        'p_value': float(p_val),
        'r2': float(r2)
    }

#########################################################################################
# OPTUNA OBJECTIVE
#########################################################################################

class OptunaObjective:
    """Objective function for Optuna optimization"""

    def __init__(self, train_loader, val_loader, n_loci,
                 max_epochs, patience, min_delta, device,
                 alpha_range, lr_range, use_pruning, reg_mode='ridge', huber_beta=1e-4):
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.reg_mode = reg_mode
        self.huber_beta = huber_beta
        self.n_loci = n_loci
        self.max_epochs = max_epochs
        self.patience = patience
        self.min_delta = min_delta
        self.device = device
        self.alpha_range = alpha_range
        self.lr_range = lr_range
        self.use_pruning = use_pruning

    def __call__(self, trial):
        # Suggest hyperparameters
        alpha = trial.suggest_float('alpha',
                                    self.alpha_range[0],
                                    self.alpha_range[1],
                                    log=True)
        learning_rate = trial.suggest_float('learning_rate',
                                           self.lr_range[0],
                                           self.lr_range[1],
                                           log=True)

        # Create model
        model = RidgeRegression(n_loci=self.n_loci, n_phen=1).to(self.device)

        # Train model
        model, val_mse, epochs_trained = train_with_pruning(
            model=model,
            train_loader=self.train_loader,
            val_loader=self.val_loader,
            reg_mode=self.reg_mode,
            huber_beta=self.huber_beta,
            alpha=alpha,
            learning_rate=learning_rate,
            max_epochs=self.max_epochs,
            min_delta=self.min_delta,
            patience=self.patience,
            device=self.device,
            trial=trial if self.use_pruning else None
        )

        # Store additional metrics in trial user attributes
        metrics = evaluate_model(model, self.val_loader, self.device)
        trial.set_user_attr('correlation', metrics['correlation'])
        trial.set_user_attr('r2', metrics['r2'])
        trial.set_user_attr('epochs_trained', epochs_trained)

        # Return validation MSE (lower is better)
        return val_mse

#########################################################################################
# RESULTS SAVING
#########################################################################################

def save_study_results(study, phenotype_name, phenotype_index, output_dir, args):
    """Save Optuna study results to CSV and JSON"""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Convert study to DataFrame
    trials_df = study.trials_dataframe()

    # Add phenotype information
    trials_df['phenotype_name'] = phenotype_name
    trials_df['phenotype_index'] = phenotype_index

    # Reorder columns for readability
    cols = ['phenotype_name', 'phenotype_index', 'number', 'value',
            'params_alpha', 'params_learning_rate',
            'user_attrs_correlation', 'user_attrs_r2', 'user_attrs_epochs_trained',
            'state', 'duration']
    remaining_cols = [c for c in trials_df.columns if c not in cols]
    trials_df = trials_df[cols + remaining_cols]

    # Save trials to CSV
    csv_file = output_dir / f"optuna_trials_pheno_{phenotype_index}_{args.reg_mode}.csv"
    trials_df.to_csv(csv_file, index=False)

    # Save best trial info to JSON
    best_trial = study.best_trial
    best_info = {
        'phenotype_name': phenotype_name,
        'phenotype_index': phenotype_index,
        'best_trial_number': best_trial.number,
        'reg_mode': args.reg_mode,
        'huber_beta': args.huber_beta if args.reg_mode == 'lasso' else None,
        'best_value': best_trial.value,
        'best_params': best_trial.params,
        'best_alpha': best_trial.params['alpha'],
        'best_learning_rate': best_trial.params['learning_rate'],
        'correlation': best_trial.user_attrs.get('correlation'),
        'r2': best_trial.user_attrs.get('r2'),
        'epochs_trained': best_trial.user_attrs.get('epochs_trained'),
        'n_trials': len(study.trials),
        'n_pruned_trials': len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED]),
        'n_complete_trials': len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]),
        'search_space': {
            'alpha': [args.alpha_min, args.alpha_max],
            'learning_rate': [args.lr_min, args.lr_max]
        },
        'timestamp': datetime.now().isoformat()
    }

    json_file = output_dir / f"optuna_best_pheno_{phenotype_index}_{args.reg_mode}.json"
    with open(json_file, 'w') as f:
        json.dump(best_info, f, indent=2)

    return csv_file, json_file, best_info

#########################################################################################
# MAIN
#########################################################################################

def main():
    args = parse_args()

    # Setup
    device = setup_device(args.device)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Using device: {device}")
    print(f"Output directory: {output_dir}")

    # ---------------------------------------------------------------------------------
    # Load feather data
    # ---------------------------------------------------------------------------------
    print("\nLoading training data...")
    train_geno_cols, train_pheno_cols, train_merged = load_feather_data(
        args.train_geno, args.train_pheno
    )

    print("Loading test data...")
    test_geno_cols, test_pheno_cols, test_merged = load_feather_data(
        args.test_geno, args.test_pheno
    )

    # Validate consistency across train/test
    if train_geno_cols != test_geno_cols:
        raise ValueError(
            f"Genotype columns differ between train ({len(train_geno_cols)} SNPs) "
            f"and test ({len(test_geno_cols)} SNPs). Ensure files share the same SNP set."
        )
    if train_pheno_cols != test_pheno_cols:
        raise ValueError(
            f"Phenotype columns differ between train and test files. "
            f"Ensure both files contain the same trait columns in the same order."
        )

    geno_cols = train_geno_cols
    pheno_cols = train_pheno_cols
    n_loci = len(geno_cols)
    n_train = len(train_merged)

    # Resolve phenotype name to index
    if args.phenotype_name not in pheno_cols:
        raise ValueError(
            f"--phenotype-name '{args.phenotype_name}' not found in phenotype file.\n"
            f"Available traits:\n" + "\n".join(f"  {c}" for c in pheno_cols)
        )
    phenotype_name = args.phenotype_name
    phenotype_index = pheno_cols.index(phenotype_name)

    print(f"\n{'='*80}")
    print(f"HYPERPARAMETER TUNING")
    print(f"{'='*80}")
    print(f"Phenotype: {phenotype_name} (index {phenotype_index})")
    print(f"Training samples: {n_train}")
    print(f"Test samples: {len(test_merged)}")
    print(f"Number of loci: {n_loci}")
    print(f"Number of trials: {args.n_trials}")
    print(f"Regularisation mode: {args.reg_mode}")
    print(f"Search space:")
    print(f"  Alpha: [{args.alpha_min}, {args.alpha_max}] (log scale)")
    print(f"  Learning rate: [{args.lr_min}, {args.lr_max}] (log scale)")
    print(f"{'='*80}\n")

    # ---------------------------------------------------------------------------------
    # Validation split (index-based, no copying of data)
    # ---------------------------------------------------------------------------------
    np.random.seed(args.val_seed)
    indices = np.arange(n_train)
    np.random.shuffle(indices)

    val_size = int(n_train * args.val_fraction)
    train_idx = indices[val_size:]
    val_idx = indices[:val_size]

    print(f"Validation split: {len(train_idx)} train, {len(val_idx)} validation\n")

    # ---------------------------------------------------------------------------------
    # Create datasets
    # ---------------------------------------------------------------------------------
    train_dataset = FeatherPhenotypeDataset(train_merged, geno_cols, phenotype_name, train_idx)
    val_dataset   = FeatherPhenotypeDataset(train_merged, geno_cols, phenotype_name, val_idx)

    train_loader = torch.utils.data.DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=3,
        pin_memory=True if device.type == 'cuda' else False
    )

    val_loader = torch.utils.data.DataLoader(
        val_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=3
    )

    # ---------------------------------------------------------------------------------
    # Setup Optuna
    # ---------------------------------------------------------------------------------
    study_name = args.study_name or f"ridge_pheno_{phenotype_name}"

    if args.db_path is None:
        db_path = output_dir / "optuna_studies.db"
    else:
        db_path = Path(args.db_path)

    storage = f"sqlite:///{db_path}"

    sampler = TPESampler(seed=42)
    pruner = MedianPruner(n_startup_trials=5, n_warmup_steps=8) if args.pruning else None

    study = optuna.create_study(
        study_name=study_name,
        storage=storage,
        sampler=sampler,
        pruner=pruner,
        direction='minimize',
        load_if_exists=True
    )

    # Starting point param combo
    study.enqueue_trial({
        'alpha': 0.1,
        'learning_rate': 0.00001
    })

    # Create objective
    objective = OptunaObjective(
        train_loader=train_loader,
        val_loader=val_loader,
        reg_mode=args.reg_mode,
        huber_beta=args.huber_beta,
        n_loci=n_loci,
        max_epochs=args.max_epochs,
        patience=args.patience,
        min_delta=args.min_delta,
        device=device,
        alpha_range=(args.alpha_min, args.alpha_max),
        lr_range=(args.lr_min, args.lr_max),
        use_pruning=args.pruning
    )

    # ---------------------------------------------------------------------------------
    # Run optimization
    # ---------------------------------------------------------------------------------
    print("Starting optimization...\n")

    if args.n_jobs > 1:
        print(f"Running {args.n_jobs} trials in parallel\n")
        study.optimize(
            objective,
            n_trials=args.n_trials,
            timeout=args.timeout,
            n_jobs=args.n_jobs,
            show_progress_bar=True
        )
    else:
        study.optimize(
            objective,
            n_trials=args.n_trials,
            timeout=args.timeout,
            show_progress_bar=True
        )

    # ---------------------------------------------------------------------------------
    # Print results
    # ---------------------------------------------------------------------------------
    print(f"\n{'='*80}")
    print("OPTIMIZATION COMPLETE")
    print(f"{'='*80}")
    print(f"Number of finished trials: {len(study.trials)}")
    print(f"Number of pruned trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED])}")
    print(f"Number of complete trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])}")

    best_trial = study.best_trial
    print(f"\nBest trial (#{best_trial.number}):")
    print(f"  Validation MSE: {best_trial.value:.6f}")
    print(f"  Alpha: {best_trial.params['alpha']:.6e}")
    print(f"  Learning rate: {best_trial.params['learning_rate']:.6e}")
    print(f"  Correlation: {best_trial.user_attrs.get('correlation', 'N/A'):.4f}")
    print(f"  R²: {best_trial.user_attrs.get('r2', 'N/A'):.4f}")
    print(f"  Epochs trained: {best_trial.user_attrs.get('epochs_trained', 'N/A')}")
    print(f"{'='*80}\n")

    # Save results
    csv_file, json_file, best_info = save_study_results(
        study, phenotype_name, phenotype_index, output_dir, args
    )

    print(f"Results saved:")
    print(f"  All trials: {csv_file}")
    print(f"  Best parameters: {json_file}")
    print(f"  Optuna database: {db_path}")

    # Summary line for parsing
    print(f"\nRESULT: pheno={phenotype_name} "
          f"best_alpha={best_trial.params['alpha']:.6e} "
          f"best_lr={best_trial.params['learning_rate']:.6e} "
          f"val_mse={best_trial.value:.6f} "
          f"corr={best_trial.user_attrs.get('correlation', 0):.4f}")

    return best_info

if __name__ == "__main__":
    main()
