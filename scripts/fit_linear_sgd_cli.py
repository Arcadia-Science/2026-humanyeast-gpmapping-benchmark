#!/usr/bin/env python3

import argparse
import json
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np
from scipy.stats import pearsonr
import pandas as pd
from pathlib import Path

#########################################################################################
# ARGUMENT PARSING
#########################################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Final ridge regression fit using tuned hyperparameters'
    )

    # Data arguments
    parser.add_argument('--train-geno', type=str, required=True,
                        help='Path to training genotype feather file (IID + SNP columns)')
    parser.add_argument('--test-geno', type=str, required=True,
                        help='Path to test genotype feather file (IID + SNP columns)')
    parser.add_argument('--train-pheno', type=str, required=True,
                        help='Path to training phenotype feather file (IID + trait columns)')
    parser.add_argument('--test-pheno', type=str, required=True,
                        help='Path to test phenotype feather file (IID + trait columns)')

    # Phenotype and tuning results
    parser.add_argument('--phenotype-name', type=str, required=True,
                        help='Name of the phenotype column to fit (must match a column in the pheno files)')
    parser.add_argument('--tuning-dir', type=str, required=True,
                        help='Directory containing Optuna JSON results '
                             '(files of the form optuna_best_pheno_{index}_ridge.json)')

    # Output
    parser.add_argument('--which-seed', type=str, required=True,
                        help='Seed identifier used in output filenames')
    parser.add_argument('--output-dir', type=str, default='final_fit_results',
                        help='Directory to save results (default: final_fit_results)')

    # Training overrides — these fall back to sensible defaults if not supplied
    parser.add_argument('--max-epochs', type=int, default=200,
                        help='Maximum training epochs (default: 200)')
    parser.add_argument('--batch-size', type=int, default=128,
                        help='Batch size (default: 128)')
    parser.add_argument('--patience', type=int, default=3,
                        help='Early stopping patience (default: 3)')
    parser.add_argument('--min-delta', type=float, default=0.0001,
                        help='Minimum improvement for early stopping (default: 0.0001)')
    parser.add_argument('--device', type=str, default='auto',
                        choices=['auto', 'cpu', 'cuda'],
                        help='Device to use (default: auto)')

    return parser.parse_args()

#########################################################################################
# DATA LOADING
#########################################################################################

def load_feather_data(geno_path, pheno_path):
    """
    Load and inner-join genotype and phenotype feather files on IID.

    Returns:
        geno_cols:  list of SNP column names
        pheno_cols: list of trait column names
        merged_df:  DataFrame with IID, all SNP columns, all trait columns
    """
    geno_df  = pd.read_feather(geno_path)
    pheno_df = pd.read_feather(pheno_path)

    geno_cols  = [c for c in geno_df.columns if c != 'IID']
    pheno_cols = [c for c in pheno_df.columns if c != 'IID']

    merged_df = pd.merge(geno_df, pheno_df, on='IID', how='inner')

    if len(merged_df) == 0:
        raise ValueError(
            f"No matching IIDs found between genotype file ({geno_path}) "
            f"and phenotype file ({pheno_path}). Check that IID values match."
        )

    return geno_cols, pheno_cols, merged_df

#########################################################################################
# DATASET
#########################################################################################

class FeatherPhenotypeDataset(torch.utils.data.Dataset):
    """
    Dataset backed by in-memory DataFrames (already merged geno + pheno).
    Materialises tensors once at construction for fast __getitem__.
    """

    def __init__(self, merged_df, geno_cols, phenotype_col):
        self.geno  = torch.tensor(merged_df[geno_cols].values,      dtype=torch.float32)
        self.pheno = torch.tensor(merged_df[[phenotype_col]].values, dtype=torch.float32)

    def __len__(self):
        return len(self.pheno)

    def __getitem__(self, idx):
        return self.pheno[idx], self.geno[idx]

#########################################################################################
# MODEL
#########################################################################################

class RidgeRegression(nn.Module):
    """Linear model for ridge regression"""
    def __init__(self, n_loci, n_phen):
        super(RidgeRegression, self).__init__()
        self.linear = nn.Linear(n_loci, n_phen)

    def forward(self, x):
        return self.linear(x)

#########################################################################################
# LOSS
#########################################################################################

def l2_penalty(model):
    """L2 penalty for ridge regression (squared L2 norm of weights, excluding bias)"""
    penalty = 0
    for name, param in model.named_parameters():
        if 'bias' not in name:
            penalty += torch.sum(param ** 2)
    return penalty

#########################################################################################
# TRAINING
#########################################################################################

def train_ridge(model, train_loader, test_loader,
                alpha, learning_rate,
                max_epochs, min_delta, patience,
                device):
    """
    Train ridge regression model using Adam with ReduceLROnPlateau scheduling
    and early stopping on test MSE.
    """

    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=3
    )

    history = {
        'train_loss':     [],
        'test_loss':      [],
        'train_mse':      [],
        'test_mse':       [],
        'epochs_trained': 0
    }

    best_loss        = float('inf')
    best_epoch       = 0
    best_model_state = None
    patience_counter = 0

    for epoch in range(max_epochs):
        model.train()
        train_loss = 0
        train_mse  = 0

        for phens, gens in train_loader:
            phens = phens.to(device)
            gens  = gens.to(device)

            output     = model(gens)
            mse_loss   = F.mse_loss(output, phens)
            l2_term    = l2_penalty(model)
            total_loss = mse_loss + alpha * l2_term

            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            train_loss += total_loss.item()
            train_mse  += mse_loss.item()

        avg_train_loss = train_loss / len(train_loader)
        avg_train_mse  = train_mse  / len(train_loader)
        history['train_loss'].append(avg_train_loss)
        history['train_mse'].append(avg_train_mse)

        if test_loader is not None:
            model.eval()
            test_loss = 0
            test_mse  = 0

            with torch.no_grad():
                for phens, gens in test_loader:
                    phens = phens.to(device)
                    gens  = gens.to(device)
                    output   = model(gens)
                    mse_loss = F.mse_loss(output, phens)
                    test_mse += mse_loss.item()

                    l2_term    = l2_penalty(model)
                    total_loss = mse_loss + alpha * l2_term
                    test_loss += total_loss.item()

            avg_test_loss = test_loss / len(test_loader)
            avg_test_mse  = test_mse  / len(test_loader)
            history['test_loss'].append(avg_test_loss)
            history['test_mse'].append(avg_test_mse)

            print(f'Epoch: {epoch+1}/{max_epochs}, '
                  f'Train Loss: {avg_train_loss:.6f}, Test Loss: {avg_test_loss:.6f}, '
                  f'Test MSE: {avg_test_mse:.6f}')

            scheduler.step(avg_test_mse)

            if avg_test_mse < (best_loss - min_delta):
                best_loss        = avg_test_mse
                best_epoch       = epoch
                patience_counter = 0
                best_model_state = {k: v.cpu().detach().clone()
                                    for k, v in model.state_dict().items()}
                print(f"New best model at epoch {epoch+1} with test MSE: {best_loss:.6f}")
            else:
                patience_counter += 1
                print(f"No improvement for {patience_counter} epochs (best: {best_loss:.6f})")

            if patience_counter >= patience:
                print(f"Early stopping triggered after {epoch+1} epochs")
                break

    history['epochs_trained'] = epoch + 1

    if best_model_state is not None:
        print(f"Restoring best model from epoch {best_epoch+1}")
        model.load_state_dict(best_model_state)

    return model, best_loss, history

#########################################################################################
# EVALUATION
#########################################################################################

def evaluate_model(model, test_loader, device):
    """Evaluate model on test set and return true and predicted values"""
    model.eval()

    true_phenotypes      = []
    predicted_phenotypes = []

    with torch.no_grad():
        for phens, gens in test_loader:
            phens = phens.to(device)
            gens  = gens.to(device)

            predictions = model(gens)
            true_phenotypes.append(phens.cpu().numpy())
            predicted_phenotypes.append(predictions.cpu().numpy())

    true_phenotypes      = np.concatenate(true_phenotypes)
    predicted_phenotypes = np.concatenate(predicted_phenotypes)

    return true_phenotypes, predicted_phenotypes

#########################################################################################
# MAIN
#########################################################################################

def main():
    args = parse_args()

    # Device
    if args.device == 'auto':
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")

    # ---------------------------------------------------------------------------------
    # Load feather data
    # ---------------------------------------------------------------------------------
    print("Loading training data...")
    train_geno_cols, train_pheno_cols, train_merged = load_feather_data(
        args.train_geno, args.train_pheno
    )

    print("Loading test data...")
    test_geno_cols, test_pheno_cols, test_merged = load_feather_data(
        args.test_geno, args.test_pheno
    )

    if train_geno_cols != test_geno_cols:
        raise ValueError(
            f"Genotype columns differ between train ({len(train_geno_cols)} SNPs) "
            f"and test ({len(test_geno_cols)} SNPs). Ensure files share the same SNP set."
        )
    if train_pheno_cols != test_pheno_cols:
        raise ValueError(
            "Phenotype columns differ between train and test files. "
            "Ensure both files contain the same trait columns in the same order."
        )

    geno_cols  = train_geno_cols
    pheno_cols = train_pheno_cols
    n_loci     = len(geno_cols)

    # ---------------------------------------------------------------------------------
    # Resolve phenotype name → index
    # ---------------------------------------------------------------------------------
    if args.phenotype_name not in pheno_cols:
        raise ValueError(
            f"--phenotype-name '{args.phenotype_name}' not found in phenotype file.\n"
            f"Available traits:\n" + "\n".join(f"  {c}" for c in pheno_cols)
        )

    phenotype_name  = args.phenotype_name
    phenotype_index = pheno_cols.index(phenotype_name)

    print(f"\nPhenotype:        {phenotype_name} (index {phenotype_index})")
    print(f"Training samples: {len(train_merged)}")
    print(f"Test samples:     {len(test_merged)}")
    print(f"Number of loci:   {n_loci}")

    # ---------------------------------------------------------------------------------
    # Load tuning results JSON
    # ---------------------------------------------------------------------------------
    json_path = Path(args.tuning_dir) / f"optuna_best_pheno_{phenotype_index}_ridge.json"

    if not json_path.exists():
        raise FileNotFoundError(
            f"Tuning results file not found: {json_path}\n"
            f"Expected pattern: optuna_best_pheno_{{index}}_ridge.json\n"
            f"Resolved phenotype '{phenotype_name}' to index {phenotype_index}."
        )

    with open(json_path) as f:
        tuning = json.load(f)

    alpha         = tuning['best_alpha']
    learning_rate = tuning['best_learning_rate']

    print(f"\nLoaded hyperparameters from: {json_path}")
    print(f"  alpha:         {alpha:.6e}")
    print(f"  learning_rate: {learning_rate:.6e}")

    # ---------------------------------------------------------------------------------
    # Build datasets and dataloaders
    # ---------------------------------------------------------------------------------
    train_dataset = FeatherPhenotypeDataset(train_merged, geno_cols, phenotype_name)
    test_dataset  = FeatherPhenotypeDataset(test_merged,  geno_cols, phenotype_name)

    train_loader = torch.utils.data.DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        num_workers=3,
        pin_memory=True if device.type == 'cuda' else False
    )
    test_loader = torch.utils.data.DataLoader(
        test_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=3
    )

    # ---------------------------------------------------------------------------------
    # Train
    # ---------------------------------------------------------------------------------
    print("\n" + "="*80)
    print("TRAINING")
    print("="*80)

    model = RidgeRegression(n_loci=n_loci, n_phen=1).to(device)

    model, best_loss, history = train_ridge(
        model=model,
        train_loader=train_loader,
        test_loader=test_loader,
        alpha=alpha,
        learning_rate=learning_rate,
        max_epochs=args.max_epochs,
        min_delta=args.min_delta,
        patience=args.patience,
        device=device
    )

    # ---------------------------------------------------------------------------------
    # Evaluate
    # ---------------------------------------------------------------------------------
    print("\n" + "="*80)
    print("EVALUATION")
    print("="*80)

    true_phenotypes, predicted_phenotypes = evaluate_model(model, test_loader, device)

    from sklearn.metrics import r2_score, mean_squared_error

    corr, p_val = pearsonr(true_phenotypes.flatten(), predicted_phenotypes.flatten())
    r2   = r2_score(true_phenotypes, predicted_phenotypes)
    mse  = mean_squared_error(true_phenotypes, predicted_phenotypes)
    rmse = np.sqrt(mse)

    print(f"\nPhenotype:           {phenotype_name}")
    print(f"Pearson correlation: {corr:.4f} (p={p_val:.2e})")
    print(f"R² score:            {r2:.4f}")
    print(f"RMSE:                {rmse:.4f}")
    print(f"Alpha:               {alpha:.6e}")
    print(f"Learning rate:       {learning_rate:.6e}")
    print(f"Epochs trained:      {history['epochs_trained']}")
    print(f"Best test MSE:       {best_loss:.6f}")
    print("="*80 + "\n")

    # ---------------------------------------------------------------------------------
    # Save outputs
    # ---------------------------------------------------------------------------------
    safe_pheno_name = phenotype_name.replace("/", "_").replace(" ", "_")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Summary metrics
    results_df = pd.DataFrame({
        'phenotype_name':      [phenotype_name],
        'phenotype_index':     [phenotype_index],
        'pearson_correlation': [corr],
        'p_value':             [p_val],
        'r2_score':            [r2],
        'mse':                 [mse],
        'rmse':                [rmse],
        'true_mean':           [np.mean(true_phenotypes)],
        'pred_mean':           [np.mean(predicted_phenotypes)],
        'true_std':            [np.std(true_phenotypes)],
        'pred_std':            [np.std(predicted_phenotypes)],
        'alpha':               [alpha],
        'learning_rate':       [learning_rate],
        'batch_size':          [args.batch_size],
        'epochs_trained':      [history['epochs_trained']],
        'final_test_mse':      [best_loss]
    })
    results_file = out / f"ridge_results_seed_{args.which_seed}_{safe_pheno_name}.csv"
    results_df.to_csv(results_file, index=False)
    print(f"Results saved to:          {results_file}")

    # Model checkpoint
    model_file = out / f"ridge_model_seed_{args.which_seed}_{safe_pheno_name}.pt"
    torch.save({
        'model_state_dict': model.state_dict(),
        'alpha':            alpha,
        'learning_rate':    learning_rate,
        'phenotype_name':   phenotype_name,
        'phenotype_index':  phenotype_index,
        'n_loci':           n_loci,
        'history':          history
    }, model_file)
    print(f"Model saved to:            {model_file}")

    # SNP weights — one row per SNP, bias appended as a final row
    weights_df = pd.concat([
        pd.DataFrame({
            'snp_id': geno_cols,
            'weight': model.linear.weight.detach().cpu().numpy().flatten(),
        }),
        pd.DataFrame({'snp_id': ['__bias__'], 'weight': [model.linear.bias.item()]})
    ], ignore_index=True)
    weights_file = out / f"ridge_weights_seed_{args.which_seed}_{safe_pheno_name}.csv"
    weights_df.to_csv(weights_file, index=False)
    print(f"Weights saved to:          {weights_file}")

    # Predictions with individual IDs
    predictions_df = pd.DataFrame({
        'IID':            test_merged['IID'].values,
        'true_phenotype': true_phenotypes.flatten(),
        'pred_phenotype': predicted_phenotypes.flatten(),
    })
    predictions_file = out / f"ridge_predictions_seed_{args.which_seed}_{safe_pheno_name}.csv"
    predictions_df.to_csv(predictions_file, index=False)
    print(f"Predictions saved to:      {predictions_file}")

    # Training history
    history_df = pd.DataFrame({
        'epoch':      range(1, history['epochs_trained'] + 1),
        'train_loss': history['train_loss'],
        'test_loss':  history['test_loss'],
        'train_mse':  history['train_mse'],
        'test_mse':   history['test_mse']
    })
    history_file = out / f"training_history_seed_{args.which_seed}_{safe_pheno_name}.csv"
    history_df.to_csv(history_file, index=False)
    print(f"Training history saved to: {history_file}")


if __name__ == "__main__":
    main()
