#!/usr/bin/env python3
"""
03_train.py - Training pipeline for Literature-Guided VAE

Two-phase training:
1. Pre-training: Unsupervised VAE on all cells
2. Fine-tuning: Multi-task learning on beta cells with labels
"""

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import json
from datetime import datetime
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Import from 02_model.py
import importlib.util
spec = importlib.util.spec_from_file_location("model", Path(__file__).parent / "02_model.py")
model_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model_module)
create_model = model_module.create_model
LiteratureGuidedVAE = model_module.LiteratureGuidedVAE
MultiTaskLoss = model_module.MultiTaskLoss

# Paths - use absolute paths for reliability
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent  # workload/code/deeplearning -> workload
RESULTS_DIR = WORKLOAD_DIR / "results" / "deep_learning"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Training configuration (adjusted for small dataset ~269 cells)
TRAINING_CONFIG = {
    # Phase 1: Pre-training
    "pretrain": {
        "epochs": 100,
        "batch_size": 32,  # Reduced for small dataset
        "learning_rate": 1e-3,
        "early_stopping": 15
    },
    # Phase 2: Fine-tuning
    "finetune": {
        "epochs": 50,
        "batch_size": 32,  # Reduced for small dataset
        "learning_rate": 1e-4,
        "early_stopping": 10
    },
    # Model architecture
    "latent_dim": 32,
    "workload_dims": 4,
    "encoder_dims": [256, 128, 64],
    "decoder_dims": [64, 128, 256],
    "dropout": 0.2,
    "n_states": 5,
    # Loss weights
    "lambda_recon": 1.0,
    "lambda_kl": 0.1,
    "lambda_condition": 1.0,
    "lambda_state": 0.5,
    "lambda_literature": 0.3,
    "lambda_cwi": 0.5
}


class EarlyStopping:
    """Early stopping to prevent overfitting."""

    def __init__(self, patience: int = 10, min_delta: float = 0.0):
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.best_loss = None
        self.should_stop = False

    def __call__(self, val_loss: float) -> bool:
        if self.best_loss is None:
            self.best_loss = val_loss
        elif val_loss > self.best_loss - self.min_delta:
            self.counter += 1
            if self.counter >= self.patience:
                self.should_stop = True
        else:
            self.best_loss = val_loss
            self.counter = 0
        return self.should_stop


def load_training_data():
    """Load prepared training data."""
    print("=" * 60)
    print("Loading training data")
    print("=" * 60)

    data_path = RESULTS_DIR / "training_data.h5ad"
    if not data_path.exists():
        raise FileNotFoundError(
            f"Training data not found at {data_path}. "
            "Run 01_prepare_data.py first."
        )

    adata = sc.read_h5ad(data_path)
    print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Get gene names
    gene_names = adata.var_names.tolist()

    # Get expression matrix (ensure dense)
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X

    # Get labels
    labels = {}

    # CWI (continuous)
    if "CWI_literature" in adata.obs.columns:
        labels["cwi"] = adata.obs["CWI_literature"].values.astype(np.float32)
        print(f"  CWI range: {labels['cwi'].min():.2f} - {labels['cwi'].max():.2f}")

    # Condition (binary)
    if "condition_binary" in adata.obs.columns:
        labels["condition"] = adata.obs["condition_binary"].values.astype(np.float32)
        print(f"  Condition: {(labels['condition'] == 1).sum()} T2D, {(labels['condition'] == 0).sum()} Normal")

    # Workload state (categorical)
    if "workload_state" in adata.obs.columns:
        state_map = {
            "S1_Resting": 0, "S2_Active": 1, "S3_Stressed": 2,
            "S4_Exhausted": 3, "S5_Failing": 4
        }
        states = adata.obs["workload_state"].map(state_map).values
        labels["state"] = states.astype(np.int64)
        print(f"  States: {np.bincount(labels['state'])}")

    return X, gene_names, labels, adata


def create_dataloaders(
    X: np.ndarray,
    labels: dict,
    batch_size: int = 128,
    val_split: float = 0.2
):
    """Create train and validation dataloaders."""
    # Convert to tensors
    X_tensor = torch.FloatTensor(X)

    tensors = [X_tensor]
    for key in ["cwi", "condition", "state"]:
        if key in labels:
            tensors.append(torch.tensor(labels[key]))
        else:
            # Placeholder if label not available
            tensors.append(torch.zeros(len(X)))

    dataset = TensorDataset(*tensors)

    # Split
    n_val = int(len(dataset) * val_split)
    n_train = len(dataset) - n_val
    train_dataset, val_dataset = random_split(
        dataset, [n_train, n_val],
        generator=torch.Generator().manual_seed(42)
    )

    train_loader = DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True, drop_last=True
    )
    val_loader = DataLoader(
        val_dataset, batch_size=batch_size, shuffle=False
    )

    print(f"Train: {len(train_dataset)}, Val: {len(val_dataset)}")

    return train_loader, val_loader


def train_epoch(
    model: LiteratureGuidedVAE,
    train_loader: DataLoader,
    criterion: MultiTaskLoss,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    phase: str = "pretrain"
) -> dict:
    """Train for one epoch."""
    model.train()
    total_loss = 0
    loss_components = {}

    for batch in train_loader:
        X = batch[0].to(device)
        cwi = batch[1].to(device) if phase == "finetune" else None
        condition = batch[2].to(device) if phase == "finetune" else None
        state = batch[3].to(device) if phase == "finetune" else None

        optimizer.zero_grad()

        outputs = model(X)

        if phase == "pretrain":
            # Only reconstruction + KL loss
            loss, losses = criterion(outputs, X)
        else:
            # Full multi-task loss
            loss, losses = criterion(
                outputs, X,
                condition_labels=condition,
                state_labels=state,
                cwi_labels=cwi
            )

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        total_loss += loss.item()
        for k, v in losses.items():
            loss_components[k] = loss_components.get(k, 0) + v

    n_batches = len(train_loader)
    avg_loss = total_loss / n_batches
    for k in loss_components:
        loss_components[k] /= n_batches

    return {"loss": avg_loss, **loss_components}


def validate(
    model: LiteratureGuidedVAE,
    val_loader: DataLoader,
    criterion: MultiTaskLoss,
    device: torch.device,
    phase: str = "pretrain"
) -> dict:
    """Validate model."""
    model.eval()
    total_loss = 0
    loss_components = {}

    with torch.no_grad():
        for batch in val_loader:
            X = batch[0].to(device)
            cwi = batch[1].to(device) if phase == "finetune" else None
            condition = batch[2].to(device) if phase == "finetune" else None
            state = batch[3].to(device) if phase == "finetune" else None

            outputs = model(X)

            if phase == "pretrain":
                loss, losses = criterion(outputs, X)
            else:
                loss, losses = criterion(
                    outputs, X,
                    condition_labels=condition,
                    state_labels=state,
                    cwi_labels=cwi
                )

            total_loss += loss.item()
            for k, v in losses.items():
                loss_components[k] = loss_components.get(k, 0) + v

    n_batches = len(val_loader)
    avg_loss = total_loss / n_batches
    for k in loss_components:
        loss_components[k] /= n_batches

    return {"loss": avg_loss, **loss_components}


def train_phase(
    model: LiteratureGuidedVAE,
    train_loader: DataLoader,
    val_loader: DataLoader,
    criterion: MultiTaskLoss,
    config: dict,
    device: torch.device,
    phase: str = "pretrain"
) -> dict:
    """Train model for one phase."""
    print(f"\n{'=' * 60}")
    print(f"Phase: {phase.upper()}")
    print(f"{'=' * 60}")

    phase_config = config[phase]
    epochs = phase_config["epochs"]
    lr = phase_config["learning_rate"]
    patience = phase_config["early_stopping"]

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=5
    )
    early_stopping = EarlyStopping(patience=patience)

    history = {"train": [], "val": []}
    best_val_loss = float('inf')
    best_state = None

    for epoch in range(epochs):
        # Train
        train_metrics = train_epoch(
            model, train_loader, criterion, optimizer, device, phase
        )
        history["train"].append(train_metrics)

        # Validate
        val_metrics = validate(model, val_loader, criterion, device, phase)
        history["val"].append(val_metrics)

        # Learning rate scheduling
        scheduler.step(val_metrics["loss"])

        # Save best model
        if val_metrics["loss"] < best_val_loss:
            best_val_loss = val_metrics["loss"]
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        # Print progress
        if (epoch + 1) % 5 == 0 or epoch == 0:
            print(f"Epoch {epoch + 1:3d}/{epochs} | "
                  f"Train: {train_metrics['loss']:.4f} | "
                  f"Val: {val_metrics['loss']:.4f} | "
                  f"LR: {optimizer.param_groups[0]['lr']:.2e}")

        # Early stopping
        if early_stopping(val_metrics["loss"]):
            print(f"Early stopping at epoch {epoch + 1}")
            break

    # Restore best model
    if best_state is not None:
        model.load_state_dict(best_state)

    print(f"Best validation loss: {best_val_loss:.4f}")

    return history


def compute_metrics(
    model: LiteratureGuidedVAE,
    val_loader: DataLoader,
    device: torch.device
) -> dict:
    """Compute evaluation metrics."""
    model.eval()

    all_cwi_pred = []
    all_cwi_true = []
    all_cond_pred = []
    all_cond_true = []
    all_state_pred = []
    all_state_true = []

    with torch.no_grad():
        for batch in val_loader:
            X = batch[0].to(device)
            cwi_true = batch[1]
            cond_true = batch[2]
            state_true = batch[3]

            outputs = model(X)

            all_cwi_pred.append(outputs["cwi_pred"].cpu())
            all_cwi_true.append(cwi_true)
            all_cond_pred.append(outputs["condition_pred"].cpu())
            all_cond_true.append(cond_true)
            all_state_pred.append(outputs["state_logits"].argmax(dim=1).cpu())
            all_state_true.append(state_true)

    # Concatenate
    cwi_pred = torch.cat(all_cwi_pred).numpy()
    cwi_true = torch.cat(all_cwi_true).numpy()
    cond_pred = torch.cat(all_cond_pred).numpy()
    cond_true = torch.cat(all_cond_true).numpy()
    state_pred = torch.cat(all_state_pred).numpy()
    state_true = torch.cat(all_state_true).numpy()

    # Compute metrics
    metrics = {}

    # CWI correlation
    cwi_corr = np.corrcoef(cwi_pred, cwi_true)[0, 1]
    cwi_mse = np.mean((cwi_pred - cwi_true) ** 2)
    metrics["cwi_correlation"] = float(cwi_corr)
    metrics["cwi_mse"] = float(cwi_mse)

    # Condition AUROC
    from sklearn.metrics import roc_auc_score, accuracy_score
    if len(np.unique(cond_true)) > 1:
        cond_auroc = roc_auc_score(cond_true, cond_pred)
        cond_acc = accuracy_score(cond_true, (cond_pred > 0.5).astype(int))
        metrics["condition_auroc"] = float(cond_auroc)
        metrics["condition_accuracy"] = float(cond_acc)

    # State accuracy
    state_acc = accuracy_score(state_true, state_pred)
    metrics["state_accuracy"] = float(state_acc)

    return metrics


def save_results(
    model: LiteratureGuidedVAE,
    adata,
    history: dict,
    metrics: dict,
    config: dict,
    device: torch.device
):
    """Save training results."""
    print("\n" + "=" * 60)
    print("Saving results")
    print("=" * 60)

    # 1. Save model weights
    model_path = RESULTS_DIR / "model_weights.pt"
    torch.save(model.state_dict(), model_path)
    print(f"Saved model weights to {model_path}")

    # 2. Extract CWI for all cells
    if hasattr(adata.X, 'toarray'):
        X = torch.FloatTensor(adata.X.toarray()).to(device)
    else:
        X = torch.FloatTensor(adata.X).to(device)

    model.eval()
    with torch.no_grad():
        cwi_scores = model.get_workload_index(X).cpu().numpy()

    cwi_df = pd.DataFrame({
        "cell_id": adata.obs.index,
        "cwi_predicted": cwi_scores,
        "cwi_literature": adata.obs.get("CWI_literature", np.nan),
        "workload_state": adata.obs.get("workload_state", "Unknown"),
        "condition": adata.obs.get("condition_binary", np.nan)
    })
    cwi_path = RESULTS_DIR / "cwi_scores.csv"
    cwi_df.to_csv(cwi_path, index=False)
    print(f"Saved CWI scores to {cwi_path}")

    # 3. Save state predictions
    with torch.no_grad():
        outputs = model(X)
        state_pred = outputs["state_logits"].argmax(dim=1).cpu().numpy()

    state_map = {0: "S1_Resting", 1: "S2_Active", 2: "S3_Stressed",
                 3: "S4_Exhausted", 4: "S5_Failing"}
    state_df = pd.DataFrame({
        "cell_id": adata.obs.index,
        "state_predicted": [state_map[s] for s in state_pred],
        "state_literature": adata.obs.get("workload_state", "Unknown")
    })
    state_path = RESULTS_DIR / "state_predictions.csv"
    state_df.to_csv(state_path, index=False)
    print(f"Saved state predictions to {state_path}")

    # 4. Save gene importance
    importance = model.get_gene_importance()
    importance_rows = []
    for module, genes in importance.items():
        for gene, weight in genes.items():
            importance_rows.append({
                "gene": gene,
                "module": module,
                "importance": weight
            })
    importance_df = pd.DataFrame(importance_rows)
    importance_df = importance_df.sort_values("importance", ascending=False)
    importance_path = RESULTS_DIR / "gene_importance.csv"
    importance_df.to_csv(importance_path, index=False)
    print(f"Saved gene importance to {importance_path}")

    # 5. Save training history
    history_path = RESULTS_DIR / "training_history.json"
    # Convert to serializable format
    history_serializable = {}
    for phase, data in history.items():
        history_serializable[phase] = {
            "train": data.get("train", []),
            "val": data.get("val", [])
        }
    with open(history_path, 'w') as f:
        json.dump(history_serializable, f, indent=2)
    print(f"Saved training history to {history_path}")

    # 6. Save validation metrics
    metrics["timestamp"] = datetime.now().isoformat()
    metrics["config"] = config
    metrics_path = RESULTS_DIR / "validation_metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved validation metrics to {metrics_path}")

    # 7. Print summary
    print("\n" + "=" * 60)
    print("TRAINING SUMMARY")
    print("=" * 60)
    print(f"CWI Correlation: {metrics.get('cwi_correlation', 'N/A'):.4f}")
    print(f"CWI MSE: {metrics.get('cwi_mse', 'N/A'):.4f}")
    print(f"Condition AUROC: {metrics.get('condition_auroc', 'N/A'):.4f}")
    print(f"Condition Accuracy: {metrics.get('condition_accuracy', 'N/A'):.4f}")
    print(f"State Accuracy: {metrics.get('state_accuracy', 'N/A'):.4f}")


def main():
    print("=" * 60)
    print("BETA-CELL WORKLOAD - DEEP LEARNING TRAINING")
    print("=" * 60)

    # Set device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Load data
    X, gene_names, labels, adata = load_training_data()

    # Create model
    print("\nCreating model...")
    model, criterion = create_model(gene_names, TRAINING_CONFIG)
    model = model.to(device)
    print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")

    # Full training history
    full_history = {}

    # Phase 1: Pre-training
    print("\n" + "=" * 60)
    print("PHASE 1: PRE-TRAINING (Unsupervised)")
    print("=" * 60)

    pretrain_loader, pretrain_val = create_dataloaders(
        X, labels,
        batch_size=TRAINING_CONFIG["pretrain"]["batch_size"]
    )

    pretrain_history = train_phase(
        model, pretrain_loader, pretrain_val, criterion,
        TRAINING_CONFIG, device, phase="pretrain"
    )
    full_history["pretrain"] = pretrain_history

    # Phase 2: Fine-tuning
    print("\n" + "=" * 60)
    print("PHASE 2: FINE-TUNING (Multi-task)")
    print("=" * 60)

    finetune_loader, finetune_val = create_dataloaders(
        X, labels,
        batch_size=TRAINING_CONFIG["finetune"]["batch_size"]
    )

    finetune_history = train_phase(
        model, finetune_loader, finetune_val, criterion,
        TRAINING_CONFIG, device, phase="finetune"
    )
    full_history["finetune"] = finetune_history

    # Compute final metrics
    print("\nComputing final metrics...")
    metrics = compute_metrics(model, finetune_val, device)

    # Save everything
    save_results(model, adata, full_history, metrics, TRAINING_CONFIG, device)

    print("\n" + "=" * 60)
    print("TRAINING COMPLETE")
    print("=" * 60)
    print(f"\nOutputs saved to: {RESULTS_DIR}")
    print("  - model_weights.pt")
    print("  - cwi_scores.csv")
    print("  - state_predictions.csv")
    print("  - gene_importance.csv")
    print("  - training_history.json")
    print("  - validation_metrics.json")


if __name__ == "__main__":
    main()
