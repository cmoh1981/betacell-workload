#!/usr/bin/env python3
"""
02_model.py - Literature-Guided VAE for Beta-Cell Workload Index

Implements the LG-VAE architecture with:
- Literature-guided feature selection with attention
- Multi-task learning (reconstruction, condition, state)
- Dedicated workload dimensions in latent space
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
import numpy as np
from typing import Dict, List, Tuple, Optional

# Literature gene modules (from architecture document)
FEATURE_MODULES = {
    "biosynthetic": [
        "INS", "IAPP", "PCSK1", "PCSK2", "CPE", "CHGB", "SCG2",
        "HSPA5", "HSP90B1", "PDIA4", "PDIA6", "CALR", "CANX",
        "SLC30A8", "SNAP25"
    ],
    "metabolic": [
        "GCK", "SLC2A2", "G6PC2", "PFKFB2",
        "TFAM", "PPARGC1A", "HADH",
        "PDX1", "MAFA", "NKX6-1", "UCN3", "NEUROD1", "NKX2-2", "PAX6",
        "PPARA", "PPARD", "HNF4A"
    ],
    "stress": [
        "XBP1", "ATF6", "ERN1", "EIF2AK3",
        "DDIT3", "ATF4", "TRIB3", "BBC3",
        "NFE2L2", "SOD1", "SOD2", "GPX1", "CAT",
        "NFKB1", "TNF", "IL1B", "CCL2"
    ],
    "dedifferentiation": [
        "ALDH1A3", "NEUROG3", "SOX9", "HES1", "GASTRIN",
        "LDHA", "HK1"
    ]
}


class LiteratureAttention(nn.Module):
    """
    Attention mechanism for literature-guided gene modules.
    Learns importance weights for each gene module.
    """

    def __init__(self, gene_names: List[str], hidden_dim: int = 64):
        super().__init__()
        self.gene_names = gene_names
        self.gene_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Create module indices
        self.module_indices = {}
        self.module_names = list(FEATURE_MODULES.keys())

        for module_name, genes in FEATURE_MODULES.items():
            indices = [self.gene_to_idx[g] for g in genes if g in self.gene_to_idx]
            if indices:
                self.module_indices[module_name] = torch.tensor(indices)

        # Attention parameters per module
        n_modules = len(self.module_indices)
        self.module_attention = nn.Sequential(
            nn.Linear(hidden_dim, n_modules),
            nn.Softmax(dim=-1)
        )

        # Gene-level attention within modules
        self.gene_attention = nn.ParameterDict({
            name: nn.Parameter(torch.ones(len(idx)) / len(idx))
            for name, idx in self.module_indices.items()
        })

        # Projection for computing attention context
        self.context_projection = nn.Linear(len(gene_names), hidden_dim)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
        """
        Apply literature-guided attention.

        Args:
            x: Gene expression tensor [batch, n_genes]

        Returns:
            attended_features: Weighted features [batch, n_modules]
            attention_weights: Dict of attention weights per module
        """
        batch_size = x.shape[0]
        device = x.device

        # Compute context for module attention
        context = self.context_projection(x)
        module_weights = self.module_attention(context)  # [batch, n_modules]

        # Compute module-level features with gene attention
        module_features = []
        attention_weights = {}

        for i, (name, idx) in enumerate(self.module_indices.items()):
            idx = idx.to(device)
            gene_expr = x[:, idx]  # [batch, n_genes_in_module]

            # Apply gene-level attention (softmax normalized)
            gene_attn = F.softmax(self.gene_attention[name], dim=0)
            attention_weights[name] = gene_attn

            # Weighted average of gene expressions
            module_feat = (gene_expr * gene_attn.unsqueeze(0)).sum(dim=1)
            module_features.append(module_feat)

        # Stack module features
        module_features = torch.stack(module_features, dim=1)  # [batch, n_modules]

        # Apply module-level attention
        attended = module_features * module_weights

        return attended, attention_weights


class Encoder(nn.Module):
    """VAE Encoder with literature-guided attention."""

    def __init__(
        self,
        input_dim: int,
        gene_names: List[str],
        hidden_dims: List[int] = [256, 128, 64],
        latent_dim: int = 32,
        dropout: float = 0.2
    ):
        super().__init__()

        self.input_dim = input_dim
        self.latent_dim = latent_dim

        # Literature attention
        self.literature_attention = LiteratureAttention(gene_names, hidden_dim=64)
        n_modules = len(self.literature_attention.module_indices)

        # Input processing: combine raw features with attended module features
        combined_dim = input_dim + n_modules

        # Encoder layers
        layers = []
        prev_dim = combined_dim
        for h_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, h_dim),
                nn.BatchNorm1d(h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            prev_dim = h_dim

        self.encoder = nn.Sequential(*layers)

        # Latent space parameters
        self.fc_mu = nn.Linear(hidden_dims[-1], latent_dim)
        self.fc_var = nn.Linear(hidden_dims[-1], latent_dim)

    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, Dict]:
        """
        Encode input to latent distribution.

        Returns:
            mu: Mean of latent distribution
            log_var: Log variance of latent distribution
            attention_info: Dictionary with attention weights
        """
        # Apply literature attention
        attended_features, attention_weights = self.literature_attention(x)

        # Combine raw and attended features
        combined = torch.cat([x, attended_features], dim=1)

        # Encode
        h = self.encoder(combined)

        # Latent parameters
        mu = self.fc_mu(h)
        log_var = self.fc_var(h)

        attention_info = {
            "module_features": attended_features,
            "gene_attention": attention_weights
        }

        return mu, log_var, attention_info


class Decoder(nn.Module):
    """VAE Decoder for gene expression reconstruction."""

    def __init__(
        self,
        latent_dim: int = 32,
        hidden_dims: List[int] = [64, 128, 256],
        output_dim: int = 500,
        dropout: float = 0.2
    ):
        super().__init__()

        layers = []
        prev_dim = latent_dim
        for h_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, h_dim),
                nn.BatchNorm1d(h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            prev_dim = h_dim

        self.decoder = nn.Sequential(*layers)
        self.output_layer = nn.Linear(hidden_dims[-1], output_dim)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent to gene expression."""
        h = self.decoder(z)
        return self.output_layer(h)


class ConditionHead(nn.Module):
    """Classification head for T2D vs Normal condition."""

    def __init__(self, latent_dim: int = 32, hidden_dim: int = 32):
        super().__init__()
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid()
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.classifier(z).squeeze(-1)


class StateHead(nn.Module):
    """Classification head for 5 workload states."""

    def __init__(self, latent_dim: int = 32, hidden_dim: int = 32, n_states: int = 5):
        super().__init__()
        self.classifier = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(hidden_dim, n_states)
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.classifier(z)


class WorkloadHead(nn.Module):
    """Regression head for continuous CWI prediction."""

    def __init__(self, workload_dims: int = 4):
        super().__init__()
        # Learnable weights for workload dimensions
        self.weights = nn.Parameter(torch.ones(workload_dims) / workload_dims)

    def forward(self, z_workload: torch.Tensor) -> torch.Tensor:
        """
        Compute CWI from workload latent dimensions.

        Args:
            z_workload: First 4 dimensions of latent space [batch, 4]
        """
        weights = F.softmax(self.weights, dim=0)
        cwi = (z_workload * weights).sum(dim=1)
        return cwi


class LiteratureGuidedVAE(nn.Module):
    """
    Literature-Guided Variational Autoencoder for Beta-Cell Workload.

    Architecture:
    - Encoder with literature-guided attention
    - Latent space with dedicated workload dimensions
    - Multi-task decoder heads
    """

    def __init__(
        self,
        input_dim: int,
        gene_names: List[str],
        latent_dim: int = 32,
        workload_dims: int = 4,
        encoder_dims: List[int] = [256, 128, 64],
        decoder_dims: List[int] = [64, 128, 256],
        dropout: float = 0.2,
        n_states: int = 5
    ):
        super().__init__()

        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.workload_dims = workload_dims
        self.gene_names = gene_names

        # Encoder
        self.encoder = Encoder(
            input_dim=input_dim,
            gene_names=gene_names,
            hidden_dims=encoder_dims,
            latent_dim=latent_dim,
            dropout=dropout
        )

        # Decoder
        self.decoder = Decoder(
            latent_dim=latent_dim,
            hidden_dims=decoder_dims,
            output_dim=input_dim,
            dropout=dropout
        )

        # Task heads
        self.condition_head = ConditionHead(latent_dim)
        self.state_head = StateHead(latent_dim, n_states=n_states)
        self.workload_head = WorkloadHead(workload_dims)

    def reparameterize(self, mu: torch.Tensor, log_var: torch.Tensor) -> torch.Tensor:
        """Reparameterization trick for VAE."""
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std

    def encode(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, Dict]:
        """Encode input to latent space."""
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent to reconstruction."""
        return self.decoder(z)

    def forward(self, x: torch.Tensor) -> Dict[str, torch.Tensor]:
        """
        Forward pass through full model.

        Returns dict with:
            - recon: Reconstructed gene expression
            - mu, log_var: Latent distribution parameters
            - z: Sampled latent vector
            - condition_pred: T2D probability
            - state_logits: Workload state logits
            - cwi_pred: Predicted CWI
            - attention_info: Attention weights
        """
        # Encode
        mu, log_var, attention_info = self.encode(x)

        # Sample latent
        z = self.reparameterize(mu, log_var)

        # Decode
        recon = self.decode(z)

        # Task predictions
        condition_pred = self.condition_head(z)
        state_logits = self.state_head(z)

        # CWI from workload dimensions
        z_workload = z[:, :self.workload_dims]
        cwi_pred = self.workload_head(z_workload)

        return {
            "recon": recon,
            "mu": mu,
            "log_var": log_var,
            "z": z,
            "condition_pred": condition_pred,
            "state_logits": state_logits,
            "cwi_pred": cwi_pred,
            "attention_info": attention_info
        }

    def get_workload_index(self, x: torch.Tensor) -> torch.Tensor:
        """Extract workload index for new samples."""
        with torch.no_grad():
            mu, _, _ = self.encode(x)
            z_workload = mu[:, :self.workload_dims]
            cwi = self.workload_head(z_workload)
            # Normalize to 0-5 range
            cwi_norm = torch.sigmoid(cwi) * 5
            return cwi_norm

    def get_gene_importance(self) -> Dict[str, np.ndarray]:
        """Extract gene importance from attention weights."""
        importance = {}

        for module_name, attn_param in self.encoder.literature_attention.gene_attention.items():
            weights = F.softmax(attn_param, dim=0).detach().cpu().numpy()
            idx = self.encoder.literature_attention.module_indices[module_name].cpu().numpy()
            genes = [self.gene_names[i] for i in idx]
            importance[module_name] = dict(zip(genes, weights))

        return importance


class MultiTaskLoss(nn.Module):
    """
    Multi-task loss for LG-VAE training.

    Loss = λ₁×Recon + λ₂×KL + λ₃×Condition + λ₄×State + λ₅×Literature
    """

    def __init__(
        self,
        lambda_recon: float = 1.0,
        lambda_kl: float = 0.1,
        lambda_condition: float = 1.0,
        lambda_state: float = 0.5,
        lambda_literature: float = 0.3,
        lambda_cwi: float = 0.5
    ):
        super().__init__()

        self.lambda_recon = lambda_recon
        self.lambda_kl = lambda_kl
        self.lambda_condition = lambda_condition
        self.lambda_state = lambda_state
        self.lambda_literature = lambda_literature
        self.lambda_cwi = lambda_cwi

        self.bce = nn.BCELoss()
        self.ce = nn.CrossEntropyLoss()
        self.mse = nn.MSELoss()

    def forward(
        self,
        outputs: Dict[str, torch.Tensor],
        x: torch.Tensor,
        condition_labels: Optional[torch.Tensor] = None,
        state_labels: Optional[torch.Tensor] = None,
        cwi_labels: Optional[torch.Tensor] = None,
        literature_targets: Optional[Dict[str, torch.Tensor]] = None
    ) -> Tuple[torch.Tensor, Dict[str, float]]:
        """
        Compute multi-task loss.

        Args:
            outputs: Model outputs dict
            x: Input gene expression
            condition_labels: Binary T2D labels
            state_labels: Workload state labels (0-4)
            cwi_labels: Continuous CWI targets
            literature_targets: Expected module scores from literature

        Returns:
            total_loss: Combined loss
            loss_dict: Individual loss components
        """
        losses = {}

        # 1. Reconstruction loss (MSE)
        recon_loss = self.mse(outputs["recon"], x)
        losses["recon"] = recon_loss.item()
        total = self.lambda_recon * recon_loss

        # 2. KL divergence
        kl_loss = -0.5 * torch.mean(
            1 + outputs["log_var"] - outputs["mu"].pow(2) - outputs["log_var"].exp()
        )
        losses["kl"] = kl_loss.item()
        total = total + self.lambda_kl * kl_loss

        # 3. Condition prediction loss (if labels available)
        if condition_labels is not None:
            condition_loss = self.bce(outputs["condition_pred"], condition_labels.float())
            losses["condition"] = condition_loss.item()
            total = total + self.lambda_condition * condition_loss

        # 4. State classification loss (if labels available)
        if state_labels is not None:
            state_loss = self.ce(outputs["state_logits"], state_labels.long())
            losses["state"] = state_loss.item()
            total = total + self.lambda_state * state_loss

        # 5. CWI regression loss (if labels available)
        if cwi_labels is not None:
            cwi_loss = self.mse(outputs["cwi_pred"], cwi_labels)
            losses["cwi"] = cwi_loss.item()
            total = total + self.lambda_cwi * cwi_loss

        # 6. Literature alignment loss
        if literature_targets is not None:
            module_features = outputs["attention_info"]["module_features"]
            lit_loss = 0
            for i, (name, target) in enumerate(literature_targets.items()):
                lit_loss = lit_loss + self.mse(module_features[:, i], target)
            lit_loss = lit_loss / len(literature_targets)
            losses["literature"] = lit_loss.item()
            total = total + self.lambda_literature * lit_loss

        losses["total"] = total.item()

        return total, losses


def create_model(
    gene_names: List[str],
    config: Optional[Dict] = None
) -> Tuple[LiteratureGuidedVAE, MultiTaskLoss]:
    """
    Factory function to create model and loss.

    Args:
        gene_names: List of gene names in input data
        config: Optional configuration dict

    Returns:
        model: LiteratureGuidedVAE instance
        criterion: MultiTaskLoss instance
    """
    if config is None:
        config = {
            "latent_dim": 32,
            "workload_dims": 4,
            "encoder_dims": [256, 128, 64],
            "decoder_dims": [64, 128, 256],
            "dropout": 0.2,
            "n_states": 5,
            "lambda_recon": 1.0,
            "lambda_kl": 0.1,
            "lambda_condition": 1.0,
            "lambda_state": 0.5,
            "lambda_literature": 0.3,
            "lambda_cwi": 0.5
        }

    input_dim = len(gene_names)

    model = LiteratureGuidedVAE(
        input_dim=input_dim,
        gene_names=gene_names,
        latent_dim=config.get("latent_dim", 32),
        workload_dims=config.get("workload_dims", 4),
        encoder_dims=config.get("encoder_dims", [256, 128, 64]),
        decoder_dims=config.get("decoder_dims", [64, 128, 256]),
        dropout=config.get("dropout", 0.2),
        n_states=config.get("n_states", 5)
    )

    criterion = MultiTaskLoss(
        lambda_recon=config.get("lambda_recon", 1.0),
        lambda_kl=config.get("lambda_kl", 0.1),
        lambda_condition=config.get("lambda_condition", 1.0),
        lambda_state=config.get("lambda_state", 0.5),
        lambda_literature=config.get("lambda_literature", 0.3),
        lambda_cwi=config.get("lambda_cwi", 0.5)
    )

    return model, criterion


if __name__ == "__main__":
    # Test model creation
    print("Testing LiteratureGuidedVAE...")

    # Simulate gene names
    test_genes = list(set(
        FEATURE_MODULES["biosynthetic"] +
        FEATURE_MODULES["metabolic"] +
        FEATURE_MODULES["stress"] +
        FEATURE_MODULES["dedifferentiation"]
    ))
    # Add some random genes to simulate full input
    test_genes.extend([f"GENE_{i}" for i in range(500 - len(test_genes))])

    print(f"Input genes: {len(test_genes)}")

    # Create model
    model, criterion = create_model(test_genes)
    print(f"Model created with {sum(p.numel() for p in model.parameters()):,} parameters")

    # Test forward pass
    batch_size = 32
    x = torch.randn(batch_size, len(test_genes))

    outputs = model(x)
    print(f"\nOutputs:")
    print(f"  Reconstruction: {outputs['recon'].shape}")
    print(f"  Latent z: {outputs['z'].shape}")
    print(f"  Condition pred: {outputs['condition_pred'].shape}")
    print(f"  State logits: {outputs['state_logits'].shape}")
    print(f"  CWI pred: {outputs['cwi_pred'].shape}")

    # Test loss computation
    condition_labels = torch.randint(0, 2, (batch_size,))
    state_labels = torch.randint(0, 5, (batch_size,))
    cwi_labels = torch.rand(batch_size) * 5

    loss, loss_dict = criterion(
        outputs, x,
        condition_labels=condition_labels,
        state_labels=state_labels,
        cwi_labels=cwi_labels
    )

    print(f"\nLosses:")
    for name, val in loss_dict.items():
        print(f"  {name}: {val:.4f}")

    # Test gene importance extraction
    importance = model.get_gene_importance()
    print(f"\nGene importance (top genes per module):")
    for module, genes in importance.items():
        top_genes = sorted(genes.items(), key=lambda x: x[1], reverse=True)[:3]
        print(f"  {module}: {top_genes}")

    print("\nModel test passed!")
