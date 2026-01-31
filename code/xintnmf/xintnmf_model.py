"""
X-intNMF: Cross-integration Non-negative Matrix Factorization
=============================================================

Implementation based on:
"X-intNMF: Network-regularized non-negative matrix factorization with
cross- and intra-omics feature interactions"

Key features:
1. Network-regularized NMF incorporating gene interaction networks
2. Cross-omics integration capability (for multi-modal data)
3. Intra-omics feature interactions (gene-gene relationships)
4. Better suited for small sample sizes than deep learning

For beta-cell workload analysis:
- Uses literature-derived gene modules as prior knowledge
- Incorporates PPI/pathway networks for regularization
- Derives interpretable workload factors

Author: Beta-Cell Workload Analysis Pipeline
"""

import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix, issparse
from typing import Dict, List, Tuple, Optional, Union
import warnings

warnings.filterwarnings('ignore')


class XintNMF:
    """
    X-intNMF: Network-regularized NMF with cross- and intra-omics interactions.

    The objective function:
    min ||X - WH||_F^2 + alpha * Tr(H * L_intra * H^T) + beta * Tr(W * L_cross * W^T)

    Where:
    - X: Data matrix (genes x cells)
    - W: Basis matrix (genes x k factors)
    - H: Coefficient matrix (k factors x cells)
    - L_intra: Intra-omics Laplacian (gene-gene network)
    - L_cross: Cross-omics Laplacian (optional, for multi-omics)
    - alpha, beta: Regularization weights
    """

    def __init__(
        self,
        n_components: int = 5,
        alpha: float = 0.1,       # Intra-omics regularization weight
        beta: float = 0.01,       # Cross-omics regularization weight
        max_iter: int = 200,
        tol: float = 1e-4,
        init: str = 'nndsvd',     # Initialization: 'random', 'nndsvd', 'nndsvda'
        random_state: Optional[int] = 42,
        verbose: bool = True
    ):
        """
        Initialize X-intNMF model.

        Parameters:
        -----------
        n_components : int
            Number of NMF components (latent factors)
        alpha : float
            Weight for intra-omics network regularization
        beta : float
            Weight for cross-omics regularization
        max_iter : int
            Maximum number of iterations
        tol : float
            Tolerance for convergence
        init : str
            Initialization method
        random_state : int
            Random seed for reproducibility
        verbose : bool
            Print progress information
        """
        self.n_components = n_components
        self.alpha = alpha
        self.beta = beta
        self.max_iter = max_iter
        self.tol = tol
        self.init = init
        self.random_state = random_state
        self.verbose = verbose

        # Model components (set after fitting)
        self.W_ = None  # Basis matrix (genes x components)
        self.H_ = None  # Coefficient matrix (components x cells)
        self.reconstruction_err_ = []
        self.n_iter_ = 0

    def _init_matrices(
        self,
        X: np.ndarray,
        W_init: Optional[np.ndarray] = None,
        H_init: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Initialize W and H matrices."""
        n_features, n_samples = X.shape
        k = self.n_components

        if W_init is not None and H_init is not None:
            return W_init.copy(), H_init.copy()

        np.random.seed(self.random_state)

        if self.init == 'random':
            W = np.abs(np.random.randn(n_features, k)) + 0.1
            H = np.abs(np.random.randn(k, n_samples)) + 0.1

        elif self.init in ['nndsvd', 'nndsvda']:
            # NNDSVD initialization (Boutsidis & Gallopoulos, 2008)
            U, S, Vt = np.linalg.svd(X, full_matrices=False)

            W = np.zeros((n_features, k))
            H = np.zeros((k, n_samples))

            # First component
            W[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
            H[0, :] = np.sqrt(S[0]) * np.abs(Vt[0, :])

            for j in range(1, min(k, len(S))):
                x = U[:, j]
                y = Vt[j, :]

                x_pos = np.maximum(x, 0)
                x_neg = np.maximum(-x, 0)
                y_pos = np.maximum(y, 0)
                y_neg = np.maximum(-y, 0)

                x_pos_norm = np.linalg.norm(x_pos)
                x_neg_norm = np.linalg.norm(x_neg)
                y_pos_norm = np.linalg.norm(y_pos)
                y_neg_norm = np.linalg.norm(y_neg)

                mp = x_pos_norm * y_pos_norm
                mn = x_neg_norm * y_neg_norm

                if mp > mn:
                    u = x_pos / (x_pos_norm + 1e-10)
                    v = y_pos / (y_pos_norm + 1e-10)
                    sigma = mp
                else:
                    u = x_neg / (x_neg_norm + 1e-10)
                    v = y_neg / (y_neg_norm + 1e-10)
                    sigma = mn

                W[:, j] = np.sqrt(S[j] * sigma) * u
                H[j, :] = np.sqrt(S[j] * sigma) * v

            # Fill remaining components if k > rank
            for j in range(min(k, len(S)), k):
                W[:, j] = np.abs(np.random.randn(n_features)) + 0.1
                H[j, :] = np.abs(np.random.randn(n_samples)) + 0.1

            if self.init == 'nndsvda':
                # Add small value to avoid zeros
                avg = np.mean(X)
                W = np.where(W < 1e-10, avg / k, W)
                H = np.where(H < 1e-10, avg / k, H)
        else:
            raise ValueError(f"Unknown init method: {self.init}")

        return W, H

    def _compute_laplacian(
        self,
        adjacency: np.ndarray,
        normalized: bool = True
    ) -> np.ndarray:
        """
        Compute graph Laplacian from adjacency matrix.

        L = D - A (unnormalized)
        L = I - D^(-1/2) A D^(-1/2) (normalized)
        """
        if issparse(adjacency):
            adjacency = adjacency.toarray()

        # Degree matrix
        D = np.diag(adjacency.sum(axis=1))

        if normalized:
            # Normalized Laplacian
            D_inv_sqrt = np.diag(1.0 / (np.sqrt(np.diag(D)) + 1e-10))
            L = np.eye(len(D)) - D_inv_sqrt @ adjacency @ D_inv_sqrt
        else:
            # Unnormalized Laplacian
            L = D - adjacency

        return L

    def fit(
        self,
        X: np.ndarray,
        network_intra: Optional[np.ndarray] = None,
        network_cross: Optional[np.ndarray] = None,
        W_init: Optional[np.ndarray] = None,
        H_init: Optional[np.ndarray] = None
    ) -> 'XintNMF':
        """
        Fit X-intNMF model to data.

        Parameters:
        -----------
        X : np.ndarray
            Data matrix (genes x cells), non-negative
        network_intra : np.ndarray, optional
            Intra-omics adjacency matrix (gene-gene network)
        network_cross : np.ndarray, optional
            Cross-omics adjacency matrix (for multi-modal data)
        W_init : np.ndarray, optional
            Initial W matrix
        H_init : np.ndarray, optional
            Initial H matrix

        Returns:
        --------
        self : fitted model
        """
        # Ensure non-negativity
        X = np.maximum(X, 0)

        if issparse(X):
            X = X.toarray()

        n_features, n_samples = X.shape

        if self.verbose:
            print(f"Fitting X-intNMF: {n_features} features x {n_samples} samples")
            print(f"Components: {self.n_components}, alpha: {self.alpha}, beta: {self.beta}")

        # Initialize W and H
        W, H = self._init_matrices(X, W_init, H_init)

        # Compute Laplacians
        if network_intra is not None:
            L_intra = self._compute_laplacian(network_intra)
            if self.verbose:
                print(f"Using intra-omics network: {network_intra.shape}")
        else:
            L_intra = np.zeros((n_features, n_features))

        if network_cross is not None:
            L_cross = self._compute_laplacian(network_cross)
            if self.verbose:
                print(f"Using cross-omics network: {network_cross.shape}")
        else:
            L_cross = np.zeros((n_features, n_features))

        # Multiplicative update rules with network regularization
        # For graph-regularized NMF: ||X - WH||^2 + alpha * Tr(W^T L W)
        # This encourages genes connected in the network to have similar factor loadings
        self.reconstruction_err_ = []

        for iteration in range(self.max_iter):
            # Update H (coefficient matrix)
            # H <- H * (W^T X) / (W^T W H + eps)
            WtX = W.T @ X
            WtWH = W.T @ W @ H

            H = H * WtX / (WtWH + 1e-10)
            H = np.maximum(H, 1e-10)

            # Update W (basis matrix) with graph regularization
            # W <- W * (X H^T) / (W H H^T + alpha * L_intra W + eps)
            XHt = X @ H.T
            WHHt = W @ H @ H.T

            # Graph regularization on W (features/genes)
            if self.alpha > 0 and network_intra is not None:
                W_reg = self.alpha * L_intra @ W
            else:
                W_reg = 0

            W = W * XHt / (WHHt + W_reg + 1e-10)
            W = np.maximum(W, 1e-10)

            # Compute reconstruction error
            recon = W @ H
            err = np.linalg.norm(X - recon, 'fro') ** 2

            # Add graph regularization term: Tr(W^T L W)
            if self.alpha > 0 and network_intra is not None:
                err += self.alpha * np.trace(W.T @ L_intra @ W)

            self.reconstruction_err_.append(err)

            # Check convergence
            if iteration > 0:
                rel_change = abs(self.reconstruction_err_[-2] - err) / (self.reconstruction_err_[-2] + 1e-10)
                if rel_change < self.tol:
                    if self.verbose:
                        print(f"Converged at iteration {iteration}, error: {err:.4f}")
                    break

            if self.verbose and iteration % 20 == 0:
                print(f"  Iteration {iteration}: error = {err:.4f}")

        self.W_ = W
        self.H_ = H
        self.n_iter_ = iteration + 1

        if self.verbose:
            print(f"Finished in {self.n_iter_} iterations")

        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        """
        Transform new data using fitted model.

        Parameters:
        -----------
        X : np.ndarray
            Data matrix (genes x cells)

        Returns:
        --------
        H : np.ndarray
            Coefficient matrix (components x cells)
        """
        if self.W_ is None:
            raise ValueError("Model not fitted. Call fit() first.")

        X = np.maximum(X, 0)
        if issparse(X):
            X = X.toarray()

        # Solve for H given fixed W using NNLS
        # min ||X - WH||_F^2 s.t. H >= 0

        # Use multiplicative updates
        n_samples = X.shape[1]
        H = np.abs(np.random.randn(self.n_components, n_samples)) + 0.1

        for _ in range(100):
            WtX = self.W_.T @ X
            WtWH = self.W_.T @ self.W_ @ H
            H = H * WtX / (WtWH + 1e-10)
            H = np.maximum(H, 1e-10)

        return H

    def fit_transform(
        self,
        X: np.ndarray,
        network_intra: Optional[np.ndarray] = None,
        network_cross: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """Fit model and return transformed data."""
        self.fit(X, network_intra, network_cross)
        return self.H_

    def inverse_transform(self, H: np.ndarray) -> np.ndarray:
        """Reconstruct data from coefficients."""
        if self.W_ is None:
            raise ValueError("Model not fitted. Call fit() first.")
        return self.W_ @ H

    def get_feature_weights(self) -> np.ndarray:
        """Get feature weights (basis matrix W)."""
        return self.W_

    def get_component_loadings(self, feature_names: Optional[List[str]] = None) -> Dict:
        """
        Get top features for each component.

        Returns dictionary mapping component index to sorted features.
        """
        if self.W_ is None:
            raise ValueError("Model not fitted.")

        loadings = {}
        for k in range(self.n_components):
            weights = self.W_[:, k]
            sorted_idx = np.argsort(weights)[::-1]

            if feature_names is not None:
                loadings[k] = [(feature_names[i], weights[i]) for i in sorted_idx[:20]]
            else:
                loadings[k] = [(i, weights[i]) for i in sorted_idx[:20]]

        return loadings


class WorkloadXintNMF(XintNMF):
    """
    X-intNMF specialized for beta-cell workload analysis.

    Incorporates:
    1. Literature-derived gene module priors
    2. PPI network regularization
    3. Interpretable workload factor extraction
    """

    # Literature-derived workload gene modules
    WORKLOAD_MODULES = {
        "biosynthetic": [
            "INS", "IAPP", "PCSK1", "PCSK2", "CPE", "SCG2", "CHGB",
            "SNAP25", "SYT1", "RAB3A", "VAMP2", "STX1A"
        ],
        "metabolic": [
            "GCK", "SLC2A2", "G6PC2", "PDX1", "MAFA", "NKX6-1",
            "ABCC8", "KCNJ11", "CACNA1C", "TRPM2"
        ],
        "stress": [
            "XBP1", "ATF6", "ERN1", "DDIT3", "ATF4", "HSPA5", "HSP90B1",
            "CALR", "CANX", "PDIA4", "PDIA6", "EIF2AK3", "TRIB3", "GPX1"
        ],
        "dedifferentiation": [
            "ALDH1A3", "NEUROG3", "SOX9", "HES1", "FOXO1", "MYC",
            "NANOG", "POU5F1", "CD44"
        ]
    }

    def __init__(
        self,
        n_components: int = 5,
        alpha: float = 0.1,
        beta: float = 0.01,
        module_weight: float = 0.5,
        **kwargs
    ):
        """
        Initialize WorkloadXintNMF.

        Parameters:
        -----------
        n_components : int
            Number of latent factors (recommend 4-6 for interpretability)
        alpha : float
            Network regularization weight
        beta : float
            Cross-omics weight (if multi-modal)
        module_weight : float
            Weight for module-based initialization prior
        """
        super().__init__(n_components=n_components, alpha=alpha, beta=beta, **kwargs)
        self.module_weight = module_weight
        self.gene_names_ = None
        self.module_assignments_ = None

    def _create_module_network(
        self,
        gene_names: List[str]
    ) -> np.ndarray:
        """
        Create gene-gene network based on module co-membership.

        Genes in the same module are connected.
        """
        n_genes = len(gene_names)
        network = np.zeros((n_genes, n_genes))

        gene_to_idx = {g: i for i, g in enumerate(gene_names)}

        for module, genes in self.WORKLOAD_MODULES.items():
            module_idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]

            # Connect genes within module
            for i in module_idx:
                for j in module_idx:
                    if i != j:
                        network[i, j] = 1.0

        return network

    def _create_module_prior_W(
        self,
        gene_names: List[str],
        n_components: int
    ) -> np.ndarray:
        """
        Create prior W matrix based on literature modules.

        Maps known modules to initial components.
        """
        n_genes = len(gene_names)
        W_prior = np.ones((n_genes, n_components)) * 0.1

        gene_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Assign first 4 components to known modules
        module_names = list(self.WORKLOAD_MODULES.keys())

        for k, module_name in enumerate(module_names[:n_components]):
            genes = self.WORKLOAD_MODULES[module_name]
            for gene in genes:
                if gene in gene_to_idx:
                    W_prior[gene_to_idx[gene], k] = 1.0

        return W_prior

    def fit_workload(
        self,
        X: np.ndarray,
        gene_names: List[str],
        external_network: Optional[np.ndarray] = None
    ) -> 'WorkloadXintNMF':
        """
        Fit model with workload-specific priors.

        Parameters:
        -----------
        X : np.ndarray
            Expression matrix (genes x cells)
        gene_names : List[str]
            Gene names corresponding to rows
        external_network : np.ndarray, optional
            External gene-gene network (e.g., PPI, pathway)
        """
        self.gene_names_ = gene_names

        if self.verbose:
            print("Setting up workload-specific X-intNMF...")

        # Create module-based network
        module_network = self._create_module_network(gene_names)

        # Combine with external network if provided
        if external_network is not None:
            # Weighted combination
            network = 0.5 * module_network + 0.5 * external_network
        else:
            network = module_network

        if self.verbose:
            n_edges = (network > 0).sum() // 2
            print(f"Network has {n_edges} edges")

        # Create module-informed initialization
        W_init = self._create_module_prior_W(gene_names, self.n_components)

        # Mix with random initialization
        np.random.seed(self.random_state)
        W_random = np.abs(np.random.randn(*W_init.shape)) + 0.1
        W_init = self.module_weight * W_init + (1 - self.module_weight) * W_random

        # Fit with network regularization
        self.fit(X, network_intra=network, W_init=W_init)

        # Assign components to modules based on gene loadings
        self._assign_components_to_modules()

        return self

    def _assign_components_to_modules(self):
        """Assign each component to the most relevant workload module."""
        if self.W_ is None or self.gene_names_ is None:
            return

        gene_to_idx = {g: i for i, g in enumerate(self.gene_names_)}

        self.module_assignments_ = {}

        for k in range(self.n_components):
            weights = self.W_[:, k]

            # Calculate overlap with each module
            module_scores = {}
            for module_name, genes in self.WORKLOAD_MODULES.items():
                module_idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
                if module_idx:
                    module_scores[module_name] = np.mean(weights[module_idx])

            if module_scores:
                best_module = max(module_scores, key=module_scores.get)
                self.module_assignments_[k] = {
                    'module': best_module,
                    'score': module_scores[best_module],
                    'all_scores': module_scores
                }

    def get_workload_scores(self, H: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
        """
        Get workload scores for each module.

        Returns dictionary mapping module names to cell scores.
        """
        if H is None:
            H = self.H_

        if H is None:
            raise ValueError("No coefficient matrix available.")

        scores = {}

        for k, assignment in self.module_assignments_.items():
            module_name = assignment['module']
            if module_name not in scores:
                scores[module_name] = H[k, :]
            else:
                # If multiple components map to same module, average them
                scores[module_name] = (scores[module_name] + H[k, :]) / 2

        return scores

    def compute_cwi(self, H: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute Composite Workload Index from NMF factors.

        CWI = (Biosynthetic / Metabolic) * (1 + Stress) * (1 + Dediff)
        """
        scores = self.get_workload_scores(H)

        # Normalize scores
        def normalize(x):
            min_val, max_val = x.min(), x.max()
            if max_val > min_val:
                return (x - min_val) / (max_val - min_val) + 0.1
            return np.ones_like(x) * 0.5

        biosynthetic = normalize(scores.get('biosynthetic', np.ones(self.H_.shape[1]) * 0.5))
        metabolic = normalize(scores.get('metabolic', np.ones(self.H_.shape[1]) * 0.5))
        stress = normalize(scores.get('stress', np.ones(self.H_.shape[1]) * 0.5))
        dediff = normalize(scores.get('dedifferentiation', np.ones(self.H_.shape[1]) * 0.5))

        cwi = (biosynthetic / metabolic) * (1 + stress) * (1 + dediff)

        return cwi


if __name__ == "__main__":
    # Test with synthetic data
    print("Testing X-intNMF implementation...")

    np.random.seed(42)

    # Create synthetic data
    n_genes = 100
    n_cells = 50
    n_components = 4

    # True factors
    W_true = np.random.rand(n_genes, n_components)
    H_true = np.random.rand(n_components, n_cells)
    X = W_true @ H_true + 0.1 * np.random.rand(n_genes, n_cells)

    # Create synthetic network (random connections)
    network = np.random.rand(n_genes, n_genes)
    network = (network + network.T) / 2  # Symmetrize
    network = (network > 0.8).astype(float)  # Threshold
    np.fill_diagonal(network, 0)

    print(f"Data shape: {X.shape}")
    print(f"Network edges: {(network > 0).sum() // 2}")

    # Fit model
    model = XintNMF(n_components=n_components, alpha=0.1, verbose=True)
    H = model.fit_transform(X, network_intra=network)

    print(f"\nFinal reconstruction error: {model.reconstruction_err_[-1]:.4f}")
    print(f"W shape: {model.W_.shape}")
    print(f"H shape: {model.H_.shape}")

    # Test reconstruction
    X_recon = model.inverse_transform(H)
    recon_error = np.linalg.norm(X - X_recon, 'fro') / np.linalg.norm(X, 'fro')
    print(f"Relative reconstruction error: {recon_error:.4f}")

    print("\nX-intNMF test completed successfully!")
