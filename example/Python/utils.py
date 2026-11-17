"""Utility functions for Bayesian clustering analysis in Python.

This module provides data generation, spatial adjacency matrix retrieval,
hyperparameter estimation for the Natarajan likelihood, and posterior
point estimation matching the BNPClust workflow.
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def generate_mixture_data(
    N: int = 100,
    K: int = 10,
    alpha: float = 10.0,
    dim: Optional[int] = None,
    radius: float = 1.0,
    sigma: float = 0.1,
    random_state: int = 44,
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate synthetic clustering mixture data on vertices of a simplex.

    Parameters
    ----------
    N : int
        Total number of observations.
    K : int
        Number of clusters.
    alpha : float
        Symmetric Dirichlet concentration parameter for cluster weights.
    dim : int, optional
        Data dimensionality (defaults to K). Must be >= K.
    radius : float
        Simplex vertex radius from origin.
    sigma : float
        Standard deviation of isotropic Gaussian noise within each cluster.
    random_state : int
        RNG seed for reproducibility.

    Returns
    -------
    points : np.ndarray of shape (N, dim)
        Generated data points sorted by cluster.
    clusts : np.ndarray of shape (N,)
        Ground-truth cluster assignments (0 to K-1).
    """
    if dim is None:
        dim = K
    if dim < K:
        raise ValueError(f"dim ({dim}) must be >= K ({K})")

    rng = np.random.default_rng(random_state)
    probs = rng.dirichlet(np.full(K, alpha))
    clusts = rng.choice(K, size=N, replace=True, p=probs)

    # Sort observations by cluster assignment
    order_idx = np.argsort(clusts)
    clusts = clusts[order_idx]

    # Cluster centres as vertices of dim-dimensional simplex
    clust_centres = np.zeros((K, dim), dtype=np.float64)
    for i in range(K):
        clust_centres[i, i] = radius

    points = np.zeros((N, dim), dtype=np.float64)
    for i in range(N):
        c = clusts[i]
        points[i] = rng.multivariate_normal(
            clust_centres[c], (sigma**2) * np.eye(dim)
        )

    return points, clusts


def retrieve_W(distance_matrix: np.ndarray, neighbours: int = 8) -> np.ndarray:
    """Construct a symmetric k-nearest neighbors adjacency matrix W.

    Parameters
    ----------
    distance_matrix : np.ndarray of shape (N, N)
        Pairwise distance matrix.
    neighbours : int
        Number of nearest neighbors to connect for each observation.

    Returns
    -------
    W : np.ndarray of shape (N, N), int32
        Symmetric binary adjacency matrix.
    """
    n = distance_matrix.shape[0]
    W = np.zeros((n, n), dtype=np.int32)
    for i in range(n):
        # Indices of nearest neighbors excluding self
        nn_indices = np.argsort(distance_matrix[i])[1 : neighbours + 1]
        W[i, nn_indices] = 1

    # Symmetrize
    W = ((W + W.T) > 0).astype(np.int32)
    return np.ascontiguousarray(W)


def classical_mds(dist_matrix: np.ndarray, n_components: int = 2) -> np.ndarray:
    """Compute classical Multidimensional Scaling (Principal Coordinate Analysis)."""
    n = dist_matrix.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (dist_matrix**2) @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    coords = eigvecs[:, :n_components] * np.sqrt(
        np.maximum(eigvals[:n_components], 0)
    )
    return coords


def set_hyperparameters(
    dist_matrix: np.ndarray,
    k_elbow: int = 3,
    random_state: int = 44,
) -> Dict[str, Any]:
    """Estimate Natarajan likelihood hyperparameters using k-means initialization.

    Fits Gamma distributions to within-cluster and between-cluster pairwise
    distances to initialize delta1, alpha, beta (cohesion) and delta2, gamma,
    zeta (repulsion).

    Parameters
    ----------
    dist_matrix : np.ndarray of shape (N, N)
        Pairwise distance matrix.
    k_elbow : int
        Initial number of clusters for k-means.
    random_state : int
        RNG seed for k-means.

    Returns
    -------
    dict containing:
        - delta1, alpha, beta : cohesion hyperparameters
        - delta2, gamma, zeta : repulsion hyperparameters
        - initial_clusters : np.ndarray of shape (N,), int32 (0-indexed)
    """
    n = dist_matrix.shape[0]

    # Step 1: Classical MDS to 2D
    mds_coords = classical_mds(dist_matrix, n_components=2)

    # Step 2: K-Means
    km = KMeans(n_clusters=k_elbow, n_init=25, random_state=random_state)
    initial_clusters = km.fit_predict(mds_coords).astype(np.int32)

    # Step 3: Split into within-cluster and between-cluster distances
    same_cluster = initial_clusters[:, None] == initial_clusters[None, :]
    i_upper, j_upper = np.triu_indices(n, k=1)

    within_mask = same_cluster[i_upper, j_upper]
    within_dists = dist_matrix[i_upper, j_upper][within_mask]
    inter_dists = dist_matrix[i_upper, j_upper][~within_mask]

    # Helper: Fit Gamma parameters by method of moments with constraints
    def _fit_gamma(
        dists: np.ndarray, is_within: bool, def_shape: float = 0.5
    ) -> Tuple[float, float, float]:
        valid = dists[dists > 1e-10]
        if len(valid) <= 1:
            return def_shape, 2.0, 2.0

        m = float(np.mean(valid))
        v = float(np.var(valid))
        shape = (m**2) / v if v > 0 else def_shape

        if is_within and shape > 1.0:
            shape = 0.9
        elif not is_within and shape < 1.0:
            shape = 1.1

        p1 = float(shape * len(valid))
        p2 = float(np.sum(valid))
        return shape, p1, p2

    delta1, alpha, beta = _fit_gamma(within_dists, is_within=True, def_shape=0.5)
    delta2, zeta, gamma = _fit_gamma(inter_dists, is_within=False, def_shape=2.0)

    return {
        "delta1": delta1,
        "alpha": alpha,
        "beta": beta,
        "delta2": delta2,
        "gamma": gamma,
        "zeta": zeta,
        "initial_clusters": initial_clusters,
    }


def binder_point_estimate(
    allocations: List[np.ndarray], burn_in: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the point cluster estimate minimizing Binder's loss to the PSM.

    Parameters
    ----------
    allocations : list of np.ndarray
        Posterior samples of allocation vectors.
    burn_in : int
        Number of initial iterations to discard as burn-in.

    Returns
    -------
    point_estimate : np.ndarray
        Best allocation vector from the posterior chain.
    psm : np.ndarray of shape (N, N)
        Posterior similarity matrix.
    """
    post_alloc = np.array(allocations[burn_in:], dtype=np.int32)
    S, N = post_alloc.shape

    # Posterior similarity matrix
    psm = np.zeros((N, N), dtype=np.float64)
    for a in post_alloc:
        psm += (a[:, None] == a[None, :])
    psm /= S

    # Find allocation minimizing distance to PSM
    best_idx = 0
    min_loss = float("inf")
    for idx, a in enumerate(post_alloc):
        adj = (a[:, None] == a[None, :]).astype(np.float64)
        loss = np.sum(np.abs(adj - psm))
        if loss < min_loss:
            min_loss = loss
            best_idx = idx

    return post_alloc[best_idx], psm


def compute_clustering_metrics(
    estimated: np.ndarray, ground_truth: np.ndarray
) -> Dict[str, Any]:
    """Compute standard clustering evaluation metrics."""
    ari = float(adjusted_rand_score(ground_truth, estimated))
    nmi = float(normalized_mutual_info_score(ground_truth, estimated))
    k_est = int(len(np.unique(estimated)))
    return {
        "K_est": k_est,
        "ARI": ari,
        "NMI": nmi,
    }
