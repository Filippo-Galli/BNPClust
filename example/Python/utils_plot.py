"""Plotting and visualization utilities for BNPClust MCMC clustering analysis in Python.

Provides functions to compute point estimates, plot posterior distributions,
traces of parameters, posterior similarity matrices, and clustering performance
metrics, matching the R workflow in example/R/utils_plot.R.
"""

import os
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend suitable for headless scripts
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.cluster import contingency_matrix

from utils import binder_point_estimate


def compute_ess(x: np.ndarray) -> float:
    """Compute Effective Sample Size (ESS) using Geyer's initial positive sequence estimator."""
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    if n <= 1:
        return 1.0
    var = np.var(x)
    if var <= 1e-12:
        return float(n)

    # Autocorrelation via FFT
    x_centered = x - np.mean(x)
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    fx = np.fft.fft(x_centered, n=n_fft)
    rx = np.fft.ifft(fx * np.conj(fx)).real[:n]
    rx /= rx[0]  # rho_k

    # Sum adjacent pairs
    pair_sums = rx[1 : n - 1 : 2] + rx[2:n:2]
    pos_mask = pair_sums > 0
    if not np.any(pos_mask):
        act = 1.0
    else:
        first_neg = np.where(~pos_mask)[0]
        cutoff = first_neg[0] if len(first_neg) > 0 else len(pair_sums)
        act = 1.0 + 2.0 * float(np.sum(pair_sums[:cutoff]))

    return float(max(1.0, n / max(act, 1.0)))


def variation_of_information(c1: np.ndarray, c2: np.ndarray) -> Tuple[float, float]:
    """Compute Variation of Information (VI) and Normalized VI (NVI).

    Returns
    -------
    vi : float
        Unnormalized Variation of Information.
    nvi : float
        Normalized Variation of Information (1 - MI / H(c1, c2)).
    """
    n = len(c1)
    cm = contingency_matrix(c1, c2).astype(np.float64)
    p = cm / n
    pi = np.sum(p, axis=1)
    pj = np.sum(p, axis=0)

    h1 = -np.sum(pi[pi > 0] * np.log(pi[pi > 0]))
    h2 = -np.sum(pj[pj > 0] * np.log(pj[pj > 0]))
    h12 = -np.sum(p[p > 0] * np.log(p[p > 0]))

    mi = h1 + h2 - h12
    vi = h1 + h2 - 2.0 * mi
    nvi = 1.0 - (mi / h12) if h12 > 0 else 0.0
    return float(vi), float(nvi)


def plot_cls_est(
    results: Dict[str, Any],
    BI: int = 0,
    save: bool = False,
    folder: str = "results/plots/",
) -> np.ndarray:
    """Compute and print cluster point estimate, cluster sizes, and ESS."""
    print("Computing point estimate using Binder's loss to PSM...")
    allocations = results["allocations"]
    pe, psm = binder_point_estimate(allocations, burn_in=BI)

    unique_cls, counts = np.unique(pe, return_counts=True)
    cluster_sizes = dict(zip(unique_cls, counts))

    k_values = np.array(results["K"][BI:], dtype=np.float64)
    ess_value = compute_ess(k_values)
    elapsed_time = results.get("elapsed_time", 1.0)
    ess_seconds = ess_value / max(elapsed_time, 1e-6)

    print("=== Clustering Results (Post Burn-in) ===")
    print("Cluster Sizes:")
    for c, cnt in cluster_sizes.items():
        print(f"  Cluster {c}: {cnt}")
    print("=== ESS ===")
    print(f"Effective Sample Size (ESS): {ess_value:.2f}")
    print(f"ESS per second: {ess_seconds:.2f}")

    if save:
        os.makedirs(folder, exist_ok=True)
        cls_file = os.path.join(folder, "cluster_estimate.txt")
        with open(cls_file, "w") as f:
            f.write("Cluster Sizes:\n")
            for c, cnt in cluster_sizes.items():
                f.write(f"Cluster {c}: {cnt}\n")
            f.write("\n=== ESS ===\n")
            f.write(f"ESS: {ess_value:.2f}\n")
            f.write(f"ESS per second: {ess_seconds:.2f}\n")

    return pe


def plot_stats(
    results: Dict[str, Any],
    ground_truth: np.ndarray,
    BI: int = 0,
    save: bool = False,
    folder: str = "results/plots/",
) -> Dict[str, float]:
    """Compute and display ARI, NMI, and NVI metrics against ground truth."""
    allocations = results["allocations"]
    pe, _ = binder_point_estimate(allocations, burn_in=BI)

    ari_val = float(adjusted_rand_score(ground_truth, pe))
    nmi_val = float(normalized_mutual_info_score(ground_truth, pe))
    vi_val, nvi_val = variation_of_information(ground_truth, pe)

    print(f"\nAdjusted Rand Index: {ari_val:.4f}")
    print(f"Normalized Mutual Information: {nmi_val:.4f}")
    print(f"Variation of Information: {vi_val:.4f} (Normalized: {nvi_val:.4f})")

    if save:
        os.makedirs(folder, exist_ok=True)
        stats_file = os.path.join(folder, "stats.txt")
        with open(stats_file, "w") as f:
            f.write(f"Adjusted Rand Index: {ari_val:.4f}\n")
            f.write(f"Normalized Mutual Information: {nmi_val:.4f}\n")
            f.write(f"Variation of Information: {vi_val:.4f}\n")
            f.write(f"Normalized Variation of Information: {nvi_val:.4f}\n")

    return {"ARI": ari_val, "NMI": nmi_val, "VI": vi_val, "NVI": nvi_val}


def plot_post_distr(
    results: Dict[str, Any],
    BI: int = 0,
    save: bool = False,
    folder: str = "results/plots/",
    ylim: Optional[Tuple[float, float]] = None,
    xlim: Optional[Tuple[float, float]] = None,
) -> None:
    """Plot posterior distribution of the number of clusters K."""
    k_values = np.array(results["K"][BI:], dtype=np.int32)
    vals, counts = np.unique(k_values, return_counts=True)
    rel_freq = counts / len(k_values)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(vals, rel_freq, width=0.8, color="#2c3e50", edgecolor="none")
    ax.set_xlabel("Number of clusters (K)", fontsize=14)
    ax.set_ylabel("Posterior Probability", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)
    if ylim:
        ax.set_ylim(ylim)
    if xlim:
        ax.set_xlim(xlim)

    if save:
        os.makedirs(folder, exist_ok=True)
        plt.savefig(os.path.join(folder, "posterior_num_clusters.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_trace_cls(
    results: Dict[str, Any],
    BI: int = 0,
    save: bool = False,
    folder: str = "results/plots/",
) -> None:
    """Plot MCMC trace of the number of clusters K."""
    k_values = np.array(results["K"][BI:], dtype=np.int32)
    iters = np.arange(1, len(k_values) + 1)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(iters, k_values, color="#1f77b4", linewidth=1.0)
    ax.set_xlabel("Iteration (post burn-in)", fontsize=14)
    ax.set_ylabel("Number of clusters", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)

    if save:
        os.makedirs(folder, exist_ok=True)
        plt.savefig(os.path.join(folder, "traceplot.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_post_sim_matrix(
    results: Dict[str, Any],
    BI: int = 0,
    save: bool = False,
    folder: str = "results/plots/",
) -> None:
    """Compute and plot the reordered posterior similarity matrix heatmap."""
    allocations = results["allocations"]
    _, psm = binder_point_estimate(allocations, burn_in=BI)
    n = psm.shape[0]

    # Hierarchical clustering reordering
    dist_psm = np.clip(1.0 - psm, 0.0, 1.0)
    np.fill_diagonal(dist_psm, 0.0)
    condensed = squareform(dist_psm, checks=False)
    z = linkage(condensed, method="average")
    order = leaves_list(z)

    psm_reordered = psm[order, :][:, order]

    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(
        psm_reordered,
        cmap="Blues",
        vmin=0.0,
        vmax=1.0,
        origin="lower",
        interpolation="nearest",
    )
    ax.set_xlabel("Data Point Index (clustered)", fontsize=13)
    ax.set_ylabel("Data Point Index (clustered)", fontsize=13)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Posterior Similarity", fontsize=12)

    if save:
        os.makedirs(folder, exist_ok=True)
        plt.savefig(os.path.join(folder, "similarity_matrix.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_trace_U(
    results: Dict[str, Any],
    BI: int = 0,
    save: bool = False,
    folder: str = "results/plots/",
) -> None:
    """Plot MCMC trace of the latent auxiliary variable U."""
    u_values = np.array(results["U"][BI:], dtype=np.float64)
    iters = np.arange(1, len(u_values) + 1)
    mean_u = float(np.mean(u_values))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(iters, u_values, color="#34495e", linewidth=1.0, label="U")
    ax.axhline(mean_u, color="red", linestyle="--", linewidth=1.5, label=f"Mean U = {mean_u:.3f}")
    ax.set_xlabel("Iteration (post burn-in)", fontsize=14)
    ax.set_ylabel("U", fontsize=14)
    ax.legend(loc="upper right", frameon=False, fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.3)

    if save:
        os.makedirs(folder, exist_ok=True)
        plt.savefig(os.path.join(folder, "U_trace.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_acf_U(
    results: Dict[str, Any],
    BI: int = 0,
    max_lag: int = 35,
    save: bool = False,
    folder: str = "results/plots/",
) -> None:
    """Plot autocorrelation function (ACF) of the latent variable U."""
    u_values = np.array(results["U"][BI:], dtype=np.float64)
    u_centered = u_values - np.mean(u_values)
    n = len(u_values)
    var = np.var(u_values)

    if var <= 1e-12:
        acf_vals = np.ones(max_lag + 1)
    else:
        acf_vals = np.correlate(u_centered, u_centered, mode="full")[n - 1 : n + max_lag]
        acf_vals /= acf_vals[0]

    lags = np.arange(len(acf_vals))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.vlines(lags, 0, acf_vals, colors="#2980b9", linewidth=2.0)
    ax.axhline(0, color="black", linewidth=0.8)
    ci = 1.96 / np.sqrt(n)
    ax.axhline(ci, color="blue", linestyle="--", linewidth=0.8)
    ax.axhline(-ci, color="blue", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Lag", fontsize=14)
    ax.set_ylabel("ACF", fontsize=14)
    ax.set_title("ACF of U over MCMC iterations", fontsize=14)
    ax.grid(True, linestyle="--", alpha=0.3)

    if save:
        os.makedirs(folder, exist_ok=True)
        plt.savefig(os.path.join(folder, "U_acf.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
