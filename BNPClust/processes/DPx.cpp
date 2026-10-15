/**
 * @file DPx.cpp
 * @brief Implementation of Dirichlet Process with module-based similarity terms
 *
 * This file contains the implementation of the DPx class, which extends the standard
 * Dirichlet Process to incorporate module-based similarity computations. Modules can
 * represent any type of similarity component (spatial, covariate-based, etc.).
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "DPx.hpp"
#include "DP.hpp"

Eigen::VectorXd DPx::gibbs_prior_existing_clusters(int obs_idx) const {
    /**
     * @brief Computes the log prior probabilities of assigning a data point to
     * all existing clusters.
     *
     * This method incorporates module-based similarity contributions to compute
     * the prior probabilities. Modules can represent any type of similarity
     * component (spatial, covariate-based, etc.).
     * @param obs_idx The index of the observation to assign.
     * @return A vector of log prior probabilities for assigning the data point to
     * each existing cluster.
     */

    // Get DP gibbs prior for existing clusters
    Eigen::VectorXd log_prior = DP::gibbs_prior_existing_clusters(obs_idx);

    // Add covariate module contributions
    for (auto &mod : modules) {
        log_prior += mod->compute_similarity_obs(obs_idx);
    }
    return log_prior;
}

double DPx::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {
    double prior = DP::gibbs_prior_existing_cluster(cls_idx, obs_idx);
    // Add covariate module contribution
    for (auto &mod : modules) {
        prior += mod->compute_similarity_obs(obs_idx, cls_idx);
    }

    return prior;
}

double DPx::gibbs_prior_new_cluster() const { return DP::gibbs_prior_new_cluster(); }

double DPx::gibbs_prior_new_cluster_obs(int obs_idx) const {
    double log_prior = DP::gibbs_prior_new_cluster();

    // add covariate module contributions
    for (auto &mod : modules) {
        log_prior += mod->compute_similarity_obs(obs_idx, -1);
    }
    return log_prior;
}

double DPx::prior_ratio_split(int ci, int cj) const {
    /**
     * @brief Computes the prior ratio for a split operation in split-merge MCMC.
     *
     * This method accounts for both the Dirichlet Process prior and module-based
     * similarity terms when computing the acceptance ratio for splitting clusters.
     * @param ci The first cluster index involved in the split.
     * @param cj The second cluster index involved in the split.
     * @return The log prior ratio for the split operation.
     */

    double log_acceptance_ratio = DP::prior_ratio_split(ci, cj);

    for (auto &mod : modules) {
        log_acceptance_ratio += mod->compute_similarity_cls(ci, false);
        log_acceptance_ratio += mod->compute_similarity_cls(cj, false);
        log_acceptance_ratio -= mod->compute_similarity_cls(ci, true);
    }

    return log_acceptance_ratio;
}

double DPx::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
    /**
     * @brief Computes the prior ratio for a merge operation in split-merge MCMC.
     *
     * This method accounts for both the Dirichlet Process prior and module-based
     * similarity terms when computing the acceptance ratio for merging clusters.
     * @param size_old_ci The size of the first cluster before the merge.
     * @param size_old_cj The size of the second cluster before the merge.
     * @return The log prior ratio for the merge operation.
     */

    // DP prior part
    double log_acceptance_ratio = DP::prior_ratio_merge(size_old_ci, size_old_cj);
    // Spatial / module part: sim(new merged ci) - sim(old ci) - sim(old cj)
    const int new_ci = data.get_allocations()[idx_i];
    for (auto &mod : modules) {
        log_acceptance_ratio += mod->compute_similarity_cls(new_ci, false);
    }

    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    for (auto &mod : modules) {
        log_acceptance_ratio -= mod->compute_similarity_cls(old_ci, true);
        log_acceptance_ratio -= mod->compute_similarity_cls(old_cj, true);
    }
    return log_acceptance_ratio;
}

double DPx::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) const {
    /**
     * @brief Computes the prior ratio for a shuffle operation in split-merge MCMC.
     *
     * This method accounts for both the Dirichlet Process prior and module-based
     * similarity terms when computing the acceptance ratio for shuffling observations
     * between clusters.
     * @param size_old_ci The size of the first cluster before the shuffle.
     * @param size_old_cj The size of the second cluster before the shuffle.
     * @param ci The first cluster index involved in the shuffle.
     * @param cj The second cluster index involved in the shuffle.
     * @return The log prior ratio for the shuffle operation.
     */

    // DP prior part
    double log_acceptance_ratio = DP::prior_ratio_shuffle(size_old_ci, size_old_cj, ci, cj);

    // Module-based similarity ratio: sim(new ci) + sim(new cj) - sim(old ci) - sim(old cj)
    const int new_ci = data.get_allocations()[idx_i];
    const int new_cj = data.get_allocations()[idx_j];
    for (auto &mod : modules) {
        log_acceptance_ratio += mod->compute_similarity_cls(new_ci, false);
        log_acceptance_ratio += mod->compute_similarity_cls(new_cj, false);
    }

    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    for (auto &mod : modules) {
        log_acceptance_ratio -= mod->compute_similarity_cls(old_ci, true);
        log_acceptance_ratio -= mod->compute_similarity_cls(old_cj, true);
    }

    return log_acceptance_ratio;
}
