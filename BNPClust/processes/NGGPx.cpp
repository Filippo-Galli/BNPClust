/**
 * @file NGGPx.cpp
 * @brief Implementation of NGGP with module-based similarity terms
 *
 * This file contains the implementation of the NGGPx class, which extends the
 * Normalized Generalized Gamma Process to incorporate module-based similarity
 * computations. Modules can represent any type of similarity component.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "./NGGPx.hpp"

double NGGPx::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {

    // NGGP gibbs prior for existing cluster
    double log_prior = NGGP::gibbs_prior_existing_cluster(cls_idx, obs_idx);

    // Add module-based similarity contribution
    for (auto &mod : modules) {
        log_prior += mod->compute_similarity_obs(obs_idx, cls_idx);
    }

    return log_prior;
}

Eigen::VectorXd NGGPx::gibbs_prior_existing_clusters(int obs_idx) const {

    // Get NGGP gibbs prior for existing clusters
    Eigen::VectorXd log_prior = NGGP::gibbs_prior_existing_clusters(obs_idx);

    // Add module-based similarity contributions
    for (auto &mod : modules) {
        log_prior += mod->compute_similarity_obs(obs_idx);
    }
    return log_prior;
}

double NGGPx::gibbs_prior_new_cluster() const { return NGGP::gibbs_prior_new_cluster(); }

double NGGPx::gibbs_prior_new_cluster_obs(int obs_idx) const {
    double log_prior = NGGP::gibbs_prior_new_cluster();

    // add module-based similarity contributions
    for (auto &mod : modules) {
        log_prior += mod->compute_similarity_obs(obs_idx, -1);
    }
    return log_prior;
}

double NGGPx::prior_ratio_split(int ci, int cj) const {
    // NGGP prior ratio for split
    double log_prior_ratio = NGGP::prior_ratio_split(ci, cj);

    // Module-based similarity ratio: sim(new ci) + sim(new cj) - sim(old merged ci)
    for (auto &mod : modules) {
        log_prior_ratio += mod->compute_similarity_cls(ci, false);
        log_prior_ratio += mod->compute_similarity_cls(cj, false);
        log_prior_ratio -= mod->compute_similarity_cls(ci, true);
    }
    return log_prior_ratio;
}

double NGGPx::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
    // NGGP prior ratio for merge
    double log_prior_ratio = NGGP::prior_ratio_merge(size_old_ci, size_old_cj);

    // Module-based similarity ratio: sim(new merged ci) - sim(old ci) - sim(old cj)
    const int new_ci = NGGP::data.get_allocations()[idx_i];
    for (auto &mod : modules) {
        log_prior_ratio += mod->compute_similarity_cls(new_ci, false);
    }

    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    for (auto &mod : modules) {
        log_prior_ratio -= mod->compute_similarity_cls(old_ci, true);
        log_prior_ratio -= mod->compute_similarity_cls(old_cj, true);
    }

    return log_prior_ratio;
}

double NGGPx::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci, int cj) const {
    // NGGP prior ratio for shuffle
    double log_prior_ratio = NGGP::prior_ratio_shuffle(size_old_ci, size_old_cj, ci, cj);

    // Module-based similarity ratio: sim(new ci) + sim(new cj) - sim(old ci) - sim(old cj)
    const int new_ci = NGGP::data.get_allocations()[idx_i];
    const int new_cj = NGGP::data.get_allocations()[idx_j];
    for (auto &mod : modules) {
        log_prior_ratio += mod->compute_similarity_cls(new_ci, false);
        log_prior_ratio += mod->compute_similarity_cls(new_cj, false);
    }

    const int old_ci = old_allocations[idx_i];
    const int old_cj = old_allocations[idx_j];
    for (auto &mod : modules) {
        log_prior_ratio -= mod->compute_similarity_cls(old_ci, true);
        log_prior_ratio -= mod->compute_similarity_cls(old_cj, true);
    }

    return log_prior_ratio;
}
