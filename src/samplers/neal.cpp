/**
 * @file neal.cpp
 * @brief Implementation of Neal's Algorithm 3 for collapsed Gibbs sampling
 *
 * This file contains the complete implementation of the Neal3 class, providing
 * efficient collapsed Gibbs sampling for Bayesian nonparametric mixture models.
 * The implementation includes sequential observation updates with automatic
 * cluster management and proper probability computations.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "neal.hpp"
#include <numeric>

// In neal.cpp, could be a private helper in Neal3 or a free function
int Neal3::sample_from_log_probs(int num_clusters) {
    // 1. Log-Sum-Exp trick for numerical stability
    double max_loglik = *std::max_element(log_likelihoods.begin(), log_likelihoods.begin() + num_clusters);

    // Reuse pre-allocated weights vector
    weights.resize(num_clusters);
    double sum_weights = 0.0;
    for (int i = 0; i < num_clusters; ++i) {
        weights[i] = exp(log_likelihoods[i] - max_loglik);
        sum_weights += weights[i];
    }

    // 2. Roulette wheel selection
    std::uniform_real_distribution<double> unif(0.0, sum_weights);
    double u = unif(gen);

    for (int i = 0; i < num_clusters; ++i) {
        u -= weights[i];
        if (u < 0.0) {
            return i;
        }
    }
    // Fallback for floating point inaccuracies
    return num_clusters - 1;
}

void Neal3::step_1_observation(int index) {
    /**
     * @brief Performs a step in the DPNeal2 sampling process.
     * @details This method is responsible for updating the allocations of the
     * data points based on the current state of the model.
     * @param index The index of the data point to update.
     */

    // Set unallocated the index
    data.set_allocation(index, -1);

    const int K = data.get_K();
    const int num_clusters = K + 1;

    // Resize vector only if needed (capacity already reserved in constructor)
    log_likelihoods.resize(num_clusters);

    // Compute combined log(likelihood * prior) for existing clusters
    // This combines likelihood computation with prior multiplication in one pass
    for (int k = 0; k < K; ++k) {
        const double log_prior = process.gibbs_prior_existing_cluster(k, index);
        log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k) + log_prior;
    }

    // Compute for new cluster
    log_likelihoods[K] = likelihood.point_loglikelihood_cond(index, K) + process.gibbs_prior_new_cluster_obs(index);

    // Sample a cluster based on the probabilities
    int sampled_cluster = sample_from_log_probs(num_clusters);

    // Set the allocation for the data point
    data.set_allocation(index, sampled_cluster);
}

void Neal3::step() {
    /**
     * @brief Performs a single step of the Neal 3 algorithm for all the
     * dataset.
     */

    for (const int idx : indices) {
        step_1_observation(idx);
    }
}
