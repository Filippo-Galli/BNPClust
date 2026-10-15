/**
 * @file binary_covariate_module_cache.cpp
 * @brief Implementation of `BinaryCovariatesModuleCache`.
 */

#include "binary_covariate_module_cache.hpp"

double BinaryCovariatesModuleCache::compute_similarity_cls(int cls_idx, bool old_allo) const {

    int counts = 0;
    int num_covariates = 0;

    // Retrieve cluster members based on allocation type
    Eigen::VectorXi cluster_members;
    if (old_allo && old_cluster_members_provider) {
        auto it = old_cluster_members_provider->find(cls_idx);
        if (it != old_cluster_members_provider->end()) {
            cluster_members = Eigen::Map<const Eigen::VectorXi>(it->second.data(), it->second.size());
        }
        // Initialize counts for binary covariates
        num_covariates = cluster_members.size();

        // Count occurrences of '1' in the binary covariates for the cluster members
        for (int i = 0; i < num_covariates; ++i) {
            int obs_idx = cluster_members(i);
            counts += cache.binary_covariates(obs_idx);
        }
    } else {
        num_covariates = cache.get_cluster_stats_ref(cls_idx).n;
        counts = cache.get_cluster_stats_ref(cls_idx).binary_sum;
    }

    // If the cluster is empty, return 0 similarity
    if (num_covariates == 0) {
        return 0.0;
    }

    // Compute log marginal likelihood using Beta-Binomial model
    double alpha = beta_prior_alpha;
    double beta = beta_prior_beta;

    double log_marginal_likelihood = lgamma(counts + alpha) + lgamma(num_covariates - counts + beta) -
                                     lgamma(num_covariates + alpha + beta) - log_beta_prior;

    return log_marginal_likelihood;
}

double BinaryCovariatesModuleCache::compute_similarity_obs(int obs_idx, int cls_idx) const {

    int num_covariates = 0;
    int counts = 0;

    // Check if checking against an existing cluster
    if (cls_idx >= 0) {
        // Initialize counts for binary covariates
        num_covariates = cache.get_cluster_stats_ref(cls_idx).n;

        // Count occurrences of '1' in the binary covariates for the cluster members
        counts = cache.get_cluster_stats_ref(cls_idx).binary_sum;

        // Check if this cluster currently contains obs_idx
        if (data.get_allocations()(obs_idx) == cls_idx) {
            num_covariates--;
            counts -= cache.binary_covariates(obs_idx);
        }
    }

    const double denominator = num_covariates + beta_prior_alpha + beta_prior_beta;
    const double numerator = (beta_prior_alpha + counts) * cache.binary_covariates(obs_idx) +
                             (beta_prior_beta + num_covariates - counts) * (1 - cache.binary_covariates(obs_idx));
    return std::log(numerator / denominator);
}

Eigen::VectorXd BinaryCovariatesModuleCache::compute_similarity_obs(int obs_idx) const {
    const int K = data.get_K();
    const Eigen::VectorXi &allocations = data.get_allocations();

    Eigen::VectorXd similarities(K);
    bool is_success = (cache.binary_covariates(obs_idx) == 1);
    double alpha = beta_prior_alpha;
    double beta = beta_prior_beta;
    int obs_cluster = allocations(obs_idx);

    // Compute similarities using gathered stats
    for (int k = 0; k < K; ++k) {
        double num;
        int cluster_counts = cache.get_cluster_stats_ref(k).binary_sum;
        int cluster_sizes = cache.get_cluster_stats_ref(k).n;

        if (k == obs_cluster) {
            cluster_sizes--;
            cluster_counts -= cache.binary_covariates(obs_idx);
        }

        if (is_success) {
            num = alpha + cluster_counts;
        } else {
            num = beta + cluster_sizes - cluster_counts;
        }

        double den = cluster_sizes + alpha + beta;
        similarities(k) = std::log(num) - std::log(den);
    }

    return similarities;
}
