/**
 * @file binary_covariate_module.cpp
 * @brief Implementation of `BinaryCovariatesModule`.
 */

#include "binary_covariate_module.hpp"

double BinaryCovariatesModule::compute_similarity_cls(int cls_idx, bool old_allo) const {
    // Retrieve cluster members based on allocation type
    Eigen::VectorXi cluster_members;
    if (old_allo && old_cluster_members_provider) {
        auto it = old_cluster_members_provider->find(cls_idx);
        if (it != old_cluster_members_provider->end()) {
            cluster_members = Eigen::Map<const Eigen::VectorXi>(it->second.data(), it->second.size());
        }
    } else {
        cluster_members = data.get_cluster_assignments(cls_idx);
    }

    // If the cluster is empty, return 0 similarity
    if (cluster_members.size() == 0) {
        return 0.0;
    }

    // Initialize counts for binary covariates
    int num_covariates = cluster_members.size();
    int counts = 0;

    // Count occurrences of '1' in the binary covariates for the cluster members
    for (int i = 0; i < num_covariates; ++i) {
        int obs_idx = cluster_members(i);
        counts += binary_covariate_data(obs_idx);
    }

    // Compute log marginal likelihood using Beta-Binomial model
    double alpha = beta_prior_alpha;
    double beta = beta_prior_beta;

    double log_marginal_likelihood = lgamma(counts + alpha) +
                                     lgamma(num_covariates - counts + beta) -
                                     lgamma(num_covariates + alpha + beta) - log_beta_prior;

    return log_marginal_likelihood;
}

double BinaryCovariatesModule::compute_similarity_obs(int obs_idx, int cls_idx) const {
    
    int num_covariates = 0;
    int counts = 0;

    // Check if checking against an existing cluster
    if (cls_idx >= 0) {
        const auto & cluster_members = data.get_cluster_assignments(cls_idx);

        // Initialize counts for binary covariates
        num_covariates = cluster_members.size();

        // Count occurrences of '1' in the binary covariates for the cluster members
        for (int i = 0; i < num_covariates; ++i) {
            int member_idx = cluster_members(i);
            counts += binary_covariate_data(member_idx);
        }
    }

    const double denominator = num_covariates + beta_prior_alpha + beta_prior_beta;
    const double numerator = (beta_prior_alpha + counts)*binary_covariate_data(obs_idx) + 
                             (beta_prior_beta + num_covariates - counts)*(1 - binary_covariate_data(obs_idx));
    return std::log(numerator / denominator);
}

Eigen::VectorXd BinaryCovariatesModule::compute_similarity_obs(int obs_idx) const{
    const int K = data.get_K();
    const Eigen::VectorXi &allocations = data.get_allocations();
    
    // Pre-allocate temporary storage for cluster statistics
    // You could make these member variables (mutable) to avoid re-allocating memory every step
    std::vector<int> cluster_counts(K, 0);
    std::vector<int> cluster_sizes(K, 0);

    // Single pass over allocations to compute stats for ALL clusters
    for (int i = 0; i < allocations.size(); ++i) {
        // Skip the observation itself to ensure predictive distribution is correct
        if (i == obs_idx) continue;

        int k = allocations[i];
        if (k >= 0 && k < K) {
            cluster_sizes[k]++;
            if (binary_covariate_data(i) == 1) {
                cluster_counts[k]++;
            }
        }
    }

    Eigen::VectorXd similarities(K);
    bool is_success = (binary_covariate_data(obs_idx) == 1);
    double alpha = beta_prior_alpha;
    double beta = beta_prior_beta;

    // Compute similarities using gathered stats
    for (int k = 0; k < K; ++k) {
        // If the cluster ended up empty (e.g. only contained obs_idx), handle gracefully
        // Usually implementation details might vary, but log(0) for empty cluster logic:
        // Use prior predictive for empty cluster logic if it's considered "active"
        
        double num;
        if (is_success) {
            num = alpha + cluster_counts[k];
        } else {
            num = beta + cluster_sizes[k] - cluster_counts[k];
        }
        
        double den = cluster_sizes[k] + alpha + beta;
        similarities(k) = std::log(num) - std::log(den);
    }
    
    return similarities;
}
