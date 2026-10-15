/**
 * @file categorical_covariate_module.cpp
 * @brief Implementation of `CategoricalCovariatesModule`.
 */

#include "categorical_covariate_module.hpp"

double CategoricalCovariatesModule::compute_similarity_cls(int cls_idx, bool old_allo) const {
    // 1. Retrieve cluster members
    Eigen::VectorXi cluster_members;
    if (old_allo && old_cluster_members_provider) {
        auto it = old_cluster_members_provider->find(cls_idx);
        if (it != old_cluster_members_provider->end()) {
            cluster_members = Eigen::Map<const Eigen::VectorXi>(it->second.data(), it->second.size());
        }
    } else {
        cluster_members = data.get_cluster_assignments(cls_idx);
    }

    if (cluster_members.size() == 0)
        return 0.0;

    // 2. Count occurrences for each category (0 to C)
    int num_categories = prior_alpha.size();
    std::vector<int> n_ji(num_categories, 0);
    for (int i = 0; i < cluster_members.size(); ++i) {
        int val = categorical_covariate_data(cluster_members(i));
        if (val >= 0 && val < num_categories) {
            n_ji[val]++;
        }
    }

    // 3. Compute Log Marginal Likelihood using your derivation:
    // log g(x) = [lgamma(alpha_0) - lgamma(alpha_0 + n_j)] + sum(lgamma(alpha_i + n_ji)) - sum(lgamma(alpha_i))
    double sum_lgamma_data = 0.0;
    for (int i = 0; i < num_categories; ++i) {
        sum_lgamma_data += std::lgamma(prior_alpha[i] + n_ji[i]);
    }

    double log_marginal_likelihood =
        lgamma_alpha_0 - std::lgamma(alpha_0 + cluster_members.size()) + sum_lgamma_data - prod_lgamma_prior;

    return log_marginal_likelihood;
}

double CategoricalCovariatesModule::compute_similarity_obs(int obs_idx, int cls_idx) const {
    int n_j = 0;                                 // Total count in cluster
    int n_jl = 0;                                // Count of category 'l' in cluster
    int l = categorical_covariate_data(obs_idx); // The category of the current observation

    if (cls_idx >= 0) {
        const auto &cluster_members = data.get_cluster_assignments(cls_idx);
        n_j = cluster_members.size();
        for (const int &member_idx : cluster_members) {
            if (categorical_covariate_data(member_idx) == l) {
                n_jl++;
            }
        }
    }

    // Using: log( (alpha_l + n_jl) / (alpha_0 + n_j) )
    return std::log(prior_alpha[l] + n_jl) - std::log(alpha_0 + n_j);
}

Eigen::VectorXd CategoricalCovariatesModule::compute_similarity_obs(int obs_idx) const {
    const int K = data.get_K();                                  // Number of clusters
    const Eigen::VectorXi &allocations = data.get_allocations(); // Current cluster assignments
    int l = categorical_covariate_data(obs_idx);                 // The category of the current observation

    std::vector<int> cluster_sizes(K, 0);
    std::vector<int> cluster_counts_l(K, 0); // Only need counts for category 'l'

    // Pass over data to get stats for category 'l' across all clusters
    for (int i = 0; i < allocations.size(); ++i) {
        if (i == obs_idx)
            continue;

        int k = allocations[i];
        if (k >= 0 && k < K) {
            cluster_sizes[k]++;
            if (categorical_covariate_data(i) == l) {
                cluster_counts_l[k]++;
            }
        }
    }

    Eigen::VectorXd similarities(K);
    for (int k = 0; k < K; ++k) {
        similarities(k) = std::log(prior_alpha[l] + cluster_counts_l[k]) - std::log(alpha_0 + cluster_sizes[k]);
    }

    return similarities;
}
