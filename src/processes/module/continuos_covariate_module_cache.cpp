/**
 * @file continuos_covariate_module_cache.cpp
 * @brief Implementation of `ContinuosCovariatesModuleCache`.
 */

#include "continuos_covariate_module_cache.hpp"

ContinuosCache::ClusterStats
ContinuosCovariatesModuleCache::compute_cluster_statistics(const Eigen::Ref<const Eigen::VectorXi> obs) const {
    ContinuosCache::ClusterStats stats;

    for (int idx = 0; idx < obs.size(); ++idx) {
        const int obs_idx = obs(idx);
        const double value = continuos_cache.continuos_covariates(obs_idx);
        stats.n += 1;
        stats.sum += value;
        stats.sumsq += value * value;
    }

    return stats;
}

double ContinuosCovariatesModuleCache::compute_similarity_cls(int cls_idx, bool old_allo) const {

    if (old_allo && old_cluster_members_provider) {
        ContinuosCache::ClusterStats stats;
        const auto old_cls_it = old_cluster_members_provider->find(cls_idx);
        if (old_cls_it == old_cluster_members_provider->end()) {
            return compute_log_marginal_likelihood(stats);
        }

        const auto &old_cls_allo = old_cls_it->second;
        stats = compute_cluster_statistics(Eigen::Map<const Eigen::VectorXi>(old_cls_allo.data(), old_cls_allo.size()));
        return compute_log_marginal_likelihood(stats);

    } else {
        const Eigen::VectorXi &obs = data.get_cluster_assignments_ref(cls_idx);
        const auto &stats = continuos_cache.get_cluster_stats_ref(cls_idx);
        return compute_log_marginal_likelihood(stats);
    }
}

double ContinuosCovariatesModuleCache::compute_similarity_obs(int obs_idx, int cls_idx) const {

    ContinuosCache::ClusterStats base_stats;

    const double covariate_val = continuos_cache.continuos_covariates(obs_idx);
    
    // Handle new cluster case
    if (cls_idx > -1 && cls_idx < data.get_K())
        base_stats = continuos_cache.get_cluster_stats_ref(cls_idx);

    return compute_log_predictive_likelihood(base_stats, covariate_val);
}

Eigen::VectorXd ContinuosCovariatesModuleCache::compute_similarity_obs(int obs_idx) const {

    const int num_clusters = data.get_K();
    // Compute log similarities for each cluster
    Eigen::VectorXd log_similarities(num_clusters);
    const double covariate_val = continuos_cache.continuos_covariates(obs_idx);

    for (int k = 0; k < num_clusters; ++k) {
        const ContinuosCache::ClusterStats &stats_ref = continuos_cache.get_cluster_stats_ref(k);
        log_similarities(k) = compute_log_predictive_likelihood(stats_ref, covariate_val);
    }

    return log_similarities;
}

double ContinuosCovariatesModuleCache::compute_log_marginal_likelihood_NNIG(const ContinuosCache::ClusterStats &stats) const {
    // if (stats.n == 0) {
    //     return 0.0;
    // }

    const double n_dbl = static_cast<double>(stats.n);
    const double inv_n = 1.0 / (n_dbl + 1e-20); // small term to avoid division by zero
    const double xbar = stats.sum * inv_n;

    // Centered sum of squares: SS = Σ(xᵢ - x̄)²
    const double ss = stats.sumsq - n_dbl * xbar * xbar;

    // Posterior parameters
    // Prior: v ~ IG(ν, S0) and μ|v ~ N(m, B v) with B interpreted as a variance multiplier.
    const double nu_n = nu + 0.5 * n_dbl;

    // Deviation from prior mean
    const double dev = xbar - m;

    // Posterior scale parameter:
    // S_n = S₀ + SS/2 + n/(2(1+nB)) (x̄-m)²
    const double one_plus_nB = 1.0 + n_dbl * B;
    const double S_n = S0 + 0.5 * ss + 0.5 * (n_dbl / one_plus_nB) * dev * dev;

    double lgamma_nu_n_temp;
    if (stats.n <= lgamma_nu_n.size() - 1) {
        lgamma_nu_n_temp = lgamma_nu_n[stats.n];
    } else {
        lgamma_nu_n_temp = std::lgamma(nu + 0.5 * n_dbl);
    }

    // log g(S) = log Γ(ν + n/2) - log Γ(ν) - n/2 log(2π)
    //            - 1/2 log(1+nB) + ν log(S₀) - (ν + n/2) log(S_n)
    return lgamma_nu_n_temp - lgamma_nu + n_dbl * const_term - 0.5 * std::log(one_plus_nB) + nu_logS0 -
           nu_n * std::log(S_n);
}

double ContinuosCovariatesModuleCache::compute_log_marginal_likelihood_NN(const ContinuosCache::ClusterStats &stats) const {

    // if (stats.n == 0) {
    //     return 0.0;
    // }

    const double n_dbl = static_cast<double>(stats.n);
    const double inv_n = 1.0 / (n_dbl + 1e-20); // small term to avoid division by zero
    const double xbar = stats.sum * inv_n;

    // Centered sum of squares: SS = Σ(x_i - x̄)²
    const double ss = stats.sumsq - n_dbl * xbar * xbar;

    // v + n_j B term from integrating μ_j ~ N(m, B)
    const double v_plus_nB = v + n_dbl * B;

    // Mean deviation from prior mean
    const double dev = xbar - m;

    // Log of posterior variance: log(τ_j) = log(B) + log(v) - log(v + n_j B)
    double log_v_plus_nB_temp;
    if (stats.n <= log_v_plus_nB.size() - 1) {
        log_v_plus_nB_temp = log_v_plus_nB[stats.n];
    } else {
        log_v_plus_nB_temp = std::log(v + n_dbl * B);
    }

    const double log_tau_j = log_B + log_v - log_v_plus_nB_temp;

    // Compute log marginal likelihood using posterior variance form:
    // log q(x_j) = -n_j/2 log(2π) - n_j/2 log(v) - 1/2 log(B) + 1/2 log(τ_j)
    //              - SS/(2v) - n_j(x̄_j - m)² / (2(v + n_j B))
    return n_dbl * const_term - 0.5 * n_dbl * log_v - 0.5 * log_B + 0.5 * log_tau_j -
           0.5 * (ss / v + n_dbl * dev * dev / v_plus_nB);
}

double ContinuosCovariatesModuleCache::compute_predictive_NN(const ContinuosCache::ClusterStats &stats, double covariate_val) const {
    double n_dbl = static_cast<double>(stats.n);
    double one_plus_nB = 1.0 + n_dbl * B;

    // Posterior mean (mu_n)
    // Formula: (m + B * sum_x) / (1 + n*B)
    double mu_n = (m + B * stats.sum) / one_plus_nB;

    // Predictive variance
    // Formula: v * (1 + (n+1)B) / (1 + nB)
    double sigma2_pred = v * (1.0 + (n_dbl + 1.0) * B) / one_plus_nB;

    // Log Normal PDF: -0.5 * log(2*pi*sigma2) - (x - mu)^2 / (2*sigma2)
    double diff = covariate_val - mu_n;
    return -0.5 * std::log(2.0 * M_PI * sigma2_pred) - 0.5 * diff * diff / sigma2_pred;
}

double ContinuosCovariatesModuleCache::compute_predictive_NNIG(const ContinuosCache::ClusterStats &stats, double covariate_val) const {
    double n_dbl = static_cast<double>(stats.n);
    double one_plus_nB = 1.0 + n_dbl * B;

    // 1. Compute Posterior Mean mu_n
    double mu_n = (m + B * stats.sum) / one_plus_nB;

    // 2. Compute S_n (Posterior Scale)
    // We compute S_n explicitly here. Note: If you cache S_n in ContinuosCache::ClusterStats, this becomes O(1).
    double S_n;
    if (stats.n == 0) {
        S_n = S0;
    } else {
        double inv_n = 1.0 / n_dbl;
        double xbar = stats.sum * inv_n;
        double ss = stats.sumsq - n_dbl * xbar * xbar; // Centered sum of squares
        double dev = xbar - m;
        S_n = S0 + 0.5 * ss + 0.5 * (n_dbl / one_plus_nB) * dev * dev;
    }

    // 3. Compute the Update Term (Delta S)
    // This is the term added to S_n when x_new is observed, without recomputing the whole sum of squares.
    // Delta S = 0.5 * (x_new - mu_n)^2 * (1+nB) / (1+(n+1)B)
    double one_plus_next_nB = 1.0 + (n_dbl + 1.0) * B;
    double diff = covariate_val - mu_n;
    double delta_S = 0.5 * diff * diff * one_plus_nB / one_plus_next_nB;

    // 4. Compute Log Probability (Student-t)
    double nu_n = nu + 0.5 * n_dbl;

    // Log of Gamma ratio: lgamma(nu_n + 0.5) - lgamma(nu_n)
    double lgamma_diff;
    if (stats.n + 1 < lgamma_nu_n.size()) {
        lgamma_diff = lgamma_nu_n[stats.n + 1] - lgamma_nu_n[stats.n];
    } else {
        lgamma_diff = std::lgamma(nu_n + 0.5) - std::lgamma(nu_n);
    }

    return lgamma_diff - 0.5 * std::log(2.0 * M_PI) -
           0.5 * std::log(one_plus_next_nB / one_plus_nB)  // Scale inflation term
           - 0.5 * std::log(S_n)                           // Normalization by current scale
           - (nu_n + 0.5) * std::log(1.0 + delta_S / S_n); // The "kernel" of the T-distribution
}