/**
 * @file continuos_covariate_module.hpp
 * @brief Covariate-related computations for clustering processes.
 */

#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Module.hpp"
#include "Eigen/Dense"
#include <cmath>
#include <vector>

/**
 * @class ContinuosCovariatesModule
 * @brief Module for covariate-related computations within clustering processes.
 *
 * This class implements the product partition model with regression on covariates
 * as described in Müller et al. (2011). It computes similarity measures based on
 * how well observations within a cluster can be explained by a common covariate
 * distribution (Normal conjugate prior).
 *
 * Reference: Müller, P., Quintana, F. (2011)
 * "A Product Partition Model With Regression on Covariates"
 */
class ContinuosCovariatesModule : public Module {
protected:
    /**
     * @name Module References
     * @{
     */

    /** @brief Reference to data object with cluster assignments */
    const Data &data;

    /** @} */

    /**
     * @name Data used
     * @{
     */
    const Eigen::VectorXd continuos_covariate_data; ///< Covariate values
    const bool fixed_v;                             ///< Whether observation variance is fixed (NN) or random (NNIG)
    const double m;                                 ///< Prior mean for covariate
    const double B;                                 ///< Prior variance for covariate
    const double v;                                 ///< Observation variance for covariate
    const double nu;                                ///< Prior shape parameter for variance (NNIG)
    const double S0;                                ///< Prior scale parameter for variance (NNIG)

    /** @} */

    /**
     * @name Helper Methods
     * @{
     */

    /**
     * @brief Sufficient statistics for covariate likelihood computations.
     */
    struct ClusterStats {
        int n = 0;
        double sum = 0.0;
        double sumsq = 0.0;
    };

    /**
     * @brief Compute cluster statistics for covariate similarity
     * @param obs Vector of observation indices in the cluster
     */
    ClusterStats compute_cluster_statistics(const Eigen::Ref<const Eigen::VectorXi> obs) const;

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal-Normal conjugate prior model:
     * - x_i ~ N(μ_j, v)  for i ∈ S_j
     * - Prior on mean μ_j: N(m, B)
     * - Observation variance v is known and fixed
     *
     * @details The marginal likelihood integrates out the cluster-specific mean μ_j.
     * With sufficient statistics:
     * - n_j = |S_j| (cluster size)
     * - x̄_j = (1/n_j) Σ_{i ∈ S_j} x_i (sample mean)
     * - SS = Σ_{i ∈ S_j} (x_i - x̄_j)² (centered sum of squares)
     *
     * The posterior distribution of μ_j is N(m̂_j, τ_j) where:
     * - τ_j = Bv / (v + n_j B)  (posterior variance)
     * - m̂_j = τ_j (n_j x̄_j / v + m / B)  (posterior mean)
     *
     * The log marginal likelihood is:
     *
     * log q(x_j) = -n_j/2 log(2π) - n_j/2 log(v) - 1/2 log(B) + 1/2 log(τ_j)
     *              - SS/(2v) - n_j(x̄_j - m)² / (2(v + n_j B))
     *
     * where log(τ_j) = log(B) + log(v) - log(v + n_j B)
     *
     *
     * @param stats Sufficient statistics for the cluster
     * @return Log marginal likelihood value
     *
     * @note This is marked as __attribute__((hot)) for performance optimization
     *       as it is called frequently in the MCMC sampling loop.
     */
    double compute_log_marginal_likelihood_NN(const ClusterStats &stats) const __attribute__((hot));

    /**
     * @brief Compute log marginal likelihood for cluster given covariates
     *
     * Implements the Normal-InverseGamma conjugate prior model:
     * x ~ N(μ, v_j)
     * - Prior on mean μ: N(m, B*v_j)
     * - Prior on variance v_j: IG(nu, S0)
     *
     * @param stats Sufficient statistics (n, sum, sum of squares)
     * @return Log marginal likelihood contribution
     *
     * @details The marginal likelihood for the NNIG model is:
     * log g(S) = log Γ(ν + n/2) - log Γ(ν) - n/2 log(2π)
     *            - 1/2 log(1 + nB) + ν log(S₀)
     *            - (ν + n/2) log(S₀ + SS/2 + n/(2(1+nB)) (x̄-m)²)
     */
    double compute_log_marginal_likelihood_NNIG(const ClusterStats &stats) const __attribute__((hot));

    /**
     * @brief Compute log predictive density for a new observation (Normal-Normal model)
     *
     * Computes the probability of observing the value at `obs_idx` given the
     * current cluster statistics, assuming the Normal-Normal conjugate prior (fixed variance).
     *
     * @param stats Sufficient statistics of the cluster (n, sum, sum of squares)
     * @param covariate_val Covariate value of the new observation
     * @return Log predictive density log p(x_new | x_cluster)
     *
     * @details The predictive distribution for the NN model is a Normal distribution:
     * x_new | x_cluster ~ N(μ_n, σ²_pred)
     *
     * Where:
     * - Posterior mean: μ_n = (m + nB x̄) / (1 + nB)
     * - Predictive variance: σ²_pred = v * (1 + (n+1)B) / (1 + nB)
     */
    double compute_predictive_NN(const ClusterStats &stats, double covariate_val) const;

    /**
     * @brief Compute log predictive density for a new observation (NNIG model)
     *
     * Computes the probability of observing the value at `obs_idx` given the
     * current cluster statistics, assuming the Normal-Normal-Inverse-Gamma conjugate prior.
     *
     * @param stats Sufficient statistics of the cluster (n, sum, sum of squares)
     * @param covariate_val Covariate value of the new observation
     * @return Log predictive density log p(x_new | x_cluster)
     *
     * @details The predictive distribution for the NNIG model is a non-standardized
     * Student-t distribution:
     * x_new | x_cluster ~ t(df=2ν_n, loc=μ_n, scale=S_n * ratio)
     *
     * Where:
     * - Degrees of freedom: 2ν_n = 2ν + n
     * - Location: μ_n = (m + nB x̄) / (1 + nB)
     * - Scale is derived from the posterior scale S_n and the variance inflation factor.
     */
    double compute_predictive_NNIG(const ClusterStats &stats, double covariate_val) const;

    /**
     * @brief Compute log marginal likelihood based on model type
     *
     * Chooses between NN and NNIG models based on covariates_data.fixed_v.
     *
     * @param stats Sufficient statistics for the cluster
     * @return Log marginal likelihood value
     */
    inline double compute_log_marginal_likelihood(const ClusterStats &stats) const __attribute__((hot)) {
        if (fixed_v) {
            return compute_log_marginal_likelihood_NN(stats);
        } else {
            return compute_log_marginal_likelihood_NNIG(stats);
        }
    }

    /**
     * @brief Compute log marginal likelihood based on model type
     *
     * Chooses between NN and NNIG models based on covariates_data.fixed_v.
     *
     * @param stats Sufficient statistics for the cluster
     * @param covariate_val Covariate value of the new observation
     * @return Log marginal likelihood value
     */
    inline double compute_log_predictive_likelihood(const ClusterStats &stats, double covariate_val) const
        __attribute__((hot, always_inline)) {
        if (fixed_v) {
            return compute_predictive_NN(stats, covariate_val);
        } else {
            return compute_predictive_NNIG(stats, covariate_val);
        }
    }

    /** @} */

    /**
     * @name Precomputed values
     * @{
     */

    const double Bv; ///< Product of prior variance and observation variance

    const double log_B; ///< Log of prior variance
    const double log_v; ///< Log of observation variance

    const double const_term; ///< Constant term in log likelihood

    const double lgamma_nu; ///< log Gamma(ν) for NNIG model (v ~ IG(ν, S₀))
    const double nu_logS0;  ///< ν log(S₀) for NNIG model (v ~ IG(ν, S₀))

    std::vector<double> log_v_plus_nB; ///< Cache for log(v_plus_nB) for NN
    std::vector<double> lgamma_nu_n;   ///< Cache for lgamma(nu_n) for NNIG

    /** @} */

public:
    /**
     * @brief Constructor for ContinuosCovariatesModule
     *
     * @param data_ Reference to Data object with cluster assignments
     * @param covariates_data_ Vector of covariate values for observations (1D array)
     * @param fixed_v_ Whether observation variance is fixed (NN) or random (NNIG)
     * @param m_ Prior mean for covariate
     * @param B_ Prior variance for covariate
     * @param v_ Observation variance for covariate
     * @param nu_ Prior shape parameter for variance (NNIG)
     * @param S0_ Prior scale parameter for variance (NNIG)
     * @param old_alloc_provider function to access old allocations
     * @param old_cluster_members_provider_ function to access old cluster members
     */
    ContinuosCovariatesModule(const Data &data_, const Eigen::VectorXd covariates_data_, bool fixed_v_, double m_ = 0,
                              double B_ = 1.0, double v_ = 1.0, double nu_ = 1.0, double S0_ = 1.0,
                              const Eigen::VectorXi *old_alloc_provider = nullptr,
                              const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : continuos_covariate_data(covariates_data_), fixed_v(fixed_v_), m(m_), B(B_), v(v_), nu(nu_), S0(S0_),
          data(data_), Module(old_alloc_provider, old_cluster_members_provider_), Bv(B * v), log_B(std::log(B)),
          log_v(std::log(v)), const_term(-0.5 * std::log(2.0 * M_PI)), lgamma_nu(std::lgamma(nu)),
          nu_logS0(nu * std::log(S0)) {
        // Precompute caches for efficiency if needed
        if (fixed_v) {
            log_v_plus_nB.reserve(data_.get_n() + 1);
            for (int n = 0; n <= data_.get_n(); ++n) {
                log_v_plus_nB.push_back(std::log(v + n * B));
            }
        } else {
            lgamma_nu_n.reserve(data_.get_n() + 1);
            for (int n = 0; n <= data_.get_n(); ++n) {
                lgamma_nu_n.push_back(std::lgamma(nu + 0.5 * static_cast<double>(n)));
            }
        }
    }

    /**
     * @name Similarity Computation Methods
     * @{
     */

    /**
     * @brief Compute covariate similarity contribution for a cluster
     *
     * Computes the log marginal likelihood of the covariates within a cluster
     * under the Normal conjugate model. Higher values indicate
     * that observations in the cluster have similar covariate values.
     *
     * @param cls_idx Index of the cluster (0 to K-1)
     * @param old_allo If true, uses old allocations from old_allocations_provider;
     *                 if false, uses current allocations from data (default: false)
     * @return Log marginal likelihood contribution (similarity score)
     *
     * @details The computation follows Müller et al. (2011):
     * 1. Compute sufficient statistics (n, sum, sum of squares)
     * 2. Update hyperparameters using conjugate update rules
     * 3. Compute log marginal likelihood using updated parameters
     *
     * This value is added to the clustering prior in split-merge moves
     * to encourage clusters with homogeneous covariate values.
     */
    double compute_similarity_cls(int cls_idx, bool old_allo = false) const override __attribute__((hot));

    /**
     * @brief Compute covariate similarity for a single observation in a cluster
     *
     * Computes the predictive contribution when adding observation obs_idx to
     * cluster cls_idx, considering the covariate values.
     *
     * @param obs_idx Index of the observation
     * @param cls_idx Index of the cluster
     * @return Log predictive density contribution
     *
     * @details Used in Gibbs sampling to compute the probability of assigning
     * an observation to a cluster based on covariate similarity.
     */
    double compute_similarity_obs(int obs_idx, int cls_idx) const override __attribute__((hot));

    /**
     * @brief Compute covariate similarity contributions for all existing clusters
     *
     * Computes the predictive contributions for adding observation obs_idx
     * to each existing cluster, considering covariate values.
     *
     * @param obs_idx Index of the observation
     * @return Vector of log predictive density contributions for each cluster
     */
    Eigen::VectorXd compute_similarity_obs(int obs_idx) const override __attribute__((hot));

    /** @} */
};