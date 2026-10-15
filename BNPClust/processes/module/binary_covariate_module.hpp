#pragma once

/**
 * @file binary_covariate_module.hpp
 * @brief Covariate-related computations for clustering processes.
 */

#include "../../utils/Data.hpp"
#include "../../utils/Module.hpp"

/**
 * @class BinaryCovariatesModule
 * @brief Module for binary covariate-related computations within clustering processes.
 *
 * This class implements covariate-related computations for clustering processes
 * where covariates are binary (0/1). It computes similarity measures based on
 * the distribution of binary covariates within clusters using a Beta-Binomial model.
 * @details alpha is associated with the number of successes (1s) and beta with the number of failures (0s).
 */

class BinaryCovariatesModule : public Module {
protected:
    /**
     * @name Module References
     * @{
     */

    const double beta_prior_alpha = 1.0; /**< Prior alpha parameter for Beta-Binomial model */
    const double beta_prior_beta = 1.0;  /**< Prior beta parameter for Beta-Binomial model */

    /** @brief Reference to data object with cluster assignments */
    const Data &data;

    const Eigen::VectorXi binary_covariate_data; /**< Binary covariate values for observations */

    /** @} */

    /**
     * @name Precomputed Values
     * @{
     */

    const double log_beta_prior =
        std::lgamma(beta_prior_alpha) + std::lgamma(beta_prior_beta) -
        std::lgamma(beta_prior_alpha + beta_prior_beta); /**< Log Beta function for prior parameters */
    /** @} */

public:
    /**
     * @brief Constructor for BinaryCovariatesModule
     * @param data_ Reference to data object with cluster assignments
     * @param binary_covariate Covariate values (1D vector of 0s and 1s)
     * @param beta_prior_alpha_ Prior alpha parameter for Beta-Binomial model
     * @param beta_prior_beta_ Prior beta parameter for Beta-Binomial model
     * @param old_alloc_provider Optional pointer to old allocations for split-merge moves
     * @param old_cluster_members_provider_ Optional pointer to old cluster members for split-merge moves
     */
    BinaryCovariatesModule(const Data &data_, const Eigen::VectorXi binary_covariate, double beta_prior_alpha_,
                           double beta_prior_beta_, const Eigen::VectorXi *old_alloc_provider = nullptr,
                           const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : beta_prior_alpha(beta_prior_alpha_), beta_prior_beta(beta_prior_beta_), data(data_),
          binary_covariate_data(binary_covariate), Module(old_alloc_provider, old_cluster_members_provider_) {}

    /**
     * @name Similarity Computation Methods
     * @{
     */

    /**
     * @brief Compute covariate similarity contribution for a cluster
     *
     * Computes the log marginal likelihood of the covariates within a cluster
     * under the Beta-Binomial conjugate model. Higher values indicate
     * that observations in the cluster have similar covariate values.
     * ...
     * @details The computation follows Müller et al. (2011):
     * 1. Count successes (1s) and failures (0s)
     * 2. Compute log marginal likelihood using log-gamma functions
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