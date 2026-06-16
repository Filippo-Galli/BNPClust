#pragma once

/**
 * @file categorical_covariate_module.hpp
 * @brief Covariate-related computations for clustering processes.
 */

#include "../../utils/Data.hpp"
#include "../../utils/Module.hpp"
#include <numeric>

/**
 * @class CategoricalCovariatesModule
 * @brief Module for categorical covariate-related computations within clustering processes.
 *
 * This class implements covariate-related computations for clustering processes
 * where covariates are categorical ({1, ..., C}). It computes similarity measures based on
 * the distribution of categorical covariates within clusters using a Dirichlet-Multinomial model.
 * @details alphas correspond to the concentration parameters for each category.
 */

class CategoricalCovariatesModule : public Module {
protected:
    /**
     * @name Module References
     * @{
     */

    const std::vector<double> prior_alpha; /**< Prior alpha parameters for Dirichlet-Multinomial model */

    /** @brief Reference to data object with cluster assignments */
    const Data &data;

    const Eigen::VectorXi categorical_covariate_data; /**< Categorical covariate values for observations */

    /** @} */

    /**
     * @name Precomputed Values
     * @{
     */

    const double alpha_0 = std::accumulate(prior_alpha.begin(), prior_alpha.end(), 0.0); /**< Sum of prior alphas */
    const double lgamma_alpha_0 = std::lgamma(alpha_0); /**< Log Gamma of sum of prior alphas */

    double prod_lgamma_prior = 0; /**< Sum of log Gamma of prior alphas  */

    /** @} */

public:
    /**
     * @brief Constructor for CategoricalCovariatesModule
     * @param data_ Reference to data object with cluster assignments
     * @param categorical_covariate Categorical covariate values for observations
     * @param prior_alpha_ Prior alpha parameters for Dirichlet-Multinomial model
     * @param old_alloc_provider Optional pointer to previous allocation provider
     * @param old_cluster_members_provider_ Optional pointer to previous cluster members provider
     */
    CategoricalCovariatesModule(
        const Data &data_, const Eigen::VectorXi categorical_covariate, std::vector<double> prior_alpha_,
        const Eigen::VectorXi *old_alloc_provider = nullptr,
        const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : prior_alpha(prior_alpha_), data(data_), categorical_covariate_data(categorical_covariate),
          Module(old_alloc_provider, old_cluster_members_provider_) {

        for (const double &alpha : prior_alpha)
            prod_lgamma_prior += std::lgamma(alpha);
    }
    /**
     * @name Similarity Computation Methods
     * @{
     */

    /**
     * @brief Compute covariate similarity contribution for a cluster
     *
     * Computes the log marginal likelihood of the covariates within a cluster
     * under the Dirichlet-Multinomial conjugate model. Higher values indicate
     * that observations in the cluster have similar covariate values.
     * ...
     * @details The computation follows Müller et al. (2011):
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
