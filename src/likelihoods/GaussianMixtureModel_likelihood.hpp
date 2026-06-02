/**
 * @file GaussianMixtureModel_likelihood.hpp
 * @brief Gaussian Mixture Model likelihood
 */

#pragma once

#include "../params/GaussianMixtureModel-params.hpp"
#include "../utils/Likelihood.hpp"
#include <Eigen/Dense>

/**
 * @class GaussianMixtureModel_likelihood
 * @brief Gaussian Mixture Model likelihood (univariate, collapsed)
 */
class GaussianMixtureModel_likelihood : public Likelihood {
private:
    const GaussianMixtureModel_params &params; ///< Hyperparameters and data
    const Eigen::VectorXd &y;                  ///< Observed values

    // Precomputed constants
    const double log_kappa0;       ///< Log of kappa0
    const double lgamma_alpha0;    ///< Log of gamma(alpha0)
    const double alpha0_log_beta0; ///< Log of alpha0 * beta0
    const double log_pi;           ///< Log of pi

    struct ClusterStats {
        int n = 0;
        double sum = 0.0;
        double sumsq = 0.0;
    };

    /**
     * @brief Compute the statistics for a given cluster assignment.
     *
     * @param cls_ass_k Cluster assignment vector
     * @return ClusterStats struct containing the computed statistics
     */
    ClusterStats compute_stats(const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const;

    /**
     * @brief Compute the log marginal likelihood for a given cluster assignment.
     *
     * @param stats ClusterStats struct containing the computed statistics
     * @return Log marginal likelihood
     */
    double log_marginal_likelihood(const ClusterStats &stats) const;

    /**
     * @brief Compute the log predictive likelihood for a given point and cluster assignment.
     *
     * @param stats ClusterStats struct containing the computed statistics
     * @param x Point value
     * @return Log predictive likelihood
     */
    double log_predictive_from_stats(const ClusterStats &stats, double x) const;

public:
    /**
     * @brief Construct a new GaussianMixtureModel_likelihood object.
     *
     * @param data Reference to Data object with distances and allocations
     * @param params GaussianMixtureModel_params struct containing model parameters
     */
    GaussianMixtureModel_likelihood(const Data &data, const GaussianMixtureModel_params &params);

    /**
     * @brief Destroy the GaussianMixtureModel_likelihood object.
     */
    ~GaussianMixtureModel_likelihood() = default;

    /**
     * @brief Compute the log likelihood of a cluster given its assignment.
     *
     * @param cluster_index Index of the cluster
     * @return Log likelihood of the cluster
     */
    double cluster_loglikelihood(int cluster_index) const override;

    /**
     * @brief Compute the log likelihood of a cluster given its assignment and cluster assignment vector.
     *
     * @param cluster_index Index of the cluster
     * @param cls_ass_k Cluster assignment vector
     * @return Log likelihood of the cluster
     */
    double cluster_loglikelihood(int cluster_index, const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const override;

    /**
     * @brief Compute the log likelihood of a point given its cluster assignment.
     *
     * @param point_index Index of the point
     * @param cluster_index Index of the cluster
     * @return Log likelihood of the point
     */
    double point_loglikelihood_cond(int point_index, int cluster_index) const override;
};
