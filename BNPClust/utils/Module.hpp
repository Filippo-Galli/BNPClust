#pragma once

/**
 * @file Module.hpp
 * @brief Base class for modules used in processes.
 */

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

/**
 * @class Module
 * @brief Abstract base class for modules used in clustering processes.
 * This class defines the interface for modules that compute similarity
 * contributions based on covariates or other factors.
 */

class Module {
protected:
    /**
     * @brief Provider function for accessing old allocation state.
     *
     * Used when computing similarity based on previous cluster assignments
     * (e.g., in split-merge algorithms).
     */
    const Eigen::VectorXi *old_allocations_provider;

    /** @brief Provider function for accessing old cluster members map */
    const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider;

public:
    Module(const Eigen::VectorXi *old_allocations_provider_ = nullptr,
           const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : old_allocations_provider(old_allocations_provider_),
          old_cluster_members_provider(old_cluster_members_provider_) {}

    void set_old_allocations_provider(const Eigen::VectorXi *provider) { old_allocations_provider = provider; }

    void set_old_cluster_members_provider(const std::unordered_map<int, std::vector<int>> *provider) {
        old_cluster_members_provider = provider;
    }

    /**
     * @name Similarity Computation Methods
     * @{
     */

    /**
     * @brief Compute similarity contribution for a cluster
     *
     *
     * @param cls_idx Index of the cluster (0 to K-1)
     * @param old_allo If true, uses old allocations from old_allocations_provider;
     *                 if false, uses current allocations from data (default: false)
     * @return Log marginal likelihood contribution (similarity score)
     *
     */
    virtual double compute_similarity_cls(int cls_idx, bool old_allo = false) const = 0;

    /**
     * @brief Compute similarity for a single observation in a cluster
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
    virtual double compute_similarity_obs(int obs_idx, int cls_idx) const = 0;

    /**
     * @brief Compute covariate similarity contributions for all existing clusters
     *
     * Computes the predictive contributions for adding observation obs_idx
     * to each existing cluster, considering covariate values.
     *
     * @param obs_idx Index of the observation
     * @return Vector of log predictive density contributions for each cluster
     */
    virtual Eigen::VectorXd compute_similarity_obs(int obs_idx) const = 0;

    virtual ~Module() = default;
};
