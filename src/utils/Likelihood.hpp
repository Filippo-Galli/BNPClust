/**
 * @file Likelihood.hpp
 * @brief Abstract base class for likelihood computation in clustering models
 *
 * This file defines the Likelihood abstract base class that provides the
 * interface for computing cluster and point-level log-likelihoods in Bayesian
 * nonparametric clustering models.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "Data.hpp"

/**
 * @class Likelihood
 * @brief Abstract base class for likelihood computation
 *
 * This class defines the interface for computing likelihoods in clustering
 * models. Derived classes must implement methods for computing both
 * cluster-level and point-level conditional log-likelihoods, which are
 * essential for Gibbs and split-merge MCMC algorithms.
 */
class Likelihood {
protected:
  const Data &data; ///< Reference to Data object with distances and allocations

public:
  Likelihood(const Data &data) : data(data) {}

  /**
   * @brief Computes the log-likelihood for a cluster
   * @param cluster_index Index of the cluster to evaluate
   * @return Total log-likelihood of the cluster
   * @note Useful for split-merge algorithms
   */

  virtual double cluster_loglikelihood(int cluster_index) const = 0;

  /**
   * @brief Computes the log-likelihood for a cluster with given assignments
   * @param cluster_index Index of the cluster to evaluate
   * @param cls_ass_k Vector of point indices in the cluster
   * @return Total log-likelihood of the cluster
   * @note Useful for split-merge algorithms
   */

  virtual double cluster_loglikelihood(
      int cluster_index,
      const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const = 0;

  /**
   * @brief Conditional log-likelihood of a point in a particular cluster
   * @param point_index Index of the point to evaluate
   * @param cluster_index Index of the cluster
   * @return Conditional log-likelihood
   * @note Useful for Gibbs sampling
   */

  virtual double point_loglikelihood_cond(int point_index,
                                          int cluster_index) const = 0;

  virtual ~Likelihood() = default;
};
