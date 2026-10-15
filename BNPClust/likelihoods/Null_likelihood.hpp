/**
 * @file Null_likelihood.hpp
 * @brief LogLikelihood equal to 1 to ignore data contribution and focus on
 * prior
 */

#pragma once

#include "../utils/Likelihood.hpp"

/**
 * @class Null_likelihood
 * @brief LogLikelihood equal to 1 to ignore data contribution and focus on
 * prior
 */

class Null_likelihood : public Likelihood {
public:
  Null_likelihood(const Data &data) : Likelihood(data) {}
  /**
   * @brief Computes the full log-likelihood for a cluster
   * @param cluster_index Index of the cluster to evaluate
   * @return Total log-likelihood (cohesion + repulsion)
   *
   * This method computes both the within-cluster cohesion and the
   * between-cluster repulsion contributions for the specified cluster.
   */
  double cluster_loglikelihood(int cluster_index) const override final {
    return 0.0;
  };

  /**
   * @brief Computes the full log-likelihood for a cluster with given
   * assignments
   * @param cluster_index Index of the cluster to evaluate
   * @param cls_ass_k Vector of point indices in the cluster
   * @return Total log-likelihood (cohesion + repulsion)
   *
   * This overload allows computing the likelihood with a custom set of
   * cluster assignments without modifying the data structure.
   */
  double cluster_loglikelihood(
      int cluster_index,
      const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const override final
      __attribute__((hot)) {
    return 0.0;
  };

  /**
   * @brief Computes the conditional log-likelihood of a point given a cluster
   * @param point_index Index of the point to evaluate
   * @param cluster_index Index of the target cluster
   * @return Conditional log-likelihood of assigning the point to the cluster
   *
   * This method evaluates how well a point fits into a specific cluster,
   * considering both its cohesion with points in that cluster and its
   * repulsion from points in other clusters.
   */
  double point_loglikelihood_cond(int point_index,
                                  int cluster_index) const override final
      __attribute__((hot)) {
    return 0.0;
  };
};
