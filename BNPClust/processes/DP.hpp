/**
 * @file DP.hpp
 * @brief Dirichlet Process implementation for Bayesian nonparametric
 * clustering.
 */

#pragma once

#include "../params/DP-params.hpp"
#include "../utils/Process.hpp"

/**
 * @class DP
 * @brief Dirichlet Process class for Bayesian nonparametric clustering.
 *
 * This class implements a Dirichlet Process (DP) that inherits from the Process
 * base class. It provides methods for Gibbs sampling and split-merge algorithms
 * used in clustering.
 */
class DP : public Process {
  const DP_params &params;
  const double log_a;

public:
  /**
   * @brief Constructor for the Dirichlet Process.
   * @param d Reference to the data object.
   * @param p Reference to the parameters object.
   */
  DP(Data &d, DP_params &p)
      : Process(d), params(p), log_a(std::log(params.a)) {};

  /**
   * @name Gibbs Sampling Methods
   * @{
   */

  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation (default: 0).
   * @return The log prior probability of assigning the data point to the
   * existing cluster.
   */
  [[nodiscard]] double
  gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) const override;

  /**
   * @brief Computes the log prior probabilities of assigning a data point to
   * all existing clusters.
   * @param obs_idx The index of the observation.
   * @return A vector of log prior probabilities for all existing clusters.
   */
  [[nodiscard]] Eigen::VectorXd
  gibbs_prior_existing_clusters(int obs_idx) const override;

  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   * @return The log prior probability of assigning the data point to a new
   * cluster.
   */
  [[nodiscard]] double gibbs_prior_new_cluster() const override;

  /** @} */

  /**
   * @name Split-Merge Algorithm Methods
   * @{
   */

  /**
   * @brief Computes the prior ratio for a split operation in a split-merge MCMC
   * algorithm.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */
  [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;

  /**
   * @brief Computes the prior ratio for a merge operation in a split-merge MCMC
   * algorithm.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */
  [[nodiscard]] double prior_ratio_merge(int size_old_ci,
                                         int size_old_cj) const override;

  /**
   * @brief Computes the prior ratio for a shuffle operation in a split-merge
   * MCMC algorithm.
   * @param size_old_ci The size of the first cluster before the shuffle.
   * @param size_old_cj The size of the second cluster before the shuffle.
   * @param ci The first cluster index involved in the shuffle.
   * @param cj The second cluster index involved in the shuffle.
   * @return The log prior ratio for the shuffle operation.
   */
  [[nodiscard]] double prior_ratio_shuffle(int size_old_ci, int size_old_cj,
                                           int ci, int cj) const override;

  /** @} */

  /**
   * @brief Updates the parameters of the Dirichlet Process.
   *
   * This is a null implementation as the Dirichlet Process has no parameters to
   * update.
   */
  void update_params() override { return; };
};
