/**
 * @file DP.cpp
 * @brief Implementation of Dirichlet Process for Bayesian nonparametric
 * clustering
 *
 * This file contains the implementation of the Dirichlet Process class, which
 * provides prior computations for clustering with the classic stick-breaking
 * construction. The DP is fundamental for nonparametric mixture models with
 * automatic cluster discovery.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "DP.hpp"

double DP::gibbs_prior_existing_cluster(int cls_idx, int obs_idx) const {
  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation (unused in this
   * implementation).
   * @return The log prior probability of assigning the data point to its
   * current cluster.
   */

  const int cluster_size = data.get_cluster_size(cls_idx);
  return (cluster_size > 0) ? log(cluster_size) : std::numeric_limits<double>::lowest();
}

Eigen::VectorXd DP::gibbs_prior_existing_clusters(int obs_idx) const {
  /**
   * @brief Computes the log prior probabilities of assigning a data point to
   * all existing clusters.
   * @param obs_idx The index of the observation (unused in this
   * implementation).
   * @return A vector of log prior probabilities for all existing clusters.
   */

  Eigen::VectorXd log_priors(data.get_K());
  for (int k = 0; k < data.get_K(); ++k) {
      int cluster_size = data.get_cluster_size(k);
      log_priors(k) = (cluster_size > 0) ? log(cluster_size) : std::numeric_limits<double>::lowest();
  }
  return log_priors;
}

double DP::gibbs_prior_new_cluster() const {
  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   * @return The log prior probability of assigning the data point to a new
   * cluster.
   */
  return log_a;
}

double DP::prior_ratio_split(int ci, int cj) const {
  /**
   * @brief Computes the prior ratio for a split operation in a split-merge MCMC
   * algorithm.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */

  const int n_ci = data.get_cluster_size(ci);
  const int n_cj = data.get_cluster_size(cj);

  return log_a - lgamma(n_ci + n_cj) + lgamma(n_ci) + lgamma(n_cj);
}

double DP::prior_ratio_merge(int size_old_ci, int size_old_cj) const {
  /**
   * @brief Computes the prior ratio for a merge operation in a split-merge MCMC
   * algorithm.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */

  const int size_merge = size_old_ci + size_old_cj;

  double log_acceptance_ratio = -log_a;
  log_acceptance_ratio += lgamma(size_merge);
  log_acceptance_ratio -= lgamma(size_old_ci);
  log_acceptance_ratio -= lgamma(size_old_cj);

  return log_acceptance_ratio;
}

double DP::prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci,
                               int cj) const {
  /**
   * @brief Computes the prior ratio for a shuffle operation in a split-merge
   * MCMC algorithm.
   * @param size_old_ci The size of the first cluster before the shuffle.
   * @param size_old_cj The size of the second cluster before the shuffle.
   * @param ci The first cluster index involved in the shuffle.
   * @param cj The second cluster index involved in the shuffle.
   * @return The log prior ratio for the shuffle operation.
   */

  const int n_ci = data.get_cluster_size(ci);
  const int n_cj = data.get_cluster_size(cj);

  double log_acceptance_ratio = 0.0;
  log_acceptance_ratio += lgamma(n_ci);
  log_acceptance_ratio += lgamma(n_cj);
  log_acceptance_ratio -= lgamma(size_old_ci);
  log_acceptance_ratio -= lgamma(size_old_cj);

  return log_acceptance_ratio;
}
