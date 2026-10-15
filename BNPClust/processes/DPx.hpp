/**
 * @file DPx.hpp
 * @brief Dirichlet Process with covariates and modules for Bayesian
 * nonparametric clustering
 *
 * This file defines the DPx class that extends the base Dirichlet Process to
 * incorporate module-based computations for covariates and other similarity
 * terms.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Module.hpp"
#include "DP.hpp"
#include <memory>

/**
 * @class DPx
 * @brief Dirichlet Process with module-based covariates for Bayesian
 * nonparametric clustering
 *
 * This class extends the Dirichlet Process to incorporate covariates and other
 * similarity components through a Module interface. It provides methods for
 * Gibbs sampling and split-merge algorithms that account for module-based
 * dependencies between observations.
 */

class DPx : public DP {

protected:
  std::vector<std::shared_ptr<Module>> modules;

public:
  /**
   * @brief Constructor for the Dirichlet Process with modules.
   * @param d Reference to the data object containing observations and cluster
   * assignments.
   * @param p Reference to the parameters object.
   * @param mods Vector of shared pointers to Module objects providing
   * similarity computations.
   */
  DPx(Data &d, DP_params &p, const std::vector<std::shared_ptr<Module>> &mods)
      : DP(d, p), modules(mods) {
    for (auto &mod : modules) {
      mod->set_old_allocations_provider(&this->old_allocations_view());
      mod->set_old_cluster_members_provider(&this->old_cluster_members_view());
    }
  }

  /**
   * @name Gibbs Sampling Methods
   * @{
   */

  /**
   * @brief Computes the log prior probability of assigning a data point to an
   * existing cluster.
   *
   * This method incorporates module-based similarity contributions when
   * computing the prior probability.
   * @param cls_idx The index of the cluster.
   * @param obs_idx The index of the observation to assign.
   * @return The log prior probability of assigning the data point to the
   * existing cluster.
   */
  [[nodiscard]] double gibbs_prior_existing_cluster(int cls_idx,
                                                    int obs_idx) const override;

  /**
   * @brief Computes the log prior probabilities of assigning a data point to
   * every existing cluster. This method is useful for Gibbs sampling over
   * existing clusters. It returns a vector of log prior probabilities for all
   * existing clusters.
   * @param obs_idx The index of the observation to assign.
   * @return A vector of log prior probabilities for assigning the data point to
   * each existing cluster.
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

  /**
   * @brief Observation-specific new-cluster prior.
   *
   * Extends the baseline new-cluster prior by including the covariate singleton
   * similarity contribution for observation `obs_idx`.
   *
   * @param obs_idx Observation index.
   * @return Log prior (up to an additive constant) for assigning observation
   * `obs_idx` to a newly created cluster.
   */
  [[nodiscard]] double gibbs_prior_new_cluster_obs(int obs_idx) const override;

  /** @} */

  /**
   * @name Split-Merge Algorithm Methods
   * @{
   */

  /**
   * @brief Computes the prior ratio for a split operation in split-merge MCMC.
   *
   * This method accounts for both the Dirichlet Process prior and module-based
   * similarity terms when computing the acceptance ratio for splitting
   * clusters.
   * @param ci The first cluster index involved in the split.
   * @param cj The second cluster index involved in the split.
   * @return The log prior ratio for the split operation.
   */
  [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;

  /**
   * @brief Computes the prior ratio for a merge operation in split-merge MCMC.
   *
   * This method accounts for both the Dirichlet Process prior and module-based
   * similarity terms when computing the acceptance ratio for merging clusters.
   * @param size_old_ci The size of the first cluster before the merge.
   * @param size_old_cj The size of the second cluster before the merge.
   * @return The log prior ratio for the merge operation.
   */
  [[nodiscard]] double prior_ratio_merge(int size_old_ci,
                                         int size_old_cj) const override;

  /**
   * @brief Computes the prior ratio for a shuffle operation in split-merge
   * MCMC.
   *
   * This method accounts for both the Dirichlet Process prior and module-based
   * similarity terms when computing the acceptance ratio for shuffling
   * observations between clusters.
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
   * @brief Updates the module-based parameters.
   *
   * This is a null implementation as DPx has no parameters to update beyond
   * those handled by the base Process class.
   */
  void update_params() override { return; };
};
