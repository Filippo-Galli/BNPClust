/**
 * @file NGGPx.hpp
 * @brief Normalized Generalized Gamma Process with module-based covariates
 *
 * This file defines the NGGPx class that extends the NGGP process to
 * incorporate module-based computations for covariates and other similarity
 * terms.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Module.hpp"
#include "NGGP.hpp"
#include <memory>
#include <vector>

/**
 * @class NGGPx
 * @brief NGGP clustering process with module-based covariates
 *
 * This class extends the Normalized Generalized Gamma Process (NGGP) to
 * incorporate module-based computations for covariates and other similarity
 * terms. It augments Gibbs and split-merge moves with module contributions
 * computed through the Module interface.
 */

class NGGPx : public NGGP {

protected:
  std::vector<std::shared_ptr<Module>> modules;

public:
  /**
   * @brief Construct an `NGGPWx` process.
   * @param d Data container.
   * @param p Model and sampler parameters.
   * @param U_sam Sampler for the latent NGGP auxiliary variable $U$.
   * @param mods Vector of shared pointers to Module objects providing
   * similarity computations.
   */
  NGGPx(Data &d, NGGP_params &p, U_sampler &U_sam,
        const std::vector<std::shared_ptr<Module>> &mods)
      : NGGP(d, p, U_sam), modules(mods) {
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
   * Combines the NGGP term for an existing cluster with spatial weighting and
   * (when applicable in the implementation) covariate similarity contributions.
   *
   * @param cls_idx Cluster index.
   * @param obs_idx Observation index (used for observation-specific
   * spatial/covariate terms).
   * @return Log prior (up to an additive constant) for assigning observation
   * `obs_idx` to existing cluster `cls_idx`.
   */
  [[nodiscard]] double
  gibbs_prior_existing_cluster(int cls_idx, int obs_idx = 0) const override;

  /**
   * @brief Computes the log prior for assigning an observation to each existing
   * cluster.
   *
   * @param obs_idx Observation index.
   * @return Vector of size $K$ containing log priors (up to an additive
   * constant), one per existing cluster.
   */
  [[nodiscard]] Eigen::VectorXd
  gibbs_prior_existing_clusters(int obs_idx) const override;

  /**
   * @brief Computes the log prior probability of assigning a data point to a
   * new cluster.
   *
   * This is the baseline new-cluster mass from the NGGP model (without
   * observation-specific covariate effects).
   *
   * @return Log prior (up to an additive constant) for creating a new cluster.
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
   * @brief Computes the prior ratio for a split operation in a spatially-aware
   * NGGP-based split-merge MCMC algorithm.
   *
   * Accounts for the NGGP prior and spatial/covariate contributions (as
   * implemented by `NGGPW` and `CovariatesModule`).
   *
   * @param ci Index of the first (proposed) split cluster.
   * @param cj Index of the second (proposed) split cluster.
   * @return Log prior ratio contribution for the split move.
   */
  [[nodiscard]] double prior_ratio_split(int ci, int cj) const override;

  /**
   * @brief Computes the prior ratio for a merge operation in a spatially-aware
   * NGGP-based split-merge MCMC algorithm.
   *
   * Accounts for the NGGP prior and spatial/covariate contributions (as
   * implemented by `NGGPW` and `CovariatesModule`).
   *
   * @param size_old_ci Size of the first cluster before the merge.
   * @param size_old_cj Size of the second cluster before the merge.
   * @return Log prior ratio contribution for the merge move.
   */
  [[nodiscard]] double prior_ratio_merge(int size_old_ci,
                                         int size_old_cj) const override;

  /**
   * @brief Computes the prior ratio for a shuffle operation in a
   * spatially-aware NGGP-based split-merge MCMC algorithm.
   *
   * Accounts for the NGGP prior and spatial/covariate contributions (as
   * implemented by `NGGPW` and `CovariatesModule`).
   *
   * @param size_old_ci Size of the first cluster before the shuffle.
   * @param size_old_cj Size of the second cluster before the shuffle.
   * @param ci Index of the first cluster involved in the shuffle.
   * @param cj Index of the second cluster involved in the shuffle.
   * @return Log prior ratio contribution for the shuffle move.
   */
  [[nodiscard]] double prior_ratio_shuffle(int size_old_ci, int size_old_cj,
                                           int ci, int cj) const override;

  /** @} */

  /**
   * @name Parameter Update Methods
   * @{
   */

  /**
   * @brief Update process parameters.
   *
   * Delegates to the configured `U_sampler` implementation to update the latent
   * auxiliary variable $U$ conditional on the current partition.
   *
   * @see U_sampler::update_U()
   */
  void update_params() override { NGGP::U_sampler_method.update_U(); };

  /** @} */
};
