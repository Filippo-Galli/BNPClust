/**
 * @file neal.hpp
 * @brief Neal's Algorithm 3 implementation for collapsed Gibbs sampling
 *
 * This file contains the implementation of Neal's Algorithm 3, a standard
 * collapsed Gibbs sampler for Bayesian nonparametric mixture models. The
 * algorithm provides efficient sequential updates for cluster assignments with
 * automatic cluster creation and deletion.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"
#include <vector>

/**
 * @brief Implementation of Neal's Algorithm 3 for collapsed Gibbs sampling
 *
 * This class implements Neal's Algorithm 3, a standard collapsed Gibbs sampler
 * for Bayesian nonparametric mixture models. The algorithm sequentially updates
 * cluster assignments by sampling from the full conditional distribution of
 * each observation, integrating out the cluster parameters.
 *
 * @details Algorithm 3 is characterized by:
 * - **Sequential updates**: One observation at a time
 * - **Collapsed sampling**: Cluster parameters are integrated out analytically
 * - **Full conditionals**: Each assignment is sampled from its exact posterior
 * - **Automatic cluster creation/deletion**: Clusters are created when needed
 * and deleted when empty
 *
 * The algorithm is particularly effective for models where the likelihood can
 * be computed in closed form after integrating out cluster-specific parameters.
 *
 * @note
 * reference Neal, R. M. (2000). "Markov Chain Sampling Methods for Dirichlet
 * Process Mixture Models"
 *
 * @see Sampler, Process, Likelihood
 */
class Neal3 : public Sampler {
private:
  // ========== Random Number Generation ==========

  /** @brief Mersenne Twister random number generator for sampling operations */
  mutable std::mt19937 gen;

  // ========== Core Algorithm Methods ==========

  /**
   * @brief Sample cluster assignment for a single observation
   *
   * @param index Index of the observation to update
   *
   * @details This method implements the core of Algorithm 3:
   * 1. Remove the observation from its current cluster
   * 2. Compute probabilities for all existing clusters plus a new cluster
   * 3. Sample new assignment from the full conditional distribution
   * 4. Update cluster assignments and clean up empty clusters
   *
   * The probabilities combine prior information from the Process with
   * likelihood information computed by integrating out cluster parameters.
   */
  void step_1_observation(int index);

  /**
   * @brief Sample an index from log-probabilities using the log-sum-exp trick
   *
   * @param num_clusters Number of clusters (including new cluster)
   * @return Index sampled according to the probabilities
   *
   * @details This method performs sampling from a discrete distribution
   * defined by log-probabilities. It uses the log-sum-exp trick for numerical
   * stability and implements roulette wheel selection.
   */
  int sample_from_log_probs(int num_clusters);

  // Pre-allocated buffers to avoid repeated allocations
  std::vector<int> indices;
  std::vector<double> log_likelihoods;
  std::vector<double> weights;
  int n_data;

public:
  // ========== Constructor ==========

  /**
   * @brief Constructor for Neal's Algorithm 3 sampler
   *
   * @param d Reference to Data object containing observations
   * @param p Reference to Params object with hyperparameters
   * @param l Reference to Likelihood object for probability computations
   * @param pr Reference to Process object (DP, NGGP, etc.) defining the prior
   *
   * @details Initializes the Gibbs sampler with all required components.
   * The random number generator is seeded from the inherited random device.
   * Pre-allocates vectors to maximum possible size to avoid resizing during
   * sampling.
   */
  Neal3(Data &d, Likelihood &l, Process &pr)
      : Sampler(d, l, pr), gen(rd()), n_data(d.get_n()) {
    // Pre-allocate to max size (all points in separate clusters + 1 new
    // cluster)
    log_likelihoods.reserve(n_data + 1);
    weights.reserve(n_data + 1);

    indices.resize(n_data);
    std::iota(indices.begin(), indices.end(), 0);
  }

  // ========== MCMC Interface ==========

  /**
   * @brief Perform one complete iteration of Neal's Algorithm 3
   *
   * @details Executes one full sweep of the collapsed Gibbs sampler by
   * sequentially updating the cluster assignment of each observation.
   * The order of updates is typically randomized to avoid systematic bias.
   *
   * After each full sweep:
   * - All observations have been considered for reassignment
   * - Cluster structure may have changed (clusters created/destroyed)
   * - The Markov chain has advanced by one step
   *
   * @see step_1_observation()
   */
  void step() override;
};
