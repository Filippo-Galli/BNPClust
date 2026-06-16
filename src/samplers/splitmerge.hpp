/**
 * @file splitmerge.hpp
 * @brief Split-Merge MCMC sampler implementation for Bayesian nonparametric
 * models
 *
 * This file contains the implementation of the split-merge algorithm, an
 * advanced MCMC method that uses joint updates through split, merge, and
 * shuffle moves. The algorithm provides better mixing than sequential methods
 * for complex clustering structures.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"

/**
 * @brief Split-Merge sampler for Bayesian nonparametric mixture models
 *
 * This class implements the split-merge algorithm, an advanced MCMC method that
 * proposes joint updates to cluster assignments through split and merge moves.
 * This approach can achieve better mixing than sequential Gibbs sampling,
 * particularly for models with strong within-cluster correlations.
 *
 * @details The algorithm alternates between three types of moves:
 * - **Split moves**: Divide a cluster into two subclusters
 * - **Merge moves**: Combine two clusters into one
 * - **Shuffle moves**: Redistribute observations between two existing clusters
 *
 * Each move uses restricted Gibbs sampling to generate proposals, followed by
 * Metropolis-Hastings acceptance/rejection based on prior and likelihood
 * ratios. The algorithm maintains detailed balance and ergodicity.
 *
 * @note
 * reference Jain, S. and Neal, R. M. (2004). "A Split-Merge Markov Chain Monte
 * Carlo Procedure for the Dirichlet Process Mixture Model"
 * reference Martinez, A. F. and Mena, R. H. (2014). "On a Nonparametric Change
 * Point Detection Model in Markovian Regimes"
 *
 * @see Sampler, SplitMerge_SAMS
 */
class SplitMerge : public Sampler {
private:
  // ========== Random Number Generation ==========

  /** @brief Mersenne Twister random number generator for sampling operations */
  mutable std::mt19937 gen;

  // ========== Move Selection Variables ==========

  /** @brief Index of first randomly chosen observation */
  int idx_i;

  /** @brief Index of second randomly chosen observation */
  int idx_j;

  /** @brief Cluster assignment of first observation */
  int ci;

  /** @brief Cluster assignment of second observation */
  int cj;

  // ========== Algorithm Configuration ==========

  /** @brief Flag to enable shuffle moves (Mena and Martinez, 2014) */
  bool shuffle_bool = false;

  // ========== State Management ==========

  /** @brief Launch state for restricted Gibbs sampling */
  Eigen::VectorXi launch_state;

  /** @brief Indices of observations in clusters ci and cj */
  Eigen::VectorXi S;

  /** @brief Original cluster assignments before move proposal */
  Eigen::VectorXi original_allocations;

  // ========== Proposal Probabilities ==========

  /** @brief Log probability of generating current state via restricted Gibbs
   * (split direction) */
  double log_split_gibbs_prob = 0;

  /** @brief Log probability of generating current state via restricted Gibbs
   * (merge direction) */
  double log_merge_gibbs_prob = 0;

  // ========== Debug variables ==========
  int accepted_split = 0;
  int accepted_merge = 0;
  int accepted_shuffle = 0;

  // ========== Move Selection Methods ==========

  /**
   * @brief Randomly select two observations for split-merge proposal
   *
   * @details Uniformly samples two distinct observation indices that will
   * be used to determine the type of move (split, merge, or shuffle).
   */
  void choose_indeces();

  /**
   * @brief Select clusters for shuffle move
   *
   * @details Determines the cluster assignments of the selected observations
   * and prepares for a shuffle move between two existing clusters.
   */
  void choose_clusters_shuffle();

  // ========== Proposal Generation ==========

  /**
   * @brief Generate proposal state via restricted Gibbs sampling
   *
   * @param iterations Number of restricted Gibbs iterations to perform
   * @param only_probabilities If true, only compute proposal probabilities
   * without updating state
   *
   * @details Performs restricted Gibbs sampling on observations in the selected
   * clusters to generate a proposal state. Also computes the probability of
   * generating this state for use in the acceptance ratio.
   */
  void restricted_gibbs(int iterations, bool only_probabilities = false);

  // ========== Split Move Implementation ==========

  /**
   * @brief Execute a split move proposal
   *
   * @details Attempts to split a cluster containing both selected observations
   * into two separate clusters using restricted Gibbs sampling to generate
   * the proposal allocation.
   */
  void split_move();

  /**
   * @brief Compute acceptance ratio for split move
   *
   * @param likelihood_old_cluster Likelihood of the original single cluster
   * @return Log acceptance ratio for the split proposal
   *
   * @details Computes the Metropolis-Hastings acceptance probability by
   * combining prior ratios, likelihood ratios, and proposal probabilities.
   */
  double compute_acceptance_ratio_split(double likelihood_old_cluster);

  // ========== Merge Move Implementation ==========

  /**
   * @brief Execute a merge move proposal
   *
   * @details Attempts to merge two clusters containing the selected
   * observations into a single cluster.
   */
  void merge_move();

  /**
   * @brief Compute acceptance ratio for merge move
   *
   * @param likelihood_old_ci Likelihood of first original cluster
   * @param likelihood_old_cj Likelihood of second original cluster
   * @return Log acceptance ratio for the merge proposal
   */
  double compute_acceptance_ratio_merge(double likelihood_old_ci,
                                        double likelihood_old_cj);

  // ========== Shuffle Move Implementation ==========

  /**
   * @brief Execute a shuffle move proposal
   *
   * @details Redistributes observations between two existing clusters without
   * changing the total number of clusters.
   */
  void shuffle();

  /**
   * @brief Compute acceptance ratio for shuffle move
   *
   * @param likelihood_old_ci Likelihood of first cluster before shuffle
   * @param likelihood_old_cj Likelihood of second cluster before shuffle
   * @param old_ci_size Size of first cluster before shuffle
   * @param old_cj_size Size of second cluster before shuffle
   * @return Log acceptance ratio for the shuffle proposal
   */
  double compute_acceptance_ratio_shuffle(double likelihood_old_ci,
                                          double likelihood_old_cj,
                                          int old_ci_size, int old_cj_size);

public:
  // ========== Constructor ==========

  /**
   * @brief Constructor for Split-Merge sampler
   *
   * @param d Reference to Data object containing observations
   * @param l Reference to Likelihood object for probability computations
   * @param pr Reference to Process object defining the prior
   * @param shuffle Flag to enable shuffle moves in addition to split-merge
   *
   * @details Initializes the split-merge sampler with the option to include
   * shuffle moves. When shuffle is enabled, the algorithm can propose
   * redistributions between existing clusters in addition to split-merge moves.
   */
  SplitMerge(Data &d, Likelihood &l, Process &pr, bool shuffle)
      : Sampler(d, l, pr), shuffle_bool(shuffle), gen(rd()) {};

  // ========== MCMC Interface ==========

  /**
   * @brief Perform one iteration of the split-merge algorithm
   *
   * @details Executes one step of the split-merge sampler:
   * 1. Randomly select two observations
   * 2. Determine move type based on current cluster assignments
   * 3. Generate proposal via restricted Gibbs sampling
   * 4. Compute acceptance ratio and accept/reject via Metropolis-Hastings
   *
   * The algorithm automatically chooses between split, merge, and shuffle moves
   * based on the cluster assignments of the selected observations.
   */
  void step() override;

  // ========== Accessor Methods ==========
  /**
   * @brief Get number of accepted split moves for diagnostics
   * @return Count of accepted split moves
   */
  int get_accepted_split() const { return accepted_split; };

  /**
   * @brief Get number of accepted merge moves for diagnostics
   * @return Count of accepted merge moves
   */
  int get_accepted_merge() const { return accepted_merge; };

  /**
   * @brief Get number of accepted shuffle moves for diagnostics
   * @return Count of accepted shuffle moves
   */
  int get_accepted_shuffle() const { return accepted_shuffle; };
};
