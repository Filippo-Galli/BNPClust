/**
 * @file splitmerge_SAMS.hpp
 * @brief Sequential Allocation Merge-Split (SAMS) sampler implementation
 *
 * This file contains the implementation of the SAMS algorithm, an optimized
 * variant of split-merge sampling that uses sequential allocation for proposal
 * generation. SAMS provides computational advantages for large datasets while
 * maintaining theoretical properties of split-merge algorithms.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"

/**
 * @brief Sequential Allocation Merge-Split (SAMS) sampler for Bayesian
 * nonparametric models
 *
 * This class implements the SAMS algorithm, a variant of the split-merge
 * sampler that uses sequential allocation instead of restricted Gibbs sampling
 * for proposal generation. SAMS can be more efficient than standard split-merge
 * for certain model configurations and provides an alternative proposal
 * mechanism.
 *
 * @details Key differences from standard Split-Merge:
 * - **Sequential Allocation**: Observations are allocated one by one in random
 * order
 * - **Simplified Proposals**: No need for full restricted Gibbs chains
 * - **Efficient Computation**: Faster proposal generation for large clusters
 * - **Maintained Ergodicity**: Preserves theoretical properties of split-merge
 *
 * The algorithm maintains the same three types of moves (split, merge, shuffle)
 * but uses a different mechanism for generating proposals within each move
 * type.
 *
 * @note
 * reference Dahl, D. B. and Newcomb, S. (2022). "Sequentially allocated
 * merge-split samplers for conjugate Bayesian nonparametric models"
 * reference Martinez, A. F. and Mena, R. H. (2014). "On a Nonparametric Change
 * Point Detection Model in Markovian Regimes"
 *
 * @see Sampler, SplitMerge
 */
class SplitMerge_SAMS : public Sampler {
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

  /** @brief Launch state for sequential allocation */
  Eigen::VectorXi launch_state;

  /** @brief Indices of observations in clusters ci and cj */
  Eigen::VectorXi S;

  /** @brief Original cluster assignments before move proposal */
  Eigen::VectorXi original_allocations;

  // ========== Proposal Probabilities ==========

  /** @brief Log probability of generating current state via sequential
   * allocation (split direction) */
  double log_split_gibbs_prob = 0;

  /** @brief Log probability of generating current state via sequential
   * allocation (merge direction) */
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
   * determine the type of move and the clusters involved.
   */
  void choose_indeces();

  /**
   * @brief Select clusters for shuffle move
   *
   * @details Determines cluster assignments and prepares for redistribution
   * between two existing clusters.
   */
  void choose_clusters_shuffle();

  // ========== SAMS-Specific Proposal Generation ==========

  /**
   * @brief Generate proposal state via sequential allocation
   *
   * @param iterations Number of allocation passes to perform
   * @param only_probabilities If true, only compute proposal probabilities
   * without updating state
   * @param sequential If true, use sequential allocation; if false, use
   * restricted Gibbs sampling
   *
   * @details Implements the core SAMS algorithm by sequentially allocating
   * observations to clusters. Unlike restricted Gibbs sampling, this approach
   * processes observations one by one in random order, making allocation
   * decisions based on current partial assignments.
   */
  void sequential_allocation(int iterations, bool only_probabilities = false,
                             bool sequential = true);

  // ========== Split Move Implementation ==========

  /**
   * @brief Execute a split move using SAMS
   *
   * @details Attempts to split a cluster using sequential allocation to
   * generate the proposal state. The two anchor observations are placed
   * in separate subclusters initially.
   */
  void split_move();

  /**
   * @brief Compute acceptance ratio for SAMS split move
   *
   * @param likelihood_old_cluster Likelihood of the original single cluster
   * @return Log acceptance ratio for the split proposal
   */
  double compute_acceptance_ratio_split(double likelihood_old_cluster);

  // ========== Merge Move Implementation ==========

  /**
   * @brief Execute a merge move using SAMS
   *
   * @details Attempts to merge two clusters into one using sequential
   * allocation to determine the final unified assignment.
   */
  void merge_move();

  /**
   * @brief Compute acceptance ratio for SAMS merge move
   *
   * @param likelihood_old_ci Likelihood of first original cluster
   * @param likelihood_old_cj Likelihood of second original cluster
   * @return Log acceptance ratio for the merge proposal
   */
  double compute_acceptance_ratio_merge(double likelihood_old_ci,
                                        double likelihood_old_cj);

  // ========== Shuffle Move Implementation ==========

  /**
   * @brief Execute a shuffle move using SAMS
   *
   * @details Redistributes observations between two clusters using
   * sequential allocation while maintaining the two-cluster structure.
   */
  void shuffle();

  /**
   * @brief Compute acceptance ratio for SAMS shuffle move
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
   * @brief Constructor for SAMS (Sequential Allocation Merge-Split) sampler
   *
   * @param d Reference to Data object containing observations
   * @param l Reference to Likelihood object for probability computations
   * @param pr Reference to Process object defining the prior
   * @param shuffle Flag to enable shuffle moves in addition to split-merge
   *
   * @details Initializes the SAMS sampler, which uses sequential allocation
   * instead of restricted Gibbs sampling for generating proposals. This can
   * provide computational advantages for certain model configurations.
   */
  SplitMerge_SAMS(Data &d, Likelihood &l, Process &pr, bool shuffle)
      : Sampler(d, l, pr), shuffle_bool(shuffle), gen(rd()) {};

  // ========== MCMC Interface ==========

  /**
   * @brief Perform one iteration of the SAMS algorithm
   *
   * @details Executes one step of the SAMS sampler:
   * 1. Randomly select two observations as anchors
   * 2. Determine move type based on their current assignments
   * 3. Generate proposal using sequential allocation
   * 4. Compute acceptance ratio and accept/reject the proposal
   *
   * The sequential allocation approach can be more efficient than restricted
   * Gibbs sampling while maintaining the theoretical properties of split-merge.
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
