/**
 * @file splitmerge_LSS_SDDS.hpp
 * @brief Locality Sensitive Sampling (LSS) with SDDS Split-Merge sampler
 * implementation
 *
 * This file contains the implementation of the LSS Split-Merge algorithm with
 * SDDS (Smart-split, Dumb-merge, Dumb-split, Smart-merge), an optimized variant
 * of split-merge sampling that uses locality sensitive sampling for selecting
 * anchor points. The SDDS strategy adaptively uses smart (with sequential
 * allocation) or dumb (simple random) proposals depending on the move type and
 * point similarity, balancing computational efficiency with mixing quality. LSS
 * provides computational advantages for large datasets by leveraging similarity
 * information while maintaining theoretical properties of split-merge
 * algorithms.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"

/**
 * @brief Locality Sensitive Sampling (LSS) with SDDS Split-Merge sampler
 *
 * This class implements the LSS Split-Merge algorithm with SDDS (Smart-split,
 * Dumb-merge, Dumb-split, Smart-merge), a variant of the split-merge sampler
 * that uses locality sensitive sampling to select anchor points based on
 * similarity/dissimilarity information. The SDDS strategy adaptively chooses
 * between smart (sequential allocation) and dumb (simple random) proposals to
 * balance computational cost with proposal quality:
 * - When dissimilar points are selected: Smart split + Dumb merge
 * - When similar points are selected: Dumb split + Smart merge
 * This improves mixing efficiency while maintaining computational tractability.
 *
 * @details Key differences from standard Split-Merge:
 * - **Locality Sensitive Sampling**: Anchor points are selected based on
 * similarity/dissimilarity derived from distances
 * - **SDDS Strategy**: Adaptive smart/dumb pairing based on point similarity
 *   * Dissimilar points → Smart split (sequential allocation) + Dumb merge
 * (direct)
 *   * Similar points → Dumb split (random) + Smart merge (sequential
 * allocation)
 * - **Sequential Allocation**: Used in "smart" moves for intelligent proposal
 * generation
 * - **Efficient Computation**: Balances computational cost with mixing quality
 * - **Maintained Ergodicity**: Preserves theoretical properties of split-merge
 *
 * The algorithm maintains the same three types of moves (split, merge, shuffle)
 * but uses locality sensitive sampling for anchor point selection and
 * adaptively applies sequential allocation based on the SDDS strategy.
 *
 * @note References:
 * - Luo, C., Shrivastava, A. (2018). "Scaling-up Split-Merge MCMC with
 *   Locality Sensitive Sampling (LSS)"
 * - Dahl, D. B. and Newcomb, S. (2022). "Sequentially allocated merge-split
 *   samplers for conjugate Bayesian nonparametric models"
 * - Martinez, A. F. and Mena, R. H. (2014). "On a Nonparametric Change Point
 *   Detection Model in Markovian Regimes"
 *
 * @note [IMPORTANT] since there is the LSS it only works with distance based
 * models
 *
 * @see Sampler, SplitMerge
 */
class SplitMerge_LSS_SDDS : public Sampler {
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

  /** @brief Size of the launch state and S vectors */
  int launch_state_size;

  /** @brief Cached reference to original allocations stored in the process */
  const Eigen::VectorXi &original_allocations;

  // ========== Proposal Probabilities ==========

  /** @brief uniform distribution between 0 and 1
   */
  std::uniform_real_distribution<> dis_real;

  /** @brief uniform integer distribution for general use */
  std::uniform_int_distribution<> dis_int;

  /**
   * @brief Log probabilities of allocating a point to clusters ci and cj
   */
  Eigen::Vector2d log_probs;

  /** @brief Probabilities of allocating a point to clusters ci and cj
   */
  Eigen::Vector2d probs;

  /** @brief Log probability of generating current state via sequential
   * allocation (split direction) */
  double log_split_gibbs_prob = 0;

  /** @brief Log probability of generating current state via sequential
   * allocation (merge direction) */
  double log_merge_gibbs_prob = 0;

  /** @brief Constant log probability for random split allocation (dumb split)
   */
  const double rand_split_prob = log(0.5);

  // ========== Debug variables ==========
  int accepted_split = 0;
  int accepted_merge = 0;
  int accepted_shuffle = 0;
  int split_moves = 0;
  int merge_moves = 0;
  int shuffle_moves = 0;

  // ========== Move Selection Methods ==========

  /**
   * @brief Randomly select two observations for split-merge proposal
   *
   * @details Samples two distinct observation indices using locality sensitive
   * sampling. When similarity=false, selects dissimilar points (for split
   * moves). When similarity=true, selects similar points (for merge moves). The
   * first point is selected uniformly, and the second is selected based on
   * distance weights from the first point.
   *
   * @param similarity If true, select similar points; if false, select
   * dissimilar points
   */
  void choose_indeces(bool similarity = false);

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
   * without modifying allocations (used for computing reverse move
   * probabilities)
   * @param sequential If true, unallocate all points before sequential
   * allocation; if false, use restricted Gibbs sampling (unallocate one at a
   * time)
   *
   * @details Implements sequential allocation by processing observations one by
   * one in random order. Each observation is allocated to cluster ci or cj
   * based on conditional likelihood and prior probabilities. This method
   * computes the log probability of the proposed allocation path, which is used
   * in the acceptance ratio.
   */
  void sequential_allocation(int iterations, bool only_probabilities = false,
                             bool sequential = true);

  // ========== Split Move Implementation ==========

  /**
   * @brief Execute a smart split move using sequential allocation
   *
   * @details Splits a cluster ci into two clusters (ci and cj) using sequential
   * allocation for intelligent proposal generation. Anchor points idx_i and
   * idx_j are placed in separate clusters, other points are randomly
   * initialized, then refined via sequential allocation. The proposal
   * probability is computed and included in the acceptance ratio. Provides
   * better proposals but at higher computational cost than dumb_split_move().
   */
  void smart_split_move();

  /**
   * @brief Execute a dumb split move with random allocation
   *
   * @details Splits a cluster ci into two clusters (ci and cj) by randomly
   * allocating points without sequential refinement. Anchor points idx_i and
   * idx_j are placed in separate clusters, and all other points are randomly
   * assigned to ci or cj with equal probability. No proposal probability is
   * computed (log_split_gibbs_prob = 0). Faster but may have lower acceptance
   * rates than smart_split_move().
   */
  void dumb_split_move();

  /**
   * @brief Compute acceptance ratio for LSS split move
   *
   * @param likelihood_old_cluster Log-likelihood of the original single cluster
   * before split
   * @return Log acceptance ratio for the split proposal
   *
   * @details Computes log acceptance ratio as:
   * log(α) = log(prior_ratio) + log(likelihood_ratio) - log(proposal_ratio)
   * where:
   * - prior_ratio accounts for cluster size changes
   * - likelihood_ratio = L(ci_new) + L(cj_new) - L(ci_old)
   * - proposal_ratio is log_split_gibbs_prob (0 for dumb split)
   */
  double compute_acceptance_ratio_split(double likelihood_old_cluster);

  // ========== Merge Move Implementation ==========

  /**
   * @brief Execute a smart merge move using sequential allocation
   *
   * @details Merges two clusters ci and cj into one (ci) using sequential
   * allocation to determine the final merged state. All points are initially
   * assigned to ci, then sequential allocation computes the probability of
   * this configuration. The proposal probability (log_merge_gibbs_prob) is
   * included in the acceptance ratio. More accurate proposals with higher
   * computational cost than dumb_merge_move().
   */
  void smart_merge_move();

  /**
   * @brief Execute a dumb merge move with direct merging
   *
   * @details Directly merges two clusters ci and cj into one (ci) without
   * sequential allocation. All points in cj are simply reassigned to ci.
   * No proposal probability is computed (log_merge_gibbs_prob = 0).
   * Faster but simpler proposal mechanism than smart_merge_move().
   */
  void dumb_merge_move();

  /**
   * @brief Compute acceptance ratio for LSS merge move
   *
   * @param likelihood_old_ci Log-likelihood of first original cluster ci before
   * merge
   * @param likelihood_old_cj Log-likelihood of second original cluster cj
   * before merge
   * @return Log acceptance ratio for the merge proposal
   *
   * @details Computes log acceptance ratio as:
   * log(α) = log(prior_ratio) + log(likelihood_ratio) + log(proposal_ratio)
   * where:
   * - prior_ratio accounts for cluster size changes
   * - likelihood_ratio = L(ci_merged) - L(ci_old) - L(cj_old)
   * - proposal_ratio is log_merge_gibbs_prob (0 for dumb merge)
   */
  double compute_acceptance_ratio_merge(double likelihood_old_ci,
                                        double likelihood_old_cj);

  // ========== Shuffle Move Implementation ==========

  /**
   * @brief Execute a shuffle move using LSS
   *
   * @details Redistributes observations between two existing clusters ci and cj
   * using sequential allocation while maintaining the two-cluster structure.
   * Computes both forward (log_merge_gibbs_prob) and reverse
   * (log_split_gibbs_prob) proposal probabilities to ensure detailed balance.
   * This move helps improve mixing between existing clusters.
   */
  void shuffle();

  /**
   * @brief Compute acceptance ratio for LSS shuffle move
   *
   * @param likelihood_old_ci Log-likelihood of first cluster ci before shuffle
   * @param likelihood_old_cj Log-likelihood of second cluster cj before shuffle
   * @param old_ci_size Size of cluster ci before shuffle
   * @param old_cj_size Size of cluster cj before shuffle
   * @return Log acceptance ratio for the shuffle proposal
   *
   * @details Computes log acceptance ratio as:
   * log(α) = log(prior_ratio) + log(likelihood_ratio) + log(proposal_ratio)
   * where:
   * - prior_ratio accounts for cluster size changes in shuffle
   * - likelihood_ratio = L(ci_new) + L(cj_new) - L(ci_old) - L(cj_old)
   * - proposal_ratio = log_merge_gibbs_prob - log_split_gibbs_prob
   */
  double compute_acceptance_ratio_shuffle(double likelihood_old_ci,
                                          double likelihood_old_cj,
                                          int old_ci_size, int old_cj_size);

public:
  // ========== Constructor ==========

  /**
   * @brief Constructor for LSS-SDDS (Locality Sensitive Sampling with
   * Smart-split, Dumb-merge, Dumb-split, Smart-merge) Split-Merge sampler
   *
   * @param d Reference to Data object containing observations
   * @param p Reference to Params object with hyperparameters (including
   * distance matrix D)
   * @param l Reference to Likelihood object for probability computations
   * @param pr Reference to Process object defining the prior
   * @param shuffle Flag to enable shuffle moves in addition to split-merge
   *
   * @details Initializes the LSS-SDDS Split-Merge sampler, which uses:
   * - Locality sensitive sampling for anchor point selection
   * - SDDS strategy (Smart-split, Dumb-merge, Dumb-split, Smart-merge) that
   *   adaptively pairs smart and dumb moves based on point similarity
   * - Sequential allocation for "smart" proposal generation
   * - Simple random allocation for "dumb" proposals
   * This sampler provides computational advantages for large datasets by
   * intelligently balancing computational cost with proposal quality.
   */
  SplitMerge_LSS_SDDS(Data &d, Params &p, Likelihood &l, Process &pr,
                      bool shuffle)
      : Sampler(d, p, l, pr), shuffle_bool(shuffle),
        original_allocations(pr.old_allocations_view()), gen(rd()),
        dis_real(0.0, 1.0), dis_int(0, 1) {};

  // ========== MCMC Interface ==========

  /**
   * @brief Perform one iteration of the LSS-SDDS Split-Merge algorithm
   *
   * @details Executes one step of the LSS-SDDS Split-Merge sampler implementing
   * the SDDS strategy (Smart-split, Dumb-merge, Dumb-split, Smart-merge):
   *
   * 1. Randomly choose between dissimilarity mode or similarity mode (50/50)
   * 2. Select two anchor observations using locality sensitive sampling:
   *    - Dissimilarity mode (move_type=0): select dissimilar points (weighted
   * by 1/distance)
   *    - Similarity mode (move_type=1): select similar points (weighted by
   * distance)
   * 3. Determine move type based on current cluster assignments:
   *    - Same cluster (ci == cj): propose split move
   *    - Different clusters (ci != cj): propose merge move
   * 4. Execute appropriate move using SDDS pairing:
   *    - Dissimilarity mode: **Smart split** if ci==cj, **Dumb merge** if
   * ci!=cj
   *    - Similarity mode: **Dumb split** if ci==cj, **Smart merge** if ci!=cj
   * 5. Optionally perform shuffle move to improve mixing
   * 6. Compute acceptance ratio and accept/reject the proposal
   *
   * The SDDS strategy balances computational efficiency with proposal quality
   * by using smart (sequential allocation) moves where they matter most and
   * dumb (simple random) moves elsewhere, while maintaining detailed balance.
   */
  void step() override final;

  // ========== Accessor Methods ==========
  /**
   * @brief Get number of accepted split moves for diagnostics
   * @return Ratio of accepted split moves
   */
  double get_accepted_split() const {
    return static_cast<double>(accepted_split) / split_moves;
  };

  /**
   * @brief Get number of accepted merge moves for diagnostics
   * @return Ratio of accepted merge moves
   */
  double get_accepted_merge() const {
    return static_cast<double>(accepted_merge) / merge_moves;
  };

  /**
   * @brief Get number of accepted shuffle moves for diagnostics
   * @return Ratio of accepted shuffle moves
   */
  double get_accepted_shuffle() const {
    return static_cast<double>(accepted_shuffle) / shuffle_moves;
  };
};
