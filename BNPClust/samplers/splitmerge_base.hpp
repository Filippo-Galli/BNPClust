/**
 * @file splitmerge_base.hpp
 * @brief Base class for Split-Merge MCMC samplers in Bayesian nonparametrics
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"
#include <cmath>
#include <random>
#include <type_traits>
#include <vector>

/**
 * @class SplitMergeBase
 * @brief Abstract base class providing shared data structures and algorithms
 * for Split-Merge MCMC samplers.
 *
 * Encapsulates anchor observation tracking, launch state extraction,
 * sequential allocation proposal computations, acceptance ratio evaluation,
 * diagnostic counters, and template method MCMC step coordination.
 */
class SplitMergeBase : public Sampler {
public:
  // ========== Mathematical Utilities ==========

  /**
   * @struct BinaryChoiceProbs
   * @brief Holds log-sum-exp normalizing constant and choice probability for binary decisions.
   */
  struct BinaryChoiceProbs {
    double log_sum; ///< log(exp(log_probs_0) + exp(log_probs_1))
    double prob_cj; ///< P(choose cj) = exp(log_probs_1) / (exp(log_probs_0) + exp(log_probs_1))

    [[nodiscard]] constexpr double prob_ci() const noexcept {
      return 1.0 - prob_cj;
    }
  };

  /**
   * @brief Numerically stable log-sum-exp and choice probability for binary decisions.
   *
   * Avoids dynamic allocations and Eigen vector overhead by evaluating the logistic function
   * with exactly one std::exp() and std::log1p() call.
   *
   * @param log_probs_0 Unnormalized log probability of cluster ci
   * @param log_probs_1 Unnormalized log probability of cluster cj
   * @return BinaryChoiceProbs containing log_sum and prob_cj
   */
  [[nodiscard]] static inline BinaryChoiceProbs
  compute_binary_choice_probs(double log_probs_0, double log_probs_1) noexcept {
    const double delta = log_probs_1 - log_probs_0;
    if (delta >= 0.0) {
      const double exp_neg = std::exp(-delta);
      return {log_probs_1 + std::log1p(exp_neg), 1.0 / (1.0 + exp_neg)};
    } else {
      const double exp_pos = std::exp(delta);
      return {log_probs_0 + std::log1p(exp_pos), exp_pos / (1.0 + exp_pos)};
    }
  }
protected:
  // ========== Random Number Generation ==========
  mutable std::mt19937 gen;
  mutable std::uniform_real_distribution<double> dis_real{0.0, 1.0};
  mutable std::uniform_int_distribution<int> dis_n;
  mutable std::uniform_int_distribution<size_t> dis_size;

  // ========== Move Selection Variables ==========
  int idx_i = 0; ///< First anchor observation
  int idx_j = 0; ///< Second anchor observation
  int ci = 0;    ///< Cluster of first anchor
  int cj = 0;    ///< Cluster of second anchor

  // ========== Algorithm Configuration ==========
  bool shuffle_bool = false; ///< Whether shuffle moves are enabled

  // ========== State Management ==========
  std::vector<int> launch_state; ///< Clusters of points in S
  std::vector<int> S; ///< Indices of non-anchor points in clusters ci and cj
  std::vector<int> original_allocations; ///< Copy of allocations before move

  // ========== Proposal Probabilities ==========
  double log_split_gibbs_prob = 0.0;
  double log_merge_gibbs_prob = 0.0;

  // ========== Diagnostic Counters ==========
  int accepted_split = 0;
  int accepted_merge = 0;
  int accepted_shuffle = 0;

  // ========== Internal Setup & Operations ==========

  /**
   * @brief Populate S and launch_state from clusters ci and cj in O(|Ci| +
   * |Cj|)
   * @param cluster_i Index of first cluster
   * @param cluster_j Index of second cluster
   */
  void prepare_launch_state(int cluster_i, int cluster_j);

  /**
   * @brief Fisher-Yates synchronized shuffle of launch_state and S in unison
   */
  void shuffle_launch_state();

  /**
   * @brief Uniform anchor pair selection
   * @param shuffle_s If true, shuffles launch_state and S after collection
   */
  void choose_indeces_uniform(bool shuffle_s = false);

  /**
   * @brief Sequential allocation proposal generator
   * @param iterations Number of sweeps
   * @param only_probabilities If true, only accumulate proposal probabilities
   * without updating allocations
   * @param sequential If true, unallocate all points before sequential
   * assignment
   */
  void sequential_allocation(int iterations, bool only_probabilities = false,
                             bool sequential = true);

  /**
   * @brief Select two distinct clusters uniformly for shuffle move
   */
  void choose_clusters_shuffle();

  /**
   * @brief Acceptance ratio computation for split move
   */
  double compute_acceptance_ratio_split(double likelihood_old_cluster);

  /**
   * @brief Acceptance ratio computation for merge move
   */
  double compute_acceptance_ratio_merge(double likelihood_old_clusters);

  /**
   * @brief Acceptance ratio computation for shuffle move
   */
  double compute_acceptance_ratio_shuffle(double likelihood_old_clusters,
                                          int old_ci_size, int old_cj_size);

  /**
   * @brief Rollback state upon rejected proposal
   */
  virtual void rollback() { data.set_allocations(original_allocations); }

  // ========== Template Method Hooks ==========
  virtual void choose_indices_for_step() = 0;
  virtual void split_move() = 0;
  virtual void merge_move() = 0;

public:
  SplitMergeBase(Data &d, Likelihood &l, Process &pr, bool shuffle)
      : Sampler(d, l, pr), shuffle_bool(shuffle), gen(rd()),
        dis_n(0, d.get_n() - 1) {}

  virtual ~SplitMergeBase() = default;

  /**
   * @brief Executes one split-merge MCMC step
   */
  void step() override;

  /**
   * @brief Default shuffle move using sequential allocation
   */
  virtual void shuffle();

  // ========== Accessors ==========
  [[nodiscard]] int get_accepted_split() const noexcept {
    return accepted_split;
  }
  [[nodiscard]] int get_accepted_merge() const noexcept {
    return accepted_merge;
  }
  [[nodiscard]] int get_accepted_shuffle() const noexcept {
    return accepted_shuffle;
  }
};
