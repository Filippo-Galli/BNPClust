/**
 * @file neal_ZDNAM.hpp
 * @brief Neal's Algorithm 3 with ZDNAM (Zero-self Downward Nested Antithetic
 * Modification)
 *
 * This file contains the implementation of Neal's Algorithm 3 enhanced with
 * ZDNAM, which eliminates self-transitions to improve MCMC mixing. The ZDNAM
 * modification is based on Neal's Algorithm 5 from "Modifying Gibbs Sampling to
 * Avoid Self Transitions".
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../utils/Sampler.hpp"

/**
 * @class Neal3ZDNAM
 * @brief Neal's Algorithm 3 with ZDNAM for improved collapsed Gibbs sampling
 *
 * @details This class implements Neal's Algorithm 3 (collapsed Gibbs sampler)
 * enhanced with the Zero-self Downward Nested Antithetic Modification (ZDNAM)
 * from Neal's Algorithm 5. The ZDNAM modification eliminates self-transitions
 * in the Gibbs sampler, which can significantly improve mixing and reduce
 * autocorrelation in the MCMC chain.
 *
 * **Key Features:**
 * - **Collapsed sampling**: Cluster parameters are integrated out analytically
 * - **Zero self-transitions**: ZDNAM ensures observations always move (when
 * possible)
 * - **Improved mixing**: Better exploration of the state space compared to
 * standard Gibbs
 * - **Automatic cluster management**: Creates and destroys clusters as needed
 *
 * **ZDNAM Algorithm:**
 * The ZDNAM modification constructs transition probabilities P(new_state |
 * current_state) such that P(current_state | current_state) = 0 when possible,
 * while maintaining the correct stationary distribution. This is achieved
 * through a careful construction involving nested antithetic coupling.
 *
 * **When to use:**
 * - High autocorrelation in standard Gibbs sampling
 * - When mixing is slow due to frequent self-transitions
 * - Models where observations tend to stay in the same cluster
 *
 * @note The algorithm falls back to standard Gibbs sampling in edge cases where
 * zero self-transition is not achievable (e.g., when one state has probability
 * ≥ 0.5).
 *
 * @note
 * reference: Neal, R. M. (2000). "Markov Chain Sampling Methods for Dirichlet
 * Process Mixture Models" \n
 * reference: Neal, R. M. (2024). "Modifying Gibbs Sampling to Avoid Self
 * Transitions"
 *
 * @see Sampler, Process, Likelihood, Neal3
 */
class Neal3ZDNAM : public Sampler {
private:
  // ========== Random Number Generation ==========

  /** @brief Mersenne Twister random number generator for sampling operations */
  mutable std::mt19937 gen;

  // ========== Core Algorithm Methods ==========

  /**
   * @brief Update cluster assignment for a single observation using ZDNAM
   *
   * @param index Index of the observation to update (0 to n-1)
   *
   * @details This method implements Neal's Algorithm 3 with ZDNAM enhancement:
   * 1. Record the current cluster assignment
   * 2. Temporarily remove the observation from its cluster
   * 3. Compute log probabilities for all existing clusters plus a new cluster
   * 4. Use ZDNAM to sample new assignment (avoiding self-transition when
   * possible)
   * 5. Assign observation to the sampled cluster
   *
   * The log probabilities combine:
   * - Prior probabilities from the Process (DP/NGGP weights)
   * - Conditional likelihood (with cluster parameters integrated out)
   *
   * The ZDNAM sampler uses the current cluster index to construct transition
   * probabilities that avoid returning to the same cluster.
   *
   * @see sample_from_log_probs_zdnam(), compute_zdnam_probabilities()
   */
  void step_1_observation(int index);

  /**
   * @brief Sample from log probabilities using standard Gibbs sampling
   *
   * @param log_probs Vector of unnormalized log probabilities
   * @return Sampled index from 0 to log_probs.size()-1
   *
   * @details Standard categorical sampling:
   * 1. Normalize log probabilities using log-sum-exp trick
   * 2. Convert to probabilities
   * 3. Sample using cumulative sum (roulette wheel)
   *
   * @note This is the fallback when ZDNAM is not applicable
   */
  int sample_from_log_probs(const std::vector<double> &log_probs);

  /**
   * @brief Compute ZDNAM transition probabilities from current state
   *
   * @param log_probs Vector of unnormalized log probabilities (target
   * distribution π)
   * @param current_k Current state/cluster index (-1 if no current state)
   * @return Vector of transition probabilities P(·|current_k) with
   * P(current_k|current_k) = 0
   *
   * @details Implements Algorithm 5 from Neal (2019). The algorithm constructs
   * transition probabilities that:
   * - Have the correct stationary distribution π
   * - Satisfy P(current_k | current_k) = 0 when possible
   * - Use nested antithetic coupling for improved mixing
   *
   * **Algorithm Steps:**
   * 1. Normalize input to get target probabilities π
   * 2. Order states by decreasing probability (σ permutation)
   * 3. Process states sequentially, allocating transition probability
   * 4. Apply special ZDNAM construction when current state is reached
   * 5. Ensure zero self-transition probability
   *
   * **Special Cases:**
   * - If π[most_probable] ≥ 0.5: Use simpler construction (downward nesting)
   * - If current state has π = 0: Fall back to standard Gibbs
   * - If m ≤ 2: Simplified construction
   *
   * **Mathematical Guarantee:**
   * The returned transition matrix satisfies:
   * - Σⱼ P(j|current_k) = 1 (proper distribution)
   * - Σⱼ P(j|i)π(i) = π(j) (detailed balance with π)
   * - P(current_k|current_k) = 0 (zero self-transition)
   *
   * @note The construction involves careful numerical balancing to avoid
   * negative probabilities while maintaining detailed balance.
   *
   * @see sample_from_log_probs_zdnam()
   */
  std::vector<double>
  compute_zdnam_probabilities(const std::vector<double> &log_probs,
                              int current_k);

  /**
   * @brief Sample from log probabilities using ZDNAM to avoid self-transitions
   *
   * @param log_probs Vector of unnormalized log probabilities (target
   * distribution)
   * @param current_k Current state/cluster index (-1 for no current state)
   * @return Sampled index from 0 to log_probs.size()-1
   *
   * @details High-level ZDNAM sampling procedure:
   * 1. If current_k is valid: compute ZDNAM transition probabilities
   * 2. If current_k is invalid: use standard normalized probabilities
   * 3. Sample from the resulting distribution using roulette wheel
   *
   * **Behavior:**
   * - When current_k ≥ 0: Returns a state different from current_k (with high
   * probability)
   * - When current_k = -1: Samples from normalized π (standard Gibbs)
   *
   * **Use Case:**
   * Called by step_1_observation() to sample new cluster assignments while
   * avoiding the tendency to stay in the current cluster.
   *
   * @see compute_zdnam_probabilities()
   */
  int sample_from_log_probs_zdnam(const std::vector<double> &log_probs,
                                  int current_k);

public:
  // ========== Constructor ==========

  /**
   * @brief Constructor for Neal's Algorithm 3 with ZDNAM sampler
   *
   * @param d Reference to Data object containing observations and cluster
   * assignments
   * @param l Reference to Likelihood object for computing conditional
   * probabilities
   * @param pr Reference to Process object (DP, NGGP, etc.) defining the prior
   * distribution
   *
   * @details Initializes the ZDNAM-enhanced Gibbs sampler with all required
   * components. The random number generator is seeded from the inherited random
   * device to ensure reproducibility when the random device is seeded
   * externally.
   *
   * @note The ZDNAM modification is automatically applied during sampling; no
   * additional configuration is needed.
   */
  Neal3ZDNAM(Data &d, Likelihood &l, Process &pr)
      : Sampler(d, l, pr), gen(rd()) {};

  // ========== MCMC Interface ==========

  /**
   * @brief Perform one complete iteration of Neal's Algorithm 3 with ZDNAM
   *
   * @details Executes one full sweep of the ZDNAM-enhanced collapsed Gibbs
   * sampler:
   * 1. Create a vector of all observation indices [0, 1, ..., n-1]
   * 2. Randomly shuffle the order to avoid systematic bias
   * 3. For each observation (in shuffled order):
   *    - Update its cluster assignment using ZDNAM
   *    - Automatically create/destroy clusters as needed
   * 4. Return with updated cluster configuration
   *
   * **After each step:**
   * - All observations have been given a chance to move clusters
   * - Zero self-transitions have been enforced (where possible)
   * - The chain has advanced by one full iteration
   * - The stationary distribution remains correct
   *
   * **Expected Behavior:**
   * Compared to standard Neal's Algorithm 3, this version should show:
   * - Reduced autocorrelation in cluster assignments
   * - Faster mixing of the MCMC chain
   * - More frequent cluster reassignments
   *
   * @note The random shuffling ensures that the order of updates does not
   * introduce systematic bias into the sampler.
   *
   * @see step_1_observation()
   */
  void step() override;
};
