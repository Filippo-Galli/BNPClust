/**
 * @file neal_ZDNAM.cpp
 * @brief Implementation of ZDNAM (Zero-self Downward Nested Antithetic
 * Modification)
 *
 * This file implements Algorithm 5 from Neal's 2024 paper "Modifying Gibbs
 * Sampling to Avoid Self Transitions". The ZDNAM technique eliminates
 * self-transitions in Gibbs sampling, improving mixing and reducing
 * autocorrelation.
 */

#include "neal_ZDNAM.hpp"
#include <algorithm>
#include <numeric>

std::vector<double>
Neal3ZDNAM::compute_zdnam_probabilities(const std::vector<double> &log_probs,
                                        int current_k) {
  // Initialize result vector with number of possible states
  int m = log_probs.size();
  std::vector<double> p(m, 0.0);

  // Step 1: Convert log probabilities to normalized probabilities (π)
  // Use log-sum-exp trick for numerical stability
  double max_loglik = *std::max_element(log_probs.begin(), log_probs.end());
  std::vector<double> pi(m);
  double sum_pi = 0.0;
  for (int i = 0; i < m; ++i) {
    pi[i] = exp(log_probs[i] - max_loglik);
    sum_pi += pi[i];
  }
  // Normalize to get target distribution π
  for (int i = 0; i < m; ++i) {
    pi[i] /= sum_pi;
  }

  // Edge case: current state has zero probability
  // Cannot construct ZDNAM, fall back to standard Gibbs
  if (current_k >= 0 && pi[current_k] <= 0.0) {
    return pi;
  }

  // Step 2: Create ordering σ by non-increasing probability
  // This is the focal value ordering used in Algorithm 5
  std::vector<int> sigma(m);
  std::iota(sigma.begin(), sigma.end(), 0);
  std::sort(sigma.begin(), sigma.end(),
            [&pi](int a, int b) { return pi[a] > pi[b]; });

  // Step 3: Handle special case where most probable state has π ≥ 1/2
  // Use simpler downward nesting construction
  if (pi[sigma[0]] >= 0.5 || m <= 2) {
    if (current_k == sigma[0]) {
      // Current state is most probable: allow reduced self-transition
      p[current_k] = (2.0 * pi[current_k] - 1.0) / pi[current_k];
      // Allocate remaining probability to other states
      for (int i = 0; i < m; ++i) {
        if (i != current_k) {
          p[i] = std::min(1.0, pi[i] / pi[current_k]);
        }
      }
    } else {
      // Current state is not most probable: always move to most probable
      p[sigma[0]] = 1.0;
    }
    return p;
  }

  // Step 4: Main ZDNAM algorithm (Algorithm 5 from Neal 2019)
  // s = sum of π values for states not yet processed (remaining probability
  // mass) f = available transition probability from current state to allocate
  double s = 1.0;
  double f = 1.0;
  int i = 0;

  // Process states in order of decreasing probability (focal value
  // construction)
  while (f > 0 && i < m && sigma[i] != current_k) {
    // Process state sigma[i] as potential focal value
    if (i + 1 < m) {
      double q = pi[sigma[i]]; // Probability of current focal value
      s -= q;                  // Remove from remaining mass

      // Look ahead to check if special ZDNAM construction is needed
      if (i + 2 < m) {
        double q2 = pi[sigma[i + 1]];      // Next focal value probability
        double s2 = std::max(0.0, s - q2); // Mass after next focal value

        // Condition for special ZDNAM construction (Theorem 4 in paper)
        // This ensures we can achieve zero self-transition
        if (q2 >= s2 && s2 > 0.0) {
          // Special ZDNAM construction coefficients
          // A is not used in this implementation
          // B and C balance probabilities to maintain detailed balance
          double A = (q + q2 - s2) / 2.0;
          double B = (q - q2 + s2) / (2.0 * s2);
          double C = (s2 + q2 - q) / (2.0 * s2);

          // Allocate transition probabilities using special construction
          p[sigma[i]] = f * B;
          p[sigma[i + 1]] = f * C;

          // Check if current state is in remaining values
          // If so, distribute probability proportionally
          for (int j = i + 2; j < m; ++j) {
            if (sigma[j] == current_k) {
              for (int k = i + 2; k < m; ++k) {
                p[sigma[k]] = f * B * pi[sigma[k]] / q;
              }
              return p;
            }
          }

          return p;
        }
      }

      // Regular DNAM step (downward nesting)
      // Allocate as much probability as possible proportional to target
      if (s > 0.0) {
        p[sigma[i]] = std::min(f, (q / s) * f);
        f -= p[sigma[i]]; // Reduce available transition probability
      }
    }
    ++i;
  }

  // Step 5: Handle when current state is reached in the ordering
  // This is where we enforce zero self-transition
  if (i < m && sigma[i] == current_k) {
    double q = pi[current_k];
    s -= q;

    if (s > 0.0) {
      // Zero self-transition (key feature of ZDNAM)
      p[current_k] = 0.0;
      // Distribute all remaining transition probability to other states
      // proportional to their target probabilities
      for (int j = i + 1; j < m; ++j) {
        p[sigma[j]] = std::min(f, (pi[sigma[j]] / s) * f);
      }
    } else {
      // Only current state remains: must self-transition
      p[current_k] = f;
    }
  }

  return p;
}

int Neal3ZDNAM::sample_from_log_probs_zdnam(
    const std::vector<double> &log_probs, int current_k) {
  std::vector<double> probs;

  // Decide which probabilities to use based on current state
  if (current_k >= 0 && current_k < log_probs.size()) {
    // Valid current state: use ZDNAM transition probabilities
    // This will avoid self-transition when possible
    probs = compute_zdnam_probabilities(log_probs, current_k);
  } else {
    // No valid current state (e.g., first allocation)
    // Use standard normalized probabilities from target distribution
    double max_loglik = *std::max_element(log_probs.begin(), log_probs.end());
    probs.resize(log_probs.size());
    double sum_probs = 0.0;

    // Log-sum-exp normalization
    for (size_t i = 0; i < log_probs.size(); ++i) {
      probs[i] = exp(log_probs[i] - max_loglik);
      sum_probs += probs[i];
    }
    for (size_t i = 0; i < log_probs.size(); ++i) {
      probs[i] /= sum_probs;
    }
  }

  // Sample from discrete distribution using inverse CDF (roulette wheel)
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double u = unif(gen);
  double cumsum = 0.0;

  for (size_t i = 0; i < probs.size(); ++i) {
    cumsum += probs[i];
    if (u < cumsum) {
      return i;
    }
  }

  // Numerical safety: if we get here due to floating point errors,
  // return the last state
  return probs.size() - 1;
}

void Neal3ZDNAM::step_1_observation(int index) {
  // Step 1: Record current cluster before removing observation
  int current_cluster = data.get_cluster_assignment(index);

  // Step 2: Temporarily remove observation from its cluster
  // This is necessary to compute the conditional probabilities correctly
  data.set_allocation(index, -1);

  // Step 3: Compute log probabilities for all possible assignments
  // Vector size = K existing clusters + 1 potential new cluster
  std::vector<double> log_likelihoods(data.get_K() + 1, 0.0);

  // Compute for existing clusters
  for (int k = 0; k < data.get_K(); ++k) {
    // Likelihood: P(x_i | allocated to cluster k, other assignments)
    log_likelihoods[k] = likelihood.point_loglikelihood_cond(index, k);
    // Prior: P(c_i = k | other assignments) from DP/NGGP
    log_likelihoods[k] += process.gibbs_prior_existing_cluster(k, index);
  }

  // Compute for potential new cluster
  log_likelihoods[data.get_K()] =
      likelihood.point_loglikelihood_cond(index, data.get_K());
  log_likelihoods[data.get_K()] += process.gibbs_prior_new_cluster_obs(index);

  // Step 4: Sample new assignment using ZDNAM
  // Pass current_cluster to enable zero self-transition
  int sampled_cluster =
      sample_from_log_probs_zdnam(log_likelihoods, current_cluster);

  // Step 5: Assign observation to sampled cluster
  // This automatically handles cluster creation/deletion
  data.set_allocation(index, sampled_cluster);
}

void Neal3ZDNAM::step() {
  // Step 1: Create vector of all observation indices [0, 1, ..., n-1]
  std::vector<int> indices(data.get_n());
  std::iota(indices.begin(), indices.end(), 0);

  // Step 2: Randomly shuffle to avoid systematic bias in update order
  std::shuffle(indices.begin(), indices.end(), gen);

  // Step 3: Update each observation in random order using ZDNAM
  for (int j = 0; j < indices.size(); ++j) {
    step_1_observation(indices[j]);
  }

  // Optional debug output (commented out for performance)
  // Rcpp::Rcout << std::endl << "Element per cluster: " << std::endl;
  // for (int k = 0; k < data.get_K(); ++k) {
  //     Rcpp::Rcout << "\tCluster " << k << ": " << data.get_cluster_size(k) <<
  //     std::endl;
  // }
  // Rcpp::Rcout <<
  // "-----------------------------------------------------------------------------"
  // << std::endl;
}