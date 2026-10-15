/**
 * @file RWMH.cpp
 * @brief Implementation of Random Walk Metropolis-Hastings sampler methods.
 */

#include "RWMH.hpp"

// Main update method: increment counter, sample, and optionally tune
void RWMH::update_U() {

  total_iterations++;

  // Choose sampling scale
  if (use_V)
    sampling_V(); // Sample on log scale
  else
    sampling_U(); // Sample on original scale

  // Adaptively tune proposal sd if enabled
  if (tuning_enabled)
    Robbins_Monro_tuning();
}

// RWMH update on U scale with truncated normal proposal
void RWMH::sampling_U() {

  double U_current = U;

  // Propose from N(U_current, proposal_sd^2) truncated to (0, ∞)
  std::normal_distribution<double> proposal(U_current, proposal_sd);
  double U_proposed;
  do {
    U_proposed = proposal(gen);
  } while (U_proposed <= 0); // Reject non-positive proposals

  // Compute asymmetric proposal ratio due to truncation
  // q(U_c|U_p) / q(U_p|U_c) where q is truncated normal
  double log_q_current_given_proposed =
      std::log(1.0 - std::erf(-U_proposed / (std::sqrt(2.0) * proposal_sd)));
  double log_q_proposed_given_current =
      std::log(1.0 - std::erf(-U_current / (std::sqrt(2.0) * proposal_sd)));

  // Evaluate target densities
  const double log_density_current = log_conditional_density_U(U_current);
  const double log_density_proposed = log_conditional_density_U(U_proposed);

  // Metropolis-Hastings acceptance ratio
  double log_acceptance_ratio =
      (log_density_proposed - log_density_current) +
      (log_q_current_given_proposed - log_q_proposed_given_current);

  // Accept or reject
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  accept = std::log(unif(gen)) < log_acceptance_ratio;
  if (accept){
    U = U_proposed;
    accepted_U++;
  }
}

// RWMH update on V = log(U) scale with symmetric normal proposal
void RWMH::sampling_V() {

  double V_current = log(U);

  // Propose from N(V_current, proposal_sd^2) - no truncation needed
  std::normal_distribution<double> proposal(V_current, proposal_sd);
  double V_proposed = proposal(gen);

  // Evaluate target densities on V scale
  const double log_density_current = log_conditional_density_V(V_current);
  const double log_density_proposed = log_conditional_density_V(V_proposed);

  // Acceptance ratio (proposal cancels due to symmetry)
  double log_acceptance_ratio = log_density_proposed - log_density_current;

  // Accept or reject
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  accept = std::log(unif(gen)) < log_acceptance_ratio;
  if (accept){
    U = exp(V_proposed);
    accepted_U++;
  }
}

// Robbins-Monro adaptive tuning to achieve optimal acceptance rate
void RWMH::Robbins_Monro_tuning() {

  const double target_acc = 0.44;                         // Optimal for 1D RWMH
  const double c = 1.0 / (target_acc * (1 - target_acc)); // Step size constant

  // Compute deviation from target acceptance rate
  const double delta = (accept ? 1.0 : 0.0) - target_acc;

  // Update log(proposal_sd) using Robbins-Monro scheme
  proposal_sd = std::exp(std::log(proposal_sd) + c * delta / total_iterations);

  // TODO: missing restarting logic

  // Optional: print diagnostics every 1000 iterations
  // if(total_iterations % 1000 == 0)
  //     Rcpp::Rcout << "[DEBUG] - Updated proposal sd: " << proposal_sd <<
  //     std::endl;
}