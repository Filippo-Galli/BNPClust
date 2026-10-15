/**
 * @file U_sampler.cpp
 * @brief Implementation of U_sampler class methods.
 */

#include "U_sampler.hpp"

// Compute log f_{U|π}(u)
double U_sampler::log_conditional_density_U(double u) const {

  const int K = data.get_K();    // Current number of clusters
  const double tau = params.tau; // Allow for future updates of tau

  // Term 1: (n - 1) * log(u)
  const double term1 = (n - 1) * log(u);

  // Term 2: -(n - sigma*K) * log(u + tau)
  const double term2 = -(n - params.sigma * K) * std::log(u + tau);

  // Term 3: -(a/sigma) * ((u + tau)^sigma - tau^sigma)
  const double term3 = -a_over_sigma * (std::pow(u + tau, params.sigma) - tau_power_sigma);

  return term1 + term2 + term3;
}

// Compute log f_{V|π}(v) where V = log(U)
double U_sampler::log_conditional_density_V(double v) const {
  const double u = std::exp(v); // Transform back to U scale
  // Apply change of variables: log f_V(v) = log f_U(u) + log|dU/dV| = log(f_U(u)) + v
  return log_conditional_density_U(u) + v;
}