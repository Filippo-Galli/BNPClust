/**
 * @file MALA.cpp
 * @brief Implementation of the Metropolis-Adjusted Langevin Algorithm (MALA)
 * for sampling the latent variable U.
 */

#include "MALA.hpp"

void MALA::update_U() {
  // Increment iteration counter for adaptation
  total_iterations++;

  // Choose sampling scale: V = log(U) for numerical stability, or U directly
  if (use_V)
    sampling_V();
  else
    sampling_U();

  // Adaptively tune step size epsilon using Robbins-Monro algorithm
  if (tuning_enabled)
    Robbins_Monro_tuning();
}

void MALA::sampling_V() {
  // Transform current U to log scale for numerical stability
  double V_current = std::log(U);

  // Compute gradient of log density with respect to V
  double grad_V_current = grad_log_conditional_density_V(V_current);

  // Generate MALA proposal: V' = V + (ε²/2)∇log p(V) + ε*Z
  std::normal_distribution<double> noise(0.0, 1.0);
  double drift_current = (epsilon * epsilon / 2.0) * grad_V_current;
  double V_proposed = V_current + drift_current + epsilon * noise(gen);

  // Transform back to U scale for density evaluation
  double U_current = U;
  double U_proposed = std::exp(V_proposed);

  // Compute log densities on V scale (log p(V) = log p(U) + log|dU/dV| = log
  // p(U) + V)
  double log_density_V_current =
      log_conditional_density_U(U_current) + V_current;
  double log_density_V_proposed =
      log_conditional_density_U(U_proposed) + V_proposed;

  // Compute gradient at proposed point for reverse proposal density
  double grad_V_proposed = grad_log_conditional_density_V(V_proposed);
  double drift_proposed = (epsilon * epsilon / 2.0) * grad_V_proposed;

  // Helper lambda for log of Gaussian density
  auto log_norm_pdf = [](double x, double mu, double sigma) {
    return -0.5 * std::log(2.0 * M_PI * sigma * sigma) -
           0.5 * (x - mu) * (x - mu) / (sigma * sigma);
  };

  // Compute forward and backward proposal log densities
  // q(V'|V) = N(V' | V + drift_current, ε²)
  // q(V|V') = N(V | V' + drift_proposed, ε²)
  double log_q_forward =
      log_norm_pdf(V_proposed, V_current + drift_current, epsilon);
  double log_q_backward =
      log_norm_pdf(V_current, V_proposed + drift_proposed, epsilon);

  // Compute log acceptance ratio: log(α) = log p(V') - log p(V) + log q(V|V') -
  // log q(V'|V)
  double log_acceptance_ratio =
      (log_density_V_proposed - log_density_V_current) +
      (log_q_backward - log_q_forward);

  // Metropolis-Hastings accept/reject step
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  accept = std::log(unif(gen)) < log_acceptance_ratio;
  if (accept) {
    U = U_proposed;
    accepted_U++;
  }
}

double MALA::grad_log_conditional_density_V(double v) const {
  // Transform V back to U scale
  double u = std::exp(v);

  // Get gradient on U scale
  double grad_log_p_U = grad_log_conditional_density_U(u);

  // Apply chain rule: ∇_V log p(V) = ∇_U log p(U) * dU/dV + 1 = ∇_U log p(U) *
  // U + 1
  return grad_log_p_U * u + 1.0;
}

void MALA::sampling_U() {
  // Store current U value
  double U_current = U;

  // Compute gradient of log density at current point
  double grad_U_current = grad_log_conditional_density_U(U_current);

  // Generate MALA proposal: U' = U + (ε²/2)∇log p(U) + ε*Z, where Z ~ N(0,1)
  std::normal_distribution<double> noise(0.0, 1.0);
  double drift_current = (epsilon * epsilon / 2.0) * grad_U_current;
  double U_proposed = U_current + drift_current + epsilon * noise(gen);

  // Evaluate log densities at current and proposed points
  double log_density_U_current = log_conditional_density_U(U_current);
  double log_density_U_proposed = log_conditional_density_U(U_proposed);

  // Compute gradient at proposed point for backward proposal density
  double grad_U_proposed = grad_log_conditional_density_U(U_proposed);
  double drift_proposed = (epsilon * epsilon / 2.0) * grad_U_proposed;

  // Compute differences for proposal densities (numerator of Gaussian exponent)
  double forward_diff = U_proposed - U_current - drift_current;
  double backward_diff = U_current - U_proposed - drift_proposed;

  // Compute log proposal transition probabilities
  // q(U'|U) ∝ exp(-0.5 * (U' - U - drift_current)² / ε²)
  // q(U|U') ∝ exp(-0.5 * (U - U' - drift_proposed)² / ε²)
  double log_q_forward =
      -0.5 * (forward_diff * forward_diff) / (epsilon * epsilon);
  double log_q_backward =
      -0.5 * (backward_diff * backward_diff) / (epsilon * epsilon);

  // Compute log acceptance ratio: log(α) = log p(U') - log p(U) + log q(U|U') -
  // log q(U'|U)
  double log_acceptance_ratio = log_density_U_proposed - log_density_U_current +
                                log_q_backward - log_q_forward;

  // Metropolis-Hastings accept/reject step
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  accept = std::log(unif(gen)) < log_acceptance_ratio;
  if (accept) {
    U = U_proposed;
    accepted_U++;
  }
}

double MALA::grad_log_conditional_density_U(double u) const {
  // First term: derivative of (n-1)log(U) gives (n-1)/U
  const double term1 = (n - 1) / u;

  // Second term: derivative of -(n - σK)log(U+τ) gives -(n - σK)/(U+τ)
  const double term2 = -(n - params.sigma * data.get_K()) / (u + params.tau);

  // Third term: derivative of -a(U+τ)^σ gives -a(U+τ)^(σ-1)
  const double term3 = -params.a * std::pow(u + params.tau, params.sigma - 1);

  return term1 + term2 + term3;
}

void MALA::Robbins_Monro_tuning() {
  // Wait for initial burn-in period before starting adaptation
  if (total_iterations < BI_adapt)
    return;

  // Optimal acceptance rate for 1D MALA (Roberts & Rosenthal, 1998)
  const double target_acc = 0.574;

  // Step size constant for Robbins-Monro algorithm
  const double c = 1.0 / (target_acc * (1 - target_acc));

  // Compute deviation from target acceptance rate
  const double delta = (accept ? 1.0 : 0.0) - target_acc;

  // Update epsilon on log scale: log(ε) ← log(ε) + c*δ/t
  epsilon = std::exp(std::log(epsilon) + c * delta / total_iterations);

  // Enforce lower bound on epsilon to prevent numerical issues
  epsilon = epsilon < 1e-1 ? 1e-1 : epsilon;

  // Restart mechanism: reset iteration counter if epsilon changes too
  // dramatically This prevents getting stuck with very small step sizes
  double factor = epsilon / old_epsilon;
  if (factor > 3 || factor < 1.0 / 3.0) {
    total_iterations = BI_adapt;
    old_epsilon = epsilon;
  }

  // Periodic logging of epsilon value for monitoring
  // if (total_iterations % 1000 == 0)
  //   Rcpp::Rcout << "[DEBUG] - Updated epsilon: " << epsilon << std::endl;
}
