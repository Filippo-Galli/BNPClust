/**
 * @file MALA.hpp
 * @brief Metropolis-Adjusted Langevin Algorithm (MALA) sampler for the latent
 * variable U.
 */

#pragma once

#include "U_sampler.hpp"

/**
 * @class MALA
 * @brief Metropolis-Adjusted Langevin Algorithm (MALA) sampler for the latent
 * variable U.
 *
 * @details This class implements the MALA algorithm, a gradient-based MCMC
 * sampler that uses first-order derivative information to propose more informed
 * moves in the parameter space. The algorithm can operate either on U directly
 * or on V = log(U) for improved numerical stability.
 *
 * MALA uses the following proposal mechanism:
 * - U' = U + (ε²/2)∇log p(U|·) + ε*Z, where Z ~ N(0,1)
 *
 * The class supports adaptive tuning of the step size ε using the Robbins-Monro
 * algorithm to achieve an optimal acceptance rate of approximately 0.574
 * (optimal for 1D problems).
 *
 * @see U_sampler
 */
class MALA : public U_sampler {

private:
  /** @brief Step size parameter for the Langevin dynamics proposal. */
  double epsilon = 1.0;

  /** @brief Previous value of epsilon, used for restart mechanism in adaptive
   * tuning. */
  double old_epsilon = 1.0;

  /**
   * @brief Computes the gradient of log conditional density with respect to U.
   *
   * @details The gradient is given by:
   * ∇_U log p(U|·) = (n-1)/U - (n - σK)/(U+τ) - a(U+τ)^(σ-1)
   *
   * @param U The current value of the latent variable U
   * @return The gradient value at U
   */
  double grad_log_conditional_density_U(double U) const;

  /**
   * @brief Computes the gradient of log conditional density with respect to V =
   * log(U).
   *
   * @details Uses the chain rule to transform the gradient from U to V scale:
   * ∇_V log p(V|·) = ∇_U log p(U|·) * U + 1
   *
   * @param v The current value of V = log(U)
   * @return The gradient value at V
   */
  double grad_log_conditional_density_V(double v) const;

  /** @brief Flag indicating whether the last proposal was accepted. */
  bool accept = false;

  /**
   * @brief Samples U using MALA on the original U scale.
   *
   * @details Implements the MALA proposal mechanism directly on U:
   * - Proposes U' = U + (ε²/2)∇log p(U) + ε*Z
   * - Computes acceptance ratio including asymmetric proposal densities
   * - Accepts/rejects according to Metropolis-Hastings criterion
   */
  void sampling_U();

  /**
   * @brief Samples U using MALA on the transformed V = log(U) scale.
   *
   * @details Implements MALA on the log-transformed space for better numerical
   * stability:
   * - Proposes V' = V + (ε²/2)∇log p(V) + ε*Z
   * - Transforms back to U scale: U' = exp(V')
   * - Includes Jacobian correction in acceptance ratio
   * - Useful when U is constrained to be positive and has large dynamic range
   */
  void sampling_V();

  /** @brief Flag to use V = log(U) transformation instead of direct U sampling.
   */
  bool use_V = false;

  /** @brief Flag to enable/disable Robbins-Monro adaptive tuning of step size.
   */
  bool tuning_enabled = false;

  /**
   * @brief Performs Robbins-Monro tuning of the proposal step size epsilon.
   *
   * @details Adjusts epsilon to target an acceptance rate of 0.574 (optimal for
   * 1D MALA), based on the Robbins-Monro stochastic approximation algorithm.
   * The update rule is:
   *
   * epsilon = exp(log(epsilon) + c * δ / t)
   *
   * where:
   * - δ = (accepted ? 1.0 : 0.0) - target_acc
   * - c = 1.0 / (target_acc * (1 - target_acc))
   * - t = total_iterations
   *
   * Includes a restart mechanism to prevent epsilon from changing too rapidly
   * (resets iteration counter if epsilon changes by more than factor of 3).
   *
   * @note Only active after BI_adapt initial iterations and when tuning_enabled
   * is true
   * @see Garthwaite et al. (2016), "Adaptive optimal scaling of
   * Metropolis-Hastings algorithms using the Robbins-Monro process"
   */
  void Robbins_Monro_tuning();

  /**
   * @brief Number of initial iterations before adaptive tuning begins.
   *
   * @details Set to 5/(target_acc * (1 - target_acc)) ≈ 20.4 iterations.
   * This allows some initial exploration before adaptation starts.
   */
  int BI_adapt = 5 / (0.574 * (1 - 0.574));

public:
  /**
   * @brief Constructor for MALA sampler.
   *
   * @param p Reference to model parameters
   * @param d Reference to observed data
   * @param use_V If true, sample on V = log(U) scale instead of U scale
   * (default: false)
   * @param eps Initial step size for Langevin dynamics (default: 1.0)
   * @param tuning If true, enable Robbins-Monro adaptive tuning (default:
   * false)
   */
  MALA(NGGP_params &p, Data &d, bool use_V = false, double eps = 1,
       bool tuning = false)
      : U_sampler(p, d), epsilon(eps), old_epsilon(epsilon),
        tuning_enabled(tuning), use_V(use_V) {};

  /**
   * @brief Performs one MALA update step for the latent variable U.
   *
   * @details Executes the following sequence:
   * 1. Increments iteration counter
   * 2. Performs MALA sampling (either on V or U scale based on use_V flag)
   * 3. Applies Robbins-Monro tuning if enabled
   *
   * @note Overrides the pure virtual function from U_sampler base class
   */
  void update_U() override;
};
