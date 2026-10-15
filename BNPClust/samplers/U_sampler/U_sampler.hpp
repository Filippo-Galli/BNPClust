/**
 * @file U_sampler.hpp
 * @brief Base class for sampling the latent variable U in NGGP mixture models.
 */

#pragma once

#include "../../params/NGGP-params.hpp"
#include "../../utils/Data.hpp"
#include <random>

/**
 * @class U_sampler
 * @brief Abstract base class for MCMC sampling of the latent variable U.
 *
 * This class provides the interface and common functionality for sampling
 * the latent variable U in Normalized Generalized Gamma Process (NGGP)
 * mixture models. The latent variable U appears in the conditional distribution
 * of the partition and must be updated via MCMC methods.
 *
 * Derived classes must implement the update_U() method using specific
 * sampling algorithms such as Random Walk Metropolis-Hastings (RWMH)
 * or Metropolis-Adjusted Langevin Algorithm (MALA).
 *
 * @note Reference: "Clustering blood donors via mixtures of product partition
 *       models with covariates" by Argiento et al., 2024
 */
class U_sampler {

protected:
  /** @brief Reference to the parameters object containing NGGP parameters. */
  const NGGP_params &params;

  /** @brief Reference to the data object containing observations and cluster
   * assignments. */
  Data &data;

  /** @brief Counter for total MCMC iterations performed. */
  int total_iterations = 0;

  /** @brief Counter of accepted U*/
  int accepted_U = 0;

  /** @brief Random device for seeding. */
  std::random_device rd;

  /** @brief Mersenne Twister random number generator. */
  mutable std::mt19937 gen;

  /**
   * @name Cached Constants
   * @brief Pre-computed values for computational efficiency.
   * @{
   */

  /** @brief Ratio a/sigma for efficient computation. */
  const double a_over_sigma = params.a / params.sigma;

  /** @brief Pre-computed tau^sigma for efficient computation. */
  const double tau_power_sigma = std::pow(params.tau, params.sigma);

  /** @brief Number of observations n. */
  const int n = data.get_n();

  /** @brief Current value of the latent variable U (initialized to 1.0). */
  double U = 1.0;

  /** @} */

  /**
   * @name Conditional Density Functions
   * @brief Methods for computing conditional densities of U and V = log(U).
   *
   * These functions evaluate the (unnormalized) conditional density of the
   * latent variable U given the current partition, which is required for
   * MCMC acceptance probabilities in Metropolis-Hastings algorithms.
   * @{
   */

  /**
   * @brief Computes the log conditional density of U given the partition.
   *
   * This method evaluates the unnormalized log conditional density of the
   * latent variable U given the current partition π. The conditional density
   * follows:
   *
   * \f[
   * f_{U|\pi}(u) \propto u^{n-1} \cdot (u + \tau)^{-(n-\sigma|\pi|)}
   * \cdot \exp\left(-\frac{a}{\sigma}\left((u+\tau)^\sigma -
   * \tau^\sigma\right)\right) \f]
   *
   * where:
   * - \f$n\f$ is the number of observations
   * - \f$|\pi|\f$ is the number of clusters (K)
   * - \f$\sigma, a, \tau\f$ are NGGP parameters
   *
   * The log density is computed as the sum of three terms:
   * 1. \f$(n-1) \log(u)\f$
   * 2. \f$-(n - \sigma K) \log(u + \tau)\f$
   * 3. \f$-\frac{a}{\sigma}\left((u+\tau)^\sigma - \tau^\sigma\right)\f$
   *
   * @param u The value of U at which to evaluate the density.
   * @return The log of the unnormalized conditional density.
   *
   * @note Reference: "Clustering blood donors via mixtures of product partition
   *       models with covariates" by Argiento et al., 2024
   *
   * @see log_conditional_density_V()
   */
  double log_conditional_density_U(double u) const;

  /**
   * @brief Computes the log conditional density of V = log(U) given the
   * partition.
   *
   * This method evaluates the unnormalized log conditional density of the
   * transformed latent variable V = log(U) given the current partition π.
   * The transformation uses the change-of-variables formula:
   *
   * \f[
   * f_{V|\pi}(v) = f_{U|\pi}(e^v) \cdot \left|\frac{dU}{dV}\right|
   * = f_{U|\pi}(u) \cdot u
   * \f]
   *
   * Taking the logarithm:
   * \f[
   * \log f_{V|\pi}(v) = \log f_{U|\pi}(e^v) + \log(e^v) = \log f_{U|\pi}(u) + v
   * \f]
   *
   * This transformed density is useful for MCMC sampling on the log scale,
   * which can improve numerical stability and sampling efficiency.
   *
   * @param v The value of V = log(U) at which to evaluate the density.
   * @return The log of the unnormalized conditional density of V.
   *
   * @see log_conditional_density_U()
   */
  double log_conditional_density_V(double v) const;

  /** @} */

public:
  /**
   * @brief Constructor for the U_sampler class.
   *
   * Initializes the sampler with references to the parameters and data objects.
   * The latent variable U is initialized to 1.0.
   *
   * @param p Reference to the parameters object containing NGGP parameters (a,
   * sigma, tau).
   * @param d Reference to the data object containing observations and cluster
   * assignments.
   */
  U_sampler(NGGP_params &p, Data &d) : params(p), data(d) {};

  /**
   * @brief Pure virtual method to update the latent variable U.
   *
   * This method must be implemented by derived classes to perform MCMC
   * updates of the latent variable U using specific sampling algorithms
   * (e.g., Random Walk Metropolis-Hastings, MALA).
   */
  virtual void update_U() = 0;

  /**
   * @brief Getter for the current value of U.
   *
   * @return The current value of the latent variable U.
   */
  double get_U() const { return U; }

  /**
   * @brief Getter for the acceptance rate of U updates.
   * @return The acceptance rate as a double.
   */
  double get_acceptance_rate() const {
    if (total_iterations == 0)
      return 0.0;
    return static_cast<double>(accepted_U) /
           static_cast<double>(total_iterations);
  }

  /**
  @brief Virtual destructor for the U_sampler class.
  */
  virtual ~U_sampler() = default;
};
