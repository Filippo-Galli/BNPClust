/**
 * @file RWMH.hpp
 * @brief Random Walk Metropolis-Hastings sampler for the latent variable U.
 */

#pragma once

#include "U_sampler.hpp"

/**
 * @class RWMH
 * @brief Random Walk Metropolis-Hastings sampler for updating the latent
 * variable U.
 *
 * This class implements a Random Walk Metropolis-Hastings (RWMH) algorithm for
 * sampling the latent variable U in NGGP mixture models. The sampler uses a
 * Gaussian proposal distribution centered at the current value.
 *
 * The sampler can operate on either:
 * - U scale directly (using a truncated normal proposal)
 * - V = log(U) scale (using an unrestricted normal proposal)
 *
 * The class also supports automatic tuning of the proposal standard deviation
 * using the Robbins-Monro stochastic approximation algorithm to achieve an
 * optimal acceptance rate of approximately 0.44 for one-dimensional proposals.
 *
 * @see U_sampler
 */
class RWMH : public U_sampler {

private:
  /**
   * @brief Standard deviation for the Gaussian proposal distribution.
   *
   * This parameter controls the step size of the random walk proposal.
   * It can be automatically tuned using Robbins-Monro adaptation if enabled.
   */
  double proposal_sd = 1.0;

  /** @brief Flag indicating whether the last proposal was accepted. */
  bool accept = false;

  /** @brief Flag to use V = log(U) scale for sampling instead of U directly. */
  bool use_V = false;

  /** @brief Flag to enable/disable Robbins-Monro automatic tuning. */
  bool tuning_enabled = false;

  /**
   * @name Sampling Methods
   * @brief Private methods implementing RWMH updates on different scales.
   * @{
   */

  /**
   * @brief Performs one RWMH update step for U on the original scale.
   *
   * Proposes a new value from a truncated normal distribution:
   * \f[
   * U' \sim \mathcal{N}(U_{\text{current}}, \sigma_{\text{prop}}^2) \quad
   * \text{truncated to } (0, \infty) \f]
   *
   * The acceptance probability accounts for the asymmetric proposal due to
   * truncation: \f[ \alpha = \min\left(1, \frac{f(U')}{f(U)} \cdot
   * \frac{q(U|U')}{q(U'|U)}\right) \f]
   *
   * where \f$q(U'|U)\f$ is the truncated normal proposal density.
   *
   * @see sampling_V()
   */
  void sampling_U();

  /**
   * @brief Performs one RWMH update step for V = log(U).
   *
   * Proposes a new value from an unrestricted normal distribution:
   * \f[
   * V' \sim \mathcal{N}(V_{\text{current}}, \sigma_{\text{prop}}^2)
   * \f]
   *
   * The acceptance probability simplifies to:
   * \f[
   * \alpha = \min\left(1, \frac{f(V')}{f(V)}\right)
   * \f]
   *
   * @note since the proposal is symmetric. Sampling on the log scale can
   * improve efficiency and avoids truncation issues.
   *
   * @see sampling_U()
   */
  void sampling_V();

  /** @} */

  /**
   * @brief Performs Robbins-Monro adaptive tuning of the proposal standard
   * deviation.
   *
   * Adjusts the proposal standard deviation to target an optimal acceptance
   * rate of 0.44 for one-dimensional random walk proposals. The tuning follows
   * the Robbins-Monro stochastic approximation scheme:
   *
   * \f[
   * \log(\sigma_{\text{prop}}^{(t+1)}) = \log(\sigma_{\text{prop}}^{(t)})
   * + \frac{c}{t} \cdot \delta_t
   * \f]
   *
   * where:
   * - \f$\delta_t = \mathbb{1}_{\{\text{accept}\}} - \alpha_{\text{target}}\f$
   * - \f$c = 1 / (\alpha_{\text{target}} \cdot (1 - \alpha_{\text{target}}))\f$
   * - \f$\alpha_{\text{target}} = 0.44\f$ (optimal for 1D random walk)
   * - \f$t\f$ is the current iteration number
   *
   * @note Reference: "Adaptive optimal scaling of Metropolis-Hastings
   * algorithms using the Robbins-Monro process" by Garthwaite et al. (2016).
   *
   * @see update_U()
   */
  void Robbins_Monro_tuning();

public:
  /**
   * @brief Constructor for the RWMH sampler.
   *
   * Initializes the Random Walk Metropolis-Hastings sampler with specified
   * configuration parameters.
   *
   * @param p Reference to the parameters object containing NGGP parameters.
   * @param d Reference to the data object containing observations and cluster
   * assignments.
   * @param use_V If true, perform sampling on V = log(U) scale; if false,
   * sample U directly.
   * @param prop_sd Initial standard deviation for the Gaussian proposal
   * distribution (default: 1.0).
   * @param tuning If true, enable Robbins-Monro adaptive tuning of proposal_sd
   * (default: false).
   *
   * @note When use_V = false, proposals are truncated to ensure U > 0.
   * @note When tuning = true, proposal_sd will be automatically adjusted during
   * sampling.
   */
  RWMH(NGGP_params &p, Data &d, bool use_V = false, double prop_sd = 1,
       bool tuning = false)
      : U_sampler(p, d), use_V(use_V), proposal_sd(prop_sd),
        tuning_enabled(tuning) {};

  /**
   * @brief Updates the latent variable U using Random Walk Metropolis-Hastings.
   *
   * This method performs one iteration of the RWMH algorithm. It:
   * 1. Increments the iteration counter
   * 2. Performs either sampling_U() or sampling_V() depending on use_V flag
   * 3. Optionally tunes the proposal standard deviation if tuning is enabled
   *
   * The choice between U and V scale affects the proposal distribution and
   * acceptance probability calculation.
   *
   * @see sampling_U(), sampling_V(), Robbins_Monro_tuning()
   */
  void update_U() override;
};
