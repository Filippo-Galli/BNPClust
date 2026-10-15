/**
 * @file GaussianMixtureModel-params.hpp
 * @brief Parameter management for Gaussian Mixture Model likelihood
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include <Eigen/Dense>

/**
 * @brief Parameters for a univariate Gaussian Mixture Model with
 * Normal-Inverse-Gamma conjugate priors.
 */
struct GaussianMixtureModel_params {
  /** @brief Prior mean for the Gaussian component mean */
  double m0;

  /** @brief Prior strength for the Gaussian component mean */
  double kappa0;

  /** @brief Shape parameter of the Inverse-Gamma prior for variance */
  double alpha0;

  /** @brief Scale parameter of the Inverse-Gamma prior for variance */
  double beta0;

  /**
   * @brief Constructor for GaussianMixtureModel_params
   * @param m0 Prior mean
   * @param kappa0 Prior strength for mean
   * @param alpha0 Shape for variance prior
   * @param beta0 Scale for variance prior
   */
  GaussianMixtureModel_params(double m0, double kappa0, double alpha0,
                              double beta0)
      : m0(m0), kappa0(kappa0), alpha0(alpha0), beta0(beta0) {}
};
