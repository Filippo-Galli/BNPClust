/**
 * @file Natarajan-params.hpp
 * @brief Parameter management for the likelihood taken from Natarajan et al.
 * (2023) "Cohesion and Repulsion in Bayesian Distance Clustering"
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

struct Natarajan_params {

  /** @brief Shape parameter for the first gamma distribution prior */
  double delta1;

  /** @brief Shape parameter for the lambda_k gamma distribution */
  double alpha;

  /** @brief Rate parameter for the lambda_k gamma distribution */
  double beta;

  /** @brief Shape parameter for the second gamma distribution prior */
  double delta2;

  /** @brief Shape parameter for the theta_kt gamma distribution */
  double gamma;

  /** @brief Rate parameter for the theta_kt gamma distribution */
  double zeta;

  /**
   * @brief Constructor for Natarajan_params
   * @param delta1 Shape parameter for the first gamma distribution prior
   * @param alpha Shape parameter for the lambda_k gamma distribution
   * @param beta Rate parameter for the lambda_k gamma distribution
   * @param delta2 Shape parameter for the second gamma distribution prior
   * @param gamma Shape parameter for the theta_kt gamma distribution
   * @param zeta Rate parameter for the theta_kt gamma distribution
   */
  Natarajan_params(double delta1, double alpha, double beta, double delta2,
                   double gamma, double zeta)
      : delta1(delta1), alpha(alpha), beta(beta), delta2(delta2), gamma(gamma),
        zeta(zeta) {}
};
