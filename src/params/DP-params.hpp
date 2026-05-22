/**
 * @file DP-params.hpp
 * @brief Parameter management for models using DP-process (Dirichlet Process)
 *
 * This file contains the parameter management for models using DP-process.
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include "Process-params.hpp"
#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

/**
 * @brief Struct for storing parameters of the DP (Dirichlet Process)
 *
 * This struct holds the parameters for the DP model, including the total mass
 * parameter `a`.
 *
 * The parameters are stored as public members for easy access.
 *
 */

struct DP_params : public Process_params {

  /** @brief Total mass parameter of the DP (controls number of clusters) */
  double a;

  /**
   * @brief Constructor for DP_params
   * @param a Total mass parameter
   */
  DP_params(double a) : Process_params("DP"), a(a) {}
};
