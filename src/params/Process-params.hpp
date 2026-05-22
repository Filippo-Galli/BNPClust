/**
 * @file Process-params.hpp
 * @brief Abstract base class for process parameters
 *
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

/**
 * @brief Abstract base class for process parameters
 *
 * This class serves as public members for easy access.
 *
 */

struct Process_params {
  const std::string name;

  /**
   * @brief Constructor for Process_params
   * @param name Name of the process
   */
  Process_params(const std::string &name) : name(name) {}
  virtual ~Process_params() {}
};
