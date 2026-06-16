/**
 * @file NGGP-params.hpp
 * @brief Parameter management for models using NGG-process (Normalized
 * Generalized Gamma Process)
 *
 * This file contains the parameter management for models using NGG-process.
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include "Process-params.hpp"
#include <Eigen/Dense>

/**
 * @brief Struct for storing parameters of the NGGP (Normalized Generalized
 * Gamma Process)
 *
 * This struct holds the parameters for the NGGP model, including the total mass
 * parameter `a`, the second parameter `sigma`, and the third parameter `tau`.
 *
 * The parameters are stored as public members for easy access.
 *
 */

struct NGGP_params : public Process_params {

    /** @brief Total mass parameter of the NGGP (controls number of clusters) */
    double a;

    /** @brief Second parameter of the NGGP (controls cluster sizes) */
    double sigma;

    /** @brief Third parameter of the NGGP (controls tail behavior) */
    double tau;

    /**
     * @brief Constructor for NGGP_params
     * @param a Total mass parameter
     * @param sigma Second parameter
     * @param tau Third parameter
     */
    NGGP_params(double a, double sigma, double tau) : Process_params("NGGP"), a(a), sigma(sigma), tau(tau) {}
};
