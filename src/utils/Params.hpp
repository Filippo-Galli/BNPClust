/**
 * @file Params.hpp
 * @brief Parameter management for Bayesian nonparametric MCMC models
 *
 * This file defines the Params structure that centralizes all hyperparameters,
 * model configuration, and MCMC settings for DP, NGGP, and spatial variants.
 * It provides a unified interface for parameter management across different
 * algorithms and process types.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include <Eigen/Dense>

/**
 * @brief Structure containing all parameters needed for the NGGP (Normalized
 * Generalized Gamma Process) and DP (Dirichlet Process) model
 *
 * This structure holds all the hyperparameters, model parameters, and
 * configuration settings required for running the DP and NGGP Bayesian
 * nonparametric model with and without spatial dependency.
 *
 * @details The parameters are organized into several categories:
 * - Distribution hyperparameters for gamma priors (delta1, alpha, beta, delta2,
 * gamma, zeta)
 * - Simulation control parameters (BI, NI)
 * - NGGP process parameters (a, sigma, tau)
 * - Spatial dependency parameters (coefficient, W)
 */
struct Params {
    // ========== Distribution Hyperparameters ==========

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

    // ========== Simulation Control Parameters ==========

    /** @brief Number of burn-in iterations to discard for chain convergence */
    int BI;

    /** @brief Number of iterations after burn-in for posterior sampling */
    int NI;

    // ========== NGGP Process Parameters ==========

    /** @brief Total mass parameter of the NGGP (controls number of clusters) */
    double a;

    /** @brief Second parameter of the NGGP (controls cluster sizes) */
    double sigma;

    /** @brief Third parameter of the NGGP (controls tail behavior) */
    double tau;

    /** @brief Distance matrix */
    Eigen::MatrixXd D;

    /** @brief Number of points */
    int n;

    /**
     * @brief Constructor with default parameter values
     *
     * @param delta1 Shape parameter for first gamma prior (default: 0.5)
     * @param alpha Shape parameter for lambda_k gamma distribution (default: 2)
     * @param beta Rate parameter for lambda_k gamma distribution (default: 2)
     * @param delta2 Shape parameter for second gamma prior (default: 2)
     * @param gamma Shape parameter for theta_kt gamma distribution (default: 2)
     * @param zeta Rate parameter for theta_kt gamma distribution (default: 2)
     * @param BI Number of burn-in iterations (default: 1000)
     * @param NI Number of post burn-in iterations (default: 10000)
     * @param a NGGP total mass parameter (default: 1.0)
     * @param sigma NGGP second parameter (default: 1.0)
     * @param tau NGGP third parameter (default: 1.0)
     * @param D Distance matrix (default: empty matrix)
     */
    Params(double delta1 = 0.5, double alpha = 2, double beta = 2, double delta2 = 2, double gamma = 2, double zeta = 2,
           int BI = 1000, int NI = 10000, double a = 1.0, double sigma = 1.0, double tau = 1.0,
           Eigen::MatrixXd D = Eigen::MatrixXd())
        : delta1(delta1), alpha(alpha), beta(beta), delta2(delta2), gamma(gamma), zeta(zeta), BI(BI), NI(NI), a(a),
          sigma(sigma), tau(tau), D(D) {
        n = D.rows();
    }
};
