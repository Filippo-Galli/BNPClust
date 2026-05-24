/**
 * @file utils-params.hpp
 * @brief Parameter management for useful parameters and data structures
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include <Eigen/Dense>

struct utils_params {
    /** @brief Number of burn-in iterations to discard for chain convergence */
    int BI;

    /** @brief Number of iterations after burn-in for posterior sampling */
    int NI;

    /** @brief Distance matrix */
    Eigen::MatrixXd D;

    /** @brief Number of points */
    int n;

    /** @brief Number of dimensions */
    int p;

    /**
     * @brief Constructor for utils_params
     * @param BI Number of burn-in iterations
     * @param NI Number of iterations after burn-in
     * @param D Distance matrix
     */
    utils_params(int BI, int NI, const Eigen::MatrixXd &D) : BI(BI), NI(NI), D(D), n(D.rows()), p(D.cols()) {}
};
