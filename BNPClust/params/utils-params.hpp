/**
 * @file utils-params.hpp
 * @brief Parameter management for useful parameters and data structures
 *
 * @author Filippo Galli
 * @date 2026
 */

#pragma once

#include <Eigen/Dense>
#include <mdspan>

/**
 * @brief Parameter struct for utility parameters and data structures.
 */

struct utils_params {
    using MatrixConstView = std::mdspan<const double, std::dextents<size_t, 2>, std::layout_left>;

    /** @brief Number of burn-in iterations to discard for chain convergence */
    int BI;

    /** @brief Number of iterations after burn-in for posterior sampling */
    int NI;

    /** @brief Data matrix */
    Eigen::MatrixXd D;

    /** @brief Number of points */
    int n;

    /** @brief Number of dimensions */
    int p;

    /**
     * @brief Constructor for utils_params
     * @param BI Number of burn-in iterations
     * @param NI Number of iterations after burn-in
     * @param D Data matrix
     */
    utils_params(int BI, int NI, const Eigen::MatrixXd &D) : BI(BI), NI(NI), D(D), n(D.rows()), p(D.cols()) {}

    /**
     * @brief Zero-overhead 2D non-owning mdspan view of data matrix D
     */
    [[nodiscard]] MatrixConstView get_D_view() const noexcept {
        return MatrixConstView(D.data(), static_cast<size_t>(n), static_cast<size_t>(p));
    }
};
