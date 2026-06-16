/**
 * @file Natarajan_likelihood.hpp
 * @brief LogLikelihood computation for clustering with cohesion and repulsion
 */

#pragma once

#include "../params/Natarajan-params.hpp"
#include "../params/utils-params.hpp"
#include "../utils/Likelihood.hpp"
#include <vector>

/**
 * @class Natarajan_likelihood
 * @brief Computes log-likelihood for clusters based on distance-based cohesion
 * and repulsion
 *
 * This class handles the computation of log-likelihoods for clustering models
 * that incorporate both cohesion (within-cluster similarity) and repulsion
 * (between-cluster dissimilarity) components. It precomputes values for
 * efficiency and provides methods to evaluate both cluster-level and
 * point-level conditional log-likelihoods.
 * reference: Natarajan et al. (2023) "Cohesion and Repulsion in Bayesian
 * Distance Clustering"
 */
class Natarajan_likelihood : public Likelihood {
private:
  const Natarajan_params &params; ///< Reference to model parameters
  const utils_params &utils;      ///< Reference to utility parameters

  // Precomputed values for efficiency
  const double lgamma_delta1; ///< Precomputed lgamma(delta1) for cohesion
  const double
      log_beta_alpha;         ///< Precomputed log(beta) * alpha - lgamma(alpha)
  const double lgamma_delta2; ///< Precomputed lgamma(delta2) for repulsion
  const double log_gamma_zeta; ///< Precomputed log(gamma) * zeta - lgamma(zeta)
  std::vector<double>
      lgamma_zeta_mt_cache; ///< Cache for lgamma(zeta_mt) values
  std::vector<double>
      lgamma_alpha_mh_cache; ///< Cache for lgamma(alpha_mh) values

  std::vector<double>
      log_D_data;   ///< Precomputed log distance matrix (flattened)
  const int D_cols; ///< Number of columns in distance matrix

  /**
   * @brief Computes the cohesion component of the log-likelihood
   * @param point_index Index of the point being evaluated
   * @param cluster_index Index of the cluster
   * @param cls_ass_k Vector of point indices in the cluster
   * @param n_k Number of points in the cluster
   * @return Cohesion log-likelihood contribution
   */
  double compute_cohesion(int point_index, int cluster_index,
                          const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k,
                          int n_k) const;

  /**
   * @brief Computes the repulsion component of the log-likelihood
   * @param point_index Index of the point being evaluated
   * @param cluster_index Index of the cluster
   * @param cls_ass_k Vector of point indices in the cluster
   * @param n_k Number of points in the cluster
   * @return Repulsion log-likelihood contribution
   */
  double compute_repulsion(int point_index, int cluster_index,
                           const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k,
                           int n_k) const;

public:
  /**
   * @brief Constructs a Likelihood object with precomputation
   * @param data Reference to Data object with distances and allocations
   * @param param Reference to model parameters
   *
   * The constructor precomputes several values for computational efficiency:
   * - Log-gamma values for delta parameters
   * - Logarithmic combinations of hyperparameters
   * - Logarithm of the entire distance matrix
   */
  Natarajan_likelihood(const Data &data, const Natarajan_params &param,
                       const utils_params &utils)
      : Likelihood(data), params(param), utils(utils),
        lgamma_delta1(lgamma(params.delta1)),
        log_beta_alpha(log(params.beta) * params.alpha - lgamma(params.alpha)),
        lgamma_delta2(lgamma(params.delta2)),
        log_gamma_zeta(log(params.gamma) * params.zeta - lgamma(params.zeta)),
        D_cols(utils.D.cols()) {
    // Precompute log distances manually
    const int n = utils.D.rows() * utils.D.cols();
    log_D_data.resize(n);

    // Log distance precomputation
    const double *D_ptr = utils.D.data();
    for (int i = 0; i < n; ++i) {
      log_D_data[i] = std::log(D_ptr[i]);
    }

    // Cache lgamma values for efficiency
    lgamma_alpha_mh_cache.resize(data.get_n() + 1, 0.0);
    lgamma_zeta_mt_cache.resize(data.get_n() + 1, 0.0);

    for (int val = 0; val <= data.get_n(); ++val) {
      lgamma_alpha_mh_cache[val] = lgamma(params.alpha + params.delta1 * val);
      lgamma_zeta_mt_cache[val] = lgamma(params.zeta + params.delta2 * val);
    }
  }

  /**
   * @brief Computes the full log-likelihood for a cluster
   * @param cluster_index Index of the cluster to evaluate
   * @return Total log-likelihood (cohesion + repulsion)
   *
   * This method computes both the within-cluster cohesion and the
   * between-cluster repulsion contributions for the specified cluster.
   */
  double cluster_loglikelihood(int cluster_index) const override final;

  /**
   * @brief Computes the full log-likelihood for a cluster with given
   * assignments
   * @param cluster_index Index of the cluster to evaluate
   * @param cls_ass_k Vector of point indices in the cluster
   * @return Total log-likelihood (cohesion + repulsion)
   *
   * This overload allows computing the likelihood with a custom set of
   * cluster assignments without modifying the data structure.
   */
  double cluster_loglikelihood(
      int cluster_index,
      const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const override final
      __attribute__((hot));

  /**
   * @brief Computes the conditional log-likelihood of a point given a cluster
   * @param point_index Index of the point to evaluate
   * @param cluster_index Index of the target cluster
   * @return Conditional log-likelihood of assigning the point to the cluster
   *
   * This method evaluates how well a point fits into a specific cluster,
   * considering both its cohesion with points in that cluster and its
   * repulsion from points in other clusters.
   */
  double point_loglikelihood_cond(int point_index,
                                  int cluster_index) const override final
      __attribute__((hot));
};
