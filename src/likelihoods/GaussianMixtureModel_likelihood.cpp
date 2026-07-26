/**
 * @file GaussianMixtureModel_likelihood.cpp
 * @brief Gaussian Mixture Model likelihood implementation
 */

#include "GaussianMixtureModel_likelihood.hpp"
#include <cmath>
#include <stdexcept>

GaussianMixtureModel_likelihood::GaussianMixtureModel_likelihood(
    const Data &data, const GaussianMixtureModel_params &params)
    : Likelihood(data), params(params), log_kappa0(std::log(params.kappa0)),
      lgamma_alpha0(std::lgamma(params.alpha0)),
      alpha0_log_beta0(params.alpha0 * std::log(params.beta0)),
      log_pi(std::log(M_PI)) {

  if (params.kappa0 <= 0.0 || params.alpha0 <= 0.0 || params.beta0 <= 0.0) {
    throw std::invalid_argument("GaussianMixtureModel_likelihood: kappa0, "
                                "alpha0, beta0 must be positive");
  }
}

GaussianMixtureModel_likelihood::ClusterStats
GaussianMixtureModel_likelihood::compute_stats(
    const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const {
  ClusterStats stats;

  for (int r = 0; r < cls_ass_k.size(); ++r) {
    const int area_idx = cls_ass_k(r);

    for (int j = 0; j < data.get_p(); ++j) {
      const double val = data.get_data(area_idx, j);
      if (std::isnan(val))
        continue;
      stats.n += 1;
      stats.sum += val;
      stats.sumsq += val * val;
    }
  }

  return stats;
}

double GaussianMixtureModel_likelihood::log_marginal_likelihood(
    const ClusterStats &stats) const {
  const int n = stats.n;
  if (n == 0) {
    return 0.0;
  }

  const double mean = stats.sum / static_cast<double>(n);
  double sse = stats.sumsq - static_cast<double>(n) * mean * mean;
  if (sse < 0.0) {
    sse = 0.0;
  }

  const double kappa_n = params.kappa0 + static_cast<double>(n);
  const double alpha_n = params.alpha0 + 0.5 * static_cast<double>(n);
  const double mean_diff = mean - params.m0;
  const double beta_n = params.beta0 + 0.5 * sse +
                        0.5 *
                            (params.kappa0 * static_cast<double>(n) / kappa_n) *
                            mean_diff * mean_diff;

  double logp = 0.0;
  logp += 0.5 * (log_kappa0 - std::log(kappa_n));
  logp += alpha0_log_beta0 - alpha_n * std::log(beta_n);
  logp += std::lgamma(alpha_n) - lgamma_alpha0;
  logp += -0.5 * static_cast<double>(n) * log_pi;

  return logp;
}

double GaussianMixtureModel_likelihood::log_predictive_from_stats(
    const ClusterStats &stats, double x) const {
  ClusterStats updated = stats;
  updated.n += 1;
  updated.sum += x;
  updated.sumsq += x * x;
  return log_marginal_likelihood(updated) - log_marginal_likelihood(stats);
}

double GaussianMixtureModel_likelihood::cluster_loglikelihood(
    int cluster_index) const {
  auto cls_ass_k = data.get_cluster_assignments_ref(cluster_index);
  return cluster_loglikelihood(cluster_index, cls_ass_k);
}

double GaussianMixtureModel_likelihood::cluster_loglikelihood(
    int cluster_index,
    const Eigen::Ref<const Eigen::VectorXi> &cls_ass_k) const {
  (void)cluster_index;
  const ClusterStats stats = compute_stats(cls_ass_k);
  return log_marginal_likelihood(stats);
}

double GaussianMixtureModel_likelihood::point_loglikelihood_cond(
    int point_index, int cluster_index) const {
  ClusterStats point_stats;

  for (int j = 0; j < data.get_p(); ++j) {
    const double x = data.get_data(point_index, j);
    point_stats.n += 1;
    point_stats.sum += x;
    point_stats.sumsq += x * x;
  }

  if (cluster_index >= data.get_K()) {
    ClusterStats empty;
    return log_marginal_likelihood(point_stats) -
           log_marginal_likelihood(empty);
  }

  auto cls_ass_k = data.get_cluster_assignments_ref(cluster_index);
  ClusterStats cluster_stats = compute_stats(cls_ass_k);

  ClusterStats updated = cluster_stats;
  updated.n += point_stats.n;
  updated.sum += point_stats.sum;
  updated.sumsq += point_stats.sumsq;

  return log_marginal_likelihood(updated) -
         log_marginal_likelihood(cluster_stats);
}
