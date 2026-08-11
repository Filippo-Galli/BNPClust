#include "GaussianMixtureModel_likelihood.hpp"
#include <cmath>
#include <stdexcept>

GaussianMixtureModel_likelihood::GaussianMixtureModel_likelihood(
    const Data &data, const GaussianMixtureModel_params &params)
    : Likelihood(data), params(params), log_pi(std::log(M_PI)),
      lgamma_alpha0(std::lgamma(params.alpha0)) {
  if (params.kappa0 <= 0.0 || params.alpha0 <= 0.0 || params.beta0 <= 0.0) {
    throw std::invalid_argument("GaussianMixtureModel_likelihood: kappa0, "
                                "alpha0, beta0 must be positive");
  }
}

GaussianMixtureModel_likelihood::ClusterStats
GaussianMixtureModel_likelihood::compute_area_stats(int area_index) const {
  ClusterStats stats;
  for (int j = 0; j < data.get_p(); ++j) {
    const double x = data.get_data(area_index, j);
    if (std::isnan(x))
      continue;
    stats.n += 1;
    stats.sum += x;
    stats.sumsq += x * x;
  }
  return stats;
}

GaussianMixtureModel_likelihood::ClusterStats
GaussianMixtureModel_likelihood::compute_stats(
    const Eigen::Ref<const Eigen::VectorXi> &cluster_areas) const {
  ClusterStats stats;
  for (int r = 0; r < cluster_areas.size(); ++r) {
    const int area_idx = cluster_areas(r);
    ClusterStats a = compute_area_stats(area_idx);
    stats.n += a.n;
    stats.sum += a.sum;
    stats.sumsq += a.sumsq;
  }
  return stats;
}

double GaussianMixtureModel_likelihood::log_marginal_likelihood(
    const ClusterStats &stats) const {
  if (stats.n == 0)
    return 0.0;

  const double n = static_cast<double>(stats.n);
  const double mean = stats.sum / n;
  double sse = stats.sumsq - n * mean * mean;
  if (sse < 0.0)
    sse = 0.0;

  const double kappa_n = params.kappa0 + n;
  const double alpha_n = params.alpha0 + 0.5 * n;
  const double mean_diff = mean - params.m0;
  const double beta_n =
      params.beta0 + 0.5 * sse +
      0.5 * (params.kappa0 * n / kappa_n) * mean_diff * mean_diff;

  double logp = 0.0;
  logp += 0.5 * (std::log(params.kappa0) - std::log(kappa_n));
  logp += params.alpha0 * std::log(params.beta0) - alpha_n * std::log(beta_n);
  logp += std::lgamma(alpha_n) - lgamma_alpha0;
  logp += -0.5 * n * log_pi;
  return logp;
}

double GaussianMixtureModel_likelihood::log_predictive_from_stats(
    const ClusterStats &stats, const ClusterStats &area_stats) const {
  ClusterStats updated = stats;
  updated.n += area_stats.n;
  updated.sum += area_stats.sum;
  updated.sumsq += area_stats.sumsq;
  return log_marginal_likelihood(updated) - log_marginal_likelihood(stats);
}

double GaussianMixtureModel_likelihood::cluster_loglikelihood(
    int cluster_index) const {
  auto cluster_areas = data.get_cluster_assignments_ref(cluster_index);
  return cluster_loglikelihood(cluster_index, cluster_areas);
}

double GaussianMixtureModel_likelihood::cluster_loglikelihood(
    int cluster_index,
    const Eigen::Ref<const Eigen::VectorXi> &cluster_areas) const {
  (void)cluster_index;
  const ClusterStats stats = compute_stats(cluster_areas);
  return log_marginal_likelihood(stats);
}

double GaussianMixtureModel_likelihood::point_loglikelihood_cond(
    int point_index, int cluster_index) const {
  ClusterStats area_stats = compute_area_stats(point_index);

  if (cluster_index >= data.get_K()) {
    ClusterStats empty;
    return log_marginal_likelihood(area_stats) - log_marginal_likelihood(empty);
  }

  auto cluster_areas = data.get_cluster_assignments_ref(cluster_index);
  const ClusterStats cluster_stats = compute_stats(cluster_areas);
  return log_predictive_from_stats(cluster_stats, area_stats);
}
