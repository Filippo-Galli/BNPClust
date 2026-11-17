/**
 * @file Natarajan_likelihood_summaryStats.cpp
 * @brief source file for Natarajan_likelihood_summaryStats
 */

#include "Natarajan_likelihood_summaryStats.hpp"
#include <cmath>

double Natarajan_likelihood_summaryStats::cluster_loglikelihood(
    int cluster_index) const {
  auto cls_ass_k = data.get_cluster_assignments_span(cluster_index);
  return cluster_loglikelihood(cluster_index, cls_ass_k);
}

double Natarajan_likelihood_summaryStats::cluster_loglikelihood(
    int cluster_index,
    std::span<const int> cls_ass_k) const {
  const int n_k = static_cast<int>(cls_ass_k.size());

  if (n_k == 0) {
    return 0;
  }

  double rep = 0;
  const int K = data.get_K();
  const auto D_view = utils.get_D_view();

  /* -------------------- Repulsion part -------------------------- */
  for (int t = 0; t < K; ++t) {
    if (t == cluster_index)
      continue;

    auto cls_ass_t = data.get_cluster_assignments_span(t);
    const int n_t = static_cast<int>(cls_ass_t.size());

    if (n_t == 0)
      continue;

    double sum = 0;

    for (int idx_i : cls_ass_k) {
      for (int idx_j : cls_ass_t) {
        sum += D_view[idx_i, idx_j];
      }
    }

    double emp_mean = sum / (n_t * n_k);

    rep += std::log(std::max(emp_mean, 1e-12)) * (params.delta2 - 1);
    rep -= lgamma_delta2;
    rep += log_gamma_zeta;
    rep += lgamma(params.delta2 + params.zeta);
    rep -= log(params.gamma + sum) * (params.delta2 + params.zeta);
  }

  /* -------------------- Cohesion part -------------------------- */
  if (n_k == 1) {
    return rep;
  }

  const int pairs = n_k * (n_k - 1) / 2;
  double log_prod = 0;
  double sum = 0;

  for (int i = 0; i < n_k; ++i) {
    const int idx_i = cls_ass_k[i];
    for (int j = i + 1; j < n_k; ++j) {
      const int idx_j = cls_ass_k[j];
      sum += D_view[idx_i, idx_j];
      log_prod += logD_view[idx_i, idx_j];
    }
  }

  double coh = 0;
  coh += log_prod * (params.delta1 - 1);
  coh -= lgamma_delta1 * pairs;
  coh += log_beta_alpha;
  coh += lgamma(pairs * params.delta1 + params.alpha);
  coh -= log(params.beta + sum) * (pairs * params.delta1 + params.alpha);

  return rep + coh;
}

double Natarajan_likelihood_summaryStats::point_loglikelihood_cond(
    int point_index, int cluster_index) const {
  auto cls_ass_k = data.get_cluster_assignments_span(cluster_index);
  const int n_k = static_cast<int>(cls_ass_k.size());

  double coeh = compute_cohesion(point_index, cluster_index, cls_ass_k, n_k);
  double rep = compute_repulsion(point_index, cluster_index, cls_ass_k, n_k);

  return coeh + rep;
}

double Natarajan_likelihood_summaryStats::compute_cohesion(
    int point_index, int cluster_index,
    std::span<const int> cls_ass_k, int n_k) const {
  if (n_k == 0) {
    return 0.0;
  }

  const auto D_view = utils.get_D_view();

  double sum_i = 0;
  double log_prod_i = 0;

  for (int idx : cls_ass_k) {
    sum_i += D_view[point_index, idx];
    log_prod_i += logD_view[point_index, idx];
  }

  const double alpha_mh = params.alpha + params.delta1 * n_k;
  const double beta_mh = params.beta + sum_i;

  double loglik = 0;
  loglik += (-n_k) * lgamma_delta1;
  loglik += (params.delta1 - 1) * log_prod_i;
  loglik += lgamma_alpha_mh_cache[n_k];
  loglik += log_beta_alpha;
  loglik -= alpha_mh * log(beta_mh);

  return loglik;
}

double Natarajan_likelihood_summaryStats::compute_repulsion(
    int point_index, int cluster_index,
    std::span<const int> cls_ass_k, int n_k) const {
  const int num_cluster = data.get_K();

  if (num_cluster <= 1) {
    return 0;
  }

  double loglik = 0;
  const auto D_view = utils.get_D_view();

  for (int t = 0; t < num_cluster; ++t) {
    if (t == cluster_index)
      continue;

    const auto cls_ass_t = data.get_cluster_assignments_span(t);
    const int n_t = static_cast<int>(cls_ass_t.size());

    if (n_t == 0)
      continue;

    double sum_i = 0;

    for (int idx : cls_ass_t) {
      sum_i += D_view[point_index, idx];
    }

    const double zeta_mt = params.zeta + params.delta2;
    const double gamma_mt = params.gamma + sum_i;
    double emp_mean = sum_i / n_t;

    loglik -= lgamma_delta2;
    loglik += (params.delta2 - 1) * std::log(std::max(emp_mean, 1e-12));
    loglik += lgamma_zeta_mt_cache;
    loglik += log_gamma_zeta;
    loglik -= zeta_mt * log(gamma_mt);
  }

  return loglik;
}

double Natarajan_likelihood_summaryStats::pairwise_repulsion(int cluster_1,
                                                             int cluster_2) const {
  if (cluster_1 == cluster_2)
    return 0.0;

  auto cls_ass_k = data.get_cluster_assignments_span(cluster_1);
  auto cls_ass_t = data.get_cluster_assignments_span(cluster_2);
  const int n_k = static_cast<int>(cls_ass_k.size());
  const int n_t = static_cast<int>(cls_ass_t.size());

  if (n_k == 0 || n_t == 0)
    return 0.0;

  const auto D_view = utils.get_D_view();

  double sum = 0.0;

  for (int idx_i : cls_ass_k) {
    for (int idx_j : cls_ass_t) {
      sum += D_view[idx_i, idx_j];
    }
  }

  double emp_mean = sum / (n_t * n_k);

  double rep = 0.0;
  rep += std::log(std::max(emp_mean, 1e-12)) * (params.delta2 - 1);
  rep -= lgamma_delta2;
  rep += log_gamma_zeta;
  rep += lgamma(params.delta2 + params.zeta);
  rep -= log(params.gamma + sum) * (params.delta2 + params.zeta);

  return rep;
}

double Natarajan_likelihood_summaryStats::clusters_loglikelihood(int c1, int c2) const {
  return cluster_loglikelihood(c1) + cluster_loglikelihood(c2) - pairwise_repulsion(c1, c2);
}

