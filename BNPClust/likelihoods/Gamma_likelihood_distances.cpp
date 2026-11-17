/**
 * @file Gamma_likelihood_distances.cpp
 * @brief Fully optimized implementation avoiding all BLAS calls
 */

#include "Gamma_likelihood_distances.hpp"
#include <cmath>

double
Gamma_likelihood_distances::cluster_loglikelihood(int cluster_index) const {
  auto cls_ass_k = data.get_cluster_assignments_span(cluster_index);
  return cluster_loglikelihood(cluster_index, cls_ass_k);
}

double Gamma_likelihood_distances::cluster_loglikelihood(
    int cluster_index,
    std::span<const int> cls_ass_k) const {
  const int n_k = static_cast<int>(cls_ass_k.size());

  if (n_k == 0) {
    return 0;
  }

  /* -------------------- Cohesion part -------------------------- */
  if (n_k == 1) {
    return 0;
  }

  const auto D_view = utils.get_D_view();
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

  return coh;
}

double
Gamma_likelihood_distances::point_loglikelihood_cond(int point_index,
                                                     int cluster_index) const {
  auto cls_ass_k = data.get_cluster_assignments_span(cluster_index);
  const int n_k = static_cast<int>(cls_ass_k.size());

  double coeh = compute_cohesion(point_index, cluster_index, cls_ass_k, n_k);

  return coeh;
}

double Gamma_likelihood_distances::compute_cohesion(
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
