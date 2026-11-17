/**
 * @file splitmerge.cpp
 * @brief Implementation of Split-Merge sampler using restricted Gibbs sampling
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge.hpp"
#include <cmath>

void SplitMerge::choose_indeces() { choose_indeces_uniform(false); }

void SplitMerge::restricted_gibbs(int iterations, bool only_probabilities) {
  for (int i = 0; i < iterations; ++i) {
    for (size_t idx = 0; idx < S.size(); ++idx) {
      int point_idx = S[idx];
      int current_cluster = data.get_cluster_assignment(point_idx);

      // Temporarily unassign point
      data.set_allocation(point_idx, -1);

      // Compute probabilities for each cluster (ci and cj)
      const double log_prob_0 =
          likelihood.point_loglikelihood_cond(point_idx, ci) +
          process.gibbs_prior_existing_cluster(ci, point_idx);
      const double log_prob_1 =
          likelihood.point_loglikelihood_cond(point_idx, cj) +
          process.gibbs_prior_existing_cluster(cj, point_idx);

      const auto [log_sum, prob_cj] =
          compute_binary_choice_probs(log_prob_0, log_prob_1);

      if (!only_probabilities) {
        const int new_cluster = (dis_real(gen) < prob_cj) ? cj : ci;
        data.set_allocation(point_idx, new_cluster);

        if (i == iterations - 1) {
          log_split_gibbs_prob +=
              (new_cluster == ci ? log_prob_0 : log_prob_1) - log_sum;
        }
      } else {
        data.set_allocation(point_idx, current_cluster);
        log_merge_gibbs_prob +=
            (current_cluster == ci ? log_prob_0 : log_prob_1) - log_sum;
      }
    }
  }
}

void SplitMerge::split_move() {
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);
  log_split_gibbs_prob = 0.0;

  data.set_allocation(idx_j, data.get_K());
  cj = data.get_cluster_assignment(idx_j);

  for (size_t idx = 0; idx < launch_state.size(); ++idx) {
    const int new_cluster = (dis_real(gen) < 0.5) ? ci : cj;
    data.set_allocation(S[idx], new_cluster);
  }

  restricted_gibbs(1);

  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  if (std::log(dis_real(gen)) > acceptance_ratio) {
    rollback();
  } else {
    accepted_split++;
  }
}

void SplitMerge::merge_move() {
  log_merge_gibbs_prob = 0.0;
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  restricted_gibbs(1, true);

  data.set_allocation(idx_j, ci);

  for (size_t idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state[idx] == cj) {
      data.set_allocation(S[idx], ci);
    }
  }

  double acceptance_ratio =
      compute_acceptance_ratio_merge(likelihood_old_clusters);

  if (std::log(dis_real(gen)) > acceptance_ratio) {
    rollback();
  } else {
    accepted_merge++;
  }
}

void SplitMerge::shuffle() {
  log_split_gibbs_prob = 0.0;
  log_merge_gibbs_prob = 0.0;

  if (data.get_K() < 2)
    return;

  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  restricted_gibbs(1, true);
  restricted_gibbs(1, false);

  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(
      likelihood_old_clusters, old_ci_size, old_cj_size);

  if (std::log(dis_real(gen)) > log_acceptance_ratio) {
    rollback();
  } else {
    accepted_shuffle++;
  }
}
