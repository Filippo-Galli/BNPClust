/**
 * @file splitmerge_base.cpp
 * @brief Implementation of shared SplitMergeBase functionality
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_base.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

void SplitMergeBase::prepare_launch_state(int cluster_i, int cluster_j) {
  ci = cluster_i;
  cj = cluster_j;

  // size of launch state retrieval
  const int size_ci = data.get_cluster_size(ci);
  const int size_cj = data.get_cluster_size(cj);
  const size_t launch_state_size =
      (ci == cj) ? static_cast<size_t>(std::max(0, size_ci - 2))
                 : static_cast<size_t>(std::max(0, size_ci + size_cj - 2));

  launch_state.clear();
  launch_state.reserve(launch_state_size);
  S.clear();
  S.reserve(launch_state_size);

  const auto init_alloc = data.get_allocations();
  original_allocations.assign(init_alloc.begin(), init_alloc.end());

  // collect non-anchor points in launch state of the two clusters
  for (int idx : data.get_cluster_assignments(ci)) {
    if (idx != idx_i && idx != idx_j) {
      S.push_back(idx);
      launch_state.push_back(ci);
    }
  }

  if (ci != cj) {
    for (int idx : data.get_cluster_assignments(cj)) {
      if (idx != idx_i && idx != idx_j) {
        S.push_back(idx);
        launch_state.push_back(cj);
      }
    }
  }
}

void SplitMergeBase::shuffle_launch_state() {
  const size_t sz = S.size();
  if (sz > 1) {
    for (size_t i = sz - 1; i > 0; --i) {
      const size_t j = dis_size(
          gen, std::uniform_int_distribution<size_t>::param_type(0, i));
      if (i != j) {
        std::swap(launch_state[i], launch_state[j]);
        std::swap(S[i], S[j]);
      }
    }
  }
}

void SplitMergeBase::choose_indeces_uniform(bool shuffle_s) {
  idx_i = dis_n(gen);
  do {
    idx_j = dis_n(gen);
  } while (idx_j == idx_i);

  ci = data.get_cluster_assignment(idx_i);
  cj = data.get_cluster_assignment(idx_j);

  prepare_launch_state(ci, cj);
  if (shuffle_s) {
    shuffle_launch_state();
  }
}

void SplitMergeBase::sequential_allocation(int iterations,
                                           bool only_probabilities,
                                           bool sequential) {
  const size_t S_size = S.size();

  for (int i = 0; i < iterations; ++i) {
    for (size_t idx = 0; idx < S_size && sequential; ++idx) {
      data.set_allocation(S[idx], -1);
    }

    for (size_t idx = 0; idx < S_size; ++idx) {
      int point_idx = S[idx];
      int current_cluster = launch_state[idx];

      if (!sequential) {
        current_cluster = data.get_cluster_assignment(point_idx);
        data.set_allocation(point_idx, -1);
      }

      const double log_probs_0 =
          likelihood.point_loglikelihood_cond(point_idx, ci) +
          process.gibbs_prior_existing_cluster(ci, point_idx);
      const double log_probs_1 =
          likelihood.point_loglikelihood_cond(point_idx, cj) +
          process.gibbs_prior_existing_cluster(cj, point_idx);

      const auto [log_sum, prob_cj] =
          compute_binary_choice_probs(log_probs_0, log_probs_1);

      if (!only_probabilities) {
        const int new_cluster = (dis_real(gen) < prob_cj) ? cj : ci;
        data.set_allocation(point_idx, new_cluster);

        if (i == iterations - 1) {
          log_split_gibbs_prob +=
              (new_cluster == ci ? log_probs_0 : log_probs_1) - log_sum;
        }
      } else {
        data.set_allocation(point_idx, current_cluster);
        log_merge_gibbs_prob +=
            (current_cluster == ci ? log_probs_0 : log_probs_1) - log_sum;
      }
    }
  }
}

void SplitMergeBase::choose_clusters_shuffle() {
  if (data.get_K() < 2)
    return;

  std::uniform_int_distribution<> dis(0, data.get_K() - 1);
  ci = dis(gen);
  do {
    cj = dis(gen);
  } while (cj == ci);

  std::uniform_int_distribution<> dis_idx_i(0, data.get_cluster_size(ci) - 1);
  idx_i = data.get_cluster_assignments(ci)[dis_idx_i(gen)];

  std::uniform_int_distribution<> dis_idx_j(0, data.get_cluster_size(cj) - 1);
  idx_j = data.get_cluster_assignments(cj)[dis_idx_j(gen)];

  prepare_launch_state(ci, cj);
}

double
SplitMergeBase::compute_acceptance_ratio_split(double likelihood_old_cluster) {
  double log_acceptance_ratio = process.prior_ratio_split(ci, cj);

  log_acceptance_ratio +=
      likelihood.clusters_loglikelihood(ci, cj) - likelihood_old_cluster;
  log_acceptance_ratio -= log_split_gibbs_prob;

  return log_acceptance_ratio;
}

double
SplitMergeBase::compute_acceptance_ratio_merge(double likelihood_old_clusters) {
  int size_old_ci = static_cast<int>(
      std::count(original_allocations.begin(), original_allocations.end(), ci));
  int size_old_cj = static_cast<int>(
      std::count(original_allocations.begin(), original_allocations.end(), cj));
  double log_acceptance_ratio =
      process.prior_ratio_merge(size_old_ci, size_old_cj);

  const int merged_ci = data.get_cluster_assignment(idx_i);
  log_acceptance_ratio +=
      likelihood.cluster_loglikelihood(merged_ci) - likelihood_old_clusters;
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

double SplitMergeBase::compute_acceptance_ratio_shuffle(
    double likelihood_old_clusters, int old_ci_size, int old_cj_size) {
  double log_acceptance_ratio =
      process.prior_ratio_shuffle(old_ci_size, old_cj_size, ci, cj);

  log_acceptance_ratio +=
      likelihood.clusters_loglikelihood(ci, cj) - likelihood_old_clusters;
  log_acceptance_ratio -= log_split_gibbs_prob;
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMergeBase::step() {
  choose_indices_for_step();
  process.set_old_allocations(data.get_allocations());
  process.set_idx_i(idx_i);
  process.set_idx_j(idx_j);

  if (ci == cj) {
    split_move();
  } else {
    merge_move();
  }

  if (shuffle_bool) {
    choose_clusters_shuffle();
    process.set_old_allocations(data.get_allocations());
    process.set_idx_i(idx_i);
    process.set_idx_j(idx_j);
    shuffle();
  }
}

void SplitMergeBase::shuffle() {
  log_split_gibbs_prob = 0.0;
  log_merge_gibbs_prob = 0.0;

  if (data.get_K() < 2)
    return;

  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  sequential_allocation(1, true, true);
  sequential_allocation(1, false, true);

  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(
      likelihood_old_clusters, old_ci_size, old_cj_size);

  if (std::log(dis_real(gen)) > log_acceptance_ratio) {
    rollback();
  } else {
    accepted_shuffle++;
  }
}
