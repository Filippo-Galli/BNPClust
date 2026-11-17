/**
 * @file splitmerge_LSS_SDDS.cpp
 * @brief Implementation of Locality Sensitive Sampling (LSS) with SDDS Split-Merge sampler
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_LSS_SDDS.hpp"
#include <algorithm>
#include <cmath>

SplitMerge_LSS_SDDS::SplitMerge_LSS_SDDS(Data &d, utils_params &p,
                                         Likelihood &l, Process &pr,
                                         bool shuffle)
    : SplitMergeBase(d, l, pr, shuffle), params(p),
      original_allocations_ref(pr.get_old_allocations()),
      probs(d.get_n()) {

  const int n = data.get_n();
  row_max_distances.resize(n);
  for (int i = 0; i < n; ++i) {
    row_max_distances[i] = params.D.row(i).maxCoeff();
  }

  // Precompute 50th percentile of pairwise distances as temperature heuristic
  std::vector<double> pairwise_dists;
  pairwise_dists.reserve(n * (n - 1) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      pairwise_dists.push_back(params.D(i, j));
    }
  }

  if (!pairwise_dists.empty()) {
    auto q_it = pairwise_dists.begin() +
                static_cast<size_t>(pairwise_dists.size() * 0.50);
    std::nth_element(pairwise_dists.begin(), q_it, pairwise_dists.end());
    temperature = *q_it;
    if (temperature < 1e-12)
      temperature = 1.0;
  } else {
    temperature = 1.0;
  }
}

void SplitMerge_LSS_SDDS::choose_indeces(bool similarity) {
  idx_i = dis_n(gen);

  const auto distances = params.D.row(idx_i);
  const double max_distance = row_max_distances[idx_i];
  const int n = data.get_n();
  double prob_sum = 0.0;

  for (int idx = 0; idx < n; ++idx) {
    if (idx == idx_i) {
      probs[idx] = 0.0;
    } else {
      if (similarity) {
        probs[idx] = 1.0 / (1.0 + distances(idx) / temperature);
      } else {
        probs[idx] = 1.0 / (1.0 + (max_distance - distances(idx)) / temperature);
      }
    }
    prob_sum += probs[idx];
  }

  double u = dis_real(gen) * prob_sum;
  idx_j = n - 1;
  for (int idx = 0; idx < n; ++idx) {
    u -= probs[idx];
    if (u <= 0.0) {
      idx_j = idx;
      break;
    }
  }

  ci = data.get_cluster_assignment(idx_i);
  cj = data.get_cluster_assignment(idx_j);

  prepare_launch_state(ci, cj);
  shuffle_launch_state();
}

void SplitMerge_LSS_SDDS::smart_split_move() {
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);
  log_split_gibbs_prob = 0.0;

  data.set_allocation(idx_j, data.get_K());
  cj = data.get_cluster_assignment(idx_j);

  for (size_t idx = 0; idx < launch_state.size(); ++idx) {
    const int new_cluster = (dis_real(gen) < 0.5) ? ci : cj;
    data.set_allocation(S[idx], new_cluster);
  }

  sequential_allocation(1);

  double acceptance_ratio = compute_acceptance_ratio_split(likelihood_old_cluster);

  if (std::log(dis_real(gen)) > acceptance_ratio)
    rollback();
  else
    accepted_split++;
}

void SplitMerge_LSS_SDDS::dumb_split_move() {
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);
  log_split_gibbs_prob = S.size() * rand_split_prob;

  data.set_allocation(idx_j, data.get_K());
  cj = data.get_cluster_assignment(idx_j);

  for (size_t idx = 0; idx < launch_state.size(); ++idx) {
    const int new_cluster = (dis_real(gen) < 0.5) ? ci : cj;
    data.set_allocation(S[idx], new_cluster);
  }

  double acceptance_ratio = compute_acceptance_ratio_split(likelihood_old_cluster);

  if (std::log(dis_real(gen)) > acceptance_ratio)
    rollback();
  else
    accepted_split++;
}

void SplitMerge_LSS_SDDS::smart_merge_move() {
  log_merge_gibbs_prob = S.size() * rand_split_prob;
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  data.set_allocation(idx_j, ci);
  for (const auto &idx : S) {
    data.set_allocation(idx, ci);
  }

  double acceptance_ratio = compute_acceptance_ratio_merge(likelihood_old_clusters);

  if (std::log(dis_real(gen)) > acceptance_ratio)
    rollback();
  else
    accepted_merge++;
}

void SplitMerge_LSS_SDDS::dumb_merge_move() {
  log_merge_gibbs_prob = 0.0;
  sequential_allocation(1, true, true);
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  data.set_allocation(idx_j, ci);
  for (size_t idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state[idx] == cj) {
      data.set_allocation(S[idx], ci);
    }
  }

  double acceptance_ratio = compute_acceptance_ratio_merge(likelihood_old_clusters);

  if (std::log(dis_real(gen)) > acceptance_ratio)
    rollback();
  else
    accepted_merge++;
}

void SplitMerge_LSS_SDDS::shuffle() {
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

  if (std::log(dis_real(gen)) > log_acceptance_ratio)
    rollback();
  else
    accepted_shuffle++;
}

void SplitMerge_LSS_SDDS::step() {
  const bool similarity = (dis_real(gen) >= 0.5);

  process.set_old_cluster_members(data.get_cluster_map());
  process.set_old_allocations(data.get_allocations());
  process.set_old_K(data.get_K());

  choose_indeces(similarity);

  process.set_idx_i(idx_i);
  process.set_idx_j(idx_j);

  log_split_gibbs_prob = 0.0;
  log_merge_gibbs_prob = 0.0;

  if (similarity) {
    if (ci != cj) {
      merge_moves++;
      smart_merge_move();
    } else {
      split_moves++;
      dumb_split_move();
    }
  } else {
    if (ci == cj) {
      split_moves++;
      smart_split_move();
    } else {
      merge_moves++;
      dumb_merge_move();
    }
  }

  log_split_gibbs_prob = 0.0;
  log_merge_gibbs_prob = 0.0;

  if (shuffle_bool) {
    shuffle_moves++;
    process.set_old_allocations(data.get_allocations());
    process.set_old_cluster_members(data.get_cluster_map());
    process.set_old_K(data.get_K());
    choose_clusters_shuffle();
    process.set_idx_i(idx_i);
    process.set_idx_j(idx_j);
    shuffle();
  }
}
