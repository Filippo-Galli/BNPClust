/**
 * @file splitmerge_SAMS.cpp
 * @brief Implementation of Sequential Allocation Merge-Split (SAMS) sampler
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_SAMS.hpp"
#include <cmath>

void SplitMerge_SAMS::choose_indeces() {
  choose_indeces_uniform(true);
}

void SplitMerge_SAMS::split_move() {
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);
  log_split_gibbs_prob = 0.0;

  data.set_allocation(idx_j, data.get_K());
  cj = data.get_cluster_assignment(idx_j);

  sequential_allocation(1);

  double acceptance_ratio = compute_acceptance_ratio_split(likelihood_old_cluster);

  if (std::log(dis_real(gen)) > acceptance_ratio) {
    rollback();
  } else {
    accepted_split++;
  }
}

void SplitMerge_SAMS::merge_move() {
  log_merge_gibbs_prob = 0.0;
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  sequential_allocation(1, true);

  data.set_allocation(idx_j, ci);

  for (size_t idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state[idx] == cj) {
      data.set_allocation(S[idx], ci);
    }
  }

  double acceptance_ratio = compute_acceptance_ratio_merge(likelihood_old_clusters);

  if (std::log(dis_real(gen)) > acceptance_ratio) {
    rollback();
  } else {
    accepted_merge++;
  }
}