/**
 * @file splitmerge_LSS.cpp
 * @brief Implementation of Locality Sensitive Sampling (LSS) Split-Merge sampler
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_LSS.hpp"

void SplitMerge_LSS::choose_indeces(bool similarity) {
  idx_i = dis_n(gen);

  const auto distances = params.D.row(idx_i);
  const int n = data.get_n();
  double distance_sum = 0.0;

  if (similarity) {
    for (int idx = 0; idx < n; ++idx) {
      probs[idx] = (idx == idx_i) ? 0.0 : distances(idx);
      distance_sum += probs[idx];
    }
  } else {
    for (int idx = 0; idx < n; ++idx) {
      probs[idx] = (idx == idx_i) ? 0.0 : (1.0 / distances(idx));
      distance_sum += probs[idx];
    }
  }

  double u = dis_real(gen) * distance_sum;
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

void SplitMerge_LSS::choose_indices_for_step() {
  choose_indeces(dis_real(gen) < 0.5);
}