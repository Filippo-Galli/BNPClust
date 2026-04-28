/**
 * @file spatial_module_continuous.cpp
 * @brief Implementation of `SpatialModuleContinuous`.
 */

#include "spatial_module_continuous.hpp"

void SpatialModuleContinuous::neighbor_cache_compute() {
  const int N = data_module.get_n();
  neighbor_cache.resize(N);

  // For each observation, store indices of its neighbors (where W(i,j) == 1)
  for (int obs_idx = 0; obs_idx < N; ++obs_idx) {
    Eigen::RowVectorXd row = W.row(obs_idx);

    // Reserve space for the number of non-zero elements in the row
    neighbor_cache[obs_idx].reserve(row.nonZeros());

    for (int j = 0; j < row.size(); ++j) {
      if (row(j) > 0) {
        neighbor_cache[obs_idx].push_back(j);
      }
    }
  }
}

double SpatialModuleContinuous::compute_similarity_obs(int obs_idx,
                                                       int cls_idx) const {
  double similarity = 0;
  const std::vector<int> &neighbors = neighbor_cache[obs_idx];

  // Iterate through cached neighbors and count those in the specified cluster
  for (size_t i = 0; i < neighbors.size(); ++i) {
    int cluster_i = data_module.get_cluster_assignment(neighbors[i]);
    if (cluster_i != -1 && cluster_i != cls_idx) {
      similarity += W(obs_idx, neighbors[i]);
    }
  }

  return -spatial_weight * similarity;
}

double SpatialModuleContinuous::compute_similarity_cls(int cls_idx,
                                                       bool old_allo) const {

  const Eigen::VectorXi &cls_idx_allocations =
      (old_allo && old_cluster_members_provider)
          ? Eigen::Map<const Eigen::VectorXi>(
                old_cluster_members_provider->at(cls_idx).data(),
                old_cluster_members_provider->at(cls_idx).size())
          : data_module.get_cluster_assignments(cls_idx);

  double total_similarity = 0;
  for (auto &&i : cls_idx_allocations) {
    // Use cached neighbor indices instead of iterating over full adjacency
    // matrix
    const std::vector<int> &neighbors = neighbor_cache[i];

    // Count neighbors in each cluster
    for (size_t j = 0; j < neighbors.size(); ++j) {
      int neighbor_idx = neighbors[j];
      int cluster_i = data_module.get_cluster_assignment(neighbor_idx);

      if (cluster_i != -1 && cluster_i != cls_idx) {
        total_similarity += W(i, neighbor_idx);
      }
    }
  }

  return -spatial_weight * total_similarity / 2; // Each edge counted twice
}

Eigen::VectorXd
SpatialModuleContinuous::compute_similarity_obs(int obs_idx) const {
  Eigen::VectorXd cluster_adjacency =
      Eigen::VectorXd::Zero(data_module.get_K());

  // Use cached neighbor indices instead of iterating over full adjacency matrix
  const std::vector<int> &neighbors = neighbor_cache[obs_idx];

  // Count neighbors in each cluster
  for (size_t i = 0; i < neighbors.size(); ++i) {
    int neighbor_idx = neighbors[i];
    int cluster_i = data_module.get_cluster_assignment(neighbor_idx);

    if (cluster_i != -1) {
      cluster_adjacency(cluster_i) += W(obs_idx, neighbor_idx);
    }
  }

  // add spatial weight
  cluster_adjacency *=
      spatial_weight * Eigen::VectorXd::Ones(data_module.get_K());

  return -cluster_adjacency;
}
