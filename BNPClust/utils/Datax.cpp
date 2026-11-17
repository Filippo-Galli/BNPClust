/**
 * @file Datax.cpp
 * @brief Implementation of Datax class
 */

#include "Datax.hpp"

void Datax::compact_cluster(int old_cluster) {
#if VERBOSITY_LEVEL >= 1
  if (old_cluster < 0 || old_cluster >= K) {
    throw std::out_of_range("old_cluster index out of bounds in compact_cluster");
  }
#endif

  const int last_cluster = K - 1;

  // If compacting the last cluster (or the only cluster), this is a pure deletion.
  if (K <= 1 || old_cluster == last_cluster) {
    cluster_members.pop_back();
    // Update cluster info to remove the old cluster
    for (auto &&ci : cluster_info) {
      ci->remove_info(old_cluster);
    }
    K--;
    return;
  }

  // Shift allocations of the last cluster to the old cluster
  auto &last_members = cluster_members[static_cast<size_t>(last_cluster)];
  if (!last_members.empty()) {
    cluster_members[static_cast<size_t>(old_cluster)] = std::move(last_members);
    auto &moved_members = cluster_members[static_cast<size_t>(old_cluster)];
    for (size_t pos = 0; pos < moved_members.size(); ++pos) {
      int point_index = moved_members[pos];
      allocations[static_cast<size_t>(point_index)] = old_cluster;
      member_pos[static_cast<size_t>(point_index)] = static_cast<int>(pos);
    }

    // Update cluster info accordingly
    for (auto &&ci : cluster_info) {
      ci->move_cluster_info(last_cluster, old_cluster);
    }
  } else {
    // Last cluster is empty, just clear the old cluster
    cluster_members[static_cast<size_t>(old_cluster)].clear();
  }

  // Remove the last cluster
  cluster_members.pop_back();
  // Update cluster info to remove the last cluster
  for (auto &&ci : cluster_info) {
    ci->remove_info(last_cluster);
  }

  K--; // Decrease the number of clusters
}

void Datax::set_allocation(int index, int cluster) {
  int old_cluster = allocations[static_cast<size_t>(index)];

  if (old_cluster == cluster) {
    return;
  }

  Data::set_allocation_wo_compaction(index, cluster);
  for (auto &&ci : cluster_info) {
    ci->set_allocation(index, cluster, old_cluster);
  }

  // Check if old cluster became empty and needs compaction
  if (old_cluster != -1 && old_cluster < K &&
      cluster_members[static_cast<size_t>(old_cluster)].empty()) {
    Datax::compact_cluster(old_cluster);
  }
}

void Datax::set_allocations(const Eigen::VectorXi &new_allocations) {
  set_allocations(std::span<const int>(
      new_allocations.data(), static_cast<size_t>(new_allocations.size())));
}

void Datax::set_allocations(std::span<const int> new_allocations) {
  Data::set_allocations(new_allocations);
  for (auto &&ci : cluster_info) {
    ci->recompute(K, get_allocations());
  }
}

void Datax::restore_state(
    const std::vector<int> &old_allocations,
    const std::unordered_map<int, std::vector<int>> &old_cluster_members, int old_K) {
  Data::restore_state(old_allocations, old_cluster_members, old_K);
  for (auto &&ci : cluster_info) {
    ci->recompute(K, get_allocations());
  }
}

void Datax::restore_state(
    Eigen::VectorXi &old_allocations,
    std::unordered_map<int, std::vector<int>> &old_cluster_members, int old_K) {
  std::vector<int> std_old_alloc(old_allocations.data(), old_allocations.data() + old_allocations.size());
  restore_state(std_old_alloc, old_cluster_members, old_K);
}