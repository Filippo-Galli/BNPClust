/**
 * @file Data.cpp
 * @brief Implementation of the Data class
 */

#include "Data.hpp"
#include <algorithm>
#include <iostream>

Data::Data(const utils_params &p, const Eigen::VectorXi &initial_allocations)
    : p(p) {
  if (initial_allocations.size() > 0) {
    allocations.assign(initial_allocations.data(),
                       initial_allocations.data() + initial_allocations.size());
  } else {
    std::cout << "[INFO] No initial allocations provided, starting with all "
                 "points in one cluster."
              << std::endl;
    allocations.assign(static_cast<size_t>(p.n), 0);
  }

  int max_c = -1;
  for (int c : allocations) {
    if (c > max_c) {
      max_c = c;
    }
  }
  K = (max_c >= 0) ? (max_c + 1) : 0;

  cluster_members.assign(static_cast<size_t>(K), std::vector<int>());
  member_pos.assign(static_cast<size_t>(p.n), -1);

  for (size_t i = 0; i < allocations.size(); ++i) {
    int cluster_id = allocations[i];
    if (cluster_id >= 0 && cluster_id < K) {
      member_pos[i] = static_cast<int>(
          cluster_members[static_cast<size_t>(cluster_id)].size());
      cluster_members[static_cast<size_t>(cluster_id)].push_back(
          static_cast<int>(i));
    }
  }
}

Data::Data(const utils_params &p, std::span<const int> initial_allocations)
    : p(p) {
  if (!initial_allocations.empty()) {
    allocations.assign(initial_allocations.begin(), initial_allocations.end());
  } else {
    allocations.assign(static_cast<size_t>(p.n), 0);
  }

  int max_c = -1;
  for (int c : allocations) {
    if (c > max_c) {
      max_c = c;
    }
  }
  K = (max_c >= 0) ? (max_c + 1) : 0;

  cluster_members.assign(static_cast<size_t>(K), std::vector<int>());
  member_pos.assign(static_cast<size_t>(p.n), -1);

  for (size_t i = 0; i < allocations.size(); ++i) {
    int cluster_id = allocations[i];
    if (cluster_id >= 0 && cluster_id < K) {
      member_pos[i] = static_cast<int>(
          cluster_members[static_cast<size_t>(cluster_id)].size());
      cluster_members[static_cast<size_t>(cluster_id)].push_back(
          static_cast<int>(i));
    }
  }
}

double Data::get_data(int i, int j) const {
#if VERBOSITY_LEVEL >= 1
  if (i < 0 || i >= p.n || j < 0 || j >= p.p) {
    throw std::out_of_range("Index out of bounds in get_data");
  }
#endif

  return p.get_D_view()[i, j];
}

std::span<const int>
Data::get_cluster_assignments(int cluster) const noexcept {
  if (cluster >= K || cluster < 0) {
    return {};
  }
  const auto &members = cluster_members[static_cast<size_t>(cluster)];
  return std::span<const int>(members.data(), members.size());
}

std::unordered_map<int, std::vector<int>> Data::get_cluster_map_copy() const {
  std::unordered_map<int, std::vector<int>> map;
  map.reserve(static_cast<size_t>(K));
  for (int k = 0; k < K; ++k) {
    if (!cluster_members[static_cast<size_t>(k)].empty()) {
      map[k] = cluster_members[static_cast<size_t>(k)];
    }
  }
  return map;
}

void Data::compact_cluster(int old_cluster) {
#if VERBOSITY_LEVEL >= 1
  if (old_cluster < 0 || old_cluster >= K) {
    throw std::out_of_range(
        "old_cluster index out of bounds in compact_cluster");
  }
#endif

  const int last_cluster = K - 1;

  // If compacting the last cluster (or the only cluster), this is a pure
  // deletion.
  if (K <= 1 || old_cluster == last_cluster) {
    cluster_members.pop_back();
    K--;
    return;
  }

  // Shift allocations and member positions of the last cluster to the old
  // cluster
  auto &last_members = cluster_members[static_cast<size_t>(last_cluster)];
  if (!last_members.empty()) {
    cluster_members[static_cast<size_t>(old_cluster)] = std::move(last_members);
    auto &moved_members = cluster_members[static_cast<size_t>(old_cluster)];
    for (size_t pos = 0; pos < moved_members.size(); ++pos) {
      int point_index = moved_members[pos];
      allocations[static_cast<size_t>(point_index)] = old_cluster;
      member_pos[static_cast<size_t>(point_index)] = static_cast<int>(pos);
    }
  } else {
    cluster_members[static_cast<size_t>(old_cluster)].clear();
  }

  // Remove the last cluster
  cluster_members.pop_back();
  K--; // Decrease the number of clusters
}

void Data::set_allocation_wo_compaction(int index, int cluster) {
  // Bounds checking for index
#if VERBOSITY_LEVEL >= 1
  if (index < 0 || index >= p.n) {
    throw std::out_of_range("Index out of bounds in set_allocation");
  }
#endif

  // Validate cluster parameter early
#if VERBOSITY_LEVEL >= 1
  if (cluster < -1 || cluster > K) {
    throw std::out_of_range("Invalid cluster index in set_allocation");
  }
#endif

  int old_cluster = allocations[static_cast<size_t>(index)];

  // Early exit if no change needed
  if (old_cluster == cluster) {
    return;
  }

  // Remove from old cluster in O(1) using swap-and-pop
  if (old_cluster != -1 && old_cluster < K) {
    int pos = member_pos[static_cast<size_t>(index)];
    auto &members = cluster_members[static_cast<size_t>(old_cluster)];

#if VERBOSITY_LEVEL >= 1
    if (pos < 0 || pos >= static_cast<int>(members.size()) ||
        members[static_cast<size_t>(pos)] != index) {
      throw std::runtime_error(
          "Inconsistent state: member_pos does not match cluster_members");
    }
#endif

    int last_point = members.back();
    members[static_cast<size_t>(pos)] = last_point;
    member_pos[static_cast<size_t>(last_point)] = pos;
    members.pop_back();
    member_pos[static_cast<size_t>(index)] = -1;
  }

  // Handle new cluster assignment in O(1)
  if (cluster == -1) {
    allocations[static_cast<size_t>(index)] = -1;
    member_pos[static_cast<size_t>(index)] = -1;
  } else if (cluster == K) {
    // New cluster creation
    allocations[static_cast<size_t>(index)] = K;
    cluster_members.emplace_back();
    member_pos[static_cast<size_t>(index)] = 0;
    cluster_members[static_cast<size_t>(K)].push_back(index);
    K++;
  } else {
    // Existing cluster assignment
    allocations[static_cast<size_t>(index)] = cluster;
    member_pos[static_cast<size_t>(index)] =
        static_cast<int>(cluster_members[static_cast<size_t>(cluster)].size());
    cluster_members[static_cast<size_t>(cluster)].push_back(index);
  }
}

void Data::set_allocation(int index, int cluster) {
  int old_cluster = allocations[static_cast<size_t>(index)];
  set_allocation_wo_compaction(index, cluster);

  // Check if old cluster became empty and needs compaction
  if (old_cluster != -1 && old_cluster < K &&
      cluster_members[static_cast<size_t>(old_cluster)].empty()) {
    compact_cluster(old_cluster);
  }
}

void Data::set_allocations(std::span<const int> new_allocations) {
#if VERBOSITY_LEVEL >= 1
  if (static_cast<int>(new_allocations.size()) != p.n) {
    throw std::invalid_argument(
        "New allocations vector must be of size n, set_allocations failed");
  }
#endif

  allocations.assign(new_allocations.begin(), new_allocations.end());
  int max_c = -1;
  for (int c : allocations) {
    if (c > max_c) {
      max_c = c;
    }
  }
  K = (max_c >= 0) ? (max_c + 1) : 0;

  cluster_members.assign(static_cast<size_t>(K), std::vector<int>());
  member_pos.assign(static_cast<size_t>(p.n), -1);

  for (size_t i = 0; i < allocations.size(); ++i) {
    int cid = allocations[i];
    if (cid >= 0 && cid < K) {
      member_pos[i] =
          static_cast<int>(cluster_members[static_cast<size_t>(cid)].size());
      cluster_members[static_cast<size_t>(cid)].push_back(static_cast<int>(i));
    }
  }
}

void Data::set_allocations(const Eigen::VectorXi &new_allocations) {
  set_allocations(std::span<const int>(
      new_allocations.data(), static_cast<size_t>(new_allocations.size())));
}

void Data::restore_state(
    const std::vector<int> &old_allocations,
    const std::unordered_map<int, std::vector<int>> &old_cluster_members, int old_K) {
#if VERBOSITY_LEVEL >= 1
  if (old_allocations.size() != allocations.size()) {
    throw std::invalid_argument(
        "Saved allocations size mismatch in restore_state");
  }
#endif

  allocations = old_allocations;
  K = old_K;

  cluster_members.assign(static_cast<size_t>(K), std::vector<int>());
  member_pos.assign(static_cast<size_t>(p.n), -1);

  for (const auto &[k, members] : old_cluster_members) {
    if (k >= 0 && k < K) {
      cluster_members[static_cast<size_t>(k)] = members;
      for (size_t pos = 0; pos < members.size(); ++pos) {
        int pt = members[pos];
        if (pt >= 0 && pt < p.n) {
          member_pos[static_cast<size_t>(pt)] = static_cast<int>(pos);
        }
      }
    }
  }
}

void Data::restore_state(
    Eigen::VectorXi &old_allocations,
    std::unordered_map<int, std::vector<int>> &old_cluster_members, int old_K) {
  std::vector<int> std_old_alloc(old_allocations.data(), old_allocations.data() + old_allocations.size());
  restore_state(std_old_alloc, old_cluster_members, old_K);
}
