/**
 * @file Data.cpp
 * @brief Implementation of the Data class
 */

#include "Data.hpp"
#include "Eigen/src/Core/Matrix.h"
#include <iostream>

Data::Data(const utils_params &p, const Eigen::VectorXi &initial_allocations)
    : p(p), allocations(initial_allocations) {

  // Initialize K (number of clusters) based on initial allocations
  if (allocations.size() == 0) {
    std::cout << "[INFO] No initial allocations provided, starting with all "
                 "points in one cluster."
              << std::endl;
    // If no initial allocations provided, start with all points in one cluster
    allocations = Eigen::VectorXi::Zero(p.n);
  }

  // Find the maximum cluster index to determine K
  K = allocations.maxCoeff() + 1;

  // Allocate cluster mapping
  for (int i = 0; i < p.n; ++i) {
    int cluster_id = allocations(i);
    if (cluster_id >= 0) {
      cluster_members[cluster_id].push_back(i);
    }
  }
}

double Data::get_distance(int i, int j) const {
#if VERBOSITY_LEVEL >= 1
  if (i < 0 || i >= p.n || j < 0 || j >= p.p) {
    throw std::out_of_range("Index out of bounds in get_distance");
  }
#endif

  return p.D(i, j);
}

Eigen::VectorXi Data::get_cluster_assignments(int cluster) const {
#if VERBOSITY_LEVEL >= 1
  if (cluster > K || cluster < 0) {
    throw std::out_of_range("Index out of bounds in get_cluster_assignments");
  }
#endif

  // Check if the cluster exists in the map before accessing it
  if (cluster >= K || cluster < 0) {
    return Eigen::VectorXi::Zero(
        0); // Return empty vector for non-existent clusters
  }

  const auto &it = cluster_members.at(cluster);
  return Eigen::Map<const Eigen::VectorXi>(
      it.data(), static_cast<Eigen::Index>(it.size()));
}

Eigen::Map<const Eigen::VectorXi>
Data::get_cluster_assignments_ref(int cluster) const {
#if VERBOSITY_LEVEL >= 1
  if (cluster > K || cluster < 0) {
    throw std::out_of_range(
        "Index out of bounds in get_cluster_assignments_ref");
  }
#endif

  // Check if the cluster exists in the map before accessing it
  if (cluster >= K || cluster < 0) {
    static const int dummy_value = 0;
    return Eigen::Map<const Eigen::VectorXi>(&dummy_value, 0);
  }
  const auto &it = cluster_members.at(cluster);

  return Eigen::Map<const Eigen::VectorXi>(
      it.data(), static_cast<Eigen::Index>(it.size()));
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
    cluster_members.erase(old_cluster);
    K--;
    return;
  }

  auto last_cluster_it = cluster_members.find(last_cluster);
  bool last_cluster_exists = (last_cluster_it != cluster_members.end());

  // Shift allocations of the last cluster to the old cluster
  if (last_cluster_exists && !last_cluster_it->second.empty()) {
    for (int point_index : last_cluster_it->second) {
      allocations(point_index) = old_cluster;
    }

    // Move members of the last cluster to the old cluster using move semantics
    cluster_members[old_cluster] = std::move(last_cluster_it->second);
  } else {
    // Last cluster is empty or doesn't exist, just clear the old cluster
    cluster_members[old_cluster].clear();
  }

  // Remove the last cluster
  cluster_members.erase(last_cluster);
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

  int old_cluster = allocations(index);

  // Early exit if no change needed
  if (old_cluster == cluster) {
    return;
  }

  auto old_cluster_it = cluster_members.end();
  bool old_cluster_exists = false;

  // Check if old_cluster exists in the map
  if (old_cluster != -1) {
    old_cluster_it = cluster_members.find(old_cluster);
    old_cluster_exists = (old_cluster_it != cluster_members.end());

#if VERBOSITY_LEVEL >= 1
    if (!old_cluster_exists) {
      throw std::runtime_error(
          "Inconsistent state: old_cluster not found in cluster_members");
    }
#endif
  }

  // Remove from old cluster first (if applicable)
  if (old_cluster_exists) {
    auto &old_members = old_cluster_it->second;
    // Use find + erase instead of remove + erase for single element (more
    // efficient)
    auto it = std::find(old_members.begin(), old_members.end(), index);
    if (it != old_members.end()) {
      old_members.erase(it);
    }
  }

  // Update allocation
  allocations(index) = cluster;

  // Handle new cluster assignment
  if (cluster == -1) {
    // Point becomes unallocated - nothing more to do after removal above
  } else if (cluster == K) {
    // New cluster creation
    cluster_members[K].push_back(index);
    K++;
  } else {
    // Existing cluster assignment - cache the target cluster iterator
    auto target_cluster_it = cluster_members.find(cluster);
    if (target_cluster_it != cluster_members.end()) {
      target_cluster_it->second.push_back(index);
    } else {
      // Create the cluster if it doesn't exist
      cluster_members[cluster].push_back(index);
    }
  }
}

void Data::set_allocation(int index, int cluster) {

  int old_cluster = allocations(index);
  auto old_cluster_it = cluster_members.find(old_cluster);

  set_allocation_wo_compaction(index, cluster);

  // Check if old cluster became empty and needs compaction
  if (old_cluster_it != cluster_members.end() &&
      old_cluster_it->second.empty()) {
    compact_cluster(old_cluster);
  }
}

void Data::set_allocations(const Eigen::VectorXi &new_allocations) {
#if VERBOSITY_LEVEL >= 1
  if (new_allocations.size() != p.n) {
    throw std::invalid_argument(
        "New allocations vector must be of size n, set_allocations failed");
  }
#endif

  allocations = new_allocations;

  // Update cluster members mapping
  cluster_members.clear();
  for (int i = 0; i < p.n; ++i)
    cluster_members[allocations(i)].push_back(i);

  // Update K based on the new allocations
  K = allocations.maxCoeff() + 1;
}

void Data::restore_state(
    Eigen::VectorXi &old_allocations,
    std::unordered_map<int, std::vector<int>> &old_cluster_members, int old_K) {
#if VERBOSITY_LEVEL >= 1
  if (old_allocations.size() != allocations.size()) {
    throw std::invalid_argument(
        "Saved allocations size mismatch in restore_state");
  }
#endif

  allocations.swap(old_allocations);
  cluster_members.swap(old_cluster_members);
  K = old_K;
}
