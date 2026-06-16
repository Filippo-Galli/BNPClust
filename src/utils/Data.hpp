/**
 * @file Data.hpp
 * @brief Data structure for managing point distances and cluster allocations
 */

#pragma once

#include "../params/utils-params.hpp"
#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

// compilation time verbosity level
#ifndef VERBOSITY_LEVEL
#define VERBOSITY_LEVEL 0
#endif

/**
 * @class Data
 * @brief Manages distance matrices and cluster allocations for points
 *
 * This class stores a distance matrix between points and tracks their
 * cluster assignments. It provides methods to query distances, manage
 * cluster allocations, and maintain cluster membership information.
 */
class Data {
protected:
  const utils_params &p; ///< Reference to model parameters

  Eigen::VectorXi allocations; ///< Cluster allocation for each point
  int K;                       ///< Current number of clusters

  /// Maps cluster indices to vectors of point indices in that cluster
  std::unordered_map<int, std::vector<int>> cluster_members;

  /**
   * @brief Removes an empty cluster and compacts cluster indices
   * @param old_cluster Index of the cluster to remove
   */
  void compact_cluster(int old_cluster);

  /**
   * @brief Assigns a point to a cluster without compaction
   * @param index Index of the point to reassign
   * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
   * @throws std::out_of_range if index or cluster is invalid
   */
  void set_allocation_wo_compaction(int index, int cluster);

public:
  /**
   * @brief Constructs a Data object with a distance matrix
   * @param p Reference to model parameters
   * @param initial_allocations Optional initial cluster assignments (default:
   * all in one cluster)
   * @throws std::invalid_argument if distance matrix is not square
   */
  Data(const utils_params &p,
       const Eigen::VectorXi &initial_allocations = Eigen::VectorXi());

  /**
   * @brief Getters
   * @{
   */

  /**
   * @brief Gets the data point at index i column j
   * @param i Index of the point
   * @param j Index of the column
   * @return Value of the data point at (i, j)
   */
  double get_data(int i, int j = 0) const;

  /**
   * @brief Gets the total number of points
   * @return Number of points
   */
  int get_n() const { return p.n; }

  /**
   * @brief Gets the current number of clusters
   * @return Number of clusters
   */
  int get_K() const { return K; }

  /**
   * @brief Gets the cluster allocations vector
   * @return Reference to the allocations vector
   */
  const Eigen::VectorXi &get_allocations() const { return allocations; }

  /**
   * @brief Gets the size of a specific cluster
   * @param cluster_index Index of the cluster
   * @return Number of points in the cluster (0 if cluster doesn't exist)
   */
  int get_cluster_size(unsigned cluster_index) const {
    return (cluster_index < K) ? cluster_members.at(cluster_index).size() : 0;
  }

  /**
   * @brief Gets the cluster assignment of a specific point
   * @param index Index of the point
   * @return Cluster index the point is assigned to
   * @throws std::out_of_range if index is out of bounds
   */
  int get_cluster_assignment(int index) const {
    if (index < 0 || index >= p.n) {
      throw std::out_of_range("Index out of bounds in get_cluster_assignment");
    }
    return allocations(index);
  }

  /**
   * @brief Gets all point indices assigned to a specific cluster
   * @param cluster Index of the cluster
   * @return Vector of point indices in the cluster
   * @throws std::out_of_range if cluster index is invalid
   */
  Eigen::VectorXi get_cluster_assignments(int cluster) const;

  /**
   * @brief Gets all point indices assigned to a specific cluster (map form)
   * @param cluster Index of the cluster
   * @return Map to vector of point indices in the cluster
   * @throws std::out_of_range if cluster index is invalid
   */
  Eigen::Map<const Eigen::VectorXi>
  get_cluster_assignments_ref(int cluster) const;

  /**
   * @brief Gets the full mapping of clusters to their member points
   * @return Unordered map of cluster indices to vectors of point indices
   */
  std::unordered_map<int, std::vector<int>> get_cluster_map_copy() const {
    return cluster_members;
  }

  /** @} */

  /**
   * @brief Setters
   * @{
   */

  /**
   * @brief Assigns a point to a cluster
   * @param index Index of the point to reassign
   * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
   * @throws std::out_of_range if index or cluster is invalid
   */
  virtual void set_allocation(int index, int cluster);

  /**
   * @brief Sets all cluster allocations at once
   * @param new_allocations Vector of cluster assignments for all points
   * @throws std::invalid_argument if vector size doesn't match number of points
   */
  virtual void set_allocations(const Eigen::VectorXi &new_allocations);

  /**
   * @brief Restores allocations, cluster memberships, and cluster count from a
   * saved state
   *
   * @param old_allocations Vector holding the saved allocations (swapped into
   * place)
   * @param old_cluster_members Map of saved cluster memberships (swapped into
   * place)
   * @param old_K Number of clusters in the saved state
   */
  virtual void
  restore_state(Eigen::VectorXi &old_allocations,
                std::unordered_map<int, std::vector<int>> &old_cluster_members,
                int old_K);

  /**
   * @brief Gets the full mapping of clusters to their member points
   * @return Reference to the vector of cluster members
   */
  const std::unordered_map<int, std::vector<int>> &get_cluster_map() const {
    return cluster_members;
  }

  /** @} */

  virtual ~Data() = default;
};
