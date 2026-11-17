/**
 * @file Data.hpp
 * @brief Data structure for managing point distances and cluster allocations
 */

#pragma once

#include "../params/utils-params.hpp"
#include <Eigen/Dense>
#include <mdspan>
#include <span>
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

  std::vector<int> allocations; ///< Cluster allocation for each point (size n)
  int K;                        ///< Current number of clusters

  /// Position of each point in its cluster's member vector: member_pos[i] =
  /// index in cluster_members[allocations[i]]
  std::vector<int> member_pos;

  /// Contiguous cluster member vectors: cluster_members[k] = list of points in
  /// cluster k
  std::vector<std::vector<int>> cluster_members;

  /**
   * @brief Removes an empty cluster and compacts cluster indices
   * @param old_cluster Index of the cluster to remove
   */
  virtual void compact_cluster(int old_cluster);

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
   * @brief Constructs a Data object with a span of initial allocations
   * @param p Reference to model parameters
   * @param initial_allocations Initial cluster assignments as std::span
   */
  Data(const utils_params &p, std::span<const int> initial_allocations);

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
   * @brief Gets the number of columns in the data
   * @return Number of columns
   */
  int get_p() const { return p.p; }

  /**
   * @brief Gets the cluster allocations vector as a std::span
   * @return Non-owning view of the allocations vector
   */
  [[nodiscard]] std::span<const int> get_allocations() const noexcept {
    return allocations;
  }

  /**
   * @brief Gets the size of a specific cluster
   * @param cluster_index Index of the cluster
   * @return Number of points in the cluster (0 if cluster doesn't exist)
   */
  int get_cluster_size(unsigned cluster_index) const {
    if (cluster_index >= static_cast<unsigned>(K)) {
      return 0;
    }
    return static_cast<int>(cluster_members[cluster_index].size());
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
    return allocations[static_cast<size_t>(index)];
  }

  /**
   * @brief Gets all point indices assigned to a specific cluster as a zero-overhead std::span
   * @param cluster Index of the cluster
   * @return Non-owning view of point indices in the cluster
   */
  [[nodiscard]] std::span<const int> get_cluster_assignments(int cluster) const noexcept;

  /**
   * @brief Direct access to contiguous cluster member lists
   * @return Const reference to vector of cluster member vectors
   */
  const std::vector<std::vector<int>> &get_cluster_members() const noexcept {
    return cluster_members;
  }

  /**
   * @brief Direct access to member position lookup array
   * @return Const reference to position lookup vector
   */
  const std::vector<int> &get_member_pos() const noexcept { return member_pos; }

  /**
   * @brief Gets the full mapping of clusters to their member points
   * @return Unordered map of cluster indices to vectors of point indices
   */
  std::unordered_map<int, std::vector<int>> get_cluster_map_copy() const;

  /**
   * @brief Gets the full mapping of clusters to their member points
   * @return Unordered map of cluster indices to vectors of point indices
   */
  std::unordered_map<int, std::vector<int>> get_cluster_map() const {
    return get_cluster_map_copy();
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
  virtual void set_allocations(std::span<const int> new_allocations);

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
  restore_state(const std::vector<int> &old_allocations,
                const std::unordered_map<int, std::vector<int>> &old_cluster_members,
                int old_K);

  virtual void
  restore_state(Eigen::VectorXi &old_allocations,
                std::unordered_map<int, std::vector<int>> &old_cluster_members,
                int old_K);

  /** @} */

  virtual ~Data() = default;
};
