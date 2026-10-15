/**
 * @file Datax.hpp
 * @brief Extended data structure integrating cluster information or caching
 * mechanisms
 *
 * This file defines the Datax class that extends the base Data class
 * with cluster information management through a vector of ClusterInfo objects.
 * It provides a unified interface for managing cluster allocations with
 * optional caching mechanisms.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "ClusterInfo.hpp"
#include "Data.hpp"
#include <memory>

/**
 * @class Datax
 * @brief Data container with integrated cluster information management
 *
 * This class extends the base Data class to include a vector of ClusterInfo
 * objects that can provide additional cluster-related computations or caching.
 * It maintains synchronization between data allocations and cluster
 * information.
 */
class Datax : public Data {

protected:
  std::vector<std::shared_ptr<ClusterInfo>> cluster_info;

  void compact_cluster(int old_cluster);

public:
  Datax(const utils_params &p, std::vector<std::shared_ptr<ClusterInfo>> ci,
        const Eigen::VectorXi &initial_allocations = Eigen::VectorXi())
      : Data(p, initial_allocations), cluster_info(std::move(ci)) {}

  /**
   * @brief Assigns a point to a cluster
   * @param index Index of the point to reassign
   * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
   * @throws std::out_of_range if index or cluster is invalid
   */
  void set_allocation(int index, int cluster) override;

  /**
   * @brief Sets all cluster allocations at once
   * @param new_allocations Vector of cluster assignments for all points
   * @throws std::invalid_argument if vector size doesn't match number of points
   */
  void set_allocations(const Eigen::VectorXi &new_allocations) override;

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
  void
  restore_state(Eigen::VectorXi &old_allocations,
                std::unordered_map<int, std::vector<int>> &old_cluster_members,
                int old_K) override;

  virtual ~Datax() = default;
};
