/**
 * @file ClusterInfo.hpp
 * @brief Abstract base class for managing cluster information and caches
 *
 * This file defines the ClusterInfo abstract interface for managing cluster-related
 * information and caches that work in conjunction with the Data class.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include <Eigen/Dense>
#include <unordered_map>
#include <vector>

/**
 * @class ClusterInfo
 * @brief Abstract interface for cluster information management
 *
 * This class provides an abstract interface for managing cluster-related information,
 * including caching mechanisms that integrate with the Data class. Derived classes
 * implement specific caching strategies for efficient cluster computations.
 */
class ClusterInfo {

public:
    ClusterInfo() = default;

    /**
     * @brief Assigns a point to a cluster
     * @param index Index of the point to reassign
     * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
     * @param old_cluster Previous cluster index of the point
     * @throws std::out_of_range if index or cluster is invalid
     */
    virtual void set_allocation(int index, int cluster, int old_cluster) = 0;

    /**
     * @brief Recomputes all cluster information from current allocations
     * @param K Current number of clusters
     * @param allocations Vector of current cluster assignments for all points
     */
    virtual void recompute(const int K, const Eigen::VectorXi &allocations) = 0;

    /**
     * @brief Moves cluster information from one cluster to another
     * @param from_cluster Index of the source cluster
     * @param to_cluster Index of the target cluster
     */
    virtual void move_cluster_info(int from_cluster, int to_cluster) = 0;

    /**
     * @brief Removes information related to a specific cluster
     * @param cluster Index of the cluster to remove
     */
    virtual void remove_info(int cluster) = 0;

    virtual ~ClusterInfo() = default;
};
