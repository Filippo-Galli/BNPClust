#pragma once

/**
 * @file binary_cache.hpp
 * @brief Cache for spatial model with binary covariates.
 */

#include "../../utils/ClusterInfo.hpp"
#include <vector>

/**
 * @class BinaryCache
 * @brief Cache for spatial model with binary covariates.
 * @details This class maintains cluster statistics for a spatial model that includes binary covariates.
 */

class BinaryCache : public ClusterInfo {
public:
    // Important to declare the struct useful.
    // HINT: If the name changed modified the file which uses it accordingly.
    /**
    * @struct ClusterStats
    * @brief Structure to hold statistics for each cluster.
    * @details Contains the sum of binary covariates and the count of observations in the cluster.
    */
    struct ClusterStats {
        int binary_sum = 0;
        int n = 0;
    };

private:
    std::vector<ClusterStats> cluster_stats;

public:
    const Eigen::VectorXi binary_covariates;

    BinaryCache(const Eigen::VectorXi &allocations_ref, const Eigen::VectorXi &binary_covariates)
        : binary_covariates(binary_covariates) {

        const int K = allocations_ref.maxCoeff() + 1;
        recompute(K > 0 ? K : 0, allocations_ref);
    }

    /**
     * @brief Assigns a point to a cluster
     * @param index Index of the point to reassign
     * @param cluster Target cluster index (K for new cluster, -1 for unallocated)
     * @param old_cluster Previous cluster index of the point
     * @throws std::out_of_range if index or cluster is invalid
     */
    void set_allocation(int index, int cluster, int old_cluster) override;

    /**
     * @brief Get cluster statistics for a specific cluster
     * @param cluster Index of the cluster
     * @return ClusterStats struct
     */
    inline ClusterStats get_cluster_stats(int cluster) const { return cluster_stats[cluster]; }

    /**
     * @brief Get cluster statistics reference for a specific cluster
     * @param cluster Index of the cluster
     * @return Const reference to ClusterStats struct
     */
    inline const ClusterStats &get_cluster_stats_ref(int cluster) const { return cluster_stats[cluster]; }

    /**
     * @brief Recomputes all cluster information from current allocations
     * @param K Current number of clusters
     * @param allocations_in Current allocations vector
     */
    void recompute(const int K, const Eigen::VectorXi &allocations_in) override;

    /**
     * @brief Moves cluster information from one cluster to another
     * @param from_cluster Index of the source cluster
     * @param to_cluster Index of the target cluster
     */
    inline void move_cluster_info(int from_cluster, int to_cluster) override {
        cluster_stats[to_cluster] = std::move(cluster_stats[from_cluster]);
    };

    /**
     * @brief Removes information related to a specific cluster
     * @param cluster Index of the cluster to remove
     */
    void remove_info(int cluster) override {
        if (cluster == static_cast<int>(cluster_stats.size()) - 1) {
            cluster_stats.pop_back();
        } else {
            cluster_stats.erase(cluster_stats.begin() + cluster);
        }
    }
};