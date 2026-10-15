/**
 * @file spatial_module_cache.hpp 
 * @brief Spatial helper methods for clustering processes (WIP).
 */

#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Module.hpp"
#include "../caches/spatial_cache.hpp"

/**
 * @class SpatialModuleCache (WIP)
 * @brief Module providing spatial methods for processes utilizing spatial
 * information.
 * @note This is a work-in-progress class, intended to replace SpatialModule but at the moment it does not work properly.
 *
 * This class offers utility functions to compute neighbor counts based on an
 * adjacency matrix W, facilitating the incorporation of spatial dependencies
 * in clustering algorithms. The neighbor relationships are cached at
 * construction time for efficient repeated access.
 */
class SpatialModuleCache : public Module {
protected:
    const Data &data_module;           ///< Reference to data object with cluster assignments
    const double spatial_weight = 1.0; ///< Weighting factor for spatial similarity

    SpatialCache &cache; ///< Spatial cache for additional optimizations

public:
    /**
     * @brief Constructs a SpatialModuleCache with parameter and data references.
     *
     * Initializes the neighbor cache by calling neighbor_cache_compute().
     *
     * @param data_ Reference to the Data object with cluster assignments.
     * @param spatial_cache_ Precomputed spatial cache.
     * @param spatial_coeff Weighting factor for spatial similarity.
     * @param old_alloc_provider function to access old allocations for
     * split-merge.
     * @param old_cluster_members_provider_ function to access old cluster members for
     * split-merge.
     */
    SpatialModuleCache(const Data &data_, SpatialCache &spatial_cache_, double spatial_coeff,
                       const Eigen::VectorXi *old_alloc_provider = nullptr,
                       const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : data_module(data_), cache(spatial_cache_), spatial_weight(spatial_coeff),
          Module(old_alloc_provider, old_cluster_members_provider_) {
        cache.set_allocation_ptr(&data_module.get_allocations());
    }

    /**
     * @name Spatial Methods
     * @{
     */

    /**
     * @brief Counts neighbors of an observation grouped by cluster membership.
     *
     * Returns a vector where element k contains the number of neighbors of
     * obs_idx that belong to cluster k.
     *
     * @param obs_idx The index of the observation (0 to N-1).
     * @return A K-dimensional vector of neighbor counts per cluster.
     */
    Eigen::VectorXd compute_similarity_obs(int obs_idx) const override;

    /**
     * @brief Counts internal edges within a cluster.
     *
     * Computes the number of edges where both endpoints belong to the specified
     * cluster. This is calculated as: (1/2) * sum_{i,j in cluster} W(i,j).
     * The division by 2 accounts for the symmetry of W (each edge counted twice).
     *
     * @param cls_idx The index of the cluster (0 to K-1).
     * @param old_allo If true, uses old allocations from
     * old_allocations_provider; if false, uses current allocations from
     * data_module (default: false).
     * @return The total number of internal edges in the cluster.
     */
    double compute_similarity_cls(int cls_idx, bool old_allo = false) const override;

    /**
     * @brief Counts neighbors of an observation within a specific cluster.
     *
     * Uses the cached neighbor list to efficiently count how many neighbors
     * of observation obs_idx belong to cluster cls_idx.
     *
     * @param obs_idx The index of the observation (0 to N-1).
     * @param cls_idx The index of the cluster to consider for neighbor counting.
     * @return The number of neighbors for the observation in the specified
     * cluster.
     */
    double compute_similarity_obs(int obs_idx, int cls_idx) const override;
    /** @} */
};
