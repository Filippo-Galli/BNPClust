/**
 * @file spatial_module.hpp
 * @brief Spatial helper methods for clustering processes.
 */

#pragma once

#include "../../utils/Data.hpp"
#include "../../utils/Module.hpp"

/**
 * @class SpatialModule
 * @brief Module providing spatial methods for processes utilizing spatial
 * information.
 *
 * This class offers utility functions to compute neighbor counts based on an
 * adjacency matrix W, facilitating the incorporation of spatial dependencies
 * in clustering algorithms. The neighbor relationships are cached at
 * construction time for efficient repeated access.
 */
class SpatialModule : public Module {
protected:
    /**
     * @brief Cache storing neighbor indices for each observation.
     *
     * neighbor_cache[i] contains a vector of indices j where W(i,j) == 1.
     * This avoids repeatedly scanning the full adjacency matrix during MCMC
     * sampling.
     */
    std::vector<std::vector<int>> neighbor_cache;

    /**
     * @brief Precomputes and stores neighbor indices for all observations.
     *
     * This method initializes neighbor_cache by extracting non-zero entries
     * from each row of the adjacency matrix W. Called once during construction
     * to enable O(neighbors) instead of O(N) lookup time.
     */
    void neighbor_cache_compute();

    const Data &data_module;           ///< Reference to data object with cluster assignments
    const Eigen::MatrixXi W;           ///< Reference to adjacency matrix W from covariates
    const double spatial_weight = 1.0; ///< Weighting factor for spatial similarity

public:
    /**
     * @brief Constructs a SpatialModule with parameter and data references.
     *
     * Initializes the neighbor cache by calling neighbor_cache_compute().
     *
     * @param data_ Reference to the Data object with cluster assignments.
     * @param W_ Reference to the adjacency matrix W.
     * matrix.
     * @param spatial_coeff Weighting factor for spatial similarity.
     * @param old_alloc_provider function to access old allocations for
     * split-merge.
     * @param old_cluster_members_provider_ function to access old cluster members for
     * split-merge.
     */
    SpatialModule(const Data &data_, const Eigen::MatrixXi W_, double spatial_coeff,
                  const Eigen::VectorXi *old_alloc_provider = nullptr,
                  const std::unordered_map<int, std::vector<int>> *old_cluster_members_provider_ = nullptr)
        : W(W_), data_module(data_), spatial_weight(spatial_coeff),
          Module(old_alloc_provider, old_cluster_members_provider_) {

        neighbor_cache_compute();
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
