/**
 * @file spatial_cache.cpp
 * @brief Implementation of `SpatialCache`.
 */

#include "spatial_cache.hpp"

void SpatialCache::neighbor_cache_compute() {
    const int N = W.rows();
    neighbor_cache.resize(N);

    // For each observation, store indices of its neighbors (where W(i,j) == 1)
    for (int obs_idx = 0; obs_idx < N; ++obs_idx) {
        Eigen::RowVectorXi row = W.row(obs_idx);
        neighbor_cache[obs_idx].reserve(row.sum());

        for (int j = 0; j < row.size(); ++j) {
            if (row(j) == 1) {
                neighbor_cache[obs_idx].push_back(j);
            }
        }
    }
}

void SpatialCache::recompute(int K, const Eigen::VectorXi &allocations) {
    cluster_stats.clear();
    cluster_stats.resize(K);

    const int N = allocations.size();
    for (int i = 0; i < N; ++i) {
        int cls_idx = allocations(i);
        if (cls_idx != -1) {
            // Sum spatial covariates from neighbors
            const std::vector<int> &row = neighbor_cache[i];
            for (size_t j = 0; j < row.size(); ++j) {
                cluster_stats[cls_idx].spatial_sum +=
                    (allocations(row[j]) == cls_idx
                         ? 1
                         : 0); // Assuming each neighbor contributes '1' if in the same cluster
            }
        }
    }
}

void SpatialCache::set_allocation(int index, int cluster, int old_cluster) {

    // Remove from old cluster
    if (old_cluster != -1) {
        ClusterStats &stats = cluster_stats[old_cluster];
        for (const auto &neighbor_idx : neighbor_cache[index]) {
            stats.spatial_sum -= (allocations_ptr->operator()(neighbor_idx) == old_cluster ? 2 : 0);
        }
    }

    // Add to new cluster
    if (cluster != -1) {

        if (cluster >= cluster_stats.size()) {
            cluster_stats.resize(cluster + 1);
        }

        ClusterStats &stats = cluster_stats[cluster];
        for (const auto &neighbor_idx : neighbor_cache[index]) {
            stats.spatial_sum += (allocations_ptr->operator()(neighbor_idx) == cluster ? 2 : 0);
        }
    }
}