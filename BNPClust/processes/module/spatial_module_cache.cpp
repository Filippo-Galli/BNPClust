/**
 * @file spatial_module_cache.cpp
 * @brief Implementation of `SpatialModuleCache`.
 */

#include "spatial_module_cache.hpp"

double SpatialModuleCache::compute_similarity_obs(int obs_idx, int cls_idx) const {
    int neighbors = 0;
    const std::vector<int> &row = cache.neighbor_cache[obs_idx];
    const int *allocations = data_module.get_allocations().data();

    for (size_t i = 0; i < row.size(); ++i) {
        const int cluster_i = allocations[row[i]];
        if (cluster_i == cls_idx && cluster_i != -1) {
            ++neighbors;
        }
    }

    return spatial_weight * neighbors;
}

double SpatialModuleCache::compute_similarity_cls(int cls_idx, bool old_allo) const {

    if (old_allo && old_cluster_members_provider) {
        const auto members_it = old_cluster_members_provider->find(cls_idx);
        if (members_it == old_cluster_members_provider->end()) {
            return 0.0;
        }
        const auto &members = members_it->second;
        const int *allocations = old_allocations_provider->data();

        double total_neighbors = 0;
        for (size_t i = 0; i < members.size(); ++i) {
            const std::vector<int> &row = cache.neighbor_cache[members[i]];

            for (size_t j = 0; j < row.size(); ++j) {
                if (allocations[row[j]] == cls_idx && allocations[row[j]] != -1) {
                    ++total_neighbors;
                }
            }
        }
        return spatial_weight * total_neighbors / 2; // Each edge counted twice
    }

    return cache.get_cluster_stats_ref(cls_idx).spatial_sum * spatial_weight / 2; // Each edge counted twice
}

Eigen::VectorXd SpatialModuleCache::compute_similarity_obs(int obs_idx) const {
    Eigen::VectorXd cluster_adjacency = Eigen::VectorXd::Zero(data_module.get_K());
    const std::vector<int> &row = cache.neighbor_cache[obs_idx];
    const int *allocations = data_module.get_allocations().data();

    for (size_t i = 0; i < row.size(); ++i) {
        const int cluster_i = allocations[row[i]];
        if (cluster_i != -1) {
            cluster_adjacency(cluster_i) += spatial_weight;
        }
    }

    return cluster_adjacency;
}