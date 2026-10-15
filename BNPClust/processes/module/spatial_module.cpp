/**
 * @file spatial_module.cpp
 * @brief Implementation of `SpatialModule`.
 */

#include "spatial_module.hpp"

double SpatialModule::compute_similarity_obs(int obs_idx, int cls_idx) const {
    int neighbors = 0;
    const std::vector<int> &row = neighbor_cache[obs_idx];

    // Iterate through cached neighbors and count those in the specified cluster
    for (size_t i = 0; i < row.size(); ++i) {
        int cluster_i = data_module.get_cluster_assignment(row[i]);
        if (cluster_i != -1 && cluster_i == cls_idx) {
            neighbors += 1;
        }
    }

    return spatial_weight * neighbors;
}

void SpatialModule::neighbor_cache_compute() {
    const int N = data_module.get_n();
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

double SpatialModule::compute_similarity_cls(int cls_idx, bool old_allo) const {

    double total_neighbors = 0;
    if (old_allo && old_cluster_members_provider) {
        const auto members_it = old_cluster_members_provider->find(cls_idx);
        if (members_it == old_cluster_members_provider->end()) {
            return 0.0;
        }

        for (const auto &i : members_it->second) {
            // Use cached neighbor indices instead of iterating over full adjacency matrix
            const std::vector<int> &row = neighbor_cache[i];

            // Count neighbors in each cluster
            for (size_t j = 0; j < row.size(); ++j) {
                int neighbor_idx = row[j];
                int cluster_i = data_module.get_cluster_assignment(neighbor_idx);

                if (cluster_i != -1 && cluster_i == cls_idx) {
                    ++total_neighbors;
                }
            }
        }
    } else {
        const Eigen::VectorXi &cls_idx_allocations =
            data_module.get_cluster_assignments(cls_idx);
        for(auto && i : cls_idx_allocations){
            // Use cached neighbor indices instead of iterating over full adjacency matrix
            const std::vector<int> &row = neighbor_cache[i];

            // Count neighbors in each cluster
            for (size_t j = 0; j < row.size(); ++j) {
                int neighbor_idx = row[j];
                int cluster_i = data_module.get_cluster_assignment(neighbor_idx);

                if (cluster_i != -1 && cluster_i == cls_idx) {
                    ++total_neighbors;
                }
            }
        }
    }

    return spatial_weight * total_neighbors / 2; // Each edge counted twice
}

Eigen::VectorXd SpatialModule::compute_similarity_obs(int obs_idx) const {
    Eigen::VectorXd cluster_adjacency = Eigen::VectorXd::Zero(data_module.get_K());

    // Use cached neighbor indices instead of iterating over full adjacency matrix
    const std::vector<int> &row = neighbor_cache[obs_idx];

    // Count neighbors in each cluster
    for (size_t i = 0; i < row.size(); ++i) {
        int neighbor_idx = row[i];
        int cluster_i = data_module.get_cluster_assignment(neighbor_idx);

        if (cluster_i != -1) {
            cluster_adjacency(cluster_i) += spatial_weight;
        }
    }

    return cluster_adjacency;
}