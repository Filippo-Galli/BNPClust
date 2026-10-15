/**
* @file continuos_cache.cpp
* @brief Implementation of `ContinuosCache`.
*/

#include "continuos_cache.hpp"

void ContinuosCache::set_allocation(int index, int cluster, int old_cluster) {

    const double value = continuos_covariates(index);
    // Remove point from old cluster stats
    if (old_cluster != -1) {
        ClusterStats &old_stats = cluster_stats[old_cluster];
        old_stats.n--;
        old_stats.sum -= value;
        old_stats.sumsq -= value * value;
    }

    // Add point to new cluster stats
    if (cluster != -1) {

        // Ensure cluster_stats vector is large enough
        if (cluster >= static_cast<int>(cluster_stats.size())) {
            cluster_stats.resize(cluster + 1);
        }

        ClusterStats &new_stats = cluster_stats[cluster];
        new_stats.n++;
        new_stats.sum += value;
        new_stats.sumsq += value * value;
    }
}

void ContinuosCache::recompute(const int K, const Eigen::VectorXi &allocations) {
    // Clear existing stats
    cluster_stats.clear();

    // Get the number of clusters from covariates allocations
    cluster_stats.resize(K);

    // Recompute stats from scratch
    for (int i = 0; i < allocations.size(); ++i) {
        int cluster = allocations(i);
        if (cluster < 0)
            continue; // Skip unallocated points

        double value = continuos_covariates(i);

        ClusterStats &stats = cluster_stats[cluster];
        stats.n++;
        stats.sum += value;
        stats.sumsq += value * value;
    }
}

void ContinuosCache::remove_info(int cluster) {
    // Remove stats for the specified cluster
    if (cluster == static_cast<int>(cluster_stats.size()) - 1) {
        cluster_stats.pop_back();
    } else {
        cluster_stats.erase(cluster_stats.begin() + cluster);
    }
}