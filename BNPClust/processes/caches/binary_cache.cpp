/**
 * @file binary_cache.cpp
 * @brief Implementation of `BinaryCache`.
 */

#include "binary_cache.hpp"

void BinaryCache::recompute(const int K, const Eigen::VectorXi &allocations_in) {
    // Clear existing stats
    cluster_stats.clear();

    // Resize to number of clusters
    cluster_stats.resize(K);

    // Compute stats
    for (int i = 0; i < allocations_in.size(); ++i) {
        int cluster = allocations_in(i);
        if (cluster < 0)
            continue;

        ClusterStats &stats = cluster_stats[cluster];
        stats.n++;
        stats.binary_sum += binary_covariates(i);
    }
}

void BinaryCache::set_allocation(int index, int cluster, int old_cluster) {

    // Remove from old cluster
    if (old_cluster != -1) {
        ClusterStats &stats = cluster_stats[old_cluster];
        stats.n--;
        stats.binary_sum -= binary_covariates(index);
    }

    // Add to new cluster
    if (cluster != -1) {

        if (cluster >= cluster_stats.size()) {
            cluster_stats.resize(cluster + 1);
        }

        ClusterStats &stats = cluster_stats[cluster];
        stats.n++;
        stats.binary_sum += binary_covariates(index);
    }
}