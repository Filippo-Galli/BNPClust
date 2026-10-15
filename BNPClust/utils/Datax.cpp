/**
 * @file Datax.cpp
 * @brief Implementation of Datax class
 */


#include "Datax.hpp"

void Datax::compact_cluster(int old_cluster) {
#if VERBOSITY_LEVEL >= 1
    if (old_cluster < 0 || old_cluster >= K) {
        throw std::out_of_range("old_cluster index out of bounds in compact_cluster");
    }
#endif

    const int last_cluster = K - 1;

    // If compacting the last cluster (or the only cluster), this is a pure deletion.
    if (K <= 1 || old_cluster == last_cluster) {
        cluster_members.erase(old_cluster);
        // Update cluster info to remove the old cluster
        for (auto && ci : cluster_info) {
            ci->remove_info(old_cluster);
        }
        K--;
        return;
    }

    auto last_cluster_it = cluster_members.find(last_cluster);
    bool last_cluster_exists = (last_cluster_it != cluster_members.end());

    // Shift allocations of the last cluster to the old cluster
    if (last_cluster_exists && !last_cluster_it->second.empty()) {
        for (int point_index : last_cluster_it->second) {
            allocations(point_index) = old_cluster;
        }

        // Move members of the last cluster to the old cluster using move semantics
        cluster_members[old_cluster] = std::move(last_cluster_it->second);
        // Update cluster info accordingly
        for (auto && ci : cluster_info) 
            ci->move_cluster_info(last_cluster, old_cluster);

    } else {
        // Last cluster is empty or doesn't exist, just clear the old cluster
        cluster_members[old_cluster].clear();
    }

    // Remove the last cluster
    cluster_members.erase(last_cluster);
    // Update cluster info to remove the last cluster
    for (auto && ci : cluster_info)
        ci->remove_info(last_cluster);

    K--; // Decrease the number of clusters
}

void Datax::set_allocation(int index, int cluster) {

    int old_cluster = allocations(index);

    if (old_cluster == cluster) {
        return;
    }

    Data::set_allocation_wo_compaction(index, cluster);
    for(auto && ci : cluster_info)
        ci->set_allocation(index, cluster, old_cluster);

    // Check if old cluster became empty and needs compaction
    auto old_cluster_it = cluster_members.find(old_cluster);
    if (old_cluster_it != cluster_members.end() && old_cluster_it->second.empty()) {
        Datax::compact_cluster(old_cluster);
    }
}

void Datax::set_allocations(const Eigen::VectorXi &new_allocations) {

    Data::set_allocations(new_allocations);
    for(auto && ci : cluster_info)
        ci->recompute(K, allocations);
}

void Datax::restore_state(Eigen::VectorXi &old_allocations,
                                      std::unordered_map<int, std::vector<int>> &old_cluster_members, int old_K) {

    Data::restore_state(old_allocations, old_cluster_members, old_K);
    for(auto && ci : cluster_info)
        ci->recompute(K, allocations);
}