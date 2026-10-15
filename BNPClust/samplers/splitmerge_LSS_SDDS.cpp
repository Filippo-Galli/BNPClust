/**
 * @file splitmerge_LSS_SDDS.cpp
 * @brief Implementation of Locality Sensitive Sampling (LSS) with SDDS
 * Split-Merge sampler
 *
 * This file contains the complete implementation of the SplitMerge_LSS_SDDS
 * class, which provides an optimized variant of split-merge sampling using:
 * - Locality sensitive sampling for anchor point selection
 * - SDDS (Smart-split-Dumb-merge, Dumb-split-Smart-merge) for adaptive move
 * selection to achieve better mixing and perfomance.
 *
 * This approach offers computational advantages for large datasets while
 * maintaining theoretical guarantees of split-merge samplers.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_LSS_SDDS.hpp"
#include <random>

void SplitMerge_LSS_SDDS::choose_indeces(bool similarity) {

  // Select first index idx_i uniformly at random
  std::uniform_int_distribution<> dis(0, data.get_n() - 1);
  idx_i = dis(gen);

  const auto distances = params.D.row(idx_i);
  const double max_distance = row_max_distances[idx_i];
  double prob_sum = 0.0;
  std::vector<double> probs(data.get_n());
  const double eps = 1e-4;

  // Compute probabilities for selecting idx_j based on distances from idx_i
  // If similarity is true, prefer closer points; otherwise, prefer distant
  // points using a Cauchy (rational) kernel scaled by temperature (median heuristic)
  for (auto idx = 0; idx < data.get_n(); ++idx) {
    if (idx == idx_i) {
      probs[idx] = 0.0;
    } else {
      if (similarity) {
        probs[idx] = 1.0 / (1.0 + distances(idx) / temperature);
      } else {
        probs[idx] = 1.0 / (1.0 + (max_distance - distances(idx)) / temperature);
      }
    }
    prob_sum += probs[idx];
  }

  // // Normalize probabilities
  // for (auto &p : probs) {
  //     p /= distance_sum;
  // }

  // Select second index idx_j based on computed probabilities
  std::discrete_distribution<> dis1(probs.begin(), probs.end());
  do {
    idx_j = dis1(gen);
  } while (idx_j == idx_i); // Ensure i and j are distinct

  // Get the clusters of the chosen indices
  ci = data.get_cluster_assignment(idx_i);
  cj = data.get_cluster_assignment(idx_j);

  // Get number of points in clusters ci and cj
  int size_ci = data.get_cluster_size(ci);
  int size_cj = data.get_cluster_size(cj);

  // Initialize launch_state with current allocations of points in clusters ci
  // and cj while S with their indices
  launch_state_size =
      ci == cj ? size_ci - 2
               : size_ci + size_cj -
                     2; // Exclude points i and j from the launch state

  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);

  // Properly collect all points from clusters ci and cj
  int s_idx = 0;
  for (int idx = 0; idx < data.get_n(); ++idx) {
    if (idx == idx_i || idx == idx_j)
      continue; // Skip points i and j

    int temp_cluster = data.get_cluster_assignment(idx);
    if (temp_cluster == ci || temp_cluster == cj) {
      S(s_idx) = idx;
      launch_state(s_idx) = temp_cluster;
      s_idx++;
    }
  }

  // Shuffle launch_state and S in unison using Fisher-Yates algorithm
  // This maintains alignment between the two vectors while shuffling in-place
  for (int i = launch_state_size - 1; i > 0; --i) {
    std::uniform_int_distribution<> dis_shuffle(0, i);
    int j = dis_shuffle(gen);
    if (i != j) {
      std::swap(launch_state(i), launch_state(j));
      std::swap(S(i), S(j));
    }
  }

// Ensure we collected the expected number of points
#if VERBOSITY_LEVEL >= 1
  if (s_idx != launch_state_size) { // since s_idx is zero-based
    // Rcpp::Rcout << "[ERROR] s_idx = " << s_idx << ", launch_state_size = " <<
    // launch_state_size << std::endl;
    throw std::runtime_error(
        "Mismatch in expected cluster sizes during split-merge initialization");
  }
#endif
}

void SplitMerge_LSS_SDDS::sequential_allocation(int iterations,
                                                bool only_probabilities,
                                                bool sequential) {

  for (int i = 0; i < iterations; ++i) {

    // Unallocate all points in S
    for (int idx = 0; idx < launch_state_size && sequential; ++idx) {
      data.set_allocation(S(idx), -1); // Unallocate point
    }

    // Allocate each point in S using RGSM or sequential allocation
    for (int idx = 0; idx < launch_state_size; ++idx) {
      int point_idx = S(idx);
      int current_cluster = launch_state(idx);

      if (!sequential) {
        current_cluster = data.get_cluster_assignment(point_idx);
        // Unallocate point only if using restricted Gibbs sampling
        data.set_allocation(point_idx, -1);
      }

      // Compute probabilities for each cluster (ci and cj)
      log_probs(0) = likelihood.point_loglikelihood_cond(point_idx, ci);
      log_probs(1) = likelihood.point_loglikelihood_cond(point_idx, cj);

      log_probs(0) += process.gibbs_prior_existing_cluster(ci, point_idx);
      log_probs(1) += process.gibbs_prior_existing_cluster(cj, point_idx);

      // Normalize log probabilities using log-sum-exp trick
      double max_log_prob = log_probs.maxCoeff();
      double log_sum =
          max_log_prob + log((log_probs.array() - max_log_prob).exp().sum());

      if (!only_probabilities) {
        // Compute normalized probabilities only for sampling
        probs = (log_probs.array() - max_log_prob).exp();

        // Sample new cluster based on computed probabilities
        std::discrete_distribution<int> dist(probs.data(),
                                             probs.data() + probs.size());
        int new_cluster_idx = dist(gen);
        int new_cluster = (new_cluster_idx == 0) ? ci : cj;

        // Assign point to the new cluster
        data.set_allocation(point_idx, new_cluster);

        // if last iteration, accumulate the log probability of the move
        if (i == iterations - 1)
          log_split_gibbs_prob += log_probs(new_cluster_idx) - log_sum;
      } else {
        // Just restore the previous allocation
        data.set_allocation(point_idx, current_cluster);

        // accumulate the log probability of the move
        int cluster_idx = (current_cluster == ci) ? 0 : 1;
        log_merge_gibbs_prob += log_probs(cluster_idx) - log_sum;
      }
    }
  }
}

double SplitMerge_LSS_SDDS::compute_acceptance_ratio_merge(
    double likelihood_old_clusters) {

  // Prior ratio
  int size_old_ci = (original_allocations.array() == ci).count();
  int size_old_cj = (original_allocations.array() == cj).count();
  double log_acceptance_ratio =
      process.prior_ratio_merge(size_old_ci, size_old_cj);

  // Likelihood ratio
  log_acceptance_ratio +=
      likelihood.cluster_loglikelihood(ci) - likelihood_old_clusters;

  // Proposal ratio of the reverse move (smart or dumb split)
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS_SDDS::smart_merge_move() {

  // Reverse dumb split proposal probability
  log_merge_gibbs_prob = S.size() * rand_split_prob;

  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  // Assign both anchor points to ci
  data.set_allocation(idx_j, ci);

  // Initially assign all points from both clusters to ci
  for (const auto &idx : S) {
    data.set_allocation(idx, ci);
  }

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_merge(likelihood_old_clusters);

  // Accept or reject the move
  if (log(dis_real(gen)) > acceptance_ratio) // move not accepted
    process.restore_state();
  else
    accepted_merge++;
}

void SplitMerge_LSS_SDDS::dumb_merge_move() {

  // Reverse smart split proposal probability
  sequential_allocation(1, true, true);

  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  // Direct merge: assign all points to ci
  data.set_allocation(idx_j, ci);
  for (int idx = 0; idx < launch_state_size; ++idx) {
    if (launch_state(idx) == cj) {
      data.set_allocation(S(idx), ci); // Merge cj into ci
    }
  }

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_merge(likelihood_old_clusters);

  // Accept or reject the move
  if (log(dis_real(gen)) > acceptance_ratio) // move not accepted
    process.restore_state();
  else
    accepted_merge++;
}

void SplitMerge_LSS_SDDS::smart_split_move() {

  // old cluster likelihood
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);

  // Create new cluster for point j
  data.set_allocation(idx_j, data.get_K()); // Create a new cluster for point j
  cj = data.get_cluster_assignment(idx_j); // Update cj to the new cluster index

  // Perform sequential allocation to determine the allocations
  sequential_allocation(1, false, true);

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  // Accept or reject the move
  if (log(dis_real(gen)) > acceptance_ratio) // move not accepted
    process.restore_state();
  else
    accepted_split++;
}

void SplitMerge_LSS_SDDS::dumb_split_move() {

  // old cluster likelihood
  double likelihood_old_cluster = likelihood.cluster_loglikelihood(ci);

  // Create new cluster for point j
  data.set_allocation(idx_j, data.get_K()); // Create a new cluster for point j
  cj = data.get_cluster_assignment(idx_j); // Update cj to the new cluster index

  // Randomly allocate points in S to either ci or cj
  for (const auto &idx : S) {
    int new_cluster = (dis_int(gen) == 0) ? ci : cj;
    data.set_allocation(idx, new_cluster);
  }

  // Random split proposal probability
  log_split_gibbs_prob =
      S.size() * rand_split_prob; // Probability of the dumb split move

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  // Accept or reject the move
  if (log(dis_real(gen)) > acceptance_ratio) // move not accepted
    process.restore_state();
  else
    accepted_split++;
}

double SplitMerge_LSS_SDDS::compute_acceptance_ratio_split(
    double likelihood_old_cluster) {

  // Prior ratio
  double log_acceptance_ratio = process.prior_ratio_split(ci, cj);

  // Likelihood ratio
  log_acceptance_ratio +=
      likelihood.clusters_loglikelihood(ci, cj) - likelihood_old_cluster;

  // Proposal ratio (dumb or smart split as forwad move)
  log_acceptance_ratio -= log_split_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS_SDDS::shuffle() {

  if (data.get_K() < 2)
    return; // No point in shuffling if there's only one cluster

  // Get number of points in clusters ci and cj and likelihoods
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  // Compute probabilities
  sequential_allocation(1, true, true);

  // Use restricted gibbs to refine the allocations
  sequential_allocation(1, false, true);

  // Compute acceptance ratio
  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(
      likelihood_old_clusters, old_ci_size, old_cj_size);

  // Accept or reject the move
  if (log(dis_real(gen)) > log_acceptance_ratio) // move not accepted
    process.restore_state();
  else
    accepted_shuffle++;
}

double SplitMerge_LSS_SDDS::compute_acceptance_ratio_shuffle(
    double likelihood_old_clusters, int old_ci_size, int old_cj_size) {

  // Prior ratio
  double log_acceptance_ratio =
      process.prior_ratio_shuffle(old_ci_size, old_cj_size, ci, cj);

  // Likelihood ratio
  log_acceptance_ratio +=
      likelihood.clusters_loglikelihood(ci, cj) - likelihood_old_clusters;

  // Proposal ratio
  log_acceptance_ratio -= log_split_gibbs_prob;
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS_SDDS::choose_clusters_shuffle() {

  if (data.get_K() < 2) {
    // std::cout << "Not enough clusters to perform shuffle." << std::endl;
    return;
  }

  // Choose two distinct clusters ci and cj uniformly at random
  std::uniform_int_distribution<> dis(0, data.get_K() - 1);
  ci = dis(gen);
  do {
    cj = dis(gen);
  } while (cj == ci); // Ensure ci and cj are distinct

  // Choose random points from clusters ci and cj
  std::uniform_int_distribution<> dis_idx_i(0, data.get_cluster_size(ci) - 1);
  idx_i = data.get_cluster_assignments(ci)[dis_idx_i(gen)];

  std::uniform_int_distribution<> dis_idx_j(0, data.get_cluster_size(cj) - 1);
  idx_j = data.get_cluster_assignments(cj)[dis_idx_j(gen)];

  // Pre-allocate launch_state and S
  const int size_ci = data.get_cluster_size(ci);
  const int size_cj = data.get_cluster_size(cj);
  launch_state_size =
      size_ci + size_cj - 2; // Exclude points i and j from the launch state
  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);

  // Create S and launch_state
  int s_idx = 0;
  for (int idx = 0; idx < data.get_n(); ++idx) {
    if (idx == idx_i || idx == idx_j)
      continue; // Skip points i and j

    int temp_cluster = data.get_cluster_assignment(idx);
    if (temp_cluster == ci || temp_cluster == cj) {
      S(s_idx) = idx;
      launch_state(s_idx) = temp_cluster;
      s_idx++;
    }
  }
}

void SplitMerge_LSS_SDDS::step() {

  std::discrete_distribution<> dist({0.5, 0.5});
  const int similarity_dist =
      dist(gen); // 0 for dissimilarity (split), 1 for similarity (merge)
  process.set_old_cluster_members(
      data.get_cluster_map()); // Update old cluster members in the process
  process.set_old_allocations(
      data.get_allocations());     // Update old allocations in the process
  process.set_old_K(data.get_K()); // Cache number of clusters for rollback
  choose_indeces(similarity_dist);
  process.set_idx_i(idx_i);
  process.set_idx_j(idx_j);

  // Reset log proposal probabilities
  log_split_gibbs_prob = 0;
  log_merge_gibbs_prob = 0;

  // Determine which move to perform based on strategy
  if (similarity_dist) {
    // smart merge - dumb split
    if (ci != cj) {
      merge_moves++;
      smart_merge_move();
    } else {
      split_moves++;
      dumb_split_move();
    }
  } else {
    // smart split - dumb merge
    if (ci == cj) {
      split_moves++;
      smart_split_move();
    } else {
      merge_moves++;
      dumb_merge_move();
    }
  }

  // Reset log proposal probabilities
  log_split_gibbs_prob = 0;
  log_merge_gibbs_prob = 0;

  if (shuffle_bool) {
    shuffle_moves++;
    process.set_old_allocations(
        data.get_allocations()); // Update old allocations in the process
    process.set_old_cluster_members(
        data.get_cluster_map()); // Update old cluster members in the process
    process.set_old_K(data.get_K()); // Cache number of clusters for rollback
    choose_clusters_shuffle();
    process.set_idx_i(idx_i);
    process.set_idx_j(idx_j);
    shuffle();
  }
}
