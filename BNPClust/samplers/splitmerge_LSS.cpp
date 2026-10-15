/**
 * @file splitmerge_LSS.cpp
 * @brief Implementation of Locality Sensitive Sampling (LSS) Split-Merge
 * sampler
 *
 * This file contains the complete implementation of the SplitMerge_LSS class,
 * which provides an optimized variant of split-merge sampling using locality
 * sensitive sampling for anchor point selection and sequential allocation for
 * proposal generation. This approach offers computational advantages for large
 * datasets.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge_LSS.hpp"
#include <random>

void SplitMerge_LSS::choose_indeces(bool similarity) {
  /**
   * @brief Randomly choose two distinct indices i and j from the data points.
   *        Identify their clusters ci and cj.
   *        Prepare the launch_state and S vectors for the split-merge
   * operation.
   */

  std::uniform_int_distribution<> dis(0, data.get_n() - 1);

  idx_i = dis(gen);

  auto distances = params.D.row(idx_i);
  double distance_sum = 0.0;
  std::vector<double> probs(data.get_n());

  if (similarity) {
    for (auto idx = 0; idx < data.get_n(); ++idx) {
      if (idx == idx_i)
        probs[idx] = 0.0;
      else
        probs[idx] = distances(idx);
      distance_sum += probs[idx];
    }

    for (auto idx = 0; idx < data.get_n(); ++idx) {
      probs[idx] /= distance_sum;
    }
  } else {
    for (auto idx = 0; idx < data.get_n(); ++idx) {
      if (idx == idx_i)
        probs[idx] = 0.0;
      else
        probs[idx] = 1 / distances(idx);
      distance_sum += probs[idx];
    }

    for (auto idx = 0; idx < data.get_n(); ++idx) {
      probs[idx] /= distance_sum;
    }
  }

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
  size_t launch_state_size =
      ci == cj ? size_ci - 2
               : size_ci + size_cj -
                     2; // Exclude points i and j from the launch state

  launch_state.resize(launch_state_size);
  S.resize(launch_state_size);
  original_allocations =
      data.get_allocations(); // Store original allocations in case of rejection

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

void SplitMerge_LSS::sequential_allocation(int iterations,
                                           bool only_probabilities,
                                           bool sequential) {
  /**
   * @brief Perform restricted Gibbs sampling on the points in S to propose new
   * allocations for iter iterations.
   * @param iterations Number of Gibbs sampling iterations to perform.
   * @param only_probabilities If true, only compute the probabilities without
   * @param sequential If true, unallocate points before sampling (default true)
   * changing allocations (used in merge move).
   */

  const int S_size = S.size();

  for (int i = 0; i < iterations; ++i) {

    // Unallocate all points in S
    for (int idx = 0; idx < S_size && sequential; ++idx) {
      data.set_allocation(S(idx), -1); // Unallocate point
    }

    for (int idx = 0; idx < S_size; ++idx) {
      int point_idx = S(idx);
      int current_cluster = launch_state(idx);

      if (!sequential) {
        current_cluster = data.get_cluster_assignment(point_idx);
        // Unallocate point only if using restricted Gibbs sampling
        data.set_allocation(point_idx, -1);
      }

      // Compute probabilities for each cluster (ci and cj)
      Eigen::Vector2d log_probs;

      log_probs(0) = likelihood.point_loglikelihood_cond(point_idx, ci);
      log_probs(1) = likelihood.point_loglikelihood_cond(point_idx, cj);

      log_probs(0) += process.gibbs_prior_existing_cluster(ci, point_idx);
      log_probs(1) += process.gibbs_prior_existing_cluster(cj, point_idx);

      // Normalize log probabilities using log-sum-exp trick
      double max_log_prob = log_probs.maxCoeff();
      double log_sum = max_log_prob + log((log_probs.array() - max_log_prob).exp().sum());

      if (!only_probabilities) {
        // Compute normalized probabilities only for sampling
        Eigen::Vector2d probs = (log_probs.array() - max_log_prob).exp();
        //probs /= probs.sum();

        // Sample new cluster based on computed probabilities
        std::discrete_distribution<int> dist(probs.data(), probs.data() + probs.size());
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

double
SplitMerge_LSS::compute_acceptance_ratio_merge(double likelihood_old_clusters) {
  /**
   * @brief Compute the log acceptance ratio for a merge move.
   * @param likelihood_old_clusters The joint log likelihood of clusters ci and cj before the merge.
   * @return The log acceptance ratio for the merge move.
   */

  // Prior ratio
  int size_old_ci = (original_allocations.array() == ci).count();
  int size_old_cj = (original_allocations.array() == cj).count();
  double log_acceptance_ratio = process.prior_ratio_merge(size_old_ci, size_old_cj);

  // Likelihood ratio
  log_acceptance_ratio += likelihood.cluster_loglikelihood(ci) - likelihood_old_clusters;

  // Proposal ratio
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS::merge_move() {
  /**
   * @brief Propose a merge move by combining clusters ci and cj.
   *        Compute the acceptance ratio and decide whether to accept or reject
   * the move.
   */

  // Reset log probabilities
  log_merge_gibbs_prob = 0;

  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  sequential_allocation(1, true); // only compute probabilities

  data.set_allocation(
      idx_j, ci); // Temporarily assign j to ci for likelihood computation

  // Propose new allocations by merging clusters ci and cj
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state(idx) == cj) {
      data.set_allocation(S(idx), ci); // Merge cj into ci
    }
  }

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_merge(likelihood_old_clusters);

  // Accept or reject the move
  std::uniform_real_distribution<> dis(0.0, 1.0);
  if (log(dis(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_merge++;
}

void SplitMerge_LSS::split_move() {
  /**
   * @brief Propose a split move by dividing cluster ci into two clusters.
   *        Compute the acceptance ratio and decide whether to accept or reject
   * the move.
   */

  double likelihood_old_cluster = likelihood.cluster_loglikelihood(
      ci); // likelihood of the merged cluster old configuration

  log_split_gibbs_prob = 0; // reset the log probability of the split move

  data.set_allocation(idx_j, data.get_K()); // Create a new cluster for point j
                                            // if they are in the same cluster
  cj = data.get_cluster_assignment(idx_j); // Update cj to the new cluster index

  // Allocate randomly points in S to either ci or cj
  std::uniform_int_distribution<> dis(0, 1);
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    int new_cluster = (dis(gen) == 0) ? ci : cj;
    data.set_allocation(S(idx), new_cluster);
  }

  // Perform restricted Gibbs sampling to refine the allocations
  sequential_allocation(1);

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  // Accept or reject the move
  std::uniform_real_distribution<> dis2(0.0, 1.0);
  if (log(dis2(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_split++;
}

double
SplitMerge_LSS::compute_acceptance_ratio_split(double likelihood_old_cluster) {
  /**
   * @brief Compute the log acceptance ratio for a split move.
   * @param likelihood_old_cluster The log likelihood of the original cluster
   * before the split.
   * @return The log acceptance ratio for the split move.
   */

  // Prior ratio
  double log_acceptance_ratio = process.prior_ratio_split(ci, cj);

  // Likelihood ratio
  log_acceptance_ratio +=
      likelihood.clusters_loglikelihood(ci, cj) - likelihood_old_cluster;

  // Proposal ratio
  log_acceptance_ratio -= log_split_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge_LSS::shuffle() {
  /**
   * @brief Perform a shuffle move to refine the allocations of points in S.
   *        This move helps to improve mixing by allowing points to switch
   * clusters.
   */

  log_split_gibbs_prob = 0; // reset the log probability of the split move - use
                            // it to store the proposal prob
  log_merge_gibbs_prob = 0;

  if (data.get_K() < 2)
    return; // No point in shuffling if there's only one cluster

  // Get number of points in clusters ci and cj and likelihoods
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  // Compute probabilities
  sequential_allocation(1, true, true); // only compute probabilities

  // Use restricted gibbs to refine the allocations
  sequential_allocation(1, false, true);

  // Compute acceptance ratio
  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(
      likelihood_old_clusters, old_ci_size, old_cj_size);

  // Accept or reject the move
  std::uniform_real_distribution<> acceptance_ratio_dis(0.0, 1.0);
  if (log(acceptance_ratio_dis(gen)) >
      log_acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_shuffle++;
}

double SplitMerge_LSS::compute_acceptance_ratio_shuffle(
    double likelihood_old_clusters, int old_ci_size, int old_cj_size) {
  /**
   * @brief Compute the log acceptance ratio for a shuffle move.
   * @return The log acceptance ratio for the shuffle move.
   */

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

void SplitMerge_LSS::choose_clusters_shuffle() {
  /**
   * @brief Randomly choose two distinct clusters ci and cj from the current
   * allocations. Update idx_i and idx_j to be random points from these
   * clusters.
   */

  if (data.get_K() < 2) {
    // std::cout << "Not enough clusters to perform shuffle." << std::endl;
    return;
  }

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

  // Store original allocations in case of rejection
  original_allocations = data.get_allocations();

  // Pre-allocate launch_state and S
  const int size_ci = data.get_cluster_size(ci);
  const int size_cj = data.get_cluster_size(cj);
  const int launch_state_size =
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

void SplitMerge_LSS::step() {
  /**
   * @brief Perform a single split-merge MCMC step.
   *        Randomly choose two indices and decide whether to propose a split or
   * merge move.
   */

  std::discrete_distribution<> dist({0.5, 0.5});
  choose_indeces(dist(gen) == 0 ? true : false);
  process.set_old_allocations(
      data.get_allocations()); // Update old allocations in the process
  process.set_idx_i(idx_i);
  process.set_idx_j(idx_j);

  if (ci == cj) {
    split_move();
  } else {
    merge_move();
  }

  if (shuffle_bool) {
    choose_clusters_shuffle();
    process.set_old_allocations(
        data.get_allocations()); // Update old allocations in the process
    process.set_idx_i(idx_i);
    process.set_idx_j(idx_j);
    shuffle();
  }
}