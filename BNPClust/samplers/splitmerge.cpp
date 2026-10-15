/**
 * @file splitmerge.cpp
 * @brief Implementation of Split-Merge MCMC sampler
 *
 * This file contains the complete implementation of the SplitMerge class,
 * providing split, merge, and shuffle moves for Bayesian nonparametric
 * clustering. The implementation includes restricted Gibbs sampling for
 * proposal generation and Metropolis-Hastings acceptance.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "splitmerge.hpp"

void SplitMerge::choose_indeces() {
  /**
   * @brief Randomly choose two distinct indices i and j from the data points.
   *        Identify their clusters ci and cj.
   *        Prepare the launch_state and S vectors for the split-merge
   * operation.
   */

  std::uniform_int_distribution<> dis(0, data.get_n() - 1);

  idx_i = dis(gen);
  do {
    idx_j = dis(gen);
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

void SplitMerge::restricted_gibbs(int iterations, bool only_probabilities) {
  /**
   * @brief Perform restricted Gibbs sampling on the points in S to propose new
   * allocations for iter iterations.
   * @param iterations Number of Gibbs sampling iterations to perform.
   * @param only_probabilities If true, only compute the probabilities without
   * changing allocations (used in merge move).
   */

  for (int i = 0; i < iterations; ++i) {
    for (int idx = 0; idx < S.size(); ++idx) {
      int point_idx = S(idx);
      int current_cluster = data.get_cluster_assignment(point_idx);

      // Remove point from its current cluster
      data.set_allocation(point_idx, -1); // Temporarily unassign the point

      // Compute probabilities for each cluster (ci and cj)
      Eigen::Vector2d log_probs;

      log_probs(0) = likelihood.point_loglikelihood_cond(point_idx, ci);
      log_probs(1) = likelihood.point_loglikelihood_cond(point_idx, cj);

      log_probs(0) += process.gibbs_prior_existing_cluster(ci, point_idx);
      log_probs(1) += process.gibbs_prior_existing_cluster(cj, point_idx);

      // Normalize to get probabilities
      double max_log_prob = log_probs.maxCoeff();
      Eigen::Vector2d probs = (log_probs.array() - max_log_prob).exp();
      probs /= probs.sum();

      if (!only_probabilities) {
        // Sample new cluster based on computed probabilities
        std::discrete_distribution<int> dist(probs.data(),
                                             probs.data() + probs.size());
        int new_cluster_idx = dist(gen);
        int new_cluster = (new_cluster_idx == 0) ? ci : cj;

        // Assign point to the new cluster
        data.set_allocation(point_idx, new_cluster);

        // if last iteration, accumulate the log probability of the move
        if (i == iterations - 1)
          log_split_gibbs_prob += log(probs(new_cluster_idx));

      } else {
        // Just restore the previous allocation
        data.set_allocation(point_idx, current_cluster);

        // accumulate the log probability of the move usually for the merge move
        log_merge_gibbs_prob += log(probs((current_cluster == ci) ? 0 : 1));
      }
    }
  }
}

double SplitMerge::compute_acceptance_ratio_merge(double likelihood_old_clusters) {
  /**
   * @brief Compute the log acceptance ratio for a merge move.
   * @param likelihood_old_clusters The joint log likelihood of clusters ci and cj before the merge.
   * @return The log acceptance ratio for the merge move.
   */

  // Prior ratio
  int size_old_ci = (original_allocations.array() == ci).count();
  int size_old_cj = (original_allocations.array() == cj).count();
  double log_prior_ratio = process.prior_ratio_merge(size_old_ci, size_old_cj);

  // Likelihood ratio
  double log_likelihood_ratio =
      likelihood.cluster_loglikelihood(ci) - likelihood_old_clusters;

  return log_prior_ratio + log_likelihood_ratio + log_merge_gibbs_prob;
}

void SplitMerge::merge_move() {
  /**
   * @brief Propose a merge move by combining clusters ci and cj.
   *        Compute the acceptance ratio and decide whether to accept or reject
   * the move.
   */

  // Reset log probabilities
  log_merge_gibbs_prob = 0;

  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);

  // CRITICAL: Compute the proposal probability BEFORE actually merging
  // This is the probability of generating the current split state from the
  // launch state
  restricted_gibbs(1, true); // Compute probability of current allocation

  data.set_allocation( idx_j, ci); // Temporarily assign j to ci for likelihood computation

  // Propose new allocations by merging clusters ci and cj
  for (int idx = 0; idx < launch_state.size(); ++idx) {
    if (launch_state(idx) == cj) {
      data.set_allocation(S(idx), ci); // Merge cj into ci
    }
  }

  // Compute acceptance ratio
  double acceptance_ratio = compute_acceptance_ratio_merge(likelihood_old_clusters);

  // Rcpp::Rcout << "[DEBUG - Merge] acceptance_ratio: " << acceptance_ratio
  //             << std::endl;

  // Accept or reject the move
  std::uniform_real_distribution<> dis(0.0, 1.0);
  if (log(dis(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_merge++;
}

void SplitMerge::split_move() {
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

  std::uniform_int_distribution<> dis(0, 1);
  for (int idx = 0; idx < launch_state.size(); ++idx) {
      int point = S(idx);
      int new_cluster = (dis(gen) == 0) ? ci : cj;  // Random
      data.set_allocation(point, new_cluster);
  }

  // Perform restricted Gibbs sampling to refine the allocations
  restricted_gibbs(5);

  // Compute acceptance ratio
  double acceptance_ratio =
      compute_acceptance_ratio_split(likelihood_old_cluster);

  // Rcpp::Rcout << "[DEBUG - Split] acceptance_ratio: " << acceptance_ratio
  //             << std::endl;

  // Accept or reject the move
  std::uniform_real_distribution<> dis2(0.0, 1.0);
  if (log(dis2(gen)) > acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_split++;
}

double
SplitMerge::compute_acceptance_ratio_split(double likelihood_old_cluster) {
  /**
   * @brief Compute the log acceptance ratio for a split move.
   * @param likelihood_old_cluster The log likelihood of the original cluster
   * before the split.
   * @return The log acceptance ratio for the split move.
   */

  // Prior ratio
  double log_prior_ratio = process.prior_ratio_split(ci, cj);

  // Likelihood ratio
  double log_likelihood_ratio =
      likelihood.clusters_loglikelihood(ci, cj) - likelihood_old_cluster;

  // Rcpp::Rcout << "\t[DEBUG - Split] log_prior_ratio: " << log_prior_ratio
  //             << ", log_likelihood_ratio: " << log_likelihood_ratio
  //             << ", log_split_gibbs_prob: " << log_split_gibbs_prob
  //             << std::endl;

  return -log_split_gibbs_prob + log_likelihood_ratio + log_prior_ratio;
}

void SplitMerge::shuffle() {
  /**
   * @brief Perform a shuffle move to refine the allocations of points in S.
   *        This move helps to improve mixing by allowing points to switch
   *        clusters.
   */

  log_split_gibbs_prob = 0; // reset the log probability of the split move - use
                            // it to store the proposal prob
  log_merge_gibbs_prob = 0;

  if (data.get_K() < 2) {
    return; // No point in shuffling if there's only one cluster
  }

  // Compute old likelihoods and sizes
  double likelihood_old_clusters = likelihood.clusters_loglikelihood(ci, cj);
  int old_ci_size = data.get_cluster_size(ci);
  int old_cj_size = data.get_cluster_size(cj);

  // Use restricted gibbs to refine the allocations
  restricted_gibbs(5);

  // Compute acceptance ratio
  double log_acceptance_ratio = compute_acceptance_ratio_shuffle(
      likelihood_old_clusters, old_ci_size, old_cj_size);

  // Accept or reject the move
  std::uniform_real_distribution<> acceptance_ratio_dis(0.0, 1.0);
  if (log(acceptance_ratio_dis(gen)) > log_acceptance_ratio) // move not accepted
    data.set_allocations(original_allocations);
  else
    accepted_shuffle++;
}

double SplitMerge::compute_acceptance_ratio_shuffle(double likelihood_old_clusters,
                                                    int old_ci_size,
                                                    int old_cj_size) {
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
  restricted_gibbs(1, true); // only compute probabilities
  log_acceptance_ratio += log_merge_gibbs_prob;

  return log_acceptance_ratio;
}

void SplitMerge::choose_clusters_shuffle() {
  /**
   * @brief Randomly choose two distinct clusters ci and cj from the current
   * allocations. Update idx_i and idx_j to be random points from these
   * clusters.
   */

  if (data.get_K() < 2){
    //std::cout << "Not enough clusters to perform shuffle." << std::endl;
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
  const int launch_state_size = size_ci + size_cj - 2; // Exclude points i and j from the launch state
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

void SplitMerge::step() {
  /**
   * @brief Perform a single split-merge MCMC step.
   *        Randomly choose two indices and decide whether to propose a split or
   * merge move.
   */

  choose_indeces();
  process.set_old_allocations(data.get_allocations()); // Update old allocations in the process
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