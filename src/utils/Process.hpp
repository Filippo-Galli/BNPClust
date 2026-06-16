/**
 * @file Process.hpp
 * @brief Abstract interface for Bayesian nonparametric processes
 *
 * This file defines the abstract Process base class that provides a unified
 * interface for different Bayesian nonparametric processes (DP, NGGP, DPW,
 * NGGPW). It defines the essential methods needed for both Gibbs sampling and
 * split-merge MCMC algorithms.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "Data.hpp"
#include <Eigen/Dense>
#include <unordered_map>

/**
 * @brief Abstract base class for Bayesian nonparametric processes
 *
 * This class provides the foundation for implementing different types of
 * Bayesian nonparametric processes including Dirichlet Process (DP), Normalized
 * Generalized Gamma Process (NGGP), and their weighted variants (DPW, NGGPW).
 * It defines the interface for prior computations needed in both Gibbs sampling
 * and split-merge MCMC algorithms.
 *
 * @details The class handles:
 * - Prior probability computations for existing and new clusters
 * - Prior ratio calculations for split-merge moves
 * - Parameter updates during MCMC sampling
 * - State management for split-merge algorithms
 *
 * Derived classes must implement the pure virtual methods to define their
 * specific process characteristics (DP, NGGP, NGGPW, DPW).
 *
 * @see DP, NGGP, NGGPW, DPW
 */
class Process {
protected:
  /** @brief Reference to the data object containing observations and
   * allocations */
  Data &data;

  /** @brief Storage for previous allocations to enable rollback in case of
   * rejection or for computation requiring it*/
  Eigen::VectorXi old_allocations;

  /** @brief Storage for previous cluster members to enable rollback in case of
   * rejection or for computation requiring it*/
  std::unordered_map<int, std::vector<int>> old_cluster_members;

  /** @brief Number of clusters associated with the stored previous state */
  int old_K = 0;

  /**
   * @brief Provides read-only access to the stored previous cluster members.
   *
   * @return Const reference to the unordered map of previous cluster members.
   */
  [[nodiscard]] const std::unordered_map<int, std::vector<int>> &
  old_cluster_members_view() const {
    return old_cluster_members;
  }

  /** @brief Index of first observation involved in split-merge move */
  int idx_i;

  /** @brief Index of second observation involved in split-merge move */
  int idx_j;

public:
  /**
   * @brief Constructor initializing process with data and parameters
   *
   * @param d Reference to Data object containing observations and current
   * allocations
   * @param p Reference to Params object containing process hyperparameters
   *
   * @details Initializes the process and stores a copy of current allocations
   * for potential rollback operations during split-merge moves.
   */
  Process(Data &d) : data(d) {
    old_allocations = data.get_allocations();
    old_cluster_members = data.get_cluster_map();
    old_K = data.get_K();
  };

  // ========== Gibbs Sampling Methods ==========

  /**
   * @brief Compute prior probability for assigning observation to existing
   * cluster
   *
   * @param cls_idx Index of the existing cluster
   * @param obs_idx Index of the observation to be assigned
   * @return Log prior probability of assignment to existing cluster
   *
   * @details This method computes the prior component of the probability for
   * assigning an observation to an existing cluster in Gibbs sampling. The
   * exact computation depends on the specific process (DP, NGGP, etc.).
   */
  virtual double gibbs_prior_existing_cluster(int cls_idx,
                                              int obs_idx) const = 0;

  /**
   * @brief Compute prior probabilities for assigning observation to all
   * existing clusters
   * @param obs_idx Index of the observation to be assigned
   * @return Vector of log prior probabilities for all existing clusters
   * @details This method computes the prior components of the probabilities for
   * assigning an observation to each existing cluster in Gibbs sampling. The
   * exact computation depends on the specific process (DP, NGGP, etc.).
   * The returned vector has a length equal to the current number of clusters.
   */

  virtual Eigen::VectorXd gibbs_prior_existing_clusters(int obs_idx) const = 0;

  /**
   * @brief Compute prior probability for creating a new cluster
   *
   * @return Log prior probability of creating a new cluster
   *
   * @details This method computes the prior probability of assigning an
   * observation to a new cluster in Gibbs sampling. The computation is
   * process-specific.
   */
  virtual double gibbs_prior_new_cluster() const = 0;

  /**
   * @brief Compute prior probability for creating a new cluster for a specific
   * observation.
   *
   * @param obs_idx Index of the observation being (re)assigned.
   * @return Log prior probability of assigning obs_idx to a new cluster.
   *
   * @details
   * Some models (e.g., PPMx covariate similarity) add an observation-dependent
   * term to the "new cluster" option. By default, this falls back to the
   * observation-agnostic prior returned by gibbs_prior_new_cluster().
   */
  virtual double gibbs_prior_new_cluster_obs(int obs_idx) const {
    (void)obs_idx;
    return gibbs_prior_new_cluster();
  }

  // ========== Split-Merge Algorithm Methods ==========

  /**
   * @brief Compute prior ratio for split move in split-merge algorithm
   *
   * @param ci Index of first resulting cluster after split
   * @param cj Index of second resulting cluster after split
   * @return Log prior ratio for the split move
   *
   * @details Computes the ratio of prior probabilities when splitting a cluster
   * into two clusters. Used in the acceptance probability of split moves.
   */
  virtual double prior_ratio_split(int ci, int cj) const = 0;

  /**
   * @brief Compute prior ratio for merge move in split-merge algorithm
   *
   * @param size_old_ci Size of first cluster before merge
   * @param size_old_cj Size of second cluster before merge
   * @return Log prior ratio for the merge move
   *
   * @details Computes the ratio of prior probabilities when merging two
   * clusters into one cluster. Used in the acceptance probability of merge
   * moves.
   */
  virtual double prior_ratio_merge(int size_old_ci, int size_old_cj) const = 0;

  /**
   * @brief Compute prior ratio for shuffle move in split-merge algorithm
   *
   * @param size_old_ci Size of first cluster before shuffle
   * @param size_old_cj Size of second cluster before shuffle
   * @param ci Index of first cluster after shuffle
   * @param cj Index of second cluster after shuffle
   * @return Log prior ratio for the shuffle move
   *
   * @details Computes the ratio of prior probabilities when redistributing
   * observations between two existing clusters. Used in restricted Gibbs
   * sampling.
   */
  virtual double prior_ratio_shuffle(int size_old_ci, int size_old_cj, int ci,
                                     int cj) const = 0;

  // ========== State Management Methods ==========

  /**
   * @brief Store current allocations for potential rollback
   *
   * @param new_allocations Current allocation vector to store
   *
   * @details Saves the current state of cluster allocations before attempting
   * a split-merge move, enabling rollback if the move is rejected.
   */
  void set_old_allocations(const Eigen::VectorXi &new_allocations) {
    old_allocations = new_allocations;
  };

  /**
   * @brief Store current cluster members for potential rollback
   *
   * @param new_cluster_members Current cluster members map to store
   *
   * @details Saves the current state of cluster members before attempting
   * a split-merge move, enabling rollback if the move is rejected.
   */
  void set_old_cluster_members(
      const std::unordered_map<int, std::vector<int>> &new_cluster_members) {
    old_cluster_members = new_cluster_members;
  };

  /**
   * @brief Store current number of clusters for potential rollback
   *
   * @param new_K Number of clusters to cache
   */
  void set_old_K(int new_K) { old_K = new_K; }

  /**
   * @brief Restores the cached state into the data object
   */
  void restore_state() {
    data.restore_state(old_allocations, old_cluster_members, old_K);
  }

  /**
   * @brief Provides read-only access to the stored previous allocations.
   *
   * @return Const reference to the vector of previous allocations.
   */
  [[nodiscard]] const Eigen::VectorXi &old_allocations_view() const {
    return old_allocations;
  }

  /**
   * @brief Set index of first observation in split-merge pair
   *
   * @param i Index of first observation
   */
  void set_idx_i(int i) { idx_i = i; };

  /**
   * @brief Set index of second observation in split-merge pair
   *
   * @param j Index of second observation
   */
  void set_idx_j(int j) { idx_j = j; };

  // ========== Parameter Update Methods ==========

  /**
   * @brief Update process parameters during MCMC sampling
   *
   * @details This method allows derived classes to update their specific
   * parameters during the MCMC chain. For example, some processes may need to
   * sample hyperparameters or update auxiliary variables. Implementation is
   * process-specific.
   */
  virtual void update_params() = 0;

  /**
   * @brief Virtual destructor for proper cleanup of derived classes
   */
  virtual ~Process() {};
};
