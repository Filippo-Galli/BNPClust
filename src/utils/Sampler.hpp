/**
 * @file Sampler.hpp
 * @brief Abstract base class for MCMC sampling algorithms
 *
 * This file defines the abstract Sampler base class that provides a unified
 * interface for different MCMC sampling strategies (Neal3, SplitMerge, SAMS).
 * It establishes the common framework for Bayesian nonparametric clustering
 * with different algorithmic approaches.
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "Data.hpp"
#include "Likelihood.hpp"
#include "Process.hpp"

#include <Eigen/Dense>
#include <random>

/**
 * @brief Abstract base class for MCMC sampler implementations
 *
 * This class provides the foundation for implementing different MCMC sampling
 * algorithms for Bayesian nonparametric models including Dirichlet Process
 * (DP), Normalized Generalized Gamma Process (NGGP), and their weighted
 * variants (DPW, NGGPW). It holds references to all necessary components and
 * provides a unified interface for different sampling strategies.
 *
 * @details The class supports various MCMC algorithms through derived
 * implementations:
 * - **Gibbs Samplers**: Standard collapsed Gibbs sampling (Neal's Algorithm 3)
 * - **Split-Merge Samplers**: Advanced algorithms for better mixing (SAMS
 * variants)
 * - **Hybrid Approaches**: Combinations of different sampling strategies
 *
 * Each sampler operates on the same core components but implements different
 * strategies for exploring the posterior distribution of cluster assignments.
 * The composition pattern allows flexible combinations of different data types,
 * likelihood models, process priors, and parameter configurations.
 *
 * @see Data, Likelihood, Process
 * @see Neal3, SplitMerge, SplitMerge_SAMS
 */
class Sampler {
protected:
    // ========== Core Components ==========

    /** @brief Reference to the data object containing observations and current
     * allocations */
    Data &data;

    /** @brief Reference to the likelihood computation object for evaluating
     * cluster assignments */
    const Likelihood &likelihood;

    /** @brief Reference to the stochastic process object (DP, NGGP, DPW, NGGPW)
     */
    Process &process;

    // ========== Random Number Generation ==========

    /** @brief Random device for generating random numbers across sampling
     * algorithms */
    std::random_device rd;

public:
    // ========== Constructor ==========

    /**
     * @brief Constructor initializing sampler with required components
     *
     * @param d Reference to Data object containing observations and current
     * cluster assignments
     * @param l Reference to Likelihood object for evaluating assignment
     * probabilities
     * @param pr Reference to Process object (DP, NGGP, DPW, or NGGPW) defining
     * the prior
     *
     * @details The constructor establishes references to all components needed
     * for MCMC sampling. The specific behavior depends on the concrete
     * implementations of each component:
     * - Different Process types (DP vs NGGP) yield different clustering behaviors
     * - Different Likelihood models handle various data types and distributions
     * - Parameter settings control burn-in, iterations, and hyperparameter values
     */
    Sampler(Data &d, const Likelihood &l, Process &pr) : data(d), likelihood(l), process(pr) {};

    // ========== Pure Virtual Interface ==========

    /**
     * @brief Pure virtual method to perform one MCMC sampling step
     *
     * This method must be implemented by derived classes to define their specific
     * sampling algorithm. Each call should update the current state of the Markov
     * chain by sampling from the appropriate conditional distributions.
     *
     * @details Different sampler implementations use different strategies:
     * - **Gibbs samplers** (e.g., Neal3): Sequential sampling of individual
     * assignments
     * - **Split-Merge samplers**: Joint updates of multiple assignments via
     * split/merge moves
     * - **Hybrid samplers**: Combinations of multiple update mechanisms
     *
     * The method should update cluster assignments in the Data object and may
     * also update auxiliary variables or hyperparameters as needed by the
     * specific algorithm.
     *
     * @note This is a pure virtual function that must be overridden by concrete
     * sampler implementations
     *
     * @see Neal3::step(), SplitMerge::step(), SplitMerge_SAMS::step()
     */
    virtual void step() = 0;

    /**
     * @brief Virtual destructor for proper cleanup of derived classes
     */
    virtual ~Sampler() = default;
};
