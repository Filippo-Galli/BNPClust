/**
 * @file splitmerge.hpp
 * @brief Split-Merge MCMC sampler implementation for Bayesian nonparametric models
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "splitmerge_base.hpp"

/**
 * @brief Split-Merge sampler for Bayesian nonparametric mixture models (Jain and Neal, 2004)
 *
 * Inherits from SplitMergeBase and uses restricted Gibbs sampling for proposal generation.
 */
class SplitMerge : public SplitMergeBase {
protected:
  void choose_indices_for_step() override {
    choose_indeces();
  }

public:
  SplitMerge(Data &d, Likelihood &l, Process &pr, bool shuffle)
      : SplitMergeBase(d, l, pr, shuffle) {}

  virtual ~SplitMerge() = default;

  void choose_indeces();

  void restricted_gibbs(int iterations, bool only_probabilities = false);

  void split_move() override;
  void merge_move() override;
  void shuffle() override;
};
