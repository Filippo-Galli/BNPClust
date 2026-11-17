/**
 * @file splitmerge_SAMS.hpp
 * @brief Sequential Allocation Merge-Split (SAMS) sampler implementation
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "splitmerge_base.hpp"

/**
 * @brief Sequential Allocation Merge-Split (SAMS) sampler for Bayesian
 * nonparametric models (Dahl and Newcomb, 2022)
 *
 * Inherits from SplitMergeBase and uses sequential allocation for proposal generation.
 */
class SplitMerge_SAMS : public SplitMergeBase {
protected:
  void choose_indices_for_step() override {
    choose_indeces();
  }

public:
  SplitMerge_SAMS(Data &d, Likelihood &l, Process &pr, bool shuffle)
      : SplitMergeBase(d, l, pr, shuffle) {}

  virtual ~SplitMerge_SAMS() = default;

  void choose_indeces();

  void split_move() override;
  void merge_move() override;
};
