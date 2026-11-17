/**
 * @file splitmerge_LSS_SDDS.hpp
 * @brief Locality Sensitive Sampling (LSS) with SDDS Split-Merge sampler implementation
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../params/utils-params.hpp"
#include "splitmerge_base.hpp"

/**
 * @brief Locality Sensitive Sampling (LSS) with SDDS Split-Merge sampler
 *
 * Inherits from SplitMergeBase and implements the Smart-split, Dumb-merge,
 * Dumb-split, Smart-merge adaptive strategy with Cauchy kernel distance weighting.
 */
class SplitMerge_LSS_SDDS : public SplitMergeBase {
private:
  const utils_params &params;

  std::vector<double> row_max_distances;
  double temperature;

  const std::vector<int> &original_allocations_ref;

  std::vector<double> probs;

  const double rand_split_prob = std::log(0.5);

  int split_moves = 0;
  int merge_moves = 0;
  int shuffle_moves = 0;

protected:
  void choose_indices_for_step() override {}
  void split_move() override {}
  void merge_move() override {}

  void rollback() override {
    process.restore_state();
  }

public:
  SplitMerge_LSS_SDDS(Data &d, utils_params &p, Likelihood &l, Process &pr,
                      bool shuffle);

  virtual ~SplitMerge_LSS_SDDS() = default;

  void choose_indeces(bool similarity = false);

  void smart_split_move();
  void dumb_split_move();
  void smart_merge_move();
  void dumb_merge_move();

  void shuffle() override;
  void step() override;

  [[nodiscard]] double get_accepted_split() const noexcept {
    return split_moves > 0 ? static_cast<double>(accepted_split) / split_moves : 0.0;
  }

  [[nodiscard]] double get_accepted_merge() const noexcept {
    return merge_moves > 0 ? static_cast<double>(accepted_merge) / merge_moves : 0.0;
  }

  [[nodiscard]] double get_accepted_shuffle() const noexcept {
    return shuffle_moves > 0 ? static_cast<double>(accepted_shuffle) / shuffle_moves : 0.0;
  }
};
