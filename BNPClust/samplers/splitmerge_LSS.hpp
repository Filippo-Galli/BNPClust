/**
 * @file splitmerge_LSS.hpp
 * @brief Locality Sensitive Sampling (LSS) Split-Merge sampler implementation
 *
 * @author Filippo Galli
 * @date 2025
 */

#pragma once

#include "../params/utils-params.hpp"
#include "splitmerge_SAMS.hpp"

/**
 * @brief Locality Sensitive Sampling (LSS) Split-Merge sampler for Bayesian
 * nonparametric models (Luo & Shrivastava, 2018)
 *
 * Inherits proposal generation and move execution from SplitMerge_SAMS,
 * overriding anchor point selection with distance-weighted locality sensitive sampling.
 */
class SplitMerge_LSS : public SplitMerge_SAMS {
protected:
  const utils_params &params;
  std::vector<double> probs;

  void choose_indices_for_step() override;

public:
  SplitMerge_LSS(Data &d, utils_params &p, Likelihood &l, Process &pr,
                 bool shuffle)
      : SplitMerge_SAMS(d, l, pr, shuffle), params(p), probs(d.get_n()) {}

  virtual ~SplitMerge_LSS() = default;

  void choose_indeces(bool similarity = false);
};
