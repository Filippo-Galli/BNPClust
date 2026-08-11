#pragma once

#include "../params/GaussianMixtureModel-params.hpp"
#include "../utils/Likelihood.hpp"

class GaussianMixtureModel_likelihood : public Likelihood {
private:
  const GaussianMixtureModel_params &params;

  const double log_pi;
  const double lgamma_alpha0;

  struct ClusterStats {
    int n = 0;
    double sum = 0.0;
    double sumsq = 0.0;
  };

  ClusterStats
  compute_stats(const Eigen::Ref<const Eigen::VectorXi> &cluster_areas) const;
  ClusterStats compute_area_stats(int area_index) const;
  double log_marginal_likelihood(const ClusterStats &stats) const;
  double log_predictive_from_stats(const ClusterStats &stats,
                                   const ClusterStats &area_stats) const;

public:
  GaussianMixtureModel_likelihood(const Data &data,
                                  const GaussianMixtureModel_params &params);
  ~GaussianMixtureModel_likelihood() = default;

  double cluster_loglikelihood(int cluster_index) const override;
  double cluster_loglikelihood(
      int cluster_index,
      const Eigen::Ref<const Eigen::VectorXi> &cluster_areas) const override;
  double point_loglikelihood_cond(int point_index,
                                  int cluster_index) const override;
};
