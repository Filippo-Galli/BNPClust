/**
 * @file bindings.cpp
 * @brief R language bindings for Bayesian nonparametric clustering models
 *
 * This file contains the R/C++ interface code that exposes the clustering
 * algorithms (DP, NGGP, and their variants) to R through Rcpp. It handles
 * object construction, method calls, and data conversion between R and C++
 * representations.
 *
 * @author Filippo Galli
 * @date 2025
 */

#include "Rcpp/XPtr.h"
#include "params/DP-params.hpp"
#include "params/GaussianMixtureModel-params.hpp"
#include "params/NGGP-params.hpp"
#include "params/Natarajan-params.hpp"
#include "params/utils-params.hpp"
#include "utils/Data.hpp"
#include "utils/Datax.hpp"
#include "utils/Params.hpp"
#include <RcppEigen.h>
#include <memory>

#include "likelihoods/Gamma_likelihood.hpp"
#include "likelihoods/GaussianMixtureModel_likelihood.hpp"
#include "likelihoods/Natarajan_likelihood.hpp"
#include "likelihoods/Natarajan_likelihood_summaryStats.hpp"
#include "likelihoods/Null_likelihood.hpp"
#include "utils/Likelihood.hpp"

#include "processes/DP.hpp"
#include "processes/DPx.hpp"
#include "processes/NGGP.hpp"
#include "processes/NGGPx.hpp"
#include "utils/Process.hpp"

#include "processes/module/binary_covariate_module.hpp"
#include "processes/module/binary_covariate_module_cache.hpp"
#include "processes/module/categorical_covariate_module.hpp"
#include "processes/module/continuos_covariate_module.hpp"
#include "processes/module/continuos_covariate_module_cache.hpp"
#include "processes/module/spatial_module.hpp"
#include "processes/module/spatial_module_cache.hpp"
#include "processes/module/spatial_module_continuous.hpp"

#include "processes/caches/binary_cache.hpp"
#include "processes/caches/continuos_cache.hpp"
#include "processes/caches/spatial_cache.hpp"
#include "utils/ClusterInfo.hpp"

#include "samplers/U_sampler/MALA.hpp"
#include "samplers/U_sampler/RWMH.hpp"
#include "samplers/U_sampler/U_sampler.hpp"

#include "samplers/neal.hpp"
#include "samplers/neal_ZDNAM.hpp"
#include "samplers/splitmerge.hpp"
#include "samplers/splitmerge_LSS.hpp"
#include "samplers/splitmerge_LSS_SDDS.hpp"
#include "samplers/splitmerge_SAMS.hpp"
#include "utils/Sampler.hpp"

// Helper functions to handle inheritance with XPtrs
Data *get_data_ptr(SEXP sexp) {
    // Try Data
    try {
        Rcpp::XPtr<Data> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    // Try Data_wClusterInfo
    try {
        Rcpp::XPtr<Datax> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    Rcpp::stop("Expected external pointer to Data or Data_wClusterInfo");
    return nullptr;
}

ClusterInfo *get_cluster_info_ptr(SEXP sexp) {
    // Try ClusterInfo
    try {
        Rcpp::XPtr<ClusterInfo> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    // Try Covariate_cache
    try {
        Rcpp::XPtr<ContinuosCache> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    // Try BinaryContinuosCache
    try {
        Rcpp::XPtr<BinaryCache> ptr(sexp);
        return ptr.get();
    } catch (...) {
    }

    Rcpp::stop("Expected external pointer to ClusterInfo or Covariate_cache");
    return nullptr;
}

// [[Rcpp::depends(RcppEigen)]]

// Factory functions for base classes
// [[Rcpp::export]]
Rcpp::XPtr<Params> create_Params(double delta1, double alpha, double beta, double delta2, double gamma, double zeta,
                                 int BI, int NI, double a, double sigma, double tau, Eigen::MatrixXd D) {
    return Rcpp::XPtr<Params>(new Params(delta1, alpha, beta, delta2, gamma, zeta, BI, NI, a, sigma, tau, D), true);
}

// Factory functions for base classes
// [[Rcpp::export]]
Rcpp::XPtr<utils_params> create_utils_params(int BI, int NI, Eigen::MatrixXd D) {
    return Rcpp::XPtr<utils_params>(new utils_params(BI, NI, D), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGP_params> create_NGGP_params(double a, double sigma, double tau) {
    return Rcpp::XPtr<NGGP_params>(new NGGP_params(a, sigma, tau), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<DP_params> create_DP_params(double a) { return Rcpp::XPtr<DP_params>(new DP_params(a), true); }

// [[Rcpp::export]]
Rcpp::XPtr<Natarajan_params> create_Natarajan_params(double delta1, double alpha, double beta, double delta2,
                                                     double gamma, double zeta) {
    return Rcpp::XPtr<Natarajan_params>(new Natarajan_params(delta1, alpha, beta, delta2, gamma, zeta), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<GaussianMixtureModel_params> create_GaussianMixtureModel_params(Eigen::VectorXd y, double m0, double kappa0,
                                                                           double alpha0, double beta0) {
    return Rcpp::XPtr<GaussianMixtureModel_params>(new GaussianMixtureModel_params(y, m0, kappa0, alpha0, beta0), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Data> create_Data(Rcpp::XPtr<utils_params> params, Eigen::VectorXi initial_allocations) {
    return Rcpp::XPtr<Data>(new Data(*params, initial_allocations), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<ContinuosCache> create_Continuos_cache(Eigen::VectorXi &initial_allocations,
                                                  Eigen::VectorXd continuos_covariates) {
    return Rcpp::XPtr<ContinuosCache>(new ContinuosCache(initial_allocations, continuos_covariates), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<BinaryCache> create_Binary_cache(Eigen::VectorXi &initial_allocations, Eigen::VectorXi binary_covariates) {
    return Rcpp::XPtr<BinaryCache>(new BinaryCache(initial_allocations, binary_covariates), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SpatialCache> create_Spatial_cache(Eigen::VectorXi &initial_allocations, Eigen::MatrixXi W) {
    return Rcpp::XPtr<SpatialCache>(new SpatialCache(initial_allocations, W), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Datax> create_Datax(Rcpp::XPtr<utils_params> params, Rcpp::List modules_list,
                               Eigen::VectorXi initial_allocations) {

    std::vector<std::shared_ptr<ClusterInfo>> modules;
    modules.reserve(modules_list.size());

    for (int i = 0; i < modules_list.size(); ++i) {
        // Get the raw ClusterInfo pointer from XPtr (could be Covariate_cache,
        // etc.)
        ClusterInfo *raw_ptr = get_cluster_info_ptr(modules_list[i]);
        // Wrap in non-owning shared_ptr (R's XPtr manages lifetime)
        modules.push_back(std::shared_ptr<ClusterInfo>(raw_ptr, [](ClusterInfo *) {}));
    }
    return Rcpp::XPtr<Datax>(new Datax(*params, std::move(modules), initial_allocations), true);
}

// Factory functions for likelihoods
// [[Rcpp::export]]
Rcpp::XPtr<Natarajan_likelihood> create_Natarajan_likelihood(SEXP data_sexp, Rcpp::XPtr<Natarajan_params> params,
                                                             Rcpp::XPtr<utils_params> utils) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Natarajan_likelihood>(new Natarajan_likelihood(*data, *params, *utils), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Natarajan_likelihood_summaryStats>
create_Natarajan_likelihood_summaryStats(SEXP data_sexp, Rcpp::XPtr<Natarajan_params> params,
                                         Rcpp::XPtr<utils_params> utils) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Natarajan_likelihood_summaryStats>(new Natarajan_likelihood_summaryStats(*data, *params, *utils),
                                                         true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Null_likelihood> create_Null_likelihood(SEXP data_sexp) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Null_likelihood>(new Null_likelihood(*data), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Gamma_likelihood> create_Gamma_likelihood(SEXP data_sexp, Rcpp::XPtr<Natarajan_params> params,
                                                     Rcpp::XPtr<utils_params> utils) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Gamma_likelihood>(new Gamma_likelihood(*data, *params, *utils), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<GaussianMixtureModel_likelihood>
create_GaussianMixtureModel_likelihood(SEXP data_sexp, Rcpp::XPtr<GaussianMixtureModel_params> params) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<GaussianMixtureModel_likelihood>(new GaussianMixtureModel_likelihood(*data, *params), true);
}

// Factory functions for U samplers
// [[Rcpp::export]]
Rcpp::XPtr<RWMH> create_RWMH(Rcpp::XPtr<NGGP_params> params, SEXP data_sexp, bool use_V, double proposal_sd,
                             bool tuning_enabled) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<RWMH>(new RWMH(*params, *data, use_V, proposal_sd, tuning_enabled), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<MALA> create_MALA(Rcpp::XPtr<NGGP_params> params, SEXP data_sexp, bool use_V, double proposal_sd,
                             bool tuning_enabled) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<MALA>(new MALA(*params, *data, use_V, proposal_sd, tuning_enabled), true);
}

// Factory functions for processes
// [[Rcpp::export]]
Rcpp::XPtr<DP> create_DP(SEXP data_sexp, Rcpp::XPtr<DP_params> params) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<DP>(new DP(*data, *params), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGP> create_NGGP(SEXP data_sexp, Rcpp::XPtr<NGGP_params> params, Rcpp::XPtr<U_sampler> u_sampler) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<NGGP>(new NGGP(*data, *params, *u_sampler), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_SpatialModule(SEXP data_sexp, Eigen::MatrixXi W,
                                                         double spatial_coefficient) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<SpatialModule>(*data, W, spatial_coefficient);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_SpatialModuleContinuous(SEXP data_sexp, Eigen::MatrixXd W,
                                                                   double spatial_coefficient) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<SpatialModuleContinuous>(*data, W, spatial_coefficient);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_SpatialModuleCache(SEXP data_sexp, Rcpp::XPtr<SpatialCache> cache,
                                                              double spatial_coefficient) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<SpatialModuleCache>(*data, *cache, spatial_coefficient, nullptr);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_ContinuosCovariatesModule(SEXP data_sexp, Eigen::VectorXd covariates,
                                                                     bool fixed_v, double m = 0, double B = 1,
                                                                     double v = 1, double nu = 1, double S0 = 1) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<ContinuosCovariatesModule>(*data, covariates, fixed_v, m, B, v, nu, S0);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_ContinuosCovariatesModuleCache(SEXP data_sexp,
                                                                          Rcpp::XPtr<ContinuosCache> cache,
                                                                          bool fixed_v, double m = 0, double B = 1,
                                                                          double v = 1, double nu = 1, double S0 = 1) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<ContinuosCovariatesModuleCache>(*data, *cache, fixed_v, m, B, v, nu, S0);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_BinaryCovariatesModule(SEXP data_sexp, Eigen::VectorXi covariates,
                                                                  double prior_a = 1.0, double prior_b = 1.0) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<BinaryCovariatesModule>(*data, covariates, prior_a, prior_b);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::shared_ptr<Module>> create_BinaryCovariatesModuleCache(SEXP data_sexp, Rcpp::XPtr<BinaryCache> cache,
                                                                       double beta_prior_alpha = 1.0,
                                                                       double beta_prior_beta = 1.0) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<BinaryCovariatesModuleCache>(*data, *cache, beta_prior_alpha, beta_prior_beta);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

Rcpp::XPtr<std::shared_ptr<Module>> create_CategoricalCovariatesModule(SEXP data_sexp,
                                                                       Eigen::VectorXi categorical_covariate,
                                                                       std::vector<double> prior_alpha) {
    Data *data = get_data_ptr(data_sexp);
    auto ptr = std::make_shared<CategoricalCovariatesModule>(*data, categorical_covariate, prior_alpha);
    return Rcpp::XPtr<std::shared_ptr<Module>>(new std::shared_ptr<Module>(ptr), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<DPx> create_DPx(SEXP data_sexp, Rcpp::XPtr<DP_params> params, Rcpp::List modules_list) {

    Data *data = get_data_ptr(data_sexp);

    std::vector<std::shared_ptr<Module>> modules;
    modules.reserve(modules_list.size());

    for (int i = 0; i < modules_list.size(); ++i) {
        // Cast the SEXP to an XPtr<shared_ptr<Module>>
        Rcpp::XPtr<std::shared_ptr<Module>> mod_ptr_wrapper(modules_list[i]);
        modules.push_back(*mod_ptr_wrapper); // Copy the shared_ptr, incrementing ref count
    }

    return Rcpp::XPtr<DPx>(new DPx(*data, *params, modules), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<NGGPx> create_NGGPx(SEXP data_sexp, Rcpp::XPtr<NGGP_params> params, Rcpp::XPtr<U_sampler> u_sampler,
                               Rcpp::List modules_list) {
    Data *data = get_data_ptr(data_sexp);

    std::vector<std::shared_ptr<Module>> modules;
    modules.reserve(modules_list.size());

    for (int i = 0; i < modules_list.size(); ++i) {
        // Cast the SEXP to an XPtr<shared_ptr<Module>>
        Rcpp::XPtr<std::shared_ptr<Module>> mod_ptr_wrapper(modules_list[i]);
        modules.push_back(*mod_ptr_wrapper); // Copy the shared_ptr, incrementing ref count
    }

    return Rcpp::XPtr<NGGPx>(new NGGPx(*data, *params, *u_sampler, modules), true);
}

// Factory functions for samplers
// [[Rcpp::export]]
Rcpp::XPtr<Neal3> create_Neal3(SEXP data_sexp, Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Neal3>(new Neal3(*data, *likelihood, *process), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<Neal3ZDNAM> create_Neal3ZDNAM(SEXP data_sexp, Rcpp::XPtr<Likelihood> likelihood,
                                         Rcpp::XPtr<Process> process) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<Neal3ZDNAM>(new Neal3ZDNAM(*data, *likelihood, *process), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge> create_SplitMerge(SEXP data_sexp, Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process,
                                         bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge>(new SplitMerge(*data, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_SAMS> create_SplitMerge_SAMS(SEXP data_sexp, Rcpp::XPtr<Likelihood> likelihood,
                                                   Rcpp::XPtr<Process> process, bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge_SAMS>(new SplitMerge_SAMS(*data, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_LSS> create_SplitMerge_LSS(SEXP data_sexp, Rcpp::XPtr<utils_params> params,
                                                 Rcpp::XPtr<Likelihood> likelihood, Rcpp::XPtr<Process> process,
                                                 bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge_LSS>(new SplitMerge_LSS(*data, *params, *likelihood, *process, shuffle), true);
}

// [[Rcpp::export]]
Rcpp::XPtr<SplitMerge_LSS_SDDS> create_SplitMerge_LSS_SDDS(SEXP data_sexp, Rcpp::XPtr<utils_params> params,
                                                           Rcpp::XPtr<Likelihood> likelihood,
                                                           Rcpp::XPtr<Process> process, bool shuffle) {
    Data *data = get_data_ptr(data_sexp);
    return Rcpp::XPtr<SplitMerge_LSS_SDDS>(new SplitMerge_LSS_SDDS(*data, *params, *likelihood, *process, shuffle),
                                           true);
}

// Wrapper functions for methods
// [[Rcpp::export]]
void process_update_params(Rcpp::XPtr<Process> process) { process->update_params(); }

// [[Rcpp::export]]
void sampler_step(Rcpp::XPtr<Sampler> sampler) { sampler->step(); }

// [[Rcpp::export]]
Eigen::VectorXi data_get_allocations(SEXP data_sexp) {
    Data *data = get_data_ptr(data_sexp);
    return data->get_allocations();
}

// [[Rcpp::export]]
int data_get_K(SEXP data_sexp) {
    Data *data = get_data_ptr(data_sexp);
    return data->get_K();
}

// [[Rcpp::export]]
double u_sampler_get_U(Rcpp::XPtr<U_sampler> u_sampler) { return u_sampler->get_U(); }

// [[Rcpp::export]]
double u_sampler_get_acceptance_rate(Rcpp::XPtr<U_sampler> u_sampler) { return u_sampler->get_acceptance_rate(); }

// [[Rcpp::export]]
int params_get_BI(Rcpp::XPtr<utils_params> params) { return params->BI; }

// [[Rcpp::export]]
int params_get_NI(Rcpp::XPtr<utils_params> params) { return params->NI; }

// [[Rcpp::export]]
double NGGP_params_get_a(Rcpp::XPtr<NGGP_params> params) { return params->a; }

// [[Rcpp::export]]
double params_get_sigma(Rcpp::XPtr<NGGP_params> params) { return params->sigma; }

// [[Rcpp::export]]
double params_get_tau(Rcpp::XPtr<NGGP_params> params) { return params->tau; }

// [[Rcpp::export]]
std::string params_get_name(Rcpp::XPtr<Process_params> params) { return params->name; }

// [[Rcpp::export]]
double DP_params_get_a(Rcpp::XPtr<DP_params> params) { return params->a; }

// [[Rcpp::export]]
void cluster_info_set_allocation(Rcpp::XPtr<ClusterInfo> cluster_info, int index, int cluster, int old_cluster) {
    cluster_info->set_allocation(index, cluster, old_cluster);
}

// [[Rcpp::export]]
void cluster_info_recompute(Rcpp::XPtr<ClusterInfo> cluster_info, int K, Eigen::VectorXi allocations) {
    cluster_info->recompute(K, allocations);
}

// [[Rcpp::export]]
void datax_set_allocation(Rcpp::XPtr<Datax> data, int index, int cluster) { data->set_allocation(index, cluster); }

// [[Rcpp::export]]
void datax_set_allocations(Rcpp::XPtr<Datax> data, Eigen::VectorXi new_allocations) {
    data->set_allocations(new_allocations);
}

// [[Rcpp::export]]
void lss_sdds_accepted_moves(Rcpp::XPtr<SplitMerge_LSS_SDDS> sampler) {

    // Display the counts
    Rcpp::Rcout << "Ratio accepted splits: " << sampler->get_accepted_split() * 100 << " %" << std::endl;
    Rcpp::Rcout << "Ratio accepted merges: " << sampler->get_accepted_merge() * 100 << " %" << std::endl;
    Rcpp::Rcout << "Ratio accepted shuffles: " << sampler->get_accepted_shuffle() * 100 << " %" << std::endl;
}
