/**
 * @file python_bindings.cpp
 * @brief Python language bindings for Bayesian nonparametric clustering models
 *
 * This file contains the Python/C++ interface code that exposes the clustering
 * algorithms (DP, NGGP, and their variants), data structures, likelihoods,
 * modules, and samplers directly to Python through nanobind.
 *
 * @author Filippo Galli
 * @date 2026
 */

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>

#include <memory>
#include <string>
#include <vector>

#include "params/Process-params.hpp"
#include "params/DP-params.hpp"
#include "params/NGGP-params.hpp"
#include "params/utils-params.hpp"
#include "params/Natarajan-params.hpp"
#include "params/GaussianMixtureModel-params.hpp"
#include "utils/Params.hpp"

#include "utils/ClusterInfo.hpp"
#include "processes/caches/continuos_cache.hpp"
#include "processes/caches/binary_cache.hpp"
#include "processes/caches/spatial_cache.hpp"

#include "utils/Data.hpp"
#include "utils/Datax.hpp"

#include "utils/Likelihood.hpp"
#include "likelihoods/Gamma_likelihood.hpp"
#include "likelihoods/Gamma_likelihood_distances.hpp"
#include "likelihoods/GaussianMixtureModel_likelihood.hpp"
#include "likelihoods/Natarajan_likelihood.hpp"
#include "likelihoods/Natarajan_likelihood_summaryStats.hpp"
#include "likelihoods/Null_likelihood.hpp"

#include "samplers/U_sampler/U_sampler.hpp"
#include "samplers/U_sampler/RWMH.hpp"
#include "samplers/U_sampler/MALA.hpp"

#include "utils/Module.hpp"
#include "processes/module/spatial_module.hpp"
#include "processes/module/spatial_module_continuous.hpp"
#include "processes/module/spatial_module_cache.hpp"
#include "processes/module/continuos_covariate_module.hpp"
#include "processes/module/continuos_covariate_module_cache.hpp"
#include "processes/module/binary_covariate_module.hpp"
#include "processes/module/binary_covariate_module_cache.hpp"
#include "processes/module/categorical_covariate_module.hpp"

#include "utils/Process.hpp"
#include "processes/DP.hpp"
#include "processes/NGGP.hpp"
#include "processes/DPx.hpp"
#include "processes/NGGPx.hpp"

#include "utils/Sampler.hpp"
#include "samplers/neal.hpp"
#include "samplers/neal_ZDNAM.hpp"
#include "samplers/splitmerge.hpp"
#include "samplers/splitmerge_SAMS.hpp"
#include "samplers/splitmerge_LSS.hpp"
#include "samplers/splitmerge_LSS_SDDS.hpp"

namespace nb = nanobind;

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

NB_MODULE(bnpclust, m) {
  m.doc() = "Python bindings for Bayesian Nonparametric Clustering (BNPClust)";

  // =========================================================================
  // 1. Parameter Classes
  // =========================================================================

  nb::class_<Process_params>(m, "Process_params",
                             "Base class for process hyperparameters.")
      .def_ro("name", &Process_params::name, "Name of the process.")
      .def("__repr__", [](const Process_params &p) {
        return "<Process_params: name='" + p.name + "'>";
      });

  nb::class_<DP_params, Process_params>(m, "DP_params",
                                        "Parameters for Dirichlet Process.")
      .def(nb::init<double>(), nb::arg("a"),
           "Construct DP parameters with total mass 'a'.")
      .def_rw("a", &DP_params::a, "Total mass parameter.")
      .def("__repr__", [](const DP_params &p) {
        return "<DP_params: a=" + std::to_string(p.a) + ">";
      });

  nb::class_<NGGP_params, Process_params>(
      m, "NGGP_params",
      "Parameters for Normalized Generalized Gamma Process.")
      .def(nb::init<double, double, double>(), nb::arg("a"), nb::arg("sigma"),
           nb::arg("tau"),
           "Construct NGGP parameters with total mass 'a', discount 'sigma', "
           "and scale 'tau'.")
      .def_rw("a", &NGGP_params::a, "Total mass parameter.")
      .def_rw("sigma", &NGGP_params::sigma, "Discount parameter.")
      .def_rw("tau", &NGGP_params::tau, "Scale parameter.")
      .def("__repr__", [](const NGGP_params &p) {
        return "<NGGP_params: a=" + std::to_string(p.a) +
               ", sigma=" + std::to_string(p.sigma) +
               ", tau=" + std::to_string(p.tau) + ">";
      });

  nb::class_<utils_params>(m, "utils_params",
                           "MCMC and dataset configuration parameters.")
      .def(nb::init<int, int, const Eigen::MatrixXd &>(), nb::arg("BI"),
           nb::arg("NI"), nb::arg("D"),
           "Construct with burn-in iterations BI, post-burn-in iterations NI, "
           "and data matrix D.")
      .def_rw("BI", &utils_params::BI, "Number of burn-in iterations.")
      .def_rw("NI", &utils_params::NI, "Number of sampling iterations.")
      .def_rw("D", &utils_params::D, "Data / distance matrix.")
      .def_ro("n", &utils_params::n, "Number of data points.")
      .def_ro("p", &utils_params::p, "Data dimension.")
      .def("__repr__", [](const utils_params &p) {
        return "<utils_params: BI=" + std::to_string(p.BI) +
               ", NI=" + std::to_string(p.NI) +
               ", n=" + std::to_string(p.n) + ", p=" + std::to_string(p.p) +
               ">";
      });

  nb::class_<Natarajan_params>(m, "Natarajan_params",
                               "Hyperparameters for Natarajan likelihood.")
      .def(nb::init<double, double, double, double, double, double>(),
           nb::arg("delta1"), nb::arg("alpha"), nb::arg("beta"),
           nb::arg("delta2"), nb::arg("gamma"), nb::arg("zeta"))
      .def_rw("delta1", &Natarajan_params::delta1)
      .def_rw("alpha", &Natarajan_params::alpha)
      .def_rw("beta", &Natarajan_params::beta)
      .def_rw("delta2", &Natarajan_params::delta2)
      .def_rw("gamma", &Natarajan_params::gamma)
      .def_rw("zeta", &Natarajan_params::zeta)
      .def("__repr__", [](const Natarajan_params &p) {
        return "<Natarajan_params: delta1=" + std::to_string(p.delta1) +
               ", alpha=" + std::to_string(p.alpha) +
               ", beta=" + std::to_string(p.beta) +
               ", delta2=" + std::to_string(p.delta2) +
               ", gamma=" + std::to_string(p.gamma) +
               ", zeta=" + std::to_string(p.zeta) + ">";
      });

  nb::class_<GaussianMixtureModel_params>(
      m, "GaussianMixtureModel_params",
      "Parameters for Gaussian Mixture Model conjugate prior.")
      .def(nb::init<double, double, double, double>(), nb::arg("m0"),
           nb::arg("kappa0"), nb::arg("alpha0"), nb::arg("beta0"))
      .def_rw("m0", &GaussianMixtureModel_params::m0)
      .def_rw("kappa0", &GaussianMixtureModel_params::kappa0)
      .def_rw("alpha0", &GaussianMixtureModel_params::alpha0)
      .def_rw("beta0", &GaussianMixtureModel_params::beta0)
      .def("__repr__", [](const GaussianMixtureModel_params &p) {
        return "<GaussianMixtureModel_params: m0=" + std::to_string(p.m0) +
               ", kappa0=" + std::to_string(p.kappa0) +
               ", alpha0=" + std::to_string(p.alpha0) +
               ", beta0=" + std::to_string(p.beta0) + ">";
      });

  nb::class_<Params>(m, "Params", "Combined parameter container.")
      .def(nb::init<double, double, double, double, double, double, int, int,
                    double, double, double, Eigen::MatrixXd>(),
           nb::arg("delta1") = 0.5, nb::arg("alpha") = 2.0,
           nb::arg("beta") = 2.0, nb::arg("delta2") = 2.0,
           nb::arg("gamma") = 2.0, nb::arg("zeta") = 2.0, nb::arg("BI") = 1000,
           nb::arg("NI") = 10000, nb::arg("a") = 1.0, nb::arg("sigma") = 1.0,
           nb::arg("tau") = 1.0, nb::arg("D") = Eigen::MatrixXd())
      .def_rw("delta1", &Params::delta1)
      .def_rw("alpha", &Params::alpha)
      .def_rw("beta", &Params::beta)
      .def_rw("delta2", &Params::delta2)
      .def_rw("gamma", &Params::gamma)
      .def_rw("zeta", &Params::zeta)
      .def_rw("BI", &Params::BI)
      .def_rw("NI", &Params::NI)
      .def_rw("a", &Params::a)
      .def_rw("sigma", &Params::sigma)
      .def_rw("tau", &Params::tau)
      .def_rw("D", &Params::D)
      .def_ro("n", &Params::n);

  // =========================================================================
  // 2. Cluster Information and Caches
  // =========================================================================

  nb::class_<ClusterInfo>(m, "ClusterInfo",
                          "Abstract base interface for cluster statistics caches.")
      .def("set_allocation", &ClusterInfo::set_allocation, nb::arg("index"),
           nb::arg("cluster"), nb::arg("old_cluster"))
      .def("recompute", &ClusterInfo::recompute, nb::arg("K"),
           nb::arg("allocations"))
      .def("move_cluster_info", &ClusterInfo::move_cluster_info,
           nb::arg("from_cluster"), nb::arg("to_cluster"))
      .def("remove_info", &ClusterInfo::remove_info, nb::arg("cluster"));

  nb::class_<ContinuosCache, ClusterInfo>(
      m, "ContinuosCache", "Cache for continuous covariate cluster statistics.")
      .def("__init__",
           [](ContinuosCache *self, const Eigen::VectorXi &initial_allocations,
              const Eigen::VectorXd &continuos_covariates) {
             Eigen::VectorXi alloc_copy = initial_allocations;
             new (self) ContinuosCache(alloc_copy, continuos_covariates);
           },
           nb::arg("initial_allocations"), nb::arg("continuos_covariates"))
      .def_ro("continuos_covariates", &ContinuosCache::continuos_covariates);

  nb::class_<BinaryCache, ClusterInfo>(
      m, "BinaryCache", "Cache for binary covariate cluster statistics.")
      .def(nb::init<const Eigen::VectorXi &, const Eigen::VectorXi &>(),
           nb::arg("initial_allocations"), nb::arg("binary_covariates"))
      .def_ro("binary_covariates", &BinaryCache::binary_covariates);

  nb::class_<SpatialCache, ClusterInfo>(
      m, "SpatialCache", "Cache for spatial adjacency cluster statistics.")
      .def(nb::init<const Eigen::VectorXi &, const Eigen::MatrixXi &>(),
           nb::arg("initial_allocations"), nb::arg("W"))
      .def_ro("W", &SpatialCache::W);

  // =========================================================================
  // 3. Data Structures
  // =========================================================================

  nb::class_<Data>(m, "Data",
                   "Point distance and cluster allocation manager.")
      .def(nb::init<const utils_params &, const Eigen::VectorXi &>(),
           nb::arg("params"),
           nb::arg("initial_allocations") = Eigen::VectorXi(),
           nb::keep_alive<1, 2>())
      .def("get_data", &Data::get_data, nb::arg("i"), nb::arg("j") = 0)
      .def("get_n", &Data::get_n, "Total number of observations.")
      .def("get_K", &Data::get_K, "Current number of clusters.")
      .def("get_p", &Data::get_p, "Data dimension.")
      .def("get_allocations",
           [](const Data &self) -> Eigen::VectorXi {
             return self.get_allocations();
           },
           "Get vector of cluster assignments.")
      .def("get_cluster_size", &Data::get_cluster_size,
           nb::arg("cluster_index"), "Size of cluster.")
      .def("get_cluster_assignment", &Data::get_cluster_assignment,
           nb::arg("index"), "Cluster assigned to observation.")
      .def("get_cluster_assignments", &Data::get_cluster_assignments,
           nb::arg("cluster"), "Observation indices in cluster.")
      .def("get_cluster_map", &Data::get_cluster_map_copy,
           "Copy of cluster-to-indices membership map.")
      .def("set_allocation", &Data::set_allocation, nb::arg("index"),
           nb::arg("cluster"), "Reassign observation to cluster.")
      .def("set_allocations", &Data::set_allocations,
           nb::arg("new_allocations"), "Set all cluster allocations.")
      .def("__repr__", [](const Data &d) {
        return "<Data: n=" + std::to_string(d.get_n()) +
               ", K=" + std::to_string(d.get_K()) + ">";
      });

  nb::class_<Datax, Data>(m, "Datax",
                          "Extended Data container integrating covariate caches.")
      .def(nb::init<const utils_params &,
                    std::vector<std::shared_ptr<ClusterInfo>>,
                    const Eigen::VectorXi &>(),
           nb::arg("params"), nb::arg("modules_list"),
           nb::arg("initial_allocations") = Eigen::VectorXi(),
           nb::keep_alive<1, 2>())
      .def("set_allocation", &Datax::set_allocation, nb::arg("index"),
           nb::arg("cluster"))
      .def("set_allocations", &Datax::set_allocations,
           nb::arg("new_allocations"))
      .def("__repr__", [](const Datax &d) {
        return "<Datax: n=" + std::to_string(d.get_n()) +
               ", K=" + std::to_string(d.get_K()) + ">";
      });

  // =========================================================================
  // 4. Likelihood Classes
  // =========================================================================

  nb::class_<Likelihood>(m, "Likelihood",
                         "Abstract base class for cluster likelihood evaluation.")
      .def("cluster_loglikelihood",
           nb::overload_cast<int>(&Likelihood::cluster_loglikelihood,
                                  nb::const_),
           nb::arg("cluster_index"))
      .def("cluster_loglikelihood",
           nb::overload_cast<int, const Eigen::Ref<const Eigen::VectorXi> &>(
               &Likelihood::cluster_loglikelihood, nb::const_),
           nb::arg("cluster_index"), nb::arg("cls_ass_k"))
      .def("clusters_loglikelihood", &Likelihood::clusters_loglikelihood,
           nb::arg("c1"), nb::arg("c2"))
      .def("point_loglikelihood_cond", &Likelihood::point_loglikelihood_cond,
           nb::arg("point_index"), nb::arg("cluster_index"));

  nb::class_<Natarajan_likelihood, Likelihood>(
      m, "Natarajan_likelihood", "Natarajan cohesion and repulsion likelihood.")
      .def(nb::init<const Data &, const Natarajan_params &,
                    const utils_params &>(),
           nb::arg("data"), nb::arg("params"), nb::arg("utils"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>());

  nb::class_<Natarajan_likelihood_summaryStats, Likelihood>(
      m, "Natarajan_likelihood_summaryStats",
      "Natarajan likelihood computed with summary statistics.")
      .def(nb::init<const Data &, const Natarajan_params &,
                    const utils_params &>(),
           nb::arg("data"), nb::arg("params"), nb::arg("utils"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>());

  nb::class_<Null_likelihood, Likelihood>(
      m, "Null_likelihood", "Null likelihood ignoring data contributions.")
      .def(nb::init<const Data &>(), nb::arg("data"), nb::keep_alive<1, 2>());

  nb::class_<Gamma_likelihood, Likelihood>(
      m, "Gamma_likelihood", "Gamma cohesion likelihood on distances.")
      .def(nb::init<const Data &, const Natarajan_params &,
                    const utils_params &>(),
           nb::arg("data"), nb::arg("params"), nb::arg("utils"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>());

  nb::class_<Gamma_likelihood_distances, Likelihood>(
      m, "Gamma_likelihood_distances",
      "Gamma likelihood on pairwise distance matrices.")
      .def(nb::init<const Data &, const Natarajan_params &,
                    const utils_params &>(),
           nb::arg("data"), nb::arg("params"), nb::arg("utils"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>());

  nb::class_<GaussianMixtureModel_likelihood, Likelihood>(
      m, "GaussianMixtureModel_likelihood",
      "Gaussian Mixture Model conjugate likelihood.")
      .def(nb::init<const Data &, const GaussianMixtureModel_params &>(),
           nb::arg("data"), nb::arg("params"), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>());

  // =========================================================================
  // 5. U Samplers
  // =========================================================================

  nb::class_<U_sampler>(m, "U_sampler",
                        "Abstract base class for NGGP latent variable U sampler.")
      .def("update_U", &U_sampler::update_U, "Sample one MCMC step for U.")
      .def("get_U", &U_sampler::get_U, "Current value of latent variable U.")
      .def("get_acceptance_rate", &U_sampler::get_acceptance_rate,
           "Acceptance rate of U updates.");

  nb::class_<RWMH, U_sampler>(
      m, "RWMH", "Random Walk Metropolis-Hastings sampler for U.")
      .def(nb::init<NGGP_params &, Data &, bool, double, bool>(),
           nb::arg("params"), nb::arg("data"), nb::arg("use_V") = true,
           nb::arg("proposal_sd") = 1.0, nb::arg("tuning_enabled") = true,
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>());

  nb::class_<MALA, U_sampler>(
      m, "MALA", "Metropolis-Adjusted Langevin Algorithm sampler for U.")
      .def(nb::init<NGGP_params &, Data &, bool, double, bool>(),
           nb::arg("params"), nb::arg("data"), nb::arg("use_V") = true,
           nb::arg("proposal_sd") = 1.0, nb::arg("tuning_enabled") = true,
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>());

  // =========================================================================
  // 6. Covariate and Spatial Modules
  // =========================================================================

  nb::class_<Module>(m, "Module",
                     "Abstract base class for covariate/spatial similarity modules.")
      .def("compute_similarity_cls", &Module::compute_similarity_cls,
           nb::arg("cls_idx"), nb::arg("old_allo") = false)
      .def("compute_similarity_obs",
           nb::overload_cast<int, int>(&Module::compute_similarity_obs,
                                       nb::const_),
           nb::arg("obs_idx"), nb::arg("cls_idx"))
      .def("compute_similarity_obs",
           nb::overload_cast<int>(&Module::compute_similarity_obs, nb::const_),
           nb::arg("obs_idx"));

  nb::class_<SpatialModule, Module>(m, "SpatialModule",
                                    "Binary adjacency spatial module.")
      .def(nb::init<const Data &, const Eigen::MatrixXi, double>(),
           nb::arg("data"), nb::arg("W"), nb::arg("spatial_coefficient") = 1.0,
           nb::keep_alive<1, 2>());

  nb::class_<SpatialModuleContinuous, Module>(
      m, "SpatialModuleContinuous", "Continuous distance spatial module.")
      .def(nb::init<const Data &, const Eigen::MatrixXd, double>(),
           nb::arg("data"), nb::arg("W"), nb::arg("spatial_coefficient") = 1.0,
           nb::keep_alive<1, 2>());

  nb::class_<SpatialModuleCache, Module>(
      m, "SpatialModuleCache", "Cached spatial adjacency module.")
      .def(nb::init<const Data &, SpatialCache &, double>(), nb::arg("data"),
           nb::arg("spatial_cache"), nb::arg("spatial_coefficient") = 1.0,
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>());

  nb::class_<ContinuosCovariatesModule, Module>(
      m, "ContinuosCovariatesModule", "Continuous covariate module.")
      .def(nb::init<const Data &, const Eigen::VectorXd, bool, double, double,
                    double, double, double>(),
           nb::arg("data"), nb::arg("covariates"), nb::arg("fixed_v"),
           nb::arg("m") = 0.0, nb::arg("B") = 1.0, nb::arg("v") = 1.0,
           nb::arg("nu") = 1.0, nb::arg("S0") = 1.0, nb::keep_alive<1, 2>());

  nb::class_<ContinuosCovariatesModuleCache, Module>(
      m, "ContinuosCovariatesModuleCache",
      "Cached continuous covariate module.")
      .def(nb::init<const Data &, const ContinuosCache &, bool, double, double,
                    double, double, double>(),
           nb::arg("data"), nb::arg("continuos_cache"), nb::arg("fixed_v"),
           nb::arg("m") = 0.0, nb::arg("B") = 1.0, nb::arg("v") = 1.0,
           nb::arg("nu") = 1.0, nb::arg("S0") = 1.0, nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>());

  nb::class_<BinaryCovariatesModule, Module>(m, "BinaryCovariatesModule",
                                             "Binary covariate module.")
      .def(nb::init<const Data &, const Eigen::VectorXi, double, double>(),
           nb::arg("data"), nb::arg("covariates"), nb::arg("prior_a") = 1.0,
           nb::arg("prior_b") = 1.0, nb::keep_alive<1, 2>());

  nb::class_<BinaryCovariatesModuleCache, Module>(
      m, "BinaryCovariatesModuleCache", "Cached binary covariate module.")
      .def(nb::init<const Data &, const BinaryCache &, double, double>(),
           nb::arg("data"), nb::arg("cache"), nb::arg("beta_prior_alpha") = 1.0,
           nb::arg("beta_prior_beta") = 1.0, nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>());

  nb::class_<CategoricalCovariatesModule, Module>(
      m, "CategoricalCovariatesModule", "Categorical covariate module.")
      .def(nb::init<const Data &, const Eigen::VectorXi, std::vector<double>>(),
           nb::arg("data"), nb::arg("categorical_covariate"),
           nb::arg("prior_alpha"), nb::keep_alive<1, 2>());

  // =========================================================================
  // 7. Nonparametric Processes
  // =========================================================================

  nb::class_<Process>(m, "Process",
                      "Abstract base class for nonparametric processes.")
      .def("update_params", &Process::update_params,
           "Update process parameters during MCMC.");

  nb::class_<DP, Process>(m, "DP", "Dirichlet Process prior.")
      .def(nb::init<Data &, DP_params &>(), nb::arg("data"), nb::arg("params"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>());

  nb::class_<NGGP, Process>(
      m, "NGGP", "Normalized Generalized Gamma Process prior.")
      .def(nb::init<Data &, NGGP_params &, U_sampler &>(), nb::arg("data"),
           nb::arg("params"), nb::arg("u_sampler"), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>());

  nb::class_<DPx, DP>(m, "DPx", "Dirichlet Process prior with covariate modules.")
      .def(nb::init<Data &, DP_params &,
                    const std::vector<std::shared_ptr<Module>> &>(),
           nb::arg("data"), nb::arg("params"), nb::arg("modules"),
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>());

  nb::class_<NGGPx, NGGP>(
      m, "NGGPx",
      "Normalized Generalized Gamma Process prior with covariate modules.")
      .def(nb::init<Data &, NGGP_params &, U_sampler &,
                    const std::vector<std::shared_ptr<Module>> &>(),
           nb::arg("data"), nb::arg("params"), nb::arg("u_sampler"),
           nb::arg("modules"), nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>());

  // =========================================================================
  // 8. MCMC Samplers
  // =========================================================================

  nb::class_<Sampler>(m, "Sampler",
                      "Abstract base class for MCMC sampling algorithms.")
      .def("step", &Sampler::step, "Execute one iteration of the MCMC step.");

  nb::class_<Neal3, Sampler>(m, "Neal3", "Neal's Algorithm 3 collapsed Gibbs sampler.")
      .def(nb::init<Data &, Likelihood &, Process &>(), nb::arg("data"),
           nb::arg("likelihood"), nb::arg("process"), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>());

  nb::class_<Neal3ZDNAM, Sampler>(m, "Neal3ZDNAM",
                                  "Neal's Algorithm 3 ZDNAM variant sampler.")
      .def(nb::init<Data &, Likelihood &, Process &>(), nb::arg("data"),
           nb::arg("likelihood"), nb::arg("process"), nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>());

  nb::class_<SplitMerge, Sampler>(m, "SplitMerge", "Split-Merge MCMC sampler.")
      .def(nb::init<Data &, Likelihood &, Process &, bool>(), nb::arg("data"),
           nb::arg("likelihood"), nb::arg("process"),
           nb::arg("shuffle") = false, nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>())
      .def("get_accepted_split", &SplitMerge::get_accepted_split)
      .def("get_accepted_merge", &SplitMerge::get_accepted_merge)
      .def("get_accepted_shuffle", &SplitMerge::get_accepted_shuffle);

  nb::class_<SplitMerge_SAMS, Sampler>(
      m, "SplitMerge_SAMS",
      "Sequentially Allocated Merge-Split (SAMS) sampler.")
      .def(nb::init<Data &, Likelihood &, Process &, bool>(), nb::arg("data"),
           nb::arg("likelihood"), nb::arg("process"),
           nb::arg("shuffle") = false, nb::keep_alive<1, 2>(),
           nb::keep_alive<1, 3>(), nb::keep_alive<1, 4>())
      .def("get_accepted_split", &SplitMerge_SAMS::get_accepted_split)
      .def("get_accepted_merge", &SplitMerge_SAMS::get_accepted_merge)
      .def("get_accepted_shuffle", &SplitMerge_SAMS::get_accepted_shuffle);

  nb::class_<SplitMerge_LSS, Sampler>(
      m, "SplitMerge_LSS",
      "Locality Sensitive Sampling (LSS) Split-Merge sampler.")
      .def(nb::init<Data &, utils_params &, Likelihood &, Process &,
                    bool>(),
           nb::arg("data"), nb::arg("params"), nb::arg("likelihood"),
           nb::arg("process"), nb::arg("shuffle") = false,
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>(), nb::keep_alive<1, 5>())
      .def("get_accepted_split", &SplitMerge_LSS::get_accepted_split)
      .def("get_accepted_merge", &SplitMerge_LSS::get_accepted_merge)
      .def("get_accepted_shuffle", &SplitMerge_LSS::get_accepted_shuffle);

  nb::class_<SplitMerge_LSS_SDDS, Sampler>(
      m, "SplitMerge_LSS_SDDS",
      "Locality Sensitive Sampling with SDDS Split-Merge sampler.")
      .def(nb::init<Data &, utils_params &, Likelihood &, Process &,
                    bool>(),
           nb::arg("data"), nb::arg("params"), nb::arg("likelihood"),
           nb::arg("process"), nb::arg("shuffle") = false,
           nb::keep_alive<1, 2>(), nb::keep_alive<1, 3>(),
           nb::keep_alive<1, 4>(), nb::keep_alive<1, 5>())
      .def("get_accepted_split", &SplitMerge_LSS_SDDS::get_accepted_split)
      .def("get_accepted_merge", &SplitMerge_LSS_SDDS::get_accepted_merge)
      .def("get_accepted_shuffle", &SplitMerge_LSS_SDDS::get_accepted_shuffle)
      .def("accepted_moves", [](const SplitMerge_LSS_SDDS &self) {
        nb::print(nb::str(("Ratio accepted splits: " +
                           std::to_string(self.get_accepted_split() * 100.0) +
                           " %").c_str()));
        nb::print(nb::str(("Ratio accepted merges: " +
                           std::to_string(self.get_accepted_merge() * 100.0) +
                           " %").c_str()));
        nb::print(nb::str(("Ratio accepted shuffles: " +
                           std::to_string(self.get_accepted_shuffle() * 100.0) +
                           " %").c_str()));
      }, "Print summary of accepted moves ratios.");
}
