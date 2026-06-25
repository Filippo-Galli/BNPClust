##############################################################################
# BNPClust — Self-Contained Example
#
# Runs a full Bayesian nonparametric clustering chain using BNPClust's
# C++ backend via Rcpp modules. The script covers:
#   1. Library / module loading
#   2. Data loading (real or simulated)
#   3. Hyperparameter configuration (Natarajan likelihood)
#   4. Parameter object initialisation
#   5. MCMC execution (SplitMerge LSS-SDDS + Neal-3 scan)
#   6. Result saving and visualisation
#
# All alternative models / likelihood / process choices are kept as
# commented-out blocks so the file doubles as a quick reference card.
##############################################################################

source("R/utils.R")
source("R/utils_plot.R")
source("R/mcmc_loop.R")

library(Rcpp)
library(RcppEigen)

dyn.load("build/libbnpclust_r.so")
bnp_mod <- Rcpp::Module("bnpclust_module", "libbnpclust_r")

## Set random seed for reproducibility
set.seed(44)

##############################################################################
# Data Loading ====
##############################################################################

## Load real data
files_folder <- "real_data/LA"
files <- list.files(files_folder)
file_chosen <- files[2]
data_matrix <- readRDS(file = paste0(files_folder, "/", file_chosen))

##############################################################################
# Spatial Adjacency Matrix ====
##############################################################################

W <- readRDS(file = paste0(files_folder, "/adj_matrix.rds"))

# Check is W is symmetric
if (!isSymmetric(W)) {
    warning("W is not symmetric!")
}

W <- matrix(as.integer(W), nrow = nrow(W), ncol = ncol(W))

##############################################################################
# Hyperparameter Configuration ====
##############################################################################

# Set hyperparameters based on distance matrix and save it for future use
hyperparams <- set_hyperparameters(
    data_matrix,
    k_elbow = 3,
    plot_clustering = FALSE,
    plot_distribution = FALSE
)

hyperparams$initial_clusters <- as.integer(hyperparams$initial_clusters - 1)

##############################################################################
# Parameter Object Initialization ====
##############################################################################

process_param <- bnp_mod$create_NGGP_params(
    1, # a
    0.1, # sigma
    1 # tau
)

utils_param <- bnp_mod$create_utils_params(5000, 15000, data_matrix)
likelihood_param <- bnp_mod$create_Natarajan_params(
    hyperparams$delta1,
    hyperparams$alpha,
    hyperparams$beta,
    hyperparams$delta2,
    hyperparams$gamma,
    hyperparams$zeta
)

##############################################################################
# Initial Cluster Allocation ====
##############################################################################

print("Initial cluster allocation:")
print(table(hyperparams$initial_clusters))

mcmc_result <- run_mcmc(
    likelihood_param,
    process_param,
    utils_param,
    hyperparams$initial_clusters,
    W,
    # continuos_covariates,
    # binary_covariates
)

# Ensure types are correct for C++
initial_allocations <- as.integer(integer(0))

# continuos_cache <- bnp_mod$create_Continuos_cache(
#     initial_allocations,
#     continuos_covariates
# )
# binary_cache <- bnp_mod$create_Binary_cache(
#     initial_allocations,
#     binary_covariates
# )
# spatial_cache <- create_Spatial_cache(initial_allocations, W)
print("Caching system instantiated")

# Instantiate Data using factory function
# data <- create_Datax(
#     utils_param,
#     list(binary_cache, continuos_cache),
#     initial_allocations
# )
data <- bnp_mod$create_Data(utils_param, initial_allocations)

print("Data instantiated")

# Instantiate Likelihood using factory function
# likelihood <- bnp_mod$create_GaussianMixtureModel_likelihood(
#     data,
#     likelihood_param
# )
likelihood <- bnp_mod$create_Natarajan_likelihood(
    data,
    likelihood_param,
    utils_param
)
# likelihood <- create_Natarajan_likelihood_summaryStats(data, params)
# likelihood <- create_Null_likelihood(data, params) # Placeholder likelihood
# likelihood <- create_Gamma_likelihood(data, params)

print("Likelihood instantiated")

# Instantiate U_sampler (RWMH) using factory function
# Constructor: Params&, Data&, bool use_V, double proposal_sd, bool tuning_enabled
u_sampler <- bnp_mod$create_RWMH(process_param, data, TRUE, 2.0, TRUE)
# u_sampler <- NULL

# Instantiate Process (NGGPx) using modules
# 1. Spatial module
#mod_spatial <- create_SpatialModule(data, W, spatial_coefficient = 0.1)
# mod_spatial <- create_SpatialModuleContinuous(data, W, spatial_coefficient = 1.0)
# mod_spatial <- create_SpatialModuleCache(data, cache = spatial_cache, spatial_coefficient = 1.0)

# 2. Covariate module (cached)
# fixed_v <- TRUE
# B <- 10 * var(continuos_covariates) # prior variance
# m <- 0 # prior mean
# v <- 0.5 * var(continuos_covariates) # known variance
# nu <- 1
# S0 <- 1.0

# mod_cont <- bnp_mod$create_ContinuosCovariatesModuleCache(
#     data,
#     continuos_cache,
#     fixed_v,
#     m,
#     B,
#     v,
#     nu,
#     S0
# )
# mod_cont <- create_ContinuosCovariatesModule(data, continuos_covariates, fixed_v = TRUE, m = m, B = B, v = v)

# 3. Binary covariate module
# mod_binary <- create_BinaryCovariatesModule(data, binary_covariates, 0.1, 0.1)
# mod_binary <- bnp_mod$create_BinaryCovariatesModuleCache(
#     data,
#     binary_cache,
#     0.1,
#     0.1
# )

# 4. Categorical covariate module
# alphas <- rep(1.0, length(unique(categorical_covariates)))
# mod_categorical <- create_CategoricalCovariatesModule(data, categorical_covariates, alphas)

print("Covariate modules instantiated")

# Combine modules into NGGPx process
# process <- create_NGGPx(data, process_param, u_sampler, list(mod_spatial, mod_cont, mod_binary))
process <- bnp_mod$create_NGGP(data, process_param, u_sampler)
# process <- bnp_mod$create_DP(
#     data,
#     process_param
# )

print("Process instantiated")

# Instantiate Sampler (SplitMerge_LSS_SDDS) using factory function
sm <- bnp_mod$create_SplitMerge_LSS_SDDS(
    data,
    utils_param,
    likelihood,
    process
)

neal3 <- bnp_mod$create_Neal3(data, likelihood, process)

print("Sampler instantiated")

# Get parameters for loop using getter functions
BI <- bnp_mod$params_get_BI(utils_param)
NI <- bnp_mod$params_get_NI(utils_param)
total_iters <- BI + NI

# Results storage
allocations_out <- vector("list", total_iters)
K_out <- integer(total_iters)
U_out <- rep(NA_real_, total_iters)

cat("Starting MCMC with", NI, "iterations after", BI, "burn-in...\n")

start_time <- Sys.time()

for (i in 1:total_iters) {
    # Update process parameters (U)
    bnp_mod$process_update_params(process)

    # MCMC Step
    bnp_mod$sampler_step(sm)

    # Neal3 Step
    if (i %% 25 == 0) {
        bnp_mod$sampler_step(neal3)
    }
    # Store results
    allocations_out[[i]] <- bnp_mod$data_get_allocations(data)
    K_out[i] <- bnp_mod$data_get_K(data)
    if (!is.null(u_sampler)) {
        U_out[i] <- bnp_mod$u_sampler_get_U(u_sampler)
    }

    # Progress
    if (i %% max(1, floor(total_iters / 20)) == 0) {
        elapsed <- as.numeric(difftime(
            Sys.time(),
            start_time,
            units = "secs"
        ))
        iter_per_sec <- i / elapsed
        eta <- (total_iters - i) / iter_per_sec
        cat(sprintf(
            "Iteration %d: Clusters: %d - iter/s: %.2f eta: %.2f\n ",
            i,
            bnp_mod$data_get_K(data),
            iter_per_sec,
            eta
        ))
    }
}

elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("MCMC completed.\n")
cat("Total time (secs):", elapsed_time, "\n")
if (!is.null(u_sampler)) {
    cat(
        "U acceptance rate:",
        bnp_mod$u_sampler_get_acceptance_rate(u_sampler) * 100,
        "%\n"
    )
}
if (exists("sampler")) {
    bnp_mod$lss_sdds_accepted_moves(sampler)
}


##############################################################################
# 7. Save Results ====
##############################################################################

file_chosen_clean <- sub("\\.rds$", "", file_chosen)
folder_clean <- gsub("/", "_", files_folder)
data_tag <- paste0(folder_clean, "_", sub("^distance_", "", file_chosen_clean))

run_process <- "NGGP" # "DP" | "NGGP" | "NGGPW" | "NGGPWx"
run_method <- "LSS_SDDS25+Gibbs1"
run_init <- "kmeans"
run_label <- "example"

output_filename <- paste(
    data_tag,
    run_process,
    run_method,
    run_init,
    run_label,
    sep = "_"
)
save_with_name(utils_param, process_param, run_init, output_filename)


##############################################################################
# 8. Visualisation ====
##############################################################################

plot_post_distr(mcmc_result, BI = mcmc_result$BI)
plot_trace_cls(mcmc_result, BI = mcmc_result$BI)
plot_post_sim_matrix(mcmc_result, BI = mcmc_result$BI)
plot_trace_U(mcmc_result, BI = mcmc_result$BI)
plot_acf_U(mcmc_result, BI = mcmc_result$BI)
plot_cls_est(mcmc_result, BI = mcmc_result$BI)
