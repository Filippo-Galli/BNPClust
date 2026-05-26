source("R/utils.R")
source("R/utils_plot.R")
source("R/mcmc_loop.R")

## Set random seed for reproducibility
set.seed(44)

##############################################################################
# Data Loading ====
##############################################################################

## Load simulated data
# sigma <- .18
# d <- 10
# file_chosen <- paste0("Natarajan_", sigma, "sigma_", d, "d")
# files_folder <- paste0("simulation_data/", file_chosen)
# all_data <- readRDS(file = paste0(files_folder, "/all_data.rds"))
# ground_truth <- readRDS(file = paste0(files_folder, "/ground_truth.rds"))
# dist_matrix <- readRDS(file = paste0(files_folder, "/dist_matrix.rds"))

## Load real data
files_folder <- "real_data/LA"
files <- list.files(files_folder)
file_chosen <- files[3]
dist_matrix <- readRDS(file = paste0(files_folder, "/", file_chosen))
puma_age_data <- readRDS(file = paste0(files_folder, "/puma_age_stats.rds"))
puma_sex_data <- readRDS(file = paste0(files_folder, "/puma_sex_stats.rds"))
# plot_distance(dist_matrix)

if (min(dist_matrix) < 0) {
    dist_matrix <- dist_matrix + abs(min(dist_matrix))
}
diag(dist_matrix) <- 0

##############################################################################
# Spatial Adjacency Matrix ====
##############################################################################

## Retrieve spatial adjacency matrix W from distance matrix
# W <- retrieve_W(dist_matrix)
W <- readRDS(file = paste0(files_folder, "/adj_matrix.rds"))

# Check is W is symmetric
if (!isSymmetric(W)) {
    warning("W is not symmetric!")
}

##############################################################################
# C++ Integration ====
##############################################################################

## C++ implementation is compiled in R/mcmc_loop.R

##############################################################################
# Hyperparameter Configuration ====
##############################################################################

# Plot k-means elbow method to help set hyperparameters
# plot_k_means(dist_matrix, max_k = 10)

# Set hyperparameters based on distance matrix and save it for future use
# hyperparams <- set_hyperparameters(dist_matrix,
#     k_elbow = 3, plot_clustering = FALSE, plot_distribution = FALSE
# )
# saveRDS(hyperparams, file = paste0(files_folder, "/hyperparameters_", sub("\\.rds$", "", file_chosen), ".rds"))
hyperparams <- readRDS(file = paste0(files_folder, "/hyperparameters_", sub("\\.rds$", "", file_chosen), ".rds"))

##############################################################################
# Parameter Object Initialization ====
##############################################################################

# Create Params using factory function instead of module constructor
process_param <- create_NGGP_params(
    1, # a
    0.1, # sigma
    1 # tau
)

utils_param <- create_utils_params(5000, 15000, dist_matrix)
likelihood_param <- create_Natarajan_params(
    hyperparams$delta1,
    hyperparams$alpha,
    hyperparams$beta,
    hyperparams$delta2,
    hyperparams$gamma,
    hyperparams$zeta
)


##############################################################################
# Covariates Object Initialization ====
##############################################################################

# Ensure W is integer matrix
# W <- matrix(as.integer(W), nrow = nrow(W), ncol = ncol(W))

# Ensure W is double matrix
W <- matrix(as.double(W), nrow = nrow(W), ncol = ncol(W))
continuos_covariates <- as.numeric(puma_age_data$AGEP_std_mean)
binary_covariates <- as.integer(puma_sex_data$SEX_mode)

##############################################################################
# Initial Cluster Allocation ====
##############################################################################

hyperparams$initial_clusters <- hyperparams$initial_clusters - 1

## Display initial cluster allocation summary
print("Initial cluster allocation:")
print(table(hyperparams$initial_clusters))

##############################################################################
# MCMC Execution ====
##############################################################################

mcmc_result <- run_mcmc(likelihood_param, process_param, utils_param, hyperparams$initial_clusters, W, continuos_covariates, binary_covariates)
elapsed_time <- mcmc_result$elapsed_time

##############################################################################
# Save Results (Optional) ====
##############################################################################
file_chosen <- sub("\\.rds$", "", file_chosen)
files_folder_clean <- gsub("/", "_", files_folder)
data_type <- paste0(files_folder_clean, "_", sub("^distance_", "", file_chosen))
process <- "NGGPWX" # Process type: "DP", "DPW", "NGGP", "NGGPW", NGGPWx
method <- "LSS_SDDS25+Gibbs1" # MCMC method used
initialization <- "kmeans" # Initialization strategy
options <- ""
filename <- paste0(data_type, "_", process, "_", method, "_", initialization, "_", options)
save_with_name(utils_param, process_param, initialization, filename)

##############################################################################
# Visualization (Optional) ====
##############################################################################

# plot_post_distr(mcmc_result, BI = mcmc_result$BI)
# plot_trace_cls(mcmc_result, BI = mcmc_result$BI)
# plot_post_sim_matrix(mcmc_result, BI = mcmc_result$BI)
# plot_trace_U(mcmc_result, BI = mcmc_result$BI)
# plot_acf_U(mcmc_result, BI = mcmc_result$BI)
# plot_cls_est(mcmc_result, BI = mcmc_result$BI)
# plot_stats(mcmc_result, ground_truth, BI = mcmc_result$BI)

# puma_ids <- sf::st_read("input/LA/counties-pumas/counties-pumas.shp", quiet = TRUE)[["PUMA"]]
# plot_map_prior_mean(unit_ids = puma_ids, puma_dir = "input/CA/counties-pumas")
# plot_map_cls(
#   results = mcmc_result,
#   BI = mcmc_result$BI,
#   unit_ids = puma_ids,
#   puma_dir = "input/LA/counties-pumas"
# )

# plot_hist_cls(
#   results = mcmc_result,
#   BI = mcmc_result$BI,
#   # input_dir = "input/CA/",
# )
