##############################################################################
# In this script we simulate spatial data and run BNPClust on it. The simulation
# data are taken from the supplementary material of
# "Gianella, M., Quintana, F. A., and Guglielmi, A. (2026b).
# Consensus Monte Carlo for Large Spatial Dataset."
#
##############################################################################

source("R/utils.R")
source("R/utils_plot.R")

library(Rcpp)
library(RcppEigen)
library(sn)

dyn.load("build/libbnpclust_r.so")
bnp_mod <- Rcpp::Module("bnpclust_module", "libbnpclust_r")

## Set random seed for reproducibility
set.seed(44)

##############################################################################
# Data Creation ====
##############################################################################

I <- 36 # 36 areal locations on a regular 6x6 grid in the unit square
N_i <- 100 # 100 observations each
grid_dim <- 6 # sqrt(I), side of the regular grid

# Area centroids on a regular 6x6 grid in [0, 1] x [0, 1]
# (row/col kept for adjacency construction, x/y are the unit-square coordinates)
area_grid <- expand.grid(row = seq_len(grid_dim), col = seq_len(grid_dim))
area_grid$x <- (area_grid$col - 0.5) / grid_dim
area_grid$y <- (area_grid$row - 0.5) / grid_dim

# Split the grid into 4 quadrants (3x3 blocks) and assign the two
# generating distributions in a checkerboard pattern across quadrants:
# bottom-left & top-right quadrants -> Student's t
# bottom-right & top-left quadrants -> Skew-Normal
is_left <- area_grid$col <= grid_dim / 2
is_bottom <- area_grid$row <= grid_dim / 2
dist_type <- ifelse(is_bottom == is_left, "t", "sn")

data_list <- vector("list", I)

for (i in seq_len(I)) {
    if (dist_type[i] == "t") {
        # Student's t, 6 df, centred at 4, scale 1.5
        data_list[[i]] <- 4 + 1.5 * rt(N_i, df = 6)
    } else {
        # Skew-Normal, location xi = 4, scale omega = 1.3, shape alpha = -3
        data_list[[i]] <- sn::rsn(N_i, xi = 4, omega = 1.3, alpha = -3)
    }
}

# Long-format table: one row per observation, tagged with its area
sim_data <- do.call(
    rbind,
    lapply(seq_len(I), function(i) {
        data.frame(area = i, y = data_list[[i]])
    })
)

# Bookkeeping used later on for saving results / output naming
files_folder <- paste0("simulated/grid", grid_dim, "x", grid_dim)
file_chosen <- paste0("distance_grid", I, ".rds")

cat(
    "Simulated",
    I,
    "areas with",
    N_i,
    "observations each (N =",
    nrow(sim_data),
    ")\n"
)


##############################################################################
# kde distances matrix ====
##############################################################################

# Perform KDE for each area and store the full density objects
density_list <- lapply(data_list, function(y_i) {
    density(y_i, n = 512, kernel = "epanechnikov")
})

data_matrix <- matrix(0, nrow = I, ncol = I)

# Create progress bar
total_iterations <- I * (I + 1) / 2
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration <- 0

# Compute Jeffreys divergence (upper triangle + copy to lower for symmetry)
for (i in seq_along(density_list)) {
    for (k in i:length(density_list)) {
        data_matrix[i, k] <- compute_kde_distances(
            density_list[[i]],
            density_list[[k]],
            type = "Jeff"
        )

        # Copy to lower triangle for symmetry (skip diagonal)
        if (i != k) {
            data_matrix[k, i] <- data_matrix[i, k]
        }

        # Update progress bar
        iteration <- iteration + 1
        setTxtProgressBar(pb, iteration)
    }
}

# Close progress bar[cite: 1]
close(pb)

rownames(data_matrix) <- colnames(data_matrix) <- paste0("area", seq_len(I))

cat(
    "\nKDE-based Jeffreys Divergence matrix between areas computed (",
    I,
    "x",
    I,
    ")\n"
)

##############################################################################
# Spatial matrix ====
##############################################################################

# Construct the true graph G_true = {(1, 2), (3, 4), (5, 6)}
W <- matrix(0L, nrow = I, ncol = I)
W[1, 2] <- W[2, 1] <- 1L
W[3, 4] <- W[4, 3] <- 1L
W[5, 6] <- W[6, 5] <- 1L

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
    plot_distribution = FALSE
)

hyperparams$initial_clusters <- as.integer(hyperparams$initial_clusters - 1)

print(data_matrix)

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

print("Caching system instantiated")
data <- bnp_mod$create_Data(
    utils_param,
    hyperparams$initial_clusters
)

print("Data instantiated")

likelihood <- bnp_mod$create_Natarajan_likelihood(
    data,
    likelihood_param,
    utils_param
)
print("Likelihood instantiated")

# Instantiate U_sampler (RWMH) using factory function
u_sampler <- bnp_mod$create_RWMH(process_param, data, TRUE, 2.0, TRUE)

# Instantiate spatial modules
mod_spatial <- bnp_mod$create_SpatialModule(data, W, spatial_coefficient = 1)
print("Covariate modules instantiated")

# Instantiate Process (NGGPx) using factory function
process <- bnp_mod$create_NGGPx(
    data,
    process_param,
    u_sampler,
    list(mod_spatial)
)

print("Process instantiated")

# Instantiate Sampler (SplitMerge_LSS_SDDS) using factory function
sm <- bnp_mod$create_SplitMerge_LSS_SDDS(
    data,
    utils_param,
    likelihood,
    process,
    TRUE
)

# Instantiate Neal3 sampler using factory function
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
    if (i %% 2 == 0) {
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

mcmc_result <- list(
    allocations = allocations_out,
    K = K_out,
    U = U_out,
    elapsed_time = elapsed_time,
    BI = BI,
    NI = NI
)

##############################################################################
# 7. Save Results ====
##############################################################################

file_chosen_clean <- sub("\\.rds$", "", file_chosen)
folder_clean <- gsub("/", "_", files_folder)
data_tag <- paste0(folder_clean, "_", sub("^distance_", "", file_chosen_clean))

run_process <- "NGGPWx" # "DP" | "NGGP" | "NGGPW" | "NGGPWx"
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
folder <- save_with_name(utils_param, process_param, run_init, output_filename)

##############################################################################
# 8. Visualisation ====
##############################################################################

plot_folder <- paste0(folder, "/plots/")
dir.create(plot_folder, recursive = TRUE)

plot_post_distr(
    mcmc_result,
    BI = mcmc_result$BI,
    save = TRUE,
    folder = plot_folder
)
plot_trace_cls(
    mcmc_result,
    BI = mcmc_result$BI,
    save = TRUE,
    folder = plot_folder
)
plot_post_sim_matrix(
    mcmc_result,
    BI = mcmc_result$BI,
    save = TRUE,
    folder = plot_folder
)
plot_trace_U(
    mcmc_result,
    BI = mcmc_result$BI,
    save = TRUE,
    folder = plot_folder
)
plot_acf_U(mcmc_result, BI = mcmc_result$BI, save = TRUE, folder = plot_folder)
plot_cls_est(
    mcmc_result,
    BI = mcmc_result$BI,
    save = TRUE,
    folder = plot_folder
)

plot_inter_intra_histograms(
    hyperparams$initial_clusters,
    data_matrix,
    save = TRUE,
    folder = plot_folder
)
