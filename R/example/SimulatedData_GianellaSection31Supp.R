##############################################################################
# In this script we simulate spatial data and run BNPClust on it. The simulation
# data are taken from the supplementary material of
# "Gianella, M., Quintana, F. A., and Guglielmi, A. (2026b).
# Consensus Monte Carlo for Large Spatial Dataset."
#
# NOTE: the original supplementary example uses I = 9 areas (3 x 3 grid).
# Here we scale the same construction up to I = 36 areas (6 x 6 grid).
##############################################################################

source("R/utils.R")
source("R/utils_plot.R")

library(Rcpp)
library(RcppEigen)

dyn.load("build/libbnpclust_r.so")
bnp_mod <- Rcpp::Module("bnpclust_module", "libbnpclust_r")

## Set random seed for reproducibility
set.seed(44)

##############################################################################
# Data Creation ====
##############################################################################

I_side <- 6
I <- I_side^2
N_i <- 100

# Centres of each area (unit square split into I_side x I_side equal cells)
area_grid <- expand.grid(
    row = 1:I_side,
    col = 1:I_side
)
area_grid$x <- (area_grid$col - 0.5) / I_side
area_grid$y <- (area_grid$row - 0.5) / I_side

x_centre <- mean(range(area_grid$x)) # grid centre coordinate (x_bar)
y_centre <- mean(range(area_grid$y)) # grid centre coordinate (y_bar)

# alr-transformed weights, eq. (8)
w_tilde_1 <- 3 * (area_grid$x - x_centre) + 3 * (area_grid$y - y_centre)
w_tilde_2 <- -3 * (area_grid$x - x_centre) - 3 * (area_grid$y - y_centre)

# Inverse-alr (logistic) transform, eq. (5), with H = 3 mixture components
denom <- 1 + exp(w_tilde_1) + exp(w_tilde_2)
w1 <- exp(w_tilde_1) / denom
w2 <- exp(w_tilde_2) / denom
w3 <- 1 / denom

area_weights <- cbind(w1, w2, w3)
stopifnot(all(abs(rowSums(area_weights) - 1) < 1e-10))

mixture_means <- c(-5, 0, 5)
mixture_sd <- 1

simulate_area <- function(w, n, means, sd) {
    comp <- sample.int(length(means), size = n, replace = TRUE, prob = w)
    rnorm(n, mean = means[comp], sd = sd)
}

data_list <- lapply(seq_len(I), function(i) {
    simulate_area(area_weights[i, ], N_i, mixture_means, mixture_sd)
})

# Long-format table: one row per observation, tagged with its area
sim_data <- do.call(
    rbind,
    lapply(seq_len(I), function(i) {
        data.frame(area = i, y = data_list[[i]])
    })
)

# Bookkeeping used later on for saving results / output naming
files_folder <- paste0("simulated/grid", I)
file_chosen <- paste0("distance_grid", I_side, "x", I_side, ".rds")

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

W <- matrix(0L, nrow = I, ncol = I)
for (i in seq_len(I)) {
    for (k in seq_len(I)) {
        if (i == k) {
            next
        }
        same_row <- area_grid$row[i] == area_grid$row[k]
        same_col <- area_grid$col[i] == area_grid$col[k]
        adjacent_col <- abs(area_grid$col[i] - area_grid$col[k]) == 1
        adjacent_row <- abs(area_grid$row[i] - area_grid$row[k]) == 1
        if ((same_row && adjacent_col) || (same_col && adjacent_row)) {
            W[i, k] <- 1L
        }
    }
}

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
    k_elbow = 5
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

    # Neal3 Step (run Gibbs every 25 moves of Split-Merge)
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
    folder = plot_folder,
    x_lim = c(0, 7.5),
    y_lim = c(0, 1.5)
)
