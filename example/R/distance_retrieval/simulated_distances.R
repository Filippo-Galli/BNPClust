# Import utils.R for additional functions
source("example/R/utils.R")

# Generate data from 4 Gaussian distributions
set.seed(44)

### ----------------------- Natarajan Data Generation ----------------- ###
sigma <- 0.25
d <- 10
# sigma <- 0.2
# d <- 50
# sigma <- 0.18
# d <- 10
data_generation <- generate_mixture_data(N = 100, sigma = sigma, d = d)

all_data <- data_generation$points
ground_truth <- data_generation$clusts

# Distance plot
# distance_plot(all_data, ground_truth)

### ----------------------- Gaussian Data Generation ----------------- ###

# n_points <- 30
# centers <- list(
#   c(10, 10),  # Cluster 1: top-right
#   c(6, 6),  # Cluster 2: bottom-right
#   c(0, 0),  # Cluster 3: bottom-left
#   c(3, 3)   # Cluster 4: top-left
# )

# # Generate data for each cluster
# data1 <- mvrnorm(n_points, mu = centers[[1]], Sigma = 3*diag(2))
# data2 <- mvrnorm(n_points, mu = centers[[2]], Sigma = 3*diag(2))
# data3 <- mvrnorm(n_points, mu = centers[[3]], Sigma = 3*diag(2))
# data4 <- mvrnorm(n_points, mu = centers[[4]], Sigma = 3*diag(2))

# # Combine all data points + labels
# ground_truth <- c(
#   rep(0, n_points),
#   rep(1, n_points),
#   rep(2, n_points),
#   rep(3, n_points)
# )
# all_data <- rbind(data1, data2, data3, data4)

### ----------------------- Gamma Data Generation ----------------- ###
# n_points <- 30
# # Gamma Generated data
# data1 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) + 8, ncol = 2)
# data2 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) + 1, ncol = 2)
# data3 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) - 1, ncol = 2)
# data4 <- matrix(rgamma(n_points * 2, shape = 0.5, rate = 1) + 4, ncol = 2)

# # Combine all data points + labels
# ground_truth <- c(
#   rep(0, n_points),
#   rep(1, n_points),
#   rep(2, n_points),
#   rep(3, n_points)
# )
# all_data <- rbind(data1, data2, data3, data4)

# Create distance matrix (n x n)
dist_matrix <- as.matrix(dist(all_data))

# Check folder existence and create if not existing
folder <- "simulation_data"
if (!dir.exists(folder)) {
  dir.create(folder)
}

# folder of the specific simulation
subfolder <- paste0("Natarajan_", sigma, "sigma_", d, "d")
if (!dir.exists(file.path(folder, subfolder))) {
  dir.create(file.path(folder, subfolder))
}

# Save distance matrix, data and ground truth in the folder RDS
saveRDS(dist_matrix, file.path(folder, subfolder, "dist_matrix.rds"))
saveRDS(all_data, file.path(folder, subfolder, "all_data.rds"))
saveRDS(ground_truth, file.path(folder, subfolder, "ground_truth.rds"))
