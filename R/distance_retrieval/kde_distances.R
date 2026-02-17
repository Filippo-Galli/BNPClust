source("R/utils.R")
source("R/utils_plot.R")

##############################################################################
# Load .csv files from 'input/' directory ====
##############################################################################
location <- "LA"
input_folder <- paste0("input/", location, "/")
data <- read.csv(paste0(input_folder, "full_dataset.csv"))
W <- as.matrix(readRDS(paste0(input_folder, "adj_matrix.rds")))

# Split data by PUMA
data <- split(data$log_income, data$COD_PUMA)

##############################################################################
# Create histogram for each PUMA ====
##############################################################################
density_list <- list()

for (i in seq_along(data)) {
  pumas_data <- data[[i]]
  # Perform KDE
  density_list[[i]] <- density(pumas_data, n = 512, kernel = "epanechnikov")
}

##############################################################################
# Distances between histograms ====
##############################################################################
distance_jeff_divergences <- matrix(0,
  nrow = length(density_list),
  ncol = length(density_list)
)
distance_wasserstein <- matrix(0,
  nrow = length(density_list),
  ncol = length(density_list)
)

# Create progress bar
# Only compute upper triangle (including diagonal) since matrices are symmetric
n <- length(density_list)
total_iterations <- n * (n + 1) / 2
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration <- 0

for (i in seq_along(density_list)) {
  for (j in i:length(density_list)) {
    distance_jeff_divergences[i, j] <-
      compute_kde_distances(density_list[[i]],
        density_list[[j]],
        type = "Jeff"
      )
    distance_wasserstein[i, j] <-
      compute_kde_distances(density_list[[i]],
        density_list[[j]],
        type = "Wasserstein"
      )

    # Copy to lower triangle for symmetry (skip diagonal)
    if (i != j) {
      distance_jeff_divergences[j, i] <- distance_jeff_divergences[i, j]
      distance_wasserstein[j, i] <- distance_wasserstein[i, j]
    }

    # Update progress bar
    iteration <- iteration + 1
    setTxtProgressBar(pb, iteration)
  }
}

# Close progress bar
close(pb)

##############################################################################
# Plot Distance ====
##############################################################################
plot_distance(distance_jeff_divergences,
  title = "Jeffreys Divergence Distance", save = TRUE,
  folder = paste0("old_results/distance_plots/", location, "/")
)
plot_distance(distance_wasserstein,
  title = "Wasserstein Distance", save = TRUE,
  folder = paste0("old_results/distance_plots/", location, "/")
)

##############################################################################
# Save Data ====
##############################################################################

# Save distance matrices
folder <- paste0("real_data/", location)
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}

saveRDS(distance_jeff_divergences,
  file = paste0(folder, "/distance_jeff_divergences.rds")
)
saveRDS(distance_wasserstein,
  file = paste0(folder, "/distance_wasserstein.rds")
)

# Save adjacency matrix
saveRDS(W, file = paste0(folder, "/adj_matrix.rds"))
