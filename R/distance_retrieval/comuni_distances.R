source("R/utils.R")
source("R/utils_plot.R")

##############################################################################
# Load data files from 'input/' directory ====
##############################################################################
cat("Loading data files...\n")

# Load income distribution data
input_name <- "municipalities"
full_dataset <- read.csv(paste0("input/", input_name, "/full_dataset.csv"), stringsAsFactors = FALSE)
cat("  - Loaded full_dataset.csv:", nrow(full_dataset), "comuni\n")

# Load geometry (if needed for spatial operations)
geometry <- readRDS(paste0("input/", input_name, "/geometry.rds"))
cat("  - Loaded geometry.rds\n")

# Load adjacency matrix
# W <- as.matrix(read.csv("input/Comuni/adj_matrix.csv", row.names = 1))
W <- as.matrix(readRDS(paste0("input/", input_name, "/adj_matrix.rds")))
cat("  - Loaded adj_matrix.rds:", nrow(W), "x", ncol(W), "\n")

if (!isSymmetric(W)) {
  warning("Adjacency matrix W is not symmetric!")
} else {
  cat("  - Adjacency matrix W is symmetric\n")
}

cat("Data loading completed successfully!\n")

##############################################################################
# Read histogram for each Comune ====
##############################################################################
cat("Extracting histograms for each Comune...\n")

# First, let's check what columns we have
cat("Available columns in full_dataset:\n")
print(head(colnames(full_dataset), 20))

# Identify income bracket columns
# R adds 'X' prefix to column names starting with numbers
all_cols <- colnames(full_dataset)

# Pattern: columns starting with X followed by digits (income brackets)
income_columns <- all_cols[grep("^X[0-9]", all_cols)]

# Exclude the first column if it's X0. (it might be a sum or total)
# Keep it if it represents the 0-10000 bracket
if (length(income_columns) > 0 && income_columns[1] == "X0.") {
  # Check if we have X0.10000, if so, remove X0.
  if ("X0.10000" %in% income_columns) {
    income_columns <- income_columns[income_columns != "X0."]
  }
}

cat("\nDetected income columns:\n")
print(income_columns)

if (length(income_columns) == 0) {
  stop("Could not automatically detect income columns. Please specify them manually.")
}

# Initialize list to store histograms
hist_list <- list()

# Extract histogram data for each Comune
n_comuni <- nrow(full_dataset)

for (i in 1:n_comuni) {
  # Extract counts for each income bracket
  income_data <- as.numeric(full_dataset[i, income_columns])

  # Define bin breaks based on column names
  # Try to extract numbers from column names
  bin_breaks_lower <- sapply(income_columns, function(col) {
    # Extract first number from column name
    nums <- as.numeric(gsub("[^0-9]", " ", col))
    nums <- nums[!is.na(nums)]
    if (length(nums) > 0) {
      return(nums[1])
    } else {
      return(NA)
    }
  })

  # If we can't parse column names, use default breaks
  if (all(is.na(bin_breaks_lower))) {
    bin_breaks <- c(0, 10000, 15000, 26000, 55000, 75000, 120000, Inf)
  } else {
    # Create breaks from column names, adding 0 at start and Inf at end
    bin_breaks <- c(0, sort(unique(bin_breaks_lower[!is.na(bin_breaks_lower)])), Inf)
  }

  # Ensure we have the right number of breaks (n_bins + 1)
  if (length(bin_breaks) != length(income_data) + 1) {
    # Fallback: create evenly spaced breaks
    bin_breaks <- seq(0, max(bin_breaks[!is.infinite(bin_breaks)]) * 1.2,
      length.out = length(income_data) + 1
    )
    bin_breaks[length(bin_breaks)] <- Inf
  }

  # Create histogram object
  hist_obj <- list(
    breaks = bin_breaks,
    counts = income_data,
    mids = (bin_breaks[-length(bin_breaks)] + bin_breaks[-1]) / 2,
    density = income_data / sum(income_data, na.rm = TRUE),
    comune = full_dataset$NAME_MUN[i],
    cod_mun = full_dataset$COD_MUN[i]
  )

  hist_list[[i]] <- hist_obj

  if (i %% 500 == 0) {
    cat("  Processed", i, "of", n_comuni, "municipalities\n")
  }
}

# Name the list elements
names(hist_list) <- full_dataset$NAME_MUN

cat("Histogram extraction completed:", length(hist_list), "histograms created\n")

# Optional: Visualize a sample histogram
if (length(hist_list) > 0) {
  cat("\nSample histogram for:", names(hist_list)[1], "\n")
  cat("  Total count:", sum(hist_list[[1]]$counts, na.rm = TRUE), "\n")
  cat("  Income brackets:", length(hist_list[[1]]$counts), "\n")
}

##############################################################################
# Distances between histograms ====
##############################################################################
cat("\nComputing distances between histograms...\n")

distance_jeff_divergences <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)
distance_cm <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)
distance_wasserstein <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)
distance_mean <- matrix(0,
  nrow = length(hist_list),
  ncol = length(hist_list)
)

n <- length(hist_list)
total_iterations <- n * (n + 1) / 2
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
iteration <- 0

for (i in seq_along(hist_list)) {
  for (j in i:length(hist_list)) {
    distance_jeff_divergences[i, j] <-
      compute_hist_distances(hist_list[[i]],
        hist_list[[j]],
        type = "Jeff"
      )
    distance_cm[i, j] <-
      compute_hist_distances(hist_list[[i]],
        hist_list[[j]],
        type = "CM"
      )
    distance_wasserstein[i, j] <-
      compute_hist_distances(hist_list[[i]],
        hist_list[[j]],
        type = "Wasserstein"
      )

    # Compute mean difference using histogram midpoints weighted by counts
    mean_i <- sum(hist_list[[i]]$mids * hist_list[[i]]$counts, na.rm = TRUE) /
      sum(hist_list[[i]]$counts, na.rm = TRUE)
    mean_j <- sum(hist_list[[j]]$mids * hist_list[[j]]$counts, na.rm = TRUE) /
      sum(hist_list[[j]]$counts, na.rm = TRUE)
    distance_mean[i, j] <- abs(mean_i - mean_j)

    if (i != j) {
      distance_jeff_divergences[j, i] <- distance_jeff_divergences[i, j]
      distance_cm[j, i] <- distance_cm[i, j]
      distance_wasserstein[j, i] <- distance_wasserstein[i, j]
      distance_mean[j, i] <- distance_mean[i, j]
    }

    # Update progress bar
    iteration <- iteration + 1
    setTxtProgressBar(pb, iteration)
  }
}

cat("Distance computation completed\n")

##############################################################################
# Plot Distance (optional) ====
##############################################################################
plot_folder <- paste0("results/distance_plots/", input_name, "/")
if(!dir.exists(plot_folder)) {
  dir.create(plot_folder, recursive = TRUE)
}

cat("\nCreating distance plots...\n")
plot_distance(distance_jeff_divergences,
title = "Jeffreys Divergence Distance", save = TRUE,
folder = plot_folder
)
plot_distance(distance_cm,
title = "Cramer-von Mises Distance", save = TRUE,
folder = plot_folder
)
plot_distance(distance_wasserstein,
title = "Wasserstein Distance", save = TRUE,
folder = plot_folder
)
plot_distance(distance_mean,
title = "Mean Difference", save = TRUE,
folder = plot_folder
)


##############################################################################
# Save Data ====
##############################################################################
cat("\nSaving results...\n")

# Create output folder if it doesn't exist
folder <- paste0("real_data/", input_name)
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}

# Save distance matrices
saveRDS(distance_jeff_divergences,
  file = paste0(folder, "/distance_jeff_divergences.rds")
)
saveRDS(distance_cm,
  file = paste0(folder, "/distance_cm.rds")
)
saveRDS(distance_wasserstein,
  file = paste0(folder, "/distance_wasserstein.rds")
)
saveRDS(distance_mean,
  file = paste0(folder, "/distance_mean.rds")
)

# Save adjacency matrix
saveRDS(W, file = paste0(folder, "/adj_matrix.rds"))

cat("All results saved to", folder, "\n")
cat("\nSummary:\n")
cat("  - Number of comuni:", length(hist_list), "\n")
cat("  - Adjacency matrix dimension:", nrow(W), "x", ncol(W), "\n")
cat("  - Number of neighbors (edges):", sum(W) / 2, "\n")
cat("  - Distance matrices saved: 4 types\n")
cat("  - Histogram data saved\n")
