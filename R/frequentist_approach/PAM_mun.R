##############################################################################
# Setup ====
##############################################################################
source("R/utils_plot.R")
suppressWarnings(library(sf))
suppressWarnings(library(cluster))

# Data loading
files_folder <- "real_data/municipalities/"
files <- list.files(files_folder)
dist_matrix <- readRDS(file.path(files_folder, files[4]))
W <- readRDS(file.path(files_folder, files[1]))

D_mat <- as.matrix(dist_matrix)
n_units <- nrow(D_mat)
rm(dist_matrix)
gc()

##############################################################################
# Load Geometry ====
##############################################################################
states <- "municipalities"
id_col <- if (states == "municipalities") "COD_MUN" else "COD_PUMA"
subfolder <- if (states == "municipalities") "geometry" else "counties-pumas"
files_geometry <- if (states == "municipalities") "municipalities.shp" else "counties-pumas.shp"

puma_sf <- sf::st_read(
    paste0("input/", states, "/", subfolder, "/", files_geometry),
    quiet = TRUE
)

##############################################################################
# PAM Clustering ====
##############################################################################

k_val <- 6

cat("\n--- PAM k =", k_val, "---\n")

# Convert distance matrix to dist object for PAM
cat("Converting distance matrix for PAM...\n")
D_dist <- as.dist(D_mat)

# Run PAM
cat("Running PAM clustering...\n")
pam_res <- pam(D_dist, k = k_val, diss = TRUE, pamonce = 5)

cluster_id <- pam_res$clustering

cat("Clusters found:", length(unique(cluster_id)), "\n")
cat("Medoids (unit indices):", pam_res$id.med, "\n")

##############################################################################
# Save Point Estimate & Results ====
##############################################################################
unit_ids <- puma_sf[[id_col]]
input_dir <- paste0("input/", states, "/")
puma_dir <- paste0("input/", states, "/", subfolder)
folder_clean <- gsub("/", "_", files_folder)
output_folder <- paste0("results/", folder_clean, "_PAM_k", k_val)

if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
}
vi_folder <- file.path(output_folder, "VI_plots")
if (!dir.exists(vi_folder)) {
    dir.create(vi_folder, recursive = TRUE)
}

saveRDS(cluster_id, file.path(vi_folder, "point_estimate.rds"))

results <- list(
    clusters  = cluster_id,
    pam       = pam_res,
    dist_file = files[4]
)
saveRDS(results, file.path(output_folder, "results.rds"))

##############################################################################
# Plots ====
##############################################################################
if(states == "municipalities") {
  plot_hist_cls_comuni(
    results = NULL,
    BI = 0,
    point_estimate = cluster_id,
    save = TRUE, folder = output_folder,
  )
}

plot_map_cls(
    results = NULL, BI = 0, point_estimate = cluster_id,
    unit_ids = unit_ids, puma_dir = puma_dir, id_col = id_col,
    save = TRUE, folder = output_folder
)

cat("\nDone. Saved in:", output_folder, "\n")