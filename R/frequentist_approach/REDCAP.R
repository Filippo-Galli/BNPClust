##############################################################################
# Setup (rgeoda) ====
##############################################################################
source("R/utils_plot.R")
library(rgeoda)
library(sf)
library(cluster)

# Data loading
files_folder <- "real_data/LA"
files <- list.files(files_folder)
dist_matrix <- readRDS(file.path(files_folder, files[4])) 
W <- readRDS(file.path(files_folder, files[1]))

D_mat <- as.matrix(dist_matrix)
n_pumas <- nrow(D_mat)

##############################################################################
# Load PUMA sf + Weights ====
##############################################################################
puma_sf <- st_read("input/LA/counties-pumas/counties-pumas.shp")

# Queen weights for spatial constraint
w_queen <- queen_weights(puma_sf)

# Data: Using the full distance matrix as features for maximum fidelity
# This ensures the frequentist model "sees" everything the Bayesian one does
data_df <- as.data.frame(D_mat)

##############################################################################
# REDCAP Loop ====
##############################################################################
puma_ids <- puma_sf$COD_PUMA
input_dir <- "input/LA/"
puma_dir <- "input/LA/counties-pumas"
k_val <- 6 # Match posterior mode

cat("\n--- REDCAP k =", k_val, "---\n")

# Regionalization with Dynamically Constrained Agglomerative Clustering
redcap_res <- redcap(
    k = k_val,
    w = w_queen,
    df = data_df,
    method = "fullorder-averagelinkage"
)

# Extract cluster IDs (ensuring it is a simple vector)
cluster_id <- as.vector(redcap_res$Clusters)

cat("Clusters found:", length(unique(cluster_id)), "\n")
print(table(cluster_id))

##############################################################################
# Save Point Estimate & Results ====
##############################################################################

# 1. Standard Results Output
folder_clean <- gsub("/", "_", files_folder)
output_folder <- paste0("results/", folder_clean, "_REDCAP_k", k_val)
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# 2. Save specifically for VI_plots
vi_folder <- paste0(output_folder, "/VI_plots") # not since we are using VI but since other functions already expect this structure
if (!dir.exists(vi_folder)) dir.create(vi_folder, recursive = TRUE)
saveRDS(cluster_id, file.path(vi_folder, "point_estimate.rds"))

# Calculate Silhouette using the original distance matrix
sil <- silhouette(cluster_id, D_mat)
results <- list(
    clusters = cluster_id,
    redcap = redcap_res,
    sil_width = mean(sil[, 3]),
    dist_file = files[4]
)

saveRDS(results, file.path(output_folder, "results.rds"))

##############################################################################
# Plots ====
##############################################################################
plot_hist_cls_pumas(
    results = NULL, BI = 0, point_estimate = cluster_id,
    save = TRUE, folder = output_folder, unit_ids = puma_ids, input_dir = input_dir
)

plot_map_cls(
    results = NULL, BI = 0, point_estimate = cluster_id, unit_ids = puma_ids,
    puma_dir = puma_dir, id_col = "COD_PUMA",
    save = TRUE, folder = output_folder
)

cat("\nDone. Point estimate saved in VI_plots/ and full results in", output_folder, "\n")