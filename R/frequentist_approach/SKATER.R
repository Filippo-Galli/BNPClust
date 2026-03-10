##############################################################################
# spdep SKATER — Italian municipalities (with islands)
#
# Inputs assumed in your environment:
#   D_mat   : numeric matrix [n x n] of pairwise distances (n = 7903)
#   W       : binary numeric matrix [n x n], 1 = neighbours, 0 otherwise
#             (sparse or dense both fine; must have rownames/colnames OR
#              rows/cols ordered identically to D_mat)
#   k_val   : integer, number of clusters (default 6)
#
# Optional (for saving/plotting):
#   puma_sf, puma_ids, puma_dir, id_col, vi_folder — same as original script
##############################################################################

library(spdep)    # skater, mstree, mat2listw
library(sf)       # only needed for map plots at the end

##############################################################################
# 0.  Parameters ====
##############################################################################
k_val <- 6        # change as needed
set.seed(42)

##############################################################################
# 1.  Load data ====
#     Adjust paths to wherever your RDS files live
##############################################################################
states        <- "municipalities"
files_folder  <- paste0("real_data/", states)
id_col        <- "COD_COM"
subfolder     <- "geometry"
files_geometry<- "municipalities.shp"
puma_dir      <- paste0("input/", states, "/", subfolder)
input_dir     <- paste0("input/", states)

files       <- list.files(files_folder)
dist_matrix <- readRDS(file.path(files_folder, files[4]))
W           <- readRDS(file.path(files_folder, files[1]))

D_mat <- as.matrix(dist_matrix)
W_mat <- as.matrix(W)          # ensure dense for manipulations below
n     <- nrow(D_mat)
cat("n areas:", n, "\n")

##############################################################################
# 2.  Detect and fix islands (isolated nodes in W) ====
#
#  spdep::mstree() requires a CONNECTED graph.  Italian municipalities include
#  proper islands (Sardinia, Sicily, Aeolian Is., etc.) that have zero
#  neighbours in a contiguity W.  We stitch each island to its single closest
#  neighbour according to D_mat so the graph becomes connected while
#  minimally distorting the spatial structure.
##############################################################################

W_aug <- W_mat   # working copy — W itself is NOT modified

row_sums    <- rowSums(W_mat)
island_idx  <- which(row_sums == 0)
cat("Islands (isolated nodes):", length(island_idx), "\n")

if (length(island_idx) > 0) {
  for (i in island_idx) {
    # Find the closest non-island (or any) node by distance
    d_row   <- D_mat[i, ]
    d_row[i] <- Inf                        # exclude self
    nearest  <- which.min(d_row)
    W_aug[i, nearest] <- 1
    W_aug[nearest, i] <- 1
    cat(sprintf("  Island %d connected to neighbour %d (dist = %.4f)\n",
                i, nearest, D_mat[i, nearest]))
  }
}

# Verify connectivity via breadth-first search on W_aug
bfs_connected <- function(adj) {
  visited <- logical(nrow(adj))
  queue   <- 1L
  visited[1] <- TRUE
  while (length(queue) > 0) {
    node    <- queue[1]; queue <- queue[-1]
    nbrs    <- which(adj[node, ] > 0)
    new_nbrs <- nbrs[!visited[nbrs]]
    visited[new_nbrs] <- TRUE
    queue   <- c(queue, new_nbrs)
  }
  all(visited)
}

if (!bfs_connected(W_aug)) {
  # More than one disconnected component remains (e.g. very small archipelagos)
  # Iterate: connect each remaining component to nearest node in the main component
  message("Graph still disconnected after single pass — running multi-pass stitching...")

  repeat {
    # Label connected components via BFS
    comp   <- integer(n)
    comp_id <- 0L
    for (start in seq_len(n)) {
      if (comp[start] != 0) next
      comp_id <- comp_id + 1L
      queue   <- start
      comp[start] <- comp_id
      while (length(queue) > 0) {
        node  <- queue[1]; queue <- queue[-1]
        nbrs  <- which(W_aug[node, ] > 0 & comp == 0)
        comp[nbrs] <- comp_id
        queue <- c(queue, nbrs)
      }
    }
    n_comp <- max(comp)
    cat("  Components remaining:", n_comp, "\n")
    if (n_comp == 1) break

    # Connect the largest component (component 1) with every other component
    main_idx <- which(comp == 1)
    for (cid in 2:n_comp) {
      other_idx <- which(comp == cid)
      # Submatrix of distances between the two components
      sub_D <- D_mat[other_idx, main_idx, drop = FALSE]
      best  <- which(sub_D == min(sub_D), arr.ind = TRUE)[1, ]
      i_node <- other_idx[best[1]]
      j_node <- main_idx[best[2]]
      W_aug[i_node, j_node] <- 1
      W_aug[j_node, i_node] <- 1
      # Merge component cid into component 1 for next iteration
      comp[comp == cid] <- 1
    }
  }
}
cat("Graph connected: TRUE\n")

##############################################################################
# 3.  Build listw from augmented W ====
##############################################################################
w_list <- mat2listw(W_aug, style = "B", zero.policy = TRUE)

##############################################################################
# 4.  Feature matrix for SKATER (MDS on distance matrix) ====
#
#  cmdscale() returns the coordinate matrix directly (no $points wrapper).
#  We use k = 5 dimensions; scale() standardises each column so all
#  features contribute equally when SKATER minimises within-cluster variance.
##############################################################################
cat("\nRunning 5-D MDS on", n, "x", n, "distance matrix...\n")
mds_coords <- cmdscale(D_mat, k = 5)     # returns n x 5 matrix directly
mst_data   <- scale(mds_coords)          # standardise columns

##############################################################################
# 5.  Minimum Spanning Tree ====
##############################################################################
cat("Building MST...\n")
mst_edges <- mstree(w_list)

# mstree() returns a matrix with columns [node_i, node_j, weight];
# skater() needs only the first two (edge endpoints)
edge_mat <- mst_edges[, 1:2, drop = FALSE]

##############################################################################
# 6.  SKATER ====
##############################################################################
cat("\n--- spdep SKATER k =", k_val, "---\n")
gc()
ptm <- proc.time()

skater_res <- skater(
  edges      = edge_mat,
  data       = mst_data,
  ncuts      = k_val - 1,
  crit       = 1,           # minimum cluster size = 1 (handles tiny islands)
  vec.crit   = rep(1, n)
)

elapsed <- (proc.time() - ptm)[3]
cat(sprintf("Runtime: %.1f minutes\n", elapsed / 60))

##############################################################################
# 7.  Extract membership ====
#
#  skater() stores the final partition in $groups when ncuts is supplied.
#  Re-label as consecutive integers 1..k.
##############################################################################
cluster_id <- as.integer(factor(skater_res$groups))
cat("Clusters found:", length(unique(cluster_id)), "\n")
print(table(cluster_id))

##############################################################################
# 8.  Save results ====
##############################################################################
folder_clean  <- gsub("/", "_", files_folder)
output_folder <- paste0("results/", folder_clean, "_SKATER_W_k", k_val)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

vi_folder <- file.path(output_folder, "VI_plots")
dir.create(vi_folder, recursive = TRUE, showWarnings = FALSE)

saveRDS(cluster_id, file.path(vi_folder, "point_estimate.rds"))

results <- list(
  clusters   = cluster_id,
  skater     = skater_res,
  dist_file  = files[4],
  n_islands  = length(island_idx),
  W_augmented = W_aug          # save augmented W for reproducibility
)
saveRDS(results, file.path(output_folder, "results.rds"))

##############################################################################
# 9.  Plots ====
##############################################################################
source("R/utils_plot.R")

puma_sf  <- st_read(file.path(puma_dir, files_geometry))
puma_ids <- puma_sf[[id_col]]

plot_hist_cls_comuni(
  results        = NULL,
  BI             = 0,
  point_estimate = cluster_id,
  save           = TRUE,
  folder         = vi_folder
)

plot_map_cls(
  results        = NULL,
  BI             = 0,
  point_estimate = cluster_id,
  unit_ids       = puma_ids,
  puma_dir       = puma_dir,
  id_col         = id_col,
  save           = TRUE,
  folder         = vi_folder
)

cat("\nDone:", output_folder, "\n")