source("R/utils_plot.R")

result_folder <- "results/"
files <- list.files(result_folder)
files

file_chosen <- files[24] 

##############################################################################
# Load Results ====
##############################################################################
folder <- paste0(result_folder, file_chosen, "/")
filename_results <- "simulation_results.rds"
filename_gt <- "ground_truth.rds"
filename_dist <- "simulation_distance_matrix.rds"
# filename_data <- "simulation_data.rds"
filename_initial_params <- "simulation_initial_params.rds"
filename_time_elapsed <- "time_taken.txt"

filename_results <- paste0(folder, filename_results)
# filename_gt <- paste0(folder, filename_gt)
filename_dist <- paste0(folder, filename_dist)
# filename_data <- paste0(folder, filename_data)
filename_initial_params <- paste0(folder, filename_initial_params)
filename_time_elapsed <- paste0(folder, filename_time_elapsed)

results <- readRDS(file = filename_results)
#ground_truth <- readRDS(file = filename_gt)
dist_matrix <- readRDS(file = filename_dist)
# all_data <- readRDS(file = filename_data)
param <- readRDS(file = filename_initial_params)
elapsed_time <- as.numeric(readLines(con = filename_time_elapsed))

# Create plot directory if it doesn't exist
folder <- paste0(folder, "VI_plots/")
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}

#################################################################
# Plot Data ====
#################################################################

# plot_data(all_data, ground_truth, save = FALSE, folder) # Uncomment to plot data
# plot_distance(dist_matrix, save = FALSE, folder) # Plot distance matrix

#################################################################
# Plot Simulation Results ====
#################################################################

cat("\nFinal hyperparameters:\n")
cat("delta1 =", param$delta1, "\n")
cat("alpha =", param$alpha, "\n")
cat("beta =", param$beta, "\n")
cat("delta2 =", param$delta2, "\n")
cat("gamma =", param$gamma, "\n")
cat("zeta =", param$zeta, "\n")
cat("Initial clusters =", param$initial_cluster, "\n")

##############################################################################
# Visualization ====
##############################################################################

parts <- strsplit(file_chosen, "_")[[1]]
BI <- as.numeric(gsub(".*_BI([0-9]+).*", "\\1", file_chosen))

# Compute point estimate once and reuse
point_estimate <- plot_cls_est(results, BI = BI, save = TRUE, folder = folder)
saveRDS(point_estimate, file = paste0(folder, "point_estimate.rds"))

point_estimate <- readRDS(file = paste0(folder, "point_estimate.rds"))

plot_post_distr(results, BI = BI, save = TRUE, folder = folder)
plot_trace_cls(results, BI = BI, save = TRUE, folder = folder)
plot_trace_U(results, BI = BI, save = TRUE, folder = folder)
plot_acf_U(results, BI = BI, save = TRUE, folder = folder)
# plot_post_sim_matrix(results, BI = BI, save = TRUE, folder = folder)

# Extract the pattern after "real_data_" and before the next "_"
states <- parts[3] # state/regions abbreviation

# Municipalities uses COD_MUN, PUMAs use COD_PUMA as the ID column
id_col <- if (states == "municipalities") "COD_MUN" else "COD_PUMA"
subfolder <- if (states == "municipalities") "geometry" else "counties-pumas"
files_geometry <- if (states == "municipalities") "municipalities.shp" else "counties-pumas.shp"
puma_ids <- sf::st_read(paste0("input/", states, "/", subfolder, "/", files_geometry), quiet = TRUE)[[id_col]]

# plot_map_cls(
#   results = results,
#   BI = BI,
#   point_estimate = point_estimate,
#   unit_ids = puma_ids,
#   puma_dir = paste0("input/", states, "/", subfolder),
#   id_col = id_col,
#   save = TRUE, folder = folder
# )

puma_dir <- paste0("input/", states, "/", subfolder)
input_dir <- paste0("input/", states, "/")

if (states == "municipalities") {
#   plot_hist_cls_comuni(
#     results = results,
#     BI = BI,
#     point_estimate = point_estimate,
#     save = TRUE, folder = folder,
#   )

#   plot_map_prior_mean_comuni(
#     unit_ids = puma_ids,
#     puma_dir = puma_dir,
#     input_dir = input_dir,
#     save = TRUE,
#     folder = folder
#   )
} else {
  plot_hist_cls_pumas(
    results = results,
    BI = BI,
    unit_ids = puma_ids,
    point_estimate = point_estimate,
    input_dir = input_dir,
    save = TRUE,
    folder = folder
  )

#   plot_map_prior_mean(
#     unit_ids = puma_ids,
#     puma_dir = puma_dir,
#     input_dir = input_dir,
#     save = TRUE,
#     folder = folder
#   )
}

# plot_stats(results, ground_truth = ground_truth, BI = BI, save = TRUE, folder = folder)

##############################################################################
# Model comparison ====
##############################################################################

source("R/utils_plot.R")

result_folder <- "results/"
files <- list.files(result_folder)
files
files_chosen <- files[c(14, 15, 16, 17)]  # Choose multiple result folders for comparison
names_corrected <- c("LCF", "NGGP", "NGGPW", "NGGPWx")  # Corrected names for display

# Load point_estimate from different model results
model_results <- list()
for (file in files_chosen) {
  folder <- paste0(result_folder, file, "/VI_plots/")
  point_estimate <- paste0(folder, "point_estimate.rds")
  results <- readRDS(file = point_estimate)
  model_results[[file]] <- results
}

# Compute pairwise ARI and NMI
gt <- as.vector(model_results[[4]])  # Assuming ground truth is the same across models (NGGPWx bin + cont)
ari_matrix <- matrix(NA, nrow = length(model_results), ncol = length(model_results))
nmi_matrix <- matrix(NA, nrow = length(model_results), ncol = length(model_results))
rownames(ari_matrix) <- colnames(ari_matrix) <- names(model_results)
rownames(nmi_matrix) <- colnames(nmi_matrix) <- names(model_results)
for(file in names(model_results)) {
  for(other_file in names(model_results)) {
    ari_matrix[file, other_file] <- arandi(model_results[[file]], model_results[[other_file]])
    nmi_matrix[file, other_file] <- NMI(as.vector(model_results[[file]]), as.vector(model_results[[other_file]]))
  }
}

colnames(ari_matrix) <- names_corrected
rownames(ari_matrix) <- names_corrected
colnames(nmi_matrix) <- names_corrected
rownames(nmi_matrix) <- names_corrected

print("Pairwise ARI Matrix:")
print(ari_matrix)

print("Pairwise NMI Matrix:")
print(nmi_matrix)
