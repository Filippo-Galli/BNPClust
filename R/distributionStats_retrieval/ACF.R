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
# Retrieve the quantiles for each PUMA ====
##############################################################################
density_list <- list()

for (i in seq_along(data)) {
    pumas_data <- data[[i]]
    # retrieve the mean for the current PUMA
    mean_pumas_data <- mean(pumas_data)
    # Retrieve the quantiles for the current PUMA
    quantiles <- quantile(pumas_data, probs = seq(0, 1, 0.1))
    density_list[[i]] <- list(mean = mean_pumas_data, quantiles = quantiles)
}

folder <- paste0("real_data/", location)
if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
}

# Save the density list
saveRDS(density_list, paste0(folder, "/mean_quantile_list.rds"))
