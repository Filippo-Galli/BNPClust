source("R/utils_plot.R")
source("R/utils.R")

# Load required packages
library(dplyr)

##############################################################################
# CONFIGURATION ====
##############################################################################

# Configuration - modify these as needed
CONFIG <- list(
    input_dir = "input/municipalities",
    output_dir = "real_data/municipalities",
    results_dir = "results",
    # Choose format: "rds" (list format, recommended) or "csv" (long format)
    covariates_format = "csv",
    # Minimum proportion of non-NA values to keep a column (0 = keep all, 1 = no NA allowed)
    na_threshold = 1
)

# Create output directory if it doesn't exist
if (!dir.exists(CONFIG$output_dir)) {
    dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
}

##############################################################################
# LOAD DATA ====
##############################################################################

cat("=== Loading Covariates Data ===\n")

# Load covariates - supports both RDS (list format, recommended) and CSV (long format)
if (CONFIG$covariates_format == "rds") {
    covariates_file <- file.path(CONFIG$input_dir, "full_dataset_covariates.rds")
    if (!file.exists(covariates_file)) {
        stop(sprintf("File not found: %s. Try setting covariates_format = 'csv' or run municipalities.R first.", covariates_file))
    }

    # Load as list (one dataframe per MUN)
    full_dataset_covariates <- readRDS(covariates_file)
    cat(sprintf("Loaded RDS format: %d MUNs\n", length(full_dataset_covariates)))

    # Combine into single dataframe for processing, preserving COD_MUN
    all_covariates <- do.call(rbind, lapply(names(full_dataset_covariates), function(MUN_id) {
        df <- full_dataset_covariates[[MUN_id]]
        df$COD_MUN <- MUN_id
        return(df)
    }))
} else {
    covariates_file <- file.path(CONFIG$input_dir, "full_dataset_covariates.csv")
    if (!file.exists(covariates_file)) {
        stop(sprintf("File not found: %s", covariates_file))
    }

    # Load long format CSV
    all_covariates <- read.csv(covariates_file, stringsAsFactors = FALSE)
    cat(sprintf("Loaded CSV format: %d rows, %d columns\n", nrow(all_covariates), ncol(all_covariates)))

    # Convert to list format for compatibility with existing code
    MUN_ids <- unique(all_covariates$COD_MUN)
    full_dataset_covariates <- lapply(MUN_ids, function(id) {
        subset(all_covariates, COD_MUN == id)
    })
    names(full_dataset_covariates) <- MUN_ids
}

# Ensure COD_MUN exists (for merging with spatial data later)
if (!"COD_MUN" %in% colnames(all_covariates)) {
    if ("STATE_MUN" %in% colnames(all_covariates)) {
        all_covariates$COD_MUN <- all_covariates$STATE_MUN
    } else {
        stop("Neither COD_MUN nor STATE_MUN column found in data")
    }
}

print("=== Checking Municipality Population Totals ===")
problem <- FALSE

# Group by COD_MUN and check totals for each municipality
mun_groups <- all_covariates %>%
    group_by(COD_MUN) %>%
    summarise(
        totale_femmine = sum(TOTALE_FEMMINE, na.rm = TRUE),
        totale_maschi = sum(TOTALE_MASCHI, na.rm = TRUE),
        totale = sum(TOTALE, na.rm = TRUE),
        .groups = "drop"
    )

# Check if totals match for each municipality
inconsistent_muns <- mun_groups %>%
    filter(totale_femmine + totale_maschi != totale)

if (nrow(inconsistent_muns) > 0) {
    cat(sprintf("Found %d municipalities with inconsistent totals:\n", nrow(inconsistent_muns)))
    print(inconsistent_muns)
    problem <- TRUE
} else {
    cat("All municipality totals are consistent.\n")
}

##############################################################################
# Inspect DATA ====
##############################################################################

cat("=== Analysing NA Covariates Data ===\n")

# Get initial column info
initial_cols <- ncol(all_covariates)
cat(sprintf("Initial columns: %d\n", initial_cols))

# 1. Handle NA values based on threshold
if (CONFIG$na_threshold > 0) {
    na_proportions <- sapply(all_covariates, function(col) mean(is.na(col)))
    cols_to_keep <- names(na_proportions)[na_proportions <= (1 - CONFIG$na_threshold)]

    removed_cols <- setdiff(names(all_covariates), cols_to_keep)
    if (length(removed_cols) > 0) {
        cat(sprintf(
            "Removing %d columns with >%.0f%% NA values\n",
            length(removed_cols), (1 - CONFIG$na_threshold) * 100
        ))
    }

    all_covariates <- all_covariates[, cols_to_keep, drop = FALSE]
}

final_cols <- ncol(all_covariates)
cat(sprintf("Final columns after NA filtering: %d (removed %d)\n", final_cols, initial_cols - final_cols))
print("Remaining columns:")
print(colnames(all_covariates))

##############################################################################
# Percentage of women (continuous) ====
##############################################################################

cat("=== Calculating Percentage of Women ===\n")

mun_sex_stats <- all_covariates %>%
    group_by(COD_MUN) %>%
    summarise(
        Women_perc = mean(TOTALE_FEMMINE / TOTALE, na.rm = TRUE),
        number_women = sum(TOTALE_FEMMINE, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(COD_MUN)

saveRDS(mun_sex_stats, file = file.path(CONFIG$output_dir, "mun_sex_stats.rds"))
write.csv(mun_sex_stats,
    file = file.path(CONFIG$output_dir, "mun_sex_stats.csv"),
    row.names = FALSE
)
cat(sprintf("  Saved sex statistics for %d municipalities\n", nrow(mun_sex_stats)))

##############################################################################
# Women Marital Status (categorical) ====
##############################################################################

cat("=== Calculating Marital Status Proportions ===\n")
mun_marital_stats <- all_covariates %>%
    group_by(COD_MUN) %>%
    summarise(
        Married_perc = mean(CONIUGATE, na.rm = TRUE),
        Single_perc = mean(NUBILI, na.rm = TRUE),
        Divorced_perc = mean(DIVORZIATE, na.rm = TRUE),
        Widowed_perc = mean(VEDOVE, na.rm = TRUE),
        number_women = sum(TOTALE_FEMMINE, na.rm = TRUE),
        majority_status = case_when(
            Married_perc >= max(Single_perc, Divorced_perc, Widowed_perc) ~ "Married",
            Single_perc >= max(Married_perc, Divorced_perc, Widowed_perc) ~ "Single",
            Divorced_perc >= max(Married_perc, Single_perc, Widowed_perc) ~ "Divorced",
            Widowed_perc >= max(Married_perc, Single_perc, Divorced_perc) ~ "Widowed",
            TRUE ~ NA_character_
        ),
        .groups = "drop"
    ) %>%
    arrange(COD_MUN)

head(mun_marital_stats)
saveRDS(mun_marital_stats, file = file.path(CONFIG$output_dir, "mun_marital_stats.rds"))
write.csv(mun_marital_stats,
    file = file.path(CONFIG$output_dir, "mun_marital_stats.csv"),
    row.names = FALSE
)
cat(sprintf("  Saved marital status statistics for %d municipalities\n", nrow(mun_marital_stats)))

##############################################################################
# COVARIATE ANALYSIS BY CLUSTER ====
##############################################################################

analyze_clusters <- function(results_dir = CONFIG$results_dir,
                             analysis_output_dir = CONFIG$output_dir) {
    cat("\n=== Analyzing Covariates by Cluster ===\n")

    # Check if results exist
    if (!dir.exists(results_dir)) {
        cat("Results directory not found. Skipping cluster analysis.\n")
        return(NULL)
    }

    files <- list.files(results_dir)
    if (length(files) == 0) {
        cat("No result files found. Skipping cluster analysis.\n")
        return(NULL)
    }

    # Find most recent or specific result file
    #file_chosen <- tail(files[grepl("NGGPWx", files)], 1)
    file_chosen <- files[17]
    if (length(file_chosen) == 0) file_chosen <- files[min(9, length(files))]

    cat(sprintf("Using results from: %s\n", file_chosen))

    # Load point estimates
    point_estimate_path <- file.path(results_dir, file_chosen, "VI_plots", "point_estimate.rds")
    if (!file.exists(point_estimate_path)) {
        cat("Point estimate file not found. Skipping.\n")
        return(NULL)
    }

    point_estimate <- readRDS(point_estimate_path)

    # Determine state and ID column
    parts <- strsplit(file_chosen, "_")[[1]]
    states <- if (any(parts == "municipalities")) "municipalities" else "LA"
    id_col <- if (states == "municipalities") "COD_MUN" else "COD_PUMA"

    # Load unit IDs from shapefile
    shp_path <- file.path(CONFIG$input_dir, "geometry", "municipalities.shp")
    if (!file.exists(shp_path)) {
        cat("Shapefile not found. Cannot match clusters to geography.\n")
        return(NULL)
    }

    shp <- sf::st_read(shp_path, quiet = TRUE)
    unit_ids <- shp[[id_col]]

    # Create cluster dataframe
    names(point_estimate) <- unit_ids
    unique_clusters <- sort(unique(point_estimate))

    cluster_df <- data.frame(
        id = unit_ids,
        cluster = factor(point_estimate, levels = unique_clusters),
        stringsAsFactors = FALSE
    )

    # Merge with covariate data
    if (exists("mun_sex_stats")) {
        sex_df <- merge(cluster_df, mun_sex_stats, by.x = "id", by.y = "COD_MUN", all.x = TRUE)

        sex_cluster_summary <- sex_df %>%
            group_by(cluster) %>%
            summarise(
                n_municipalities = n(),
                avg_perc_female = mean(Women_perc, na.rm = TRUE),
                sd_perc_female = sd(Women_perc, na.rm = TRUE),
                .groups = "drop"
            )

        cat("\nSex Distribution by Cluster:\n")
        print(sex_cluster_summary)

        write.csv(sex_cluster_summary,
            file = file.path(analysis_output_dir, "cluster_sex_summary.csv"),
            row.names = FALSE
        )
    }

    return(list(
        sex = if (exists("sex_cluster_summary")) sex_cluster_summary else NULL
    ))
}

cluster_analysis_results <- analyze_clusters()
