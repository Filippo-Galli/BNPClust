## @file utils.R
## @brief Utility functions for Bayesian clustering analysis and MCMC visualization
## @author Filippo Galli
## @date 2025

suppressMessages(library(Rcpp))
suppressMessages(library(ggplot2))
suppressMessages(library(MASS)) # For fitdistr function
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(spam)) # For comp.psm
suppressMessages(library(fields)) # For minVI
suppressMessages(library(viridisLite)) # For color scales
suppressMessages(library(RColorBrewer)) # For color palettes
suppressMessages(library(pheatmap)) # For heatmaps
suppressMessages(library(mcclust.ext)) # For MCMC clustering functions
suppressMessages(library(mvtnorm))
suppressMessages(library(gtools))
suppressMessages(library(salso))
suppressMessages(library(aricode))
suppressMessages(library(reshape2))
suppressMessages(library(cluster))
suppressMessages(library(label.switching))
suppressMessages(library(coda))

## @brief Generate consistent colors for clusters
## @param n_clusters Number of clusters
## @return Named vector of colors
get_cluster_colors <- function(n_clusters) {
    # Use viridis turbo for many clusters
    colors <- viridisLite::turbo(n_clusters)

    # Return named vector with cluster numbers as names
    names(colors) <- as.character(seq_len(n_clusters))
    return(colors)
}

plot_distance <- function(dist_matrix, cls = NULL,
                          save = FALSE, folder = "results/plots/",
                          title = "Distance Histogram",
                          normalize = FALSE,
                          breaks = 30) {
    # Calculate distance matrix
    dist_matrix <- as.matrix(dist_matrix)

    cls_exist <- !is.null(cls)

    if (!cls_exist) {
        cls <- rep(1, nrow(dist_matrix)) # All in one cluster
    }

    # Get upper triangular indices to avoid duplicates
    upper_tri_indices <- which(upper.tri(dist_matrix), arr.ind = TRUE)

    # Extract distances and determine if they are intra or inter-cluster
    distances <- dist_matrix[upper_tri_indices]
    cluster_pairs <- cbind(
        cls[upper_tri_indices[, 1]],
        cls[upper_tri_indices[, 2]]
    )

    # Classify distances as intra-cluster or inter-cluster
    intra_cluster <- distances[cluster_pairs[, 1] == cluster_pairs[, 2]]
    inter_cluster <- distances[cluster_pairs[, 1] != cluster_pairs[, 2]]

    # Pre-compute histograms for ylim calculation
    h_intra <- hist(intra_cluster, breaks = breaks, plot = FALSE)
    h_inter <- if (cls_exist) hist(inter_cluster, breaks = breaks, plot = FALSE) else NULL

    # Determine y-axis label and limits based on normalization
    ylab <- if (normalize) "Density" else "Frequency"

    ylim_max <- if (normalize) {
        max(h_intra$density, if (cls_exist) h_inter$density else 0) * 1.1
    } else {
        max(h_intra$counts, if (cls_exist) h_inter$counts else 0) * 1.1
    }

    # Create histogram with overlaid distributions
    hist(intra_cluster,
        breaks = breaks,
        col = rgb(1, 0.5, 0, 0.7),
        main = "",
        xlab = "Distance",
        ylab = ylab,
        probability = normalize,
        xlim = range(c(intra_cluster, inter_cluster)),
        ylim = c(0, ylim_max)
    )

    if (cls_exist) {
        hist(inter_cluster,
            breaks = breaks,
            col = rgb(0, 0, 1, 0.7),
            probability = normalize,
            add = TRUE
        )
    }

    # Add legend
    legend("topright",
        legend = if (cls_exist) c("Intra-cluster", "Inter-cluster") else "Distance",
        fill = if (cls_exist) c(rgb(1, 0.5, 0, 0.7), rgb(0, 0, 1, 0.7)) else rgb(1, 0.5, 0, 0.7),
        bty = "n"
    )
    # title(main = title)

    # Save plot if needed
    if (save) {
        if (!dir.exists(folder)) {
            dir.create(folder, recursive = TRUE)
        }
        words_title <- gsub(" ", "_", title)
        dev.copy(png,
            filename = paste0(folder, words_title, ".png"),
            width = 2400, height = 1800, res = 300
        )
        dev.off()
    }
}


plot_post_distr <- function(results, BI, save = FALSE, folder = "results/plots/") {
    k_values <- unlist(results$K)

    if (BI > 0 && length(k_values) > BI) {
        k_values <- k_values[(BI + 1):length(k_values)]
    }

    post_k <- table(k_values) / length(k_values)
    df <- data.frame(
        cluster_found = as.numeric(names(post_k)),
        rel_freq = as.numeric(post_k)
    )

    # Adaptive: if few unique clusters, show all; otherwise use pretty breaks
    n_unique_clusters <- length(unique(df$cluster_found))

    p1 <- ggplot(data = df, aes(x = factor(cluster_found), y = rel_freq)) +
        geom_col(width = 0.8) +
        labs(
            x = "Cluster Found",
            y = "Relative Frequency",
            # title = paste("Posterior distribution number of clusters")
            # title = paste("Prior distribution number of clusters ")
        ) +
        theme(
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 15),
            text = element_text(size = 15),
            panel.background = element_blank(),
            panel.grid.major = element_line(color = "grey95"),
            panel.grid.minor = element_line(color = "grey95")
        )

    # Adaptive scaling
    if (n_unique_clusters <= 20) {
        # Show all cluster values when there are few
        p1 <- p1 + scale_x_discrete(drop = FALSE)
    } else {
        # Show fewer labels (every nth) when many clusters to avoid overcrowding
        step <- max(1, floor(n_unique_clusters / 10))
        labels_to_show <- seq(1, n_unique_clusters, by = step)
        p1 <- p1 + scale_x_discrete(labels = function(x) ifelse(as.numeric(x) %in% labels_to_show, x, ""))
    }

    print(p1)

    if (save) {
        ggsave(
            filename = paste0(folder, "posterior_num_clusters.png"),
            plot = p1,
            width = 8, height = 6
        )
    }
}

plot_trace_cls <- function(results, BI, save = FALSE, folder = "results/plots/") {
    k_values <- unlist(results$K)

    # Apply burn-in period
    if (BI > 0 && length(k_values) > BI) {
        k_values <- k_values[(BI + 1):length(k_values)]
    }

    ### Second plot - Trace of number of clusters
    k_df <- data.frame(
        Iteration = seq_along(k_values),
        NumClusters = k_values
    )

    p2 <- ggplot(k_df, aes(x = Iteration, y = NumClusters)) +
        geom_line() +
        labs(
            x = "Iteration",
            y = "Number of clusters",
            # title = paste("Trace Plot (Burn-in:", BI, "iterations)")
        ) +
        theme(
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            text = element_text(size = 15),
            panel.background = element_blank(),
            panel.grid.major = element_line(color = "grey95"),
            panel.grid.minor = element_line(color = "grey95"),
            legend.position = "top"
        )

    print(p2)
    if (save) {
        ggsave(
            filename = paste0(folder, "traceplot.png"),
            plot = p2, width = 8, height = 6
        )
    }
}

plot_post_sim_matrix <- function(results, BI, save = FALSE, folder = "results/plots/") {
    #### Apply burn-in to allocations
    allocations_post_burnin <- results$allocations
    if (BI > 0 && length(allocations_post_burnin) > BI) {
        allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
    }

    #### Compute posterior similarity matrix using salso::psm
    n <- length(allocations_post_burnin[[1]])
    n_iter <- length(allocations_post_burnin)

    cat("Data dimensions: n =", n, ", n_iter =", n_iter, "\n")

    # Convert allocations to matrix format (each row is one iteration)
    alloc_matrix <- matrix(unlist(allocations_post_burnin),
        nrow = n_iter,
        ncol = n,
        byrow = TRUE
    )

    cat("Computing posterior similarity matrix using salso::psm...\n")
    similarity_matrix <- salso::psm(alloc_matrix, nCores = 10)

    # Hierarchical clustering
    cat("Performing hierarchical clustering...\n")
    hc <- hclust(as.dist(1 - similarity_matrix))
    ordered_idx <- hc$order
    sim_reordered <- similarity_matrix[ordered_idx, ordered_idx]

    cat("Creating plot (full", n, "x", n, "matrix)...\n")

    # Convert reordered similarity matrix to long format for ggplot2
    # Use indices directly to avoid large melt operation
    sim_df <- data.frame(
        i = rep(1:n, each = n),
        j = rep(1:n, times = n),
        similarity = as.vector(sim_reordered)
    )

    # Create ggplot with geom_raster (more efficient than geom_tile)
    p <- ggplot(sim_df, aes(x = j, y = i, fill = similarity)) +
        geom_raster() +
        scale_fill_gradient(
            low = "white", high = "darkblue",
            limits = c(0, 1),
            name = "Posterior\nSimilarity"
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(
            x = "Data Point Index (clustered)",
            y = "Data Point Index (clustered)",
            # title = "Posterior Similarity Matrix"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid = element_blank(),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10)
        ) +
        coord_fixed()

    print(p)

    if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        ggsave(
            filename = paste0(folder, "similarity_matrix.png"),
            plot = p,
            width = 10, height = 10, dpi = 400
        )
        cat("Plot saved to:", paste0(folder, "similarity_matrix.png"), "\n")
    }
}

plot_stats <- function(results, ground_truth, BI, save = FALSE, folder = "results/plots/") {
    #### Apply burn-in to allocations
    allocations_post_burnin <- results$allocations
    if (BI > 0 && length(allocations_post_burnin) > BI) {
        allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
    }

    #### Convert allocations to matrix format for SALSO
    C <- matrix(unlist(lapply(allocations_post_burnin, function(x) x + 1)),
        nrow = length(allocations_post_burnin),
        ncol = length(ground_truth),
        byrow = TRUE
    )

    #### Get point estimate using Variation of Information (VI) loss
    point_estimate <- salso::salso(C,
        loss = "VI",
        maxNClusters = 200,
        maxZealousAttempts = 1000
    )

    #### Print results
    cat("=== SALSO Clustering Results (Post Burn-in) ===\n")
    cat("Cluster Sizes:\n")
    print(table(point_estimate))

    point_labels <- as.vector(point_estimate)

    #### Adjusted Rand Index (ARI)
    cat(
        "\nAdjusted Rand Index:",
        arandi(point_estimate, ground_truth), "\n"
    )

    #### NMI
    cat(
        "Normalized Mutual Information:",
        NMI(point_labels, ground_truth), "\n"
    )

    #### VI
    cat("Variation of Information:", NVI(point_labels, ground_truth), "\n")

    if (save) {
        stats_file <- paste0(folder, "salso_stats.txt")
        write(paste("Adjusted Rand Index:", arandi(point_estimate, ground_truth)), file = stats_file)
        write(paste("Normalized Mutual Information:", NMI(point_labels, ground_truth)), file = stats_file, append = TRUE)
        write(paste("Variation of Information:", NVI(point_labels, ground_truth)), file = stats_file, append = TRUE)
    }
}

plot_trace_U <- function(results, BI, save = FALSE, folder = "results/plots/") {
    U_after_burnin <- results$U
    if (BI > 0 && length(U_after_burnin) > BI) {
        U_after_burnin <- U_after_burnin[(BI + 1):length(U_after_burnin)]
    }

    draw_plot <- function() {
        plot(U_after_burnin, type = "l", xlab = "Iteration", ylab = "U")
        abline(h = mean(U_after_burnin), col = "red", lty = 2)
        legend("topright",
            legend = sprintf("Mean U = %.3f", mean(U_after_burnin)),
            col = "red", lty = 2, bty = "n"
        )
        # title(main = sprintf("Trace of U over MCMC iterations\n(mean U = %.3f)", mean(U_after_burnin)))
    }

    draw_plot()

    if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        file <- file.path(folder, "U_trace.png")
        if (requireNamespace("ragg", quietly = TRUE)) {
            ragg::agg_png(file, width = 2400, height = 1800, res = 300)
        } else {
            png(file, width = 2400, height = 1800, res = 300)
        }
        draw_plot()
        dev.off()
    }
}

plot_acf_U <- function(results, BI, save = FALSE, folder = "results/plots/") {
    U_after_burnin <- results$U
    if (BI > 0 && length(U_after_burnin) > BI) {
        U_after_burnin <- U_after_burnin[(BI + 1):length(U_after_burnin)]
    }

    draw_plot <- function() {
        acf(U_after_burnin, main = "ACF of U over MCMC iterations")
    }

    draw_plot()

    if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        file <- file.path(folder, "U_acf.png")
        if (requireNamespace("ragg", quietly = TRUE)) {
            ragg::agg_png(file, width = 2400, height = 1800, res = 300)
        } else {
            png(file, width = 2400, height = 1800, res = 300)
        }
        draw_plot()
        dev.off()
    }
}

plot_k_means <- function(dist_matrix, max_k = 10) {
    # Ensure dist_matrix is in the right format
    if (!inherits(dist_matrix, "dist")) {
        dist_matrix <- as.dist(dist_matrix)
    }

    mds_result <- cmdscale(dist_matrix, k = 2)

    # Step 1: Compute elbow method for K selection using kmeans
    wss <- numeric(max_k)

    for (k in 1:max_k) {
        kmeans_result <- kmeans(mds_result, centers = k, nstart = 25)
        wss[k] <- kmeans_result$tot.withinss
    }

    # Plot the elbow curve
    elbow_data <- data.frame(K = 1:max_k, WSS = wss)
    elbow_plot <- ggplot(elbow_data, aes(x = K, y = WSS)) +
        geom_line() +
        geom_point() +
        labs(
            title = "Elbow Method for Optimal K (k-means)",
            x = "Number of Clusters (K)",
            y = "Total Dissimilarity"
        ) +
        theme_minimal()

    print(elbow_plot)

    # Return the WSS values for further analysis
    return(invisible(elbow_data))
}

plot_data <- function(all_data, cluster_labels, save = FALSE, folder = "results/plots/") {
    ground_truth <- as.factor(cluster_labels)

    # Create data frame for plotting
    plot_data <- data.frame(
        x = all_data[, 1],
        y = all_data[, 2],
        cluster = ground_truth
    )

    # Plot the clusters
    p <- ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
        geom_point(size = 3) +
        labs(
            # title = "Clusters",
            x = "X coordinate",
            y = "Y coordinate"
        ) +
        theme_minimal()

    print(p)

    # Save the plot
    if (save) {
        ggsave(
            filename = paste0(folder, "data_clusters.png"),
            plot = p3, width = 8, height = 6
        )
    }
}

plot_cls_est <- function(results, BI, save = FALSE, start_time, end_time, folder = "results/plots/") {
    cat("Computing point estimate using SALSO...\n")
    #### Apply burn-in to allocations
    allocations_post_burnin <- results$allocations
    if (BI > 0 && length(allocations_post_burnin) > BI) {
        allocations_post_burnin <- allocations_post_burnin[(BI + 1):length(allocations_post_burnin)]
    }

    #### Convert allocations to matrix format for SALSO
    C <- matrix(unlist(lapply(allocations_post_burnin, function(x) x + 1)),
        nrow = length(allocations_post_burnin),
        ncol = length(allocations_post_burnin[[1]]),
        byrow = TRUE
    )

    #### Get point estimate using Variation of Information (VI) loss
    point_estimate <- salso::salso(C,
        loss = "VI",
        maxNClusters = 200,
        maxZealousAttempts = 1000
    )

    #### Print results
    cat("=== SALSO Clustering Results (Post Burn-in) ===\n")
    cat("Cluster Sizes:\n")
    print(table(point_estimate))

    cat("=== ESS ===\n")
    ess_value <- effectiveSize(as.mcmc(unlist(results$K)))
    ess_seconds <- ess_value / elapsed_time
    cat("Effective Sample Size (ESS):", ess_value, "\n")
    cat("ESS per second:", ess_seconds, "\n")

    if (save) {
        cls_file <- paste0(folder, "salso_cluster_estimate.txt")
        write("Cluster Sizes:", file = cls_file)
        write(capture.output(table(point_estimate)), file = cls_file, append = TRUE)
        write("\n=== ESS ===", file = cls_file, append = TRUE)
        write(paste("ESS:", ess_value), file = cls_file, append = TRUE)
        write(paste("ESS per second:", ess_seconds), file = cls_file, append = TRUE)
    }
    cat("Computed point estimate using SALSO...\n")
    return(invisible(point_estimate))
}

plot_map_cls <- function(results, BI, point_estimate = NULL, save = FALSE,
                         folder = "results/plots/",
                         puma_dir = "input/LA/counties-pumas",
                         id_col = "COD_PUMA",
                         unit_ids = NULL,
                         simplify_geom = TRUE,
                         dTolerance = 100) {
    cat("Computing map of cluster assignments...\n")

    if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Package 'sf' is required for plot_map_cls().")
    }

    # Only compute point_estimate if not provided
    if (is.null(point_estimate)) {
        point_estimate <- plot_cls_est(results, BI, save = FALSE)
    }

    shp <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp) == 0) {
        stop("No .shp file found in '", puma_dir, "'.")
    }

    # Assign names to point_estimate if missing
    if (is.null(names(point_estimate))) {
        candidate_ids <- unit_ids %||% results$unit_ids %||% results$puma_ids
        if (is.null(candidate_ids) || length(candidate_ids) != length(point_estimate)) {
            stop("Provide unit_ids or store results$unit_ids/results$puma_ids matching the PUMAs.")
        }
        names(point_estimate) <- candidate_ids
    }

    # Pre-compute cluster factor with proper levels (avoids recomputation)
    unique_clusters <- sort(unique(point_estimate))
    n_clusters <- length(unique_clusters)
    cluster_colors <- get_cluster_colors(n_clusters)

    # Create cluster data frame with pre-factored clusters
    cluster_df <- data.frame(
        id = names(point_estimate),
        cluster = factor(point_estimate, levels = unique_clusters),
        stringsAsFactors = FALSE
    )
    names(cluster_df)[1] <- id_col

    # Read shapefile - only columns we need
    geom <- sf::st_read(shp[1], quiet = TRUE)
    geom <- geom[, c(id_col, attr(geom, "sf_column"))]
    geom[[id_col]] <- as.character(geom[[id_col]])

    # Use match-based join (faster than merge for simple cases)
    idx <- match(geom[[id_col]], cluster_df[[id_col]])
    keep <- !is.na(idx)
    geom <- geom[keep, ]
    geom$cluster <- cluster_df$cluster[idx[keep]]

    # Simplify geometry for faster rendering (optional)
    if (simplify_geom && nrow(geom) > 50) {
        geom <- sf::st_simplify(geom, preserveTopology = TRUE, dTolerance = dTolerance)
    }

    p <- ggplot2::ggplot(geom) +
        ggplot2::geom_sf(ggplot2::aes(fill = cluster), color = "grey60", linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = cluster_colors, drop = FALSE, na.value = "lightgrey") +
        ggplot2::labs(
            # title = "PUMAs by Cluster Assignment",
            fill = "Cluster"
        ) +
        ggplot2::theme_minimal()

    print(p)

    if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        ggplot2::ggsave(file.path(folder, "puma_clusters.png"), p, width = 10, height = 8)
    }

    invisible(p)
}

plot_map_prior_mean <- function(save = FALSE, folder = "results/plots/",
                                puma_dir = "input/counties-pumas",
                                input_dir = "input/CA/",
                                id_col = "COD_PUMA",
                                unit_ids = NULL,
                                simplify_tol = 100) {
    cat("Computing map of mean income ...\n")

    # --- Check Dependencies ---
    if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required.")
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

    # --- 1. Load Data ---
    data <- NULL
    csv_file <- file.path(input_dir, "full_dataset.csv")
    rds_file <- file.path(input_dir, "full_dataset.rds")

    if (file.exists(rds_file)) {
        data <- readRDS(rds_file)
    } else if (file.exists(csv_file)) {
        data <- read.csv(csv_file)
    } else {
        stop("Could not find full_dataset in ", input_dir)
    }

    # --- 2. Extract Numeric Data and Compute Means ---
    # Remove NAME_PUMA 
    data <- data %>% select(-matches("NAME_PUMA", ignore.case = TRUE))

    cat("Calculating means grouped by PUMA...\n")
    # Assuming COD_PUMA is the identifier for PUMAs, group by it and calculate means of log_income
    prior_means <- data %>%
        group_by_at("COD_PUMA") %>%
        summarise(prior_mean = mean(exp(log_income), na.rm = TRUE)) %>%
        ungroup()

    # --- 3. Synchronize IDs ---
    # If unit_ids weren't provided, they MUST be the names of the numeric columns
    if (is.null(unit_ids)) {
        unit_ids <- names(prior_means)
    }

    # Double check length match
    if (length(prior_means$prior_mean) != length(unit_ids)) {
        stop(sprintf(
            "Length mismatch: IDs (%d) vs Data Elements (%d).",
            length(unit_ids), length(prior_means$prior_mean)
        ))
    }
    names(prior_means$prior_mean) <- unit_ids

    # --- 4. Load Map Geometry ---
    cat("Loading map geometry...\n")
    shp_files <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) stop("No .shp file found.")

    geom <- sf::st_read(shp_files[1], quiet = TRUE)

    if (!id_col %in% names(geom)) {
        stop("Column '", id_col, "' not found in shapefile.")
    }

    # --- 5. Merge and Plot ---
    prior_df <- data.frame(
        id_temp = as.character(unit_ids),
        prior_mean = as.numeric(prior_means$prior_mean),
        stringsAsFactors = FALSE
    )
    names(prior_df)[1] <- id_col

    geom[[id_col]] <- as.character(geom[[id_col]])
    map_data <- dplyr::left_join(geom, prior_df, by = id_col)

    cat("Generating plot...\n")
    lbl_fun <- if (requireNamespace("scales", quietly = TRUE)) scales::comma else waiver()

    p <- ggplot2::ggplot(map_data) +
        ggplot2::geom_sf(ggplot2::aes(fill = prior_mean), color = "white", linewidth = 0.05) +
        ggplot2::scale_fill_viridis_c(
            option = "viridis",
            na.value = "lightgrey",
            labels = lbl_fun,
            name = "Mean Income"
        ) +
        ggplot2::theme_light() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
            panel.border = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(color = "grey60", size = 8)
        )

    if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        ggplot2::ggsave(file.path(folder, "puma_prior_means.png"), p, width = 10, height = 8, bg = "white")
        cat("Saved to:", file.path(folder, "puma_prior_means.png"), "\n")
    }

    invisible(p)
}

plot_map_prior_mean_comuni <- function(save = FALSE, folder = "results/plots/",
                                       puma_dir = "input/counties-pumas",
                                       input_dir = "input/Comuni/",
                                       id_col = "COD_MUN",
                                       unit_ids = NULL) {
    cat("Computing map of mean income per comune...\n")

    if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Package 'sf' is required.")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Package 'dplyr' is required.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required.")
    }

    shp <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp) == 0) {
        stop("No .shp file found in '", puma_dir, "'.")
    }

    # Load data (unchanged)
    data <- NULL
    rds_file <- paste0(input_dir, "full_dataset.rds")
    if (file.exists(rds_file)) {
        data <- readRDS(rds_file)
    }
    if (is.null(data)) {
        csv_file <- paste0(input_dir, "full_dataset.csv")
        if (file.exists(csv_file)) {
            data <- read.csv(csv_file)
        }
    }
    if (is.null(data)) {
        dat_file <- paste0(input_dir, "full_dataset.dat")
        if (file.exists(dat_file)) {
            data <- tryCatch(
                readRDS(dat_file),
                error = function(e) {
                    env <- new.env()
                    load(dat_file, envir = env)
                    get(ls(env)[1], envir = env)
                }
            )
        }
    }
    if (is.null(data)) {
        stop("Could not load full_dataset from ", input_dir)
    }

    if (is.null(unit_ids)) {
        unit_ids <- data[[id_col]] # Assume first matching column or specify
    }

    # Define income bins and midpoints (€) for IRPEF brackets [memory:5]
    brackets <- c(0, 0, 10000, 15000, 26000, 55000, 75000, 120000, Inf)
    midpoints <- c(0, 5000, 12500, 20500, 40750, 65250, 97500, 150000)

    # Identify frequency columns (X* pattern from head(data))
    income_cols <- grep("^X", names(data), value = TRUE)
    if (length(income_cols) != length(midpoints)) {
        warning(
            "Number of income columns (", length(income_cols),
            ") doesn't match expected bins (", length(midpoints),
            "). Check column names."
        )
    }

    # Compute per-municipality mean income: sum(count * midpoint) / sum(count)
    n_bins <- min(length(income_cols), length(midpoints))
    total_income <- rowSums(sapply(1:n_bins, function(i) {
        as.numeric(data[[income_cols[i]]]) * midpoints[i]
    }))
    total_people <- rowSums(as.matrix(data[income_cols[1:n_bins]]))

    prior_means <- ifelse(total_people > 0, total_income / total_people, 0)
    names(prior_means) <- as.character(unit_ids)

    cat("Mean income range:", round(range(prior_means, na.rm = TRUE), 0), "€\n")
    cat("Global mean:", round(mean(prior_means, na.rm = TRUE), 0), "€\n")

    # Spatial join (unchanged)
    geom <- sf::st_read(shp[1], quiet = TRUE)
    prior_df <- tibble::tibble(
        !!id_col := as.character(unit_ids),
        prior_mean = prior_means
    )
    geom[[id_col]] <- as.character(geom[[id_col]])
    geom <- dplyr::left_join(geom, prior_df, by = id_col)

    p <- ggplot2::ggplot(geom) +
        ggplot2::geom_sf(aes(fill = prior_mean), color = "grey60", size = 0.2) +
        ggplot2::scale_fill_viridis_c(
            option = "viridis", na.value = "lightgrey",
            labels = scales::comma_format(accuracy = 1e3)
        ) +
        ggplot2::labs(
            # title = "Italian Comuni by Mean Income",
            fill = "Mean Income (€)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    print(p)

    if (save) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        ggplot2::ggsave(file.path(folder, "comuni_prior_mean_income.png"),
            p,
            width = 12, height = 10, dpi = 300
        )
    }

    invisible(p)
}

plot_hist_cls_pumas <- function(results, BI, input_dir = "input/LA/", point_estimate = NULL,
                                save_bool = FALSE, folder = "results/plots/", log_bool = TRUE, unit_ids = NULL) {
    cat("Computing histograms of cluster assignments (Long Format)...\n")

    # 1. Load the CSV
    csv_file <- file.path(input_dir, "full_dataset.csv")
    raw_data <- read.csv(csv_file)

    # 2. TRANSFORM & SPLIT (The critical change)
    # Convert log-income to original scale first
    if (log_bool) {
        raw_data$income_val <- exp(raw_data$log_income)
    } else {
        raw_data$income_val <- raw_data$log_income
    }

    # This creates a list where each element is named after a PUMA ID (COD_PUMA)
    # and contains a vector of all income values for that PUMA.
    data_list <- split(raw_data$income_val, raw_data$COD_PUMA)

    cat("Successfully split data into", length(data_list), "PUMAs.\n")

    # 3. Handle Cluster Estimates
    if (is.null(point_estimate)) {
        point_estimate <- plot_cls_est(results, BI, save = FALSE)
    }

    unique_clusters <- sort(unique(point_estimate))
    n_clusters <- length(unique_clusters)
    cluster_colors <- get_cluster_colors(n_clusters)
    names(point_estimate) <- unit_ids

    # 4. Group by Cluster using Names
    data_split <- lapply(unique_clusters, function(cl) {
        # Get the PUMA IDs belonging to this cluster
        pumas_in_cl <- names(point_estimate)[point_estimate == cl]

        # Grab those PUMAs from our split list
        subset_data <- data_list[pumas_in_cl]

        # Remove any that didn't have data in the CSV
        subset_data <- subset_data[!sapply(subset_data, is.null)]
        return(subset_data)
    })
    names(data_split) <- as.character(unique_clusters)

    # 5. Plotting
    render_plots <- function() {
        n_rows <- ceiling(sqrt(n_clusters))
        n_cols <- ceiling(n_clusters / n_rows)
        graphics::par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

        for (cl in names(data_split)) {
            combined_vals <- unlist(data_split[[cl]], use.names = FALSE)

            if (length(combined_vals) == 0) {
                plot.new() # title(main = paste("Cluster", cl, "(Empty)"))
                next
            }

            # Calculate stats for the legend
            cl_mean <- mean(combined_vals, na.rm = TRUE)
            cl_sd <- sd(combined_vals, na.rm = TRUE)
            cl_n <- length(combined_vals)

            cl_color <- adjustcolor(cluster_colors[cl], alpha.f = 0.6)

            # Draw Histogram
            hist(combined_vals,
                breaks = 50, probability = TRUE, col = cl_color, border = "white",
                main = paste("Cluster", cl, "\n(", length(data_split[[cl]]), "PUMAs)"),
                xlab = "Income"
            )

            # Add Density Line
            if (length(unique(combined_vals)) > 1) lines(density(combined_vals, na.rm = TRUE), lwd = 2)

            # Add Statistics Legend on the right
            legend("topright",
                legend = c(
                    paste0("Mean: ", formatC(cl_mean, format = "f", big.mark = ",", digits = 0)),
                    paste0("SD:   ", formatC(cl_sd, format = "f", big.mark = ",", digits = 0)),
                    paste0("N:    ", formatC(cl_n, format = "d", big.mark = ","))
                ),
                bty = "n", # No box around legend
                cex = 0.8, # Slightly smaller text
                text.font = 2, # Bold text
                adj = c(0, 0.5)
            )
        }
    }

    render_plots()

    if (save_bool) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        png(file.path(folder, "cluster_histograms.png"), width = 2400, height = 1800, res = 300)
        render_plots()
        dev.off()
    }

    invisible(point_estimate)
}

plot_hist_cls_comuni <- function(results, BI, input_dir = "input/Comuni/", point_estimate = NULL, save = FALSE, folder = "results/plots/") {
    cat("Computing aggregated histograms of cluster assignments for Comuni...\n")

    # Load full_dataset.csv which contains histogram data (income brackets as columns)
    csv_file <- paste0(input_dir, "full_dataset.csv")
    if (!file.exists(csv_file)) {
        stop("Could not find full_dataset.csv in ", input_dir)
    }
    data <- read.csv(csv_file, stringsAsFactors = FALSE)
    cat("Loaded", nrow(data), "comuni\n")

    # Identify income bracket columns (columns starting with X followed by digits)
    all_cols <- colnames(data)
    income_columns <- all_cols[grep("^X[0-9]", all_cols)]

    # Remove X0. if X0.10000 exists (avoid double counting)
    if ("X0." %in% income_columns && "X0.10000" %in% income_columns) {
        income_columns <- income_columns[income_columns != "X0."]
    }

    if (length(income_columns) == 0) {
        stop("Could not detect income bracket columns in the data")
    }
    n_brackets <- length(income_columns)
    cat("Detected", n_brackets, "income brackets:", paste(income_columns, collapse = ", "), "\n")

    # Extract bin breaks from column names (lower bound of each bracket)
    bin_breaks_lower <- sapply(income_columns, function(col) {
        nums <- as.numeric(unlist(regmatches(col, gregexpr("[0-9]+", col))))
        if (length(nums) > 0) {
            return(nums[1])
        } else {
            return(NA)
        }
    })

    # Create midpoints for each bracket (one per income column)
    bin_mids <- numeric(n_brackets)
    for (i in seq_len(n_brackets)) {
        lower <- bin_breaks_lower[i]
        if (i < n_brackets) {
            upper <- bin_breaks_lower[i + 1]
        } else {
            # Last bracket: use 1.5x the lower bound as approximate upper
            upper <- lower * 1.5
        }
        bin_mids[i] <- (lower + upper) / 2
    }

    # Create better labels for x-axis with k suffix for thousands
    bracket_labels <- sapply(seq_along(bin_breaks_lower), function(i) {
        lower <- bin_breaks_lower[i]
        if (is.na(lower)) {
            return("?")
        }

        # Format with k for thousands
        if (lower >= 1000) {
            lower_label <- paste0(round(lower / 1000), "k")
        } else {
            lower_label <- as.character(lower)
        }

        # Add upper bound for clarity
        if (i < length(bin_breaks_lower)) {
            upper <- bin_breaks_lower[i + 1]
            if (upper >= 1000) {
                upper_label <- paste0(round(upper / 1000), "k")
            } else {
                upper_label <- as.character(upper)
            }
            return(paste0(lower_label, "-", upper_label))
        } else {
            return(paste0(lower_label, "+"))
        }
    })

    if (is.null(point_estimate)) {
        point_estimate <- plot_cls_est(results, BI, save = FALSE)
    }

    unique_clusters <- sort(unique(point_estimate))
    n_clusters <- length(unique_clusters)

    cat("\n=== Cluster Diagnostics ===\n")
    cat("Number of clusters found:", n_clusters, "\n")
    cat("Cluster sizes:\n")
    print(table(point_estimate))

    # Extract income counts matrix
    income_matrix <- as.matrix(data[, income_columns])

    draw_histograms <- function() {
        op <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(op), add = TRUE)

        # Get consistent colors for clusters
        cluster_colors <- get_cluster_colors(n_clusters)

        # Split into batches of max 9 clusters per plot
        max_plots_per_page <- 9
        n_batches <- ceiling(n_clusters / max_plots_per_page)

        for (batch_idx in seq_len(n_batches)) {
            # Determine which clusters go in this batch
            start_idx <- (batch_idx - 1) * max_plots_per_page + 1
            end_idx <- min(batch_idx * max_plots_per_page, n_clusters)
            batch_clusters <- unique_clusters[start_idx:end_idx]
            n_plots_in_batch <- length(batch_clusters)

            # Calculate layout for this batch
            n_rows <- ceiling(sqrt(n_plots_in_batch))
            n_cols <- ceiling(n_plots_in_batch / n_rows)
            # Increase bottom margin for better x-axis label spacing
            graphics::par(mfrow = c(n_rows, n_cols), mar = c(6, 4, 3, 1))

            for (cl in batch_clusters) {
                # Get comuni indices in this cluster
                cluster_idx <- which(point_estimate == cl)
                n_comuni <- length(cluster_idx)

                # Aggregate counts across all comuni in this cluster
                aggregated_counts <- colSums(income_matrix[cluster_idx, , drop = FALSE], na.rm = TRUE)
                total_count <- sum(aggregated_counts)

                # Compute density (normalized counts)
                density_vals <- aggregated_counts / total_count

                # Compute weighted mean income for this cluster
                mean_income <- sum(bin_mids * aggregated_counts) / total_count

                # Use consistent cluster colors with transparency
                cl_color <- adjustcolor(cluster_colors[as.character(cl)], alpha.f = 0.7)

                # Create barplot (histogram-style)
                bp <- barplot(density_vals,
                    names.arg = NULL,
                    main = paste("Cluster", cl, "\n(n =", n_comuni, "comuni)"),
                    xlab = "",
                    ylab = "Density",
                    col = cl_color,
                    border = "white",
                    space = 0,
                    xaxt = "n" # Suppress default x-axis
                )

                # Add x-axis labels with better spacing
                # Show every label if <= 8 brackets, otherwise show every other label
                if (n_brackets <= 8) {
                    label_indices <- seq_along(bracket_labels)
                } else {
                    label_indices <- seq(1, n_brackets, by = 2)
                }

                axis(1,
                    at = bp[label_indices], labels = bracket_labels[label_indices],
                    las = 2, cex.axis = 0.65, padj = 0.5
                )
                mtext("Income Bracket (€)", side = 1, line = 4.5, cex = 0.7)

                cat("\nCluster", cl, "statistics:\n")
                cat("  N comuni:", n_comuni, "\n")
                cat("  Total taxpayers (N):", total_count, "\n")
                cat("  Weighted mean income:", round(mean_income, 0), "\n")

                legend("topright",
                    legend = c(
                        paste("Mean:", format(round(mean_income, 0), big.mark = ",")),
                        paste("N:", format(total_count, big.mark = ","))
                    ),
                    bty = "n",
                    cex = 0.8
                )
            }

            # If saving and multiple batches, save each batch separately
            if (save && n_batches > 1) {
                if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
                out_file <- file.path(folder, paste0("cluster_histograms_batch_", batch_idx, ".png"))
                if (requireNamespace("ragg", quietly = TRUE)) {
                    ragg::agg_png(filename = out_file, width = 2400, height = 1800, res = 300)
                } else {
                    png_type <- if (capabilities("cairo")) "cairo" else if (capabilities("X11")) "Xlib" else "cairo"
                    grDevices::png(filename = out_file, width = 2400, height = 1800, res = 300, type = png_type)
                }

                # Redraw the plots for this batch
                graphics::par(mfrow = c(n_rows, n_cols), mar = c(6, 4, 3, 1))
                for (cl in batch_clusters) {
                    cluster_idx <- which(point_estimate == cl)
                    n_comuni <- length(cluster_idx)
                    aggregated_counts <- colSums(income_matrix[cluster_idx, , drop = FALSE], na.rm = TRUE)
                    total_count <- sum(aggregated_counts)
                    density_vals <- aggregated_counts / total_count
                    mean_income <- sum(bin_mids * aggregated_counts) / total_count
                    cl_color <- adjustcolor(cluster_colors[as.character(cl)], alpha.f = 0.7)

                    bp <- barplot(density_vals,
                        names.arg = NULL,
                        main = paste("Cluster", cl, "\n(n =", n_comuni, "comuni)"),
                        xlab = "",
                        ylab = "Density",
                        col = cl_color,
                        border = "white",
                        space = 0,
                        xaxt = "n"
                    )

                    # Show labels based on number of brackets
                    if (n_brackets <= 8) {
                        label_indices <- seq_along(bracket_labels)
                    } else {
                        label_indices <- seq(1, n_brackets, by = 2)
                    }

                    axis(1,
                        at = bp[label_indices], labels = bracket_labels[label_indices],
                        las = 2, cex.axis = 0.65, padj = 0.5
                    )
                    mtext("Income Bracket (€)", side = 1, line = 4.5, cex = 0.7)

                    legend("topright",
                        legend = c(
                            paste("Mean:", format(round(mean_income, 0), big.mark = ",")),
                            paste("N:", format(total_count, big.mark = ","))
                        ),
                        bty = "n",
                        cex = 0.8
                    )
                }

                grDevices::dev.off()
                cat("\nSaved batch", batch_idx, "of", n_batches, "to", out_file, "\n")
            }
        }
    }

    draw_histograms()

    # Only save here if we have <= 9 clusters (single file case)
    device_opened <- FALSE
    if (save && n_clusters <= 9) {
        if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
        out_file <- file.path(folder, "cluster_histograms.png")
        if (requireNamespace("ragg", quietly = TRUE)) {
            ragg::agg_png(filename = out_file, width = 2400, height = 1800, res = 300)
        } else {
            png_type <- if (capabilities("cairo")) "cairo" else if (capabilities("X11")) "Xlib" else "cairo"
            grDevices::png(filename = out_file, width = 2400, height = 1800, res = 300, type = png_type)
        }
        device_opened <- TRUE
        on.exit(
            {
                if (device_opened) grDevices::dev.off()
            },
            add = TRUE
        )

        draw_histograms()

        grDevices::dev.off()
        device_opened <- FALSE
    }

    invisible(point_estimate)
}
