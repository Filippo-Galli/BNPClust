setwd("~/Documents/BNPClust/R")

source("utils.R")
source("utils_plot.R")

##############################################
# Cardinality of the PUMA ====
##############################################

plot_map_data_cardinality <- function(
    save = FALSE,
    folder = "../results/plots/",
    puma_dir = "../input/LA/counties-pumas",
    input_dir = "../input/LA/",
    id_col = "COD_PUMA",
    simplify_tol = NULL,
    rotation_deg = 0,
    zoom_to_data = TRUE,
    margin_frac = 0.05
) {
    cat("Computing cardinality of PUMAS ...\n")

    # --- Check Dependencies ---
    pkgs <- c("sf", "dplyr", "ggplot2", "scales")
    missing_pkgs <- pkgs[
        !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ",
            paste(missing_pkgs, collapse = ", ")
        )
    }

    # --- 1. Load Data ---
    csv_file <- file.path(input_dir, "full_dataset.csv")
    rds_file <- file.path(input_dir, "full_dataset.rds")

    if (file.exists(rds_file)) {
        data <- readRDS(rds_file)
    } else if (file.exists(csv_file)) {
        data <- read.csv(csv_file, stringsAsFactors = FALSE)
    } else {
        stop("Could not find full_dataset in ", input_dir)
    }

    if (!id_col %in% names(data)) {
        stop("Column '", id_col, "' not found in dataset.")
    }

    # Drop NAME_PUMA if present
    data <- data %>%
        dplyr::select(-dplyr::matches("NAME_PUMA", ignore.case = TRUE))

    # --- 2. Compute Cardinalities ---
    cat("Calculating cardinalities grouped by PUMA...\n")

    cardinalities <- data %>%
        dplyr::mutate(!!id_col := as.character(.data[[id_col]])) %>%
        dplyr::group_by(.data[[id_col]]) %>%
        dplyr::summarise(
            cardinality = mean(exp(log_income), na.rm = TRUE), # dplyr::n(),
            .groups = "drop"
        )

    prior_df <- cardinalities %>%
        dplyr::transmute(
            !!id_col := as.character(.data[[id_col]]),
            cardinality = cardinality
        )

    # --- 3. Load Geometry ---
    cat("Loading map geometry...\n")
    shp_files <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) {
        stop("No .shp file found in ", puma_dir)
    }

    geom <- sf::st_read(shp_files[1], quiet = TRUE)

    if (!id_col %in% names(geom)) {
        stop("Column '", id_col, "' not found in shapefile.")
    }

    geom[[id_col]] <- as.character(geom[[id_col]])

    if (
        !is.null(simplify_tol) && is.numeric(simplify_tol) && simplify_tol > 0
    ) {
        geom <- sf::st_simplify(geom, dTolerance = simplify_tol)
    }

    # --- 4. Optional rotation ---
    # This rotates the plotted geometry, which is what you need if you want
    # the top edge to appear visually horizontal.
    if (is.numeric(rotation_deg) && rotation_deg != 0) {
        theta <- rotation_deg * pi / 180
        rot <- matrix(
            c(cos(theta), sin(theta), -sin(theta), cos(theta)),
            nrow = 2,
            byrow = TRUE
        )

        geom_rot <- geom %>%
            dplyr::mutate(geometry_rot = sf::st_geometry(.) * rot) %>%
            sf::st_drop_geometry() %>%
            dplyr::rename(geometry = geometry_rot) %>%
            sf::st_as_sf()

        sf::st_crs(geom_rot) <- sf::st_crs(geom)
        geom <- geom_rot
    }

    # --- 5. Join data to geometry ---
    map_data <- dplyr::left_join(geom, prior_df, by = id_col)

    # Use only observed units for framing if requested
    map_data_main <- if (isTRUE(zoom_to_data)) {
        dplyr::filter(map_data, !is.na(cardinality))
    } else {
        map_data
    }

    if (nrow(map_data_main) == 0) {
        stop(
            "No matching map polygons after join. Check ID formats in data and shapefile."
        )
    }

    # --- 6. Bounding box for centered framing ---
    bbox <- sf::st_bbox(map_data_main)
    x_margin <- margin_frac * (bbox["xmax"] - bbox["xmin"])
    y_margin <- margin_frac * (bbox["ymax"] - bbox["ymin"])

    lbl_fun <- scales::label_comma()

    # --- 7. Plot ---
    cat("Generating plot...\n")

    p <- ggplot2::ggplot(map_data_main) +
        ggplot2::geom_sf(
            ggplot2::aes(fill = cardinality),
            color = "white",
            linewidth = 0.05
        ) +
        ggplot2::scale_fill_viridis_c(
            option = "viridis",
            na.value = "lightgrey",
            labels = lbl_fun,
            name = "Income ($)"
        ) +
        ggplot2::coord_sf(
            xlim = c(bbox["xmin"] - x_margin, bbox["xmax"] + x_margin),
            ylim = c(bbox["ymin"] - y_margin, bbox["ymax"] + y_margin),
            expand = FALSE,
        ) +
        ggplot2::theme_light() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            plot.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            legend.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            legend.box.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            panel.border = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(color = "grey60", size = 8)
        )

    if (save) {
        if (!dir.exists(folder)) {
            dir.create(folder, recursive = TRUE)
        }

        out_file <- file.path(folder, "puma_prior_mean.png")

        ggplot2::ggsave(
            filename = out_file,
            plot = p,
            width = 10,
            height = 8,
            bg = "transparent"
        )

        cat("Saved to:", out_file, "\n")
    }

    invisible(p)
}

plot_map_prior_mean_mun <- function(
    save = FALSE,
    folder = "../results/plots/",
    puma_dir = "../input/municipalities/geometry",
    input_dir = "../input/municipalities/",
    id_col = "COD_MUN",
    unit_ids = NULL,
    plot_cardinality = FALSE,
    cap_quantile = 0.98,
    simplify_tol = NULL,
    rotation_deg = 0,
    zoom_to_data = TRUE,
    margin_frac = 0.05
) {
    cat(
        if (plot_cardinality) {
            "Computing cardinality of municipalities ...\n"
        } else {
            "Computing map of mean income per municipality...\n"
        }
    )

    # --- Check Dependencies ---
    pkgs <- c("sf", "dplyr", "ggplot2", "scales")
    missing_pkgs <- pkgs[
        !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ",
            paste(missing_pkgs, collapse = ", ")
        )
    }

    # --- 1. Load Data ---
    rds_file <- file.path(input_dir, "full_dataset.rds")
    csv_file <- file.path(input_dir, "full_dataset.csv")
    dat_file <- file.path(input_dir, "full_dataset.dat")

    if (file.exists(rds_file)) {
        data <- readRDS(rds_file)
    } else if (file.exists(csv_file)) {
        data <- read.csv(csv_file, stringsAsFactors = FALSE)
    } else if (file.exists(dat_file)) {
        data <- tryCatch(
            readRDS(dat_file),
            error = function(e) {
                env <- new.env()
                load(dat_file, envir = env)
                get(ls(env)[1], envir = env)
            }
        )
    } else {
        stop("Could not find full_dataset in ", input_dir)
    }

    if (!id_col %in% names(data) && is.null(unit_ids)) {
        stop("Column '", id_col, "' not found in dataset.")
    }

    if (is.null(unit_ids)) {
        unit_ids <- data[[id_col]]
    }

    # --- 2. Compute Fill Variable ---
    income_cols <- grep(
        "^X0|^X10000|^X15000|^X26000|^X55000|^X75000|^X120000",
        names(data),
        value = TRUE
    )
    # Safer/more general: just reuse the same bracket-detecting logic
    # used below, since both branches need these columns.
    income_cols <- grep("^X", names(data), value = TRUE)

    if (plot_cardinality) {
        cat("Calculating cardinalities per municipality...\n")

        if (length(income_cols) == 0) {
            stop("No bracket columns (matching '^X') found in dataset.")
        }

        cardinality <- rowSums(
            sapply(income_cols, function(col) as.numeric(data[[col]]))
        )

        cat("Cardinality range:", range(cardinality, na.rm = TRUE), "\n")
        cat("Total records:", sum(cardinality, na.rm = TRUE), "\n")

        cap_value <- NULL
        if (!is.null(cap_quantile)) {
            cap_value <- stats::quantile(
                cardinality,
                cap_quantile,
                na.rm = TRUE
            )
            n_capped <- sum(cardinality > cap_value, na.rm = TRUE)
            if (n_capped > 0) {
                cat(
                    "Capping",
                    n_capped,
                    "municipalities above",
                    scales::label_comma()(round(cap_value)),
                    "(",
                    scales::label_percent()(1 - cap_quantile),
                    "of comuni)\n"
                )
            }
            cardinality_display <- pmin(cardinality, cap_value)
        } else {
            cardinality_display <- cardinality
        }

        prior_df <- data.frame(
            id = as.character(data[[id_col]]),
            value = cardinality_display,
            stringsAsFactors = FALSE
        )
        names(prior_df)[1] <- id_col

        legend_name <- "Numerosities"
        lbl_fun <- if (!is.null(cap_value)) {
            function(x) {
                labs <- scales::label_comma()(x)
                labs[x >= cap_value] <- paste0(
                    scales::label_comma()(round(cap_value)),
                    "+"
                )
                labs
            }
        } else {
            scales::label_comma()
        }
        out_name <- "comuni_cardinality.png"
    } else {
        cat("Calculating mean income per municipality...\n")

        # Income bins and midpoints (€) for IRPEF brackets
        midpoints <- c(0, 5000, 12500, 20500, 40750, 65250, 97500, 150000)

        if (length(income_cols) != length(midpoints)) {
            warning(
                "Number of income columns (",
                length(income_cols),
                ") doesn't match expected bins (",
                length(midpoints),
                "). Check column names."
            )
        }

        n_bins <- min(length(income_cols), length(midpoints))
        total_income <- rowSums(sapply(1:n_bins, function(i) {
            as.numeric(data[[income_cols[i]]]) * midpoints[i]
        }))
        total_people <- rowSums(as.matrix(data[income_cols[1:n_bins]]))

        prior_means <- ifelse(total_people > 0, total_income / total_people, 0)

        cat(
            "Mean income range:",
            round(range(prior_means, na.rm = TRUE), 0),
            "€\n"
        )
        cat("Global mean:", round(mean(prior_means, na.rm = TRUE), 0), "€\n")

        prior_df <- data.frame(
            id = as.character(
                if (is.null(unit_ids)) data[[id_col]] else unit_ids
            ),
            value = prior_means,
            stringsAsFactors = FALSE
        )
        names(prior_df)[1] <- id_col

        legend_name <- "Income (€)"
        lbl_fun <- scales::label_comma(accuracy = 1e3)
        out_name <- "comuni_prior_mean_income.png"
    }

    # --- 3. Load Geometry ---
    cat("Loading map geometry...\n")
    shp_files <- list.files(puma_dir, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) {
        stop("No .shp file found in ", puma_dir)
    }

    geom <- sf::st_read(shp_files[1], quiet = TRUE)

    if (!id_col %in% names(geom)) {
        stop("Column '", id_col, "' not found in shapefile.")
    }

    geom[[id_col]] <- as.character(geom[[id_col]])

    if (
        !is.null(simplify_tol) && is.numeric(simplify_tol) && simplify_tol > 0
    ) {
        geom <- sf::st_simplify(geom, dTolerance = simplify_tol)
    }

    # --- 4. Optional rotation ---
    if (is.numeric(rotation_deg) && rotation_deg != 0) {
        theta <- rotation_deg * pi / 180
        rot <- matrix(
            c(cos(theta), sin(theta), -sin(theta), cos(theta)),
            nrow = 2,
            byrow = TRUE
        )

        geom_rot <- geom %>%
            dplyr::mutate(geometry_rot = sf::st_geometry(.) * rot) %>%
            sf::st_drop_geometry() %>%
            dplyr::rename(geometry = geometry_rot) %>%
            sf::st_as_sf()

        sf::st_crs(geom_rot) <- sf::st_crs(geom)
        geom <- geom_rot
    }

    # --- 5. Join data to geometry ---
    map_data <- dplyr::left_join(geom, prior_df, by = id_col)

    map_data_main <- if (isTRUE(zoom_to_data)) {
        dplyr::filter(map_data, !is.na(value))
    } else {
        map_data
    }

    if (nrow(map_data_main) == 0) {
        stop(
            "No matching map polygons after join. Check ID formats in data and shapefile."
        )
    }

    # --- 6. Bounding box for centered framing ---
    bbox <- sf::st_bbox(map_data_main)
    x_margin <- margin_frac * (bbox["xmax"] - bbox["xmin"])
    y_margin <- margin_frac * (bbox["ymax"] - bbox["ymin"])

    # --- 7. Plot ---
    cat("Generating plot...\n")

    p <- ggplot2::ggplot(map_data_main) +
        ggplot2::geom_sf(
            ggplot2::aes(fill = value),
            color = "white",
            linewidth = 0.05
        ) +
        ggplot2::scale_fill_viridis_c(
            option = "viridis",
            na.value = "lightgrey",
            labels = lbl_fun,
            name = legend_name
        ) +
        ggplot2::coord_sf(
            xlim = c(bbox["xmin"] - x_margin, bbox["xmax"] + x_margin),
            ylim = c(bbox["ymin"] - y_margin, bbox["ymax"] + y_margin),
            expand = FALSE,
        ) +
        ggplot2::theme_light() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            plot.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            legend.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            legend.box.background = ggplot2::element_rect(
                fill = "transparent",
                color = NA
            ),
            panel.border = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(color = "grey60", size = 8)
        )

    if (save) {
        if (!dir.exists(folder)) {
            dir.create(folder, recursive = TRUE)
        }

        out_file <- file.path(folder, out_name)

        ggplot2::ggsave(
            filename = out_file,
            plot = p,
            width = 10,
            height = 8,
            bg = "transparent"
        )

        cat("Saved to:", out_file, "\n")
    }

    invisible(p)
}

plot_map_data_cardinality(save = TRUE, rotation_deg = 13) # rotation_deg = 13 for LA, 5 for USA
# plot_map_prior_mean_mun(save = TRUE, plot_cardinality = TRUE)

national_bin_frequency <- function(
    input_dir = "../input/municipalities/",
    bin_labels = c(
        "0-",
        "0-10000",
        "10000-15000",
        "15000-26000",
        "26000-55000",
        "55000-75000",
        "75000-120000",
        "120000+"
    )
) {
    # --- Load Data ---
    rds_file <- file.path(input_dir, "full_dataset.rds")
    csv_file <- file.path(input_dir, "full_dataset.csv")
    dat_file <- file.path(input_dir, "full_dataset.dat")

    if (file.exists(rds_file)) {
        data <- readRDS(rds_file)
    } else if (file.exists(csv_file)) {
        data <- read.csv(csv_file, stringsAsFactors = FALSE)
    } else if (file.exists(dat_file)) {
        data <- tryCatch(
            readRDS(dat_file),
            error = function(e) {
                env <- new.env()
                load(dat_file, envir = env)
                get(ls(env)[1], envir = env)
            }
        )
    } else {
        stop("Could not find full_dataset in ", input_dir)
    }

    # --- Identify bracket columns ---
    income_cols <- grep("^X", names(data), value = TRUE)
    if (length(income_cols) == 0) {
        stop("No bracket columns (matching '^X') found in dataset.")
    }

    if (length(income_cols) != length(bin_labels)) {
        warning(
            "Number of income columns (",
            length(income_cols),
            ") doesn't match number of bin_labels (",
            length(bin_labels),
            "). Check column order/names."
        )
    }

    n_bins <- min(length(income_cols), length(bin_labels))

    # --- Sum each bracket column nationally ---
    national_totals <- sapply(income_cols[1:n_bins], function(col) {
        sum(as.numeric(data[[col]]), na.rm = TRUE)
    })

    total_all <- sum(national_totals)

    result <- data.frame(
        bin = bin_labels[1:n_bins],
        column = income_cols[1:n_bins],
        count = as.numeric(national_totals),
        share = as.numeric(national_totals) / total_all,
        stringsAsFactors = FALSE
    )

    cat("National total records:", scales::label_comma()(total_all), "\n")
    print(
        data.frame(
            bin = result$bin,
            count = scales::label_comma()(result$count),
            share = scales::label_percent(accuracy = 0.1)(result$share)
        )
    )

    invisible(result)
}
# freq <- national_bin_frequency()

# print(freq)
