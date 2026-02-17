source("R/utils_plot.R")
source("R/utils.R")

# Load required packages
library(dplyr)

##############################################################################
# CONFIGURATION ====
##############################################################################

# Configuration - modify these as needed
CONFIG <- list(
    input_dir = "input/LA",
    output_dir = "real_data/LA",
    results_dir = "results/",
    # Choose format: "rds" (list format, recommended) or "csv" (long format)
    covariates_format = "csv",
    # Minimum proportion of non-NA values to keep a column (0 = keep all, 1 = no NA allowed)
    na_threshold = 0.0
)

# Create output directory if it doesn't exist
if (!dir.exists(CONFIG$output_dir)) {
    dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
}

##############################################################################
# LOAD DATA ====
##############################################################################

cat("=== Loading Covariates Data ===\n")

# Load covariates - supports both RDS (list) and CSV (long) formats
if (CONFIG$covariates_format == "rds") {
    covariates_file <- file.path(CONFIG$input_dir, "full_dataset_covariates.rds")
    if (!file.exists(covariates_file)) {
        stop(sprintf("File not found: %s. Try setting covariates_format = 'csv' or run ACF.R first.", covariates_file))
    }

    # Load as list (one dataframe per PUMA)
    full_dataset_covariates <- readRDS(covariates_file)
    cat(sprintf("Loaded RDS format: %d PUMAs\n", length(full_dataset_covariates)))

    # Combine into single dataframe for processing, preserving COD_PUMA
    all_covariates <- do.call(rbind, lapply(names(full_dataset_covariates), function(puma_id) {
        df <- full_dataset_covariates[[puma_id]]
        df$COD_PUMA <- puma_id
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
    puma_ids <- unique(all_covariates$COD_PUMA)
    full_dataset_covariates <- lapply(puma_ids, function(id) {
        subset(all_covariates, COD_PUMA == id)
    })
    names(full_dataset_covariates) <- puma_ids
}

# Ensure COD_PUMA exists (for merging with spatial data later)
if (!"COD_PUMA" %in% colnames(all_covariates)) {
    if ("STATE_PUMA" %in% colnames(all_covariates)) {
        all_covariates$COD_PUMA <- all_covariates$STATE_PUMA
    } else {
        stop("Neither COD_PUMA nor STATE_PUMA column found in data")
    }
}

##############################################################################
# PUMS VARIABLES DEFINITION ====
##############################################################################

pums_variables <- list(
    # Continuous variables (numeric amounts, counts, time, years)
    continuous = c(
        "ADJHSG", "ADJINC", "FINCP", "HINCP", "INTP", "OIP", "PAP", "PERNP",
        "PINCP", "RETP", "SEMP", "SSIP", "SSP", "WAGP", "CONP", "ELEP", "FULP",
        "GASP", "GRNTP", "INSP", "MHP", "MRGP", "RNTP", "SMP", "SMOCP", "VALP",
        "WATP", "JWMNP", "WKHP", "AGEP", "CITWP", "YOEP", "BDSP", "NOC", "NRC",
        "NPF", "RMSP", "RACNUM", "LPINCP" # Include transformed income
    ),

    # Categorical nominal (unordered categories)
    categorical_nominal = c(
        "RT", "ST", "PUMA", "DIVISION", "REGION", "POBP", "POWPUMA", "POWSP",
        "MIGPUMA", "MIGSP", "WAOB", "HHT", "HHL", "RELP", "SFN", "SFR", "MSP",
        "PAOC", "QTRBIR", "TYPE", "HFL", "RESMODE", "TEN", "VACS", "SEX", "MAR",
        "CIT", "MIG", "NATIVITY", "NOP", "ANC", "ANC1P", "ANC2P", "HISP", "RAC1P",
        "RAC2P", "RAC3P", "LANP", "SCH", "COW", "ESR", "INDP", "NAICSP", "OCCP",
        "SOCP", "FOD1P", "FOD2P", "JWTR", "JWAP", "JWDP", "MIL", "VPS", "FHINS3C",
        "FHINS4C", "FHINS5C", "COD_PUMA" # PUMA identifier
    ),

    # Categorical ordinal (ordered levels)
    categorical_ordinal = c(
        "ACR", "BLD", "YBL", "MV", "VEH", "TAXP", "AGS", "GRPIP", "OCPIP",
        "POVPIP", "SCHG", "SCHL", "ENG", "NP", "R18", "R60", "R65", "FPARC",
        "HUPAC", "HUPAOC", "HUPARC", "GCM", "MARHT", "WKW", "WKL", "WKEXREL",
        "WORKSTAT", "WIF", "JWRIP", "DRIVESP", "DECADE", "DRAT"
    ),

    # Binary (0/1 or Yes/No)
    binary = c(
        "BATH", "BUS", "KIT", "PLM", "REFR", "RWAT", "SINK", "STOV", "TEL", "TOIL",
        "HUGCL", "LNGI", "MULTG", "NPP", "NR", "PSF", "PARTNER", "RNTM", "SRNT",
        "SVAL", "DDRS", "DEAR", "DEYE", "DOUT", "DPHY", "DREM", "DIS", "DRATX",
        "HICOV", "PRIVCOV", "PUBCOV", paste0("HINS", 1:7), "FER", "LANX", "GCL",
        "GCR", "RACAIAN", "RACASN", "RACBLK", "RACNHPI", "RACSOR", "RACWHT",
        "RC", "OC", "FS", "MRGI", "MRGT", "MRGX", "SMX", "SCIENGP", "SCIENGRLP",
        "WRK", "NWAB", "NWAV", "NWLA", "NWLK", "NWRE", "MARHD", "MARHM", "MARHW",
        paste0("MLP", c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"))
    ),

    # Allocation flags (binary indicators of imputation)
    allocation_flags = c(
        paste0("F", c(
            "ACRP", "AGSP", "BATHP", "BDSP", "BLDP", "BUSP", "CONP",
            "ELEP", "FSP", "FULP", "GASP", "HFLP", "INSP", "KITP", "MHP",
            "MRGIP", "MRGP", "MRGTP", "MRGXP", "MVP", "PLMP", "REFRP",
            "RMSP", "RNTMP", "RNTP", "RWATP", "SINKP", "SMP", "SMXHP",
            "SMXSP", "STOVP", "TAXP", "TELP", "TENP", "TOILP", "VACSP",
            "VALP", "VEHP", "WATP", "YBLP"
        )),
        paste0("F", c(
            "AGEP", "ANCP", "CITP", "CITWP", "COWP", "DDRSP", "DEARP",
            "DEYEP", "DOUTP", "DPHYP", "DRATP", "DRATXP", "DREMP", "ENGP",
            "ESRP", "FERP", "FODP", "GCLP", "GCMP", "GCRP", paste0("HINS", 1:7, "P"),
            "HISP", "INDP", "INTP", "JWDP", "JWMNP", "JWRIP", "JWTRP",
            "LANP", "LANXP", "MARP", "MARHDP", "MARHMP", "MARHTP", "MARHWP",
            "MARHYP", "MIGP", "MIGSP", "MILPP", "MILSP", "OCCP", "OIP",
            "PAP", "POBP", "POWSP", "RACP", "RELP", "RETP", "SCHGP",
            "SCHLP", "SCHP", "SEMP", "SEXP", "SSIP", "SSP", "WAGP",
            "WKHP", "WKLP", "WKWP", "WRKP", "YOEP"
        ))
    ),

    # Weights (special handling)
    weights = c("WGTP", "PWGTP", paste0("WGTP", 1:80), paste0("PWGTP", 1:80))
)

# Flatten all variable categories for reference
all_pums_vars <- unique(unlist(pums_variables))
cat(sprintf("Defined %d PUMS variables\n", length(all_pums_vars)))

##############################################################################
# COVARIATE PREPROCESSING ====
##############################################################################

cat("\n=== Preprocessing Covariates ===\n")

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

# 2. Convert categorical variables to factors
convert_to_factor <- function(df, var_list) {
    for (var in var_list) {
        if (var %in% colnames(df)) {
            # Handle both character and numeric categoricals
            df[[var]] <- as.factor(as.character(df[[var]]))
        }
    }
    return(df)
}

cat("Converting categorical variables to factors...\n")
all_covariates <- convert_to_factor(all_covariates, pums_variables$categorical_nominal)
all_covariates <- convert_to_factor(all_covariates, pums_variables$categorical_ordinal)
all_covariates <- convert_to_factor(all_covariates, pums_variables$binary)
all_covariates <- convert_to_factor(all_covariates, pums_variables$allocation_flags)

# 3. Ensure continuous variables are numeric
cat("Ensuring continuous variables are numeric...\n")
for (var in pums_variables$continuous) {
    if (var %in% colnames(all_covariates)) {
        all_covariates[[var]] <- as.numeric(as.character(all_covariates[[var]]))
    }
}

# 4. Remove columns with zero variance (constants)
variances <- sapply(all_covariates, function(x) {
    if (is.numeric(x)) {
        var(x, na.rm = TRUE)
    } else if (is.factor(x)) {
        length(levels(x)) - 1 # 0 if only 1 level, >0 otherwise
    } else {
        length(unique(x)) - 1
    }
})

constant_cols <- names(variances)[variances == 0]
if (length(constant_cols) > 0) {
    cat(sprintf(
        "Removing %d constant columns: %s\n",
        length(constant_cols), paste(constant_cols, collapse = ", ")
    ))
    all_covariates <- all_covariates[, !colnames(all_covariates) %in% constant_cols]
}

covariate_names <- colnames(all_covariates)
cat(sprintf("Final columns after cleaning: %d\n", length(covariate_names)))

##############################################################################
# GAUSSIANITY TRANSFORMATION FOR INCOME ====
##############################################################################

# if ("PINCP" %in% colnames(all_covariates)) {
#     cat("\n=== Transforming Income Variable ===\n")

#     lpincp_values <- all_covariates$PINCP

#     # Box-Cox transformation to approximate normality
#     library(MASS)

#     # Remove zeros and negatives for Box-Cox
#     positive_income <- lpincp_values[lpincp_values > 0 & !is.na(lpincp_values)]

#     if (length(positive_income) > 100) {
#         boxcox_result <- boxcox(positive_income ~ 1, lambda = seq(-2, 2, 0.1), plotit = FALSE)
#         best_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
#         cat(sprintf("Best Box-Cox lambda for PINCP: %.3f\n", best_lambda))

#         # Apply transformation
#         if (abs(best_lambda) < 0.01) {
#             # Approximate log transform
#             all_covariates$LPINCP <- log(pmax(all_covariates$PINCP, 1))
#             cat("Applied log transformation (lambda ≈ 0)\n")
#         } else {
#             all_covariates$LPINCP <- (all_covariates$PINCP^best_lambda - 1) / best_lambda
#             cat(sprintf("Applied Box-Cox transformation with λ = %.3f\n", best_lambda))
#         }

#         # Update continuous variables list if not already there
#         if (!"LPINCP" %in% pums_variables$continuous) {
#             pums_variables$continuous <- c(pums_variables$continuous, "LPINCP")
#         }
#     }
# } else {
#     cat("\nNote: PINCP column not found, skipping transformation\n")
# }

##############################################################################
# HELPER FUNCTIONS FOR COVARIATE EXTRACTION ====
##############################################################################

#' Aggregate covariates at PUMA level
#'
#' @param covariates_df Dataframe with individual-level data
#' @param puma_col Name of the PUMA ID column
#' @param vars Variables to aggregate (default: all numeric)
#' @return Dataframe with one row per PUMA
aggregate_puma_covariates <- function(covariates_df, puma_col = "COD_PUMA",
                                      vars = NULL) {
    if (is.null(vars)) {
        # Default to all numeric columns except weights and IDs
        vars <- setdiff(
            names(covariates_df),
            c(puma_col, "SERIALNO", "SPORDER", pums_variables$weights)
        )
        vars <- vars[sapply(covariates_df[, vars], is.numeric)]
    }

    # Calculate aggregations
    puma_summary <- covariates_df %>%
        group_by(across(all_of(puma_col))) %>%
        summarise(
            # Count
            n = n(),

            # For each numeric variable, calculate mean and sd
            across(
                all_of(vars),
                list(
                    mean = ~ mean(., na.rm = TRUE),
                    sd = ~ sd(., na.rm = TRUE),
                    median = ~ median(., na.rm = TRUE)
                ),
                .names = "{.col}_{.fn}"
            ),
            .groups = "drop"
        )

    return(puma_summary)
}

#' Extract modal category for categorical variables by PUMA
#'
#' @param covariates_df Dataframe with individual-level data
#' @param puma_col Name of the PUMA ID column
#' @param cat_vars Categorical variables to get modes for
#' @return Dataframe with modal categories
get_puma_modes <- function(covariates_df, puma_col = "COD_PUMA", cat_vars = NULL) {
    if (is.null(cat_vars)) {
        cat_vars <- intersect(
            names(covariates_df),
            c(
                pums_variables$categorical_nominal,
                pums_variables$categorical_ordinal,
                pums_variables$binary
            )
        )
    }

    get_mode <- function(x) {
        if (all(is.na(x))) {
            return(NA)
        }
        ux <- unique(x[!is.na(x)])
        ux[which.max(tabulate(match(x[!is.na(x)], ux)))]
    }

    mode_df <- covariates_df %>%
        group_by(across(all_of(puma_col))) %>%
        summarise(
            across(all_of(cat_vars), get_mode, .names = "{.col}_mode"),
            .groups = "drop"
        )

    return(mode_df)
}

#' Get proportion for binary/categorical variables by PUMA
#'
#' @param covariates_df Dataframe with individual-level data
#' @param var Variable to calculate proportion for (should be binary or have reference_level)
#' @param reference_level Level to calculate proportion of (for factors)
#' @param puma_col Name of the PUMA ID column
get_puma_proportion <- function(covariates_df, var, reference_level = NULL,
                                puma_col = "COD_PUMA") {
    if (!var %in% names(covariates_df)) {
        stop(sprintf("Variable %s not found in data", var))
    }

    if (is.factor(covariates_df[[var]]) || is.character(covariates_df[[var]])) {
        if (is.null(reference_level)) {
            # Use first level or most common
            reference_level <- ifelse(is.factor(covariates_df[[var]]),
                levels(covariates_df[[var]])[1],
                sort(unique(covariates_df[[var]]))[1]
            )
        }

        prop_df <- covariates_df %>%
            group_by(across(all_of(puma_col))) %>%
            summarise(
                !!paste0(var, "_prop_", reference_level) := mean(.data[[var]] == reference_level, na.rm = TRUE),
                .groups = "drop"
            )
    } else {
        # For numeric (assumed binary 0/1)
        prop_df <- covariates_df %>%
            group_by(across(all_of(puma_col))) %>%
            summarise(
                !!paste0(var, "_prop") := mean(.data[[var]], na.rm = TRUE),
                .groups = "drop"
            )
    }

    return(prop_df)
}

#' Get summary with also standard deviation
summ_with_sd <- function(x) {
    c(summary(x), sd = sd(x, na.rm = TRUE))
}

##############################################################################
# EXTRACT KEY COVARIATES (AGEP & SEX) ====
##############################################################################

cat("\n=== Extracting Key Covariates by PUMA ===\n")

# AGEP Analysis
if ("AGEP" %in% colnames(all_covariates)) {
    cat("Processing AGEP (Age)...\n")

    # Standardize AGEP
    agep_stats <- list(
        mean = mean(all_covariates$AGEP, na.rm = TRUE),
        sd = sd(all_covariates$AGEP, na.rm = TRUE)
    )

    all_covariates$AGEP_std <- scale(all_covariates$AGEP,
        center = agep_stats$mean,
        scale = agep_stats$sd
    )[, 1]

    # Aggregate by PUMA
    puma_age_stats <- all_covariates %>%
        group_by(COD_PUMA) %>%
        summarise(
            AGEP_mean = mean(AGEP, na.rm = TRUE),
            AGEP_sd = sd(AGEP, na.rm = TRUE),
            AGEP_std_mean = mean(AGEP_std, na.rm = TRUE),
            AGEP_std_sd = sd(AGEP_std, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(COD_PUMA)

    # Check for ties
    agep_ties <- puma_age_stats %>%
        group_by(AGEP_mean) %>%
        filter(n() > 1)

    if (nrow(agep_ties) > 0) {
        cat(sprintf("  Warning: %d ties found in mean AGEP values\n", nrow(agep_ties)))
    }

    # Summary statistics
    cat("Summary of AGEP by PUMA:\n")
    print(summ_with_sd(puma_age_stats$AGEP_std_mean))
    print(summ_with_sd(puma_age_stats$AGEP_std_sd))

    # Save
    saveRDS(puma_age_stats, file = file.path(CONFIG$output_dir, "puma_age_stats.rds"))
    write.csv(puma_age_stats,
        file = file.path(CONFIG$output_dir, "puma_age_stats.csv"),
        row.names = FALSE
    )
    cat(sprintf("  Saved age statistics for %d PUMAs\n", nrow(puma_age_stats)))

    # Store for later use
    puma_agep_std_mean <- puma_age_stats
}

# SEX Analysis
if ("SEX" %in% colnames(all_covariates)) {
    cat("Processing SEX...\n")

    if (!is.factor(all_covariates$SEX)) {
        # Handle both numeric 1/2 and character formats
        all_covariates$SEX <- factor(as.numeric(as.character(all_covariates$SEX)),
            levels = c(1, 2),
            labels = c("Male", "Female")
        )
    }

    # Verify conversion
    cat(sprintf(
        "  SEX levels after conversion: %s\n",
        paste(levels(all_covariates$SEX), collapse = ", ")
    ))
    cat(sprintf(
        "  Distribution: %s\n",
        paste(table(all_covariates$SEX, useNA = "ifany"), collapse = ", ")
    ))

    get_mode <- function(x) {
        if (all(is.na(x))) {
            return(NA_character_)
        }
        ux <- unique(x[!is.na(x)])
        ux[which.max(tabulate(match(x[!is.na(x)], ux)))]
    }

    puma_sex_stats <- all_covariates %>%
        group_by(COD_PUMA) %>%
        summarise(
            n = n(),
            SEX_mode = get_mode(SEX),
            SEX_mode_explained = ifelse(get_mode(SEX) == 1, "Male", "Female"),
            Perc_Female = table(SEX, useNA = "ifany")[[2]] / length(SEX),
            .groups = "drop"
        ) %>%
        arrange(COD_PUMA)

    # Convert to numeric coding: 0 = Male, 1 = Female
    puma_sex_stats$SEX_mode <- ifelse(puma_sex_stats$SEX_mode_explained == "Female", 1, 0)

    # Validation warning
    if (all(puma_sex_stats$Perc_Female == 0, na.rm = TRUE)) {
        cat("  WARNING: All female proportions are zero! Check SEX encoding.\n")
    }

    # Save outputs
    saveRDS(puma_sex_stats, file = file.path(CONFIG$output_dir, "puma_sex_stats.rds"))
    write.csv(puma_sex_stats,
        file = file.path(CONFIG$output_dir, "puma_sex_stats.csv"),
        row.names = FALSE
    )
    cat(sprintf("  Saved sex statistics for %d PUMAs\n", nrow(puma_sex_stats)))

    # Preview first few rows
    print(head(puma_sex_stats[, c("COD_PUMA", "SEX_mode", "SEX_mode_explained", "Perc_Female")], 10))

    sex_puma_data <- puma_sex_stats

    print(summ_with_sd(puma_sex_stats$Perc_Female))
}

##############################################################################
# ADDITIONAL COVARIATE EXTRACTION ====
##############################################################################

cat("\n=== Extracting Additional Demographics ===\n")

# Extract race/ethnicity proportions if available
race_vars <- intersect(
    c("RAC1P", "HISP", "RACWHT", "RACBLK", "RACASN", "RACAIAN"),
    colnames(all_covariates)
)

if (length(race_vars) > 0) {
    cat(sprintf("Processing race variables: %s\n", paste(race_vars, collapse = ", ")))

    # For binary race indicators (RACWHT, RACBLK, etc.)
    binary_race <- intersect(race_vars, pums_variables$binary)
    if (length(binary_race) > 0) {
        race_props_list <- lapply(binary_race, function(var) {
            get_puma_proportion(all_covariates, var, reference_level = "1", puma_col = "COD_PUMA")
        })
        race_props <- Reduce(function(x, y) merge(x, y, by = "COD_PUMA", all = TRUE), race_props_list)

        saveRDS(race_props, file = file.path(CONFIG$output_dir, "puma_race_proportions.rds"))
        write.csv(race_props,
            file = file.path(CONFIG$output_dir, "puma_race_proportions.csv"),
            row.names = FALSE
        )
        cat(sprintf("  Saved race proportions for %d PUMAs\n", nrow(race_props)))

        sapply(race_props[, grep("^RAC", names(race_props))], summ_with_sd)
    }
}

# Education level
if ("SCHL" %in% colnames(all_covariates)) {
    cat("Processing education levels (SCHL)...\n")

    schl_map <- c(
        "bb" = NA_real_, # N/A (<3yo) → exclude from analysis
        "01" = 1, # No schooling
        "02" = 2, # Nursery/preschool
        "03" = 3, # Kindergarten
        "04" = 4, "05" = 5, "06" = 6, "07" = 7, "08" = 8, "09" = 9,
        "10" = 10, "11" = 11, "12" = 12, "13" = 13, "14" = 14, "15" = 15,
        "16" = 16, # Regular HS diploma
        "17" = 17, # GED/alternative
        "18" = 18, # Some college <1yr
        "19" = 19, # 1+yr college, no degree
        "20" = 20, # Associate's
        "21" = 21, # Bachelor's
        "22" = 22, # Master's
        "23" = 23, # Professional
        "24" = 24 # Doctorate
    )
    # Convert SCHL to numeric (safe for factor/character)
    all_covariates$SCHL_num <- case_when(
        is.na(all_covariates$SCHL) ~ NA_real_,
        TRUE ~ as.numeric(schl_map[as.character(all_covariates$SCHL)])
    )
    
    # PUMA-level aggregation (identical to race code structure)
    edu_stats <- all_covariates %>%
        filter(!is.na(SCHL_num)) %>%  # Drop invalid
        group_by(COD_PUMA) %>%
        summarise(
            SCHL_median = median(SCHL_num, na.rm = TRUE),
            SCHL_mean = mean(SCHL_num, na.rm = TRUE),
            HS_plus_prop = mean(SCHL_num >= 16, na.rm = TRUE),      # HS diploma+
            College_prop = mean(SCHL_num >= 18, na.rm = TRUE),      # Some college+
            Bachelor_plus_prop = mean(SCHL_num >= 21, na.rm = TRUE), # Bachelor's+
            .groups = "drop"
        ) %>%
        mutate(across(-COD_PUMA, ~ round(.x, 4)))  # Round for stability
    
    saveRDS(edu_stats, file = file.path(CONFIG$output_dir, "puma_education_stats.rds"))
    write.csv(edu_stats, file = file.path(CONFIG$output_dir, "puma_education_stats.csv"), row.names = FALSE)
    cat(sprintf("  Saved PUMA education stats (%d PUMAs): median/mean SCHL, HS+/College+/Bachelor+ props\n", nrow(edu_stats)))
    
    # Print summary like race code
    sapply(edu_stats[, -1], summ_with_sd)
}

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
    file_chosen <- files[4]
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
    states <- if (any(parts == "Comuni")) "Comuni" else "LA"
    id_col <- if (states == "Comuni") "COD_MUN" else "COD_PUMA"

    # Load unit IDs from shapefile
    shp_path <- file.path(CONFIG$input_dir, "counties-pumas", "counties-pumas.shp")
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

    # Add state prefix for LA PUMAs if needed
    if (states != "Comuni" && !grepl("^06_", cluster_df$id[1])) {
        cluster_df$id <- paste0("06_", cluster_df$id)
    }

    # Merge with covariate data
    if (exists("sex_puma_data")) {
        sex_df <- merge(cluster_df, sex_puma_data, by.x = "id", by.y = "COD_PUMA", all.x = TRUE)

        mode_func <- function(x) {
        if (length(x) == 0) return(NA)
        tbl <- table(x)
        names(tbl)[which.max(tbl)]
        }

        sex_cluster_summary <- sex_df %>%
            group_by(cluster) %>%
            summarise(
                n_pumas = n(),
                mode_sex = mode_func(SEX_mode),
                avg_perc_female = mean(Perc_Female, na.rm = TRUE),
                sd_perc_female = sd(Perc_Female, na.rm = TRUE),
                .groups = "drop"
            )

        cat("\nSex Distribution by Cluster:\n")
        print(sex_cluster_summary)

        write.csv(sex_cluster_summary,
            file = file.path(analysis_output_dir, "cluster_sex_summary.csv"),
            row.names = FALSE
        )
    }

    if (exists("puma_agep_std_mean")) {
        age_df <- merge(cluster_df, puma_agep_std_mean, by.x = "id", by.y = "COD_PUMA", all.x = TRUE)

        age_cluster_summary <- age_df %>%
            group_by(cluster) %>%
            summarise(
                n_pumas = n(),
                mean_age = mean(AGEP_mean, na.rm = TRUE),
                sd_between_pumas = sd(AGEP_mean, na.rm = TRUE),
                avg_age_sd_within = mean(AGEP_sd, na.rm = TRUE),
                .groups = "drop"
            )

        cat("\nAge Distribution by Cluster:\n")
        print(age_cluster_summary)

        write.csv(age_cluster_summary,
            file = file.path(analysis_output_dir, "cluster_age_summary.csv"),
            row.names = FALSE
        )
    }

    return(list(
        sex = if (exists("sex_cluster_summary")) sex_cluster_summary else NULL,
        age = if (exists("age_cluster_summary")) age_cluster_summary else NULL
    ))
}

# Run cluster analysis only if results exist
if (dir.exists(CONFIG$results_dir) && length(list.files(CONFIG$results_dir)) > 0) {
    cluster_results <- analyze_clusters()
} else {
    cat("\nSkipping cluster analysis (no results directory)\n")
}

##############################################################################
# SUMMARY AND EXPORT ====
##############################################################################

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("COVARIATE PROCESSING COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat("Output files created:\n")
output_files <- list.files(CONFIG$output_dir, full.names = TRUE)
for (f in output_files) {
    cat(sprintf(" - %s (%.2f MB)\n", basename(f), file.size(f) / 1024^2))
}

cat("\n")
cat("Variables available in all_covariates:\n")
cat(sprintf(" - Total: %d\n", length(covariate_names)))
cat(sprintf(" - Continuous: %d\n", sum(covariate_names %in% pums_variables$continuous)))
cat(sprintf(" - Categorical: %d\n", sum(covariate_names %in% c(
    pums_variables$categorical_nominal,
    pums_variables$categorical_ordinal
))))
cat(sprintf(" - Binary: %d\n", sum(covariate_names %in% pums_variables$binary)))
cat("\nPUMAs processed:", length(unique(all_covariates$COD_PUMA)), "\n")

cat("\nExample usage:\n")
cat(" # Access individual level data:\n")
cat(" head(all_covariates[, c('COD_PUMA', 'AGEP', 'SEX', 'PINCP')])\n")
cat(" \n")
cat(" # Get PUMA-level aggregates:\n")
cat(" puma_means <- aggregate_puma_covariates(all_covariates, vars = c('AGEP', 'PINCP'))\n")
cat(" \n")
cat(" # Get modes for categorical variables:\n")
cat(" puma_modes <- get_puma_modes(all_covariates, cat_vars = c('SEX', 'RAC1P'))\n")
cat(" \n")
cat(" # Merge with cluster results:\n")
cat(" cluster_data <- merge(cluster_df, puma_age_stats, by.x = 'id', by.y = 'COD_PUMA')\n")
