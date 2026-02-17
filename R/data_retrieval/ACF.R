# Script to download and process ACF data with complete covariates
# Simplified version of the original mechanism in SPMIX-applications by TeoGiane
# Author: Filippo Galli

# ========== PUMS VARIABLES DEFINITION (COMPLETE SCHEMA) ==========

pums_variables_definition <- list(
    # ========== CONTINUOUS VARIABLES ==========
    continuous = c(
        "ADJHSG", "ADJINC", # Adjustment factors
        "FINCP", "HINCP", "INTP", "OIP", "PAP", "PERNP", "PINCP",
        "RETP", "SEMP", "SSIP", "SSP", "WAGP", # Income
        "CONP", "ELEP", "FULP", "GASP", "GRNTP", "INSP", "MHP",
        "MRGP", "RNTP", "SMP", "SMOCP", "VALP", "WATP", # Housing costs
        "JWMNP", "WKHP", "AGEP", "CITWP", "YOEP", # Time/Years
        "BDSP", "NOC", "NRC", "NPF", "RMSP", "RACNUM" # Counts
    ),

    # ========== CATEGORICAL NOMINAL ==========
    categorical_nominal = c(
        "RT", "ST", "PUMA", "DIVISION", "REGION", "POBP", "POWPUMA", "POWSP",
        "MIGPUMA", "MIGSP", "WAOB", # Geography
        "HHT", "HHL", "RELP", "SFN", "SFR", "MSP", "PAOC", "QTRBIR", # Household
        "TYPE", "HFL", "RESMODE", "TEN", "VACS", # Housing
        "SEX", "MAR", "CIT", "MIG", "NATIVITY", "NOP", "ANC", "ANC1P", "ANC2P",
        "HISP", "RAC1P", "RAC2P", "RAC3P", "LANP", # Demographics
        "SCH", "COW", "ESR", "INDP", "NAICSP", "OCCP", "SOCP", "FOD1P", "FOD2P",
        "JWTR", "JWAP", "JWDP", "MIL", "VPS", # Employment
        "FHINS3C", "FHINS4C", "FHINS5C" # Insurance
    ),

    # ========== CATEGORICAL ORDINAL ==========
    categorical_ordinal = c(
        "ACR", "BLD", "YBL", "MV", "VEH", "TAXP", # Housing
        "AGS", "GRPIP", "OCPIP", "POVPIP", # Income ratios
        "SCHG", "SCHL", "ENG", "NP", "R18", "R60", "R65", "FPARC",
        "HUPAC", "HUPAOC", "HUPARC", "GCM", "MARHT", # Household
        "WKW", "WKL", "WKEXREL", "WORKSTAT", "WIF", "JWRIP", "DRIVESP",
        "DECADE", "DRAT" # Work/Disability
    ),

    # ========== BINARY (YES/NO) ==========
    binary = c(
        "BATH", "BUS", "KIT", "PLM", "REFR", "RWAT", "SINK", "STOV", "TEL", "TOIL", # Housing
        "HUGCL", "LNGI", "MULTG", "NPP", "NR", "PSF", "PARTNER", "RNTM", "SRNT", "SVAL", # Household
        "DDRS", "DEAR", "DEYE", "DOUT", "DPHY", "DREM", "DIS", "DRATX", # Health
        "HICOV", "PRIVCOV", "PUBCOV",
        paste0("HINS", 1:7), # Insurance types
        "FER", "LANX", "GCL", "GCR", # Demographics
        "RACAIAN", "RACASN", "RACBLK", "RACNHPI", "RACSOR", "RACWHT", "RC", "OC", # Race/ethnicity
        "FS", "MRGI", "MRGT", "MRGX", "SMX", # Mortgage/Financial
        "SCIENGP", "SCIENGRLP", "WRK", "NWAB", "NWAV", "NWLA", "NWLK", "NWRE", # Work
        "MARHD", "MARHM", "MARHW", # Marital
        paste0("MLPA", "MLPB", "MLPC", "MLPD", "MLPE", "MLPF", "MLPG", "MLPH", "MLPI", "MLPJ", "MLPK") # Military
    ),

    # ========== ALLOCATION FLAGS ==========
    allocation_flags = c(
        "FACRP", "FAGSP", "FBATHP", "FBDSP", "FBLDP", "FBUSP",
        "FCONP", "FELEP", "FFSP", "FFULP", "FGASP", "FHFLP",
        "FINSP", "FKITP", "FMHP", "FMRGIP", "FMRGP", "FMRGTP",
        "FMRGXP", "FMVP", "FPLMP", "FREFRP", "FRMSP", "FRNTMP",
        "FRNTP", "FRWATP", "FSINKP", "FSMP", "FSMXHP", "FSMXSP",
        "FSTOVP", "FTAXP", "FTELP", "FTENP", "FTOILP", "FVACSP",
        "FVALP", "FVEHP", "FWATP", "FYBLP",
        # Person flags
        "FAGEP", "FANCP", "FCITP", "FCITWP", "FCOWP", "FDDRSP",
        "FDEARP", "FDEYEP", "FDOUTP", "FDPHYP", "FDRATP", "FDRATXP",
        "FDREMP", "FENGP", "FESRP", "FFERP", "FFODP", "FGCLP",
        "FGCMP", "FGCRP", "FHINS1P", "FHINS2P", "FHINS3P", "FHINS4P",
        "FHINS5P", "FHINS6P", "FHINS7P", "FHISP", "FINDP", "FINTP",
        "FJWDP", "FJWMNP", "FJWRIP", "FJWTRP", "FLANP", "FLANXP",
        "FMARP", "FMARHDP", "FMARHMP", "FMARHTP", "FMARHWP", "FMARHYP",
        "FMIGP", "FMIGSP", "FMILPP", "FMILSP", "FOCCP", "FOIP",
        "FPAP", "FPOBP", "FPOWSP", "FRACP", "FRELP", "FRETP",
        "FSCHGP", "FSCHLP", "FSCHP", "FSEMP", "FSEXP", "FSSIP",
        "FSSP", "FWAGP", "FWKHP", "FWKLP", "FWKWP", "FWRKP",
        "FYOEP"
    ),

    # ========== WEIGHTS ==========
    weights = c(
        paste0("WGTP", 1:80), # Housing weight replicates
        paste0("PWGTP", 1:80), # Person weight replicates
        "WGTP", "PWGTP" # Main weights
    )
)

# ========== ALL PUMS VARIABLES (FLATTENED) ==========
all_pums_vars <- unique(unlist(pums_variables_definition))
cat("Total PUMS variables defined:", length(all_pums_vars), "\n")

# ========== Configuration ==========

library(argparser)

parser <- arg_parser("Simplified Census Data Workflow - Enhanced with Full Covariates")

parser <- add_argument(parser, "--states",
    default = "LA",
    help = "Comma-separated states or regions (West, Midwest, South, Northeast, USA, LA)"
)

parser <- add_argument(parser, "--num-datasets",
    type = "integer", default = 1,
    help = "Number of subsampled datasets to generate"
)

parser <- add_argument(parser, "--subsample-size",
    type = "integer", default = -1,
    help = "Number of observations per PUMA in each subsample (use -1 for all observations)"
)

parser <- add_argument(parser, "--cores",
    type = "integer", default = 1,
    help = "Number of cores for parallel download (default: 1)"
)

parser <- add_argument(parser, "--output-dir",
    default = "input/LA/",
    help = "Output directory for processed data"
)

parser <- add_argument(parser, "--keep-raw",
    flag = TRUE,
    help = "Keep raw downloaded files (default: delete after processing)"
)

args <- parse_args(parser)

# ========== Load packages ==========

cat("\n=== Loading packages ===\n")

suppressMessages({
    library(sf)
    library(dplyr)
    library(spdep)
    library(parallel)
})

# ========== State/Region mappings ==========

STATE_FIPS <- c(
    # West
    "California" = "06", "Nevada" = "32", "Oregon" = "41", "Washington" = "53",
    "Arizona" = "04", "Idaho" = "16", "Utah" = "49", "Montana" = "30",
    "Wyoming" = "56", "Colorado" = "08", "New Mexico" = "35", "Alaska" = "02", "Hawaii" = "15",
    # Midwest
    "North Dakota" = "38", "South Dakota" = "46", "Nebraska" = "31", "Kansas" = "20",
    "Minnesota" = "27", "Iowa" = "19", "Missouri" = "29", "Wisconsin" = "55",
    "Illinois" = "17", "Michigan" = "26", "Indiana" = "18", "Ohio" = "39",
    # South
    "Texas" = "48", "Oklahoma" = "40", "Arkansas" = "05", "Louisiana" = "22",
    "Mississippi" = "28", "Alabama" = "01", "Tennessee" = "47", "Kentucky" = "21",
    "West Virginia" = "54", "Virginia" = "51", "North Carolina" = "37", "South Carolina" = "45",
    "Georgia" = "13", "Florida" = "12",
    # Northeast
    "Maryland" = "24", "Delaware" = "10", "District of Columbia" = "11", "Pennsylvania" = "42",
    "New Jersey" = "34", "New York" = "36", "Connecticut" = "09", "Rhode Island" = "44",
    "Massachusetts" = "25", "Vermont" = "50", "New Hampshire" = "33", "Maine" = "23"
)

STATE_ABBREV <- c(
    "06" = "ca", "32" = "nv", "41" = "or", "53" = "wa", "04" = "az", "16" = "id",
    "49" = "ut", "30" = "mt", "56" = "wy", "08" = "co", "35" = "nm", "02" = "ak", "15" = "hi",
    "38" = "nd", "46" = "sd", "31" = "ne", "20" = "ks", "27" = "mn", "19" = "ia",
    "29" = "mo", "55" = "wi", "17" = "il", "26" = "mi", "18" = "in", "39" = "oh",
    "48" = "tx", "40" = "ok", "05" = "ar", "22" = "la", "28" = "ms", "01" = "al",
    "47" = "tn", "21" = "ky", "54" = "wv", "51" = "va", "37" = "nc", "45" = "sc",
    "13" = "ga", "12" = "fl", "24" = "md", "10" = "de", "11" = "dc", "42" = "pa",
    "34" = "nj", "36" = "ny", "09" = "ct", "44" = "ri", "25" = "ma", "50" = "vt",
    "33" = "nh", "23" = "me"
)

REGIONS <- list(
    West = c(
        "California", "Nevada", "Oregon", "Washington", "Arizona", "Idaho",
        "Utah", "Montana", "Wyoming", "Colorado", "New Mexico", "Alaska", "Hawaii"
    ),
    Midwest = c(
        "North Dakota", "South Dakota", "Nebraska", "Kansas", "Minnesota",
        "Iowa", "Missouri", "Wisconsin", "Illinois", "Michigan", "Indiana", "Ohio"
    ),
    South = c(
        "Texas", "Oklahoma", "Arkansas", "Louisiana", "Mississippi", "Alabama",
        "Tennessee", "Kentucky", "West Virginia", "Virginia", "North Carolina",
        "South Carolina", "Georgia", "Florida"
    ),
    Northeast = c(
        "Maryland", "Delaware", "District of Columbia", "Pennsylvania",
        "New Jersey", "New York", "Connecticut", "Rhode Island", "Massachusetts",
        "Vermont", "New Hampshire", "Maine"
    ),
    LA = c("California"),
    USA = names(STATE_FIPS)
)

# ========== Parse states argument ==========

parse_states <- function(states_string) {
    items <- trimws(strsplit(states_string, ",")[[1]])
    selected <- c()
    filter_la <- FALSE

    for (item in items) {
        region_match <- which(toupper(names(REGIONS)) == toupper(item))
        if (length(region_match) > 0) {
            region_name <- names(REGIONS)[region_match[1]]
            selected <- c(selected, REGIONS[[region_name]])
            cat(sprintf(" Region '%s' → %d states\n", region_name, length(REGIONS[[region_name]])))
            if (toupper(region_name) == "LA") {
                filter_la <- TRUE
            }
        } else {
            selected <- c(selected, item)
        }
    }
    list(states = unique(selected), filter_la_bay_area = filter_la)
}

parsed_result <- parse_states(args$states)
states_to_process <- parsed_result$states
filter_la_bay_area <- parsed_result$filter_la_bay_area

cat("\n=== Selected States ===\n")
cat(sprintf("Total: %d states\n", length(states_to_process)))
cat(paste(states_to_process, collapse = ", "), "\n")
if (filter_la_bay_area) {
    cat("Note: Will filter to LA Bay Area (93 PUMAs in LA, Orange, Ventura counties)\n")
}

# Validate states
invalid_states <- states_to_process[!states_to_process %in% names(STATE_FIPS)]
if (length(invalid_states) > 0) {
    stop(sprintf("Invalid states: %s", paste(invalid_states, collapse = ", ")))
}

# ========== Create directories ==========

dir.create("raw", showWarnings = FALSE, recursive = TRUE)
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# ========== Download data ==========
cat("\n=== Downloading Data ===\n")

download_with_progress <- function(url, destfile, quiet = FALSE) {
    if (!quiet) cat(sprintf(" Downloading: %s\n", basename(url)))

    # Try aria2c if available (faster split downloading)
    if (nzchar(Sys.which("aria2c"))) {
        dir_name <- dirname(destfile)
        file_name <- basename(destfile)
        
        # -x 4: max 4 connections
        # -s 4: split to 4 parts
        # --file-allocation=none: faster startup
        cmd <- sprintf("aria2c -x 4 -s 4 --allow-overwrite=true --file-allocation=none -d '%s' -o '%s' '%s' %s", 
                       dir_name, file_name, url, if(quiet) "-q" else "")
        
        status <- system(cmd, ignore.stdout = quiet, ignore.stderr = quiet)
        if (status == 0 && file.exists(destfile)) {
             if (!quiet) cat(sprintf(" ✓ Downloaded (aria2c): %.2f MB\n", file.size(destfile) / 1024^2))
             return(TRUE)
        }
    }

    # Fallback to standard download.file
    result <- tryCatch(
        {
            download.file(url, destfile, mode = "wb", quiet = quiet)
            TRUE
        },
        error = function(e) {
            if (!quiet) cat(sprintf(" ✗ Failed: %s\n", e$message))
            FALSE
        }
    )

    if (result && file.exists(destfile)) {
        if (!quiet) cat(sprintf(" ✓ Downloaded: %.2f MB\n", file.size(destfile) / 1024^2))
    }
    return(result)
}

# Download shapefile
shapefile_url <- "https://usa.ipums.org/usa/resources/volii/shapefiles/ipums_puma_2010_tl20.zip"
shapefile_zip <- "raw/ipums_puma_2010_tl20.zip"

cat("\n[1/2] Shapefile\n")
if (download_with_progress(shapefile_url, shapefile_zip)) {
    unzip(shapefile_zip, exdir = "raw")
    cat(" ✓ Extracted\n")
}
# Download state data
cat(sprintf("\n[2/2] Census data for %d state(s)\n", length(states_to_process)))

downloaded_states <- c()
failed_states <- c()

if (args$cores > 1) {
    cat(sprintf("Parallel download enabled: %d cores\n", args$cores))
}

process_state_download <- function(i) {
    state <- states_to_process[i]
    fips <- STATE_FIPS[state]
    abbrev <- STATE_ABBREV[fips]
    url <- sprintf("https://www2.census.gov/programs-surveys/acs/experimental/2020/data/pums/1-Year/csv_p%s.zip", abbrev)
    destfile <- sprintf("raw/csv_p%s.zip", abbrev)

    # For parallel, we suppress individual progress bars and start/end messages to avoid clutter
    is_parallel <- args$cores > 1
    
    if (!is_parallel) {
        cat(sprintf("\n [%d/%d] %s\n", i, length(states_to_process), state))
    }

    if (download_with_progress(url, destfile, quiet = is_parallel)) {
        unzip(destfile, exdir = "raw")
        if (!is_parallel) cat(" ✓ Extracted\n")
        return(list(state = state, success = TRUE))
    } else {
        return(list(state = state, success = FALSE))
    }
}

if (args$cores > 1) {
    results <- mclapply(seq_along(states_to_process), process_state_download, mc.cores = args$cores)
    
    for (res in results) {
        if (!is.null(res) && !inherits(res, "try-error")) {
            if (res$success) {
                downloaded_states <- c(downloaded_states, res$state)
                cat(sprintf(" ✓ Processed: %s\n", res$state))
            } else {
                failed_states <- c(failed_states, res$state)
                cat(sprintf(" ✗ Failed: %s\n", res$state))
            }
        } else {
            cat(" ✗ Error in parallel worker\n")
        }
    }
} else {
    for (i in seq_along(states_to_process)) {
        res <- process_state_download(i)
        if (res$success) {
            downloaded_states <- c(downloaded_states, res$state)
        } else {
            failed_states <- c(failed_states, res$state)
        }
    }
}

if (length(failed_states) > 0) {
    warning(sprintf("Failed to download: %s", paste(failed_states, collapse = ", ")))
    states_to_process <- downloaded_states
}

if (length(downloaded_states) == 0) {
    stop("No states were successfully downloaded. Cannot proceed.")
}

cat(sprintf(
    "\n✓ Successfully downloaded %d/%d states\n",
    length(downloaded_states), length(states_to_process) + length(failed_states)
))

# ========== Load and process data ==========

cat("\n=== Processing Data ===\n")

# Load shapefile
cat("Loading shapefile... ")
sf_usa <- read_sf("raw/ipums_puma_2010_tl20.shp")
cat(sprintf("✓ (%d PUMAs)\n", nrow(sf_usa)))

# Load census data
cat("Loading census data...\n")

raw_data_list <- list()

for (i in seq_along(states_to_process)) {
    state <- states_to_process[i]
    fips <- STATE_FIPS[state]
    csv_file <- sprintf("raw/psam_p%s.csv", fips)

    if (!file.exists(csv_file)) {
        warning(sprintf("Data file not found: %s", csv_file))
        next
    }

    cat(sprintf(" [%d/%d] %s... ", i, length(states_to_process), state))
    data <- read.csv(csv_file)
    cat(sprintf("%d records\n", nrow(data)))
    raw_data_list[[i]] <- data
}

if (length(raw_data_list) == 0) {
    stop("No data files could be loaded. Cannot proceed.")
}

# Combine and clean data
cat("Combining and cleaning... ")

census_data <- bind_rows(raw_data_list) %>%
    mutate(
        ST = sprintf("%02d", ST),
        PUMA = sprintf("%05d", PUMA),
        STATE_PUMA = paste0(ST, "_", PUMA)
    ) %>%
    filter(PINCP > 0) %>%
    mutate(LPINCP = log(PINCP))

cat(sprintf("✓ (%d observations)\n", nrow(census_data)))

# ========== ENHANCED: Prepare census data with ALL variables ==========

cat("\n=== Enriching Census Data with Complete Variable Set ===\n")

# Identify which variables are present in census_data
present_vars <- colnames(census_data)
cat(sprintf("Variables present in raw data: %d\n", length(present_vars)))

# Identify missing variables
missing_vars <- setdiff(all_pums_vars, present_vars)
cat(sprintf("Variables missing from raw data: %d\n", length(missing_vars)))

if (length(missing_vars) > 0) {
    cat(sprintf("  Adding %d NA columns for missing covariates...\n", length(missing_vars)))
    for (var in missing_vars) {
        census_data[[var]] <- NA
    }
}

# Reorder columns to have all PUMS variables first, then any extras
col_order <- c(
    intersect(all_pums_vars, colnames(census_data)),
    setdiff(colnames(census_data), all_pums_vars)
)

census_data <- census_data[, col_order]

cat(sprintf("✓ Census data now has %d columns\n", ncol(census_data)))

# Prepare geographic data
cat("Preparing geographic boundaries... ")

selected_fips <- STATE_FIPS[states_to_process]

sf_counties <- sf_usa %>%
    mutate(
        ST = if ("STATEFP" %in% colnames(.)) {
            sprintf("%02d", as.numeric(STATEFP))
        } else if ("STATEFP10" %in% colnames(.)) {
            sprintf("%02d", as.numeric(STATEFP10))
        } else {
            STATE_FIPS[State]
        },
        PUMA = sprintf("%05d", as.numeric(PUMA)),
        STATE_PUMA = paste0(ST, "_", PUMA)
    ) %>%
    filter(!is.na(ST), State %in% states_to_process) %>%
    select(ST, PUMA, STATE_PUMA, State, Name, geometry)

cat(sprintf("✓ (%d PUMAs)\n", nrow(sf_counties)))

# Filter for LA Bay Area if requested
if (filter_la_bay_area) {
    cat("Filtering for LA Bay Area (Los Angeles, Orange, Ventura counties)... ")

    la_pumas <- sprintf("06_%05d", 3701:3769)
    orange_pumas <- sprintf("06_%05d", 5901:5918)
    ventura_pumas <- sprintf("06_%05d", 11101:11106)
    la_bay_area_pumas <- c(la_pumas, orange_pumas, ventura_pumas)

    sf_counties <- sf_counties %>%
        filter(STATE_PUMA %in% la_bay_area_pumas)

    census_data <- census_data %>%
        filter(STATE_PUMA %in% la_bay_area_pumas)

    cat(sprintf("✓ (%d LA Bay Area PUMAs)\n", length(la_bay_area_pumas)))
}

# Match data and geography
cat("Matching data and geography... ")

common_pumas <- intersect(sf_counties$STATE_PUMA, unique(census_data$STATE_PUMA))

if (length(common_pumas) == 0) {
    stop("No matching PUMAs found between shapefile and data!")
}

sf_counties <- sf_counties %>%
    filter(STATE_PUMA %in% common_pumas) %>%
    arrange(ST, PUMA)

census_data <- census_data %>%
    filter(STATE_PUMA %in% common_pumas)

cat(sprintf("✓ (%d matched PUMAs)\n", length(common_pumas)))

# Split data by PUMA
data_by_puma <- census_data %>%
    group_by(STATE_PUMA) %>%
    group_split() %>%
    setNames(sf_counties$STATE_PUMA)

# ========== Compute adjacency matrix ==========

cat("Computing spatial adjacency matrix... ")

adj_list <- poly2nb(sf_counties, queen = FALSE, snap = 0.01)
isolated <- which(card(adj_list) == 0)

if (length(isolated) > 0) {
    cat(sprintf("\n Warning: %d isolated PUMAs (no neighbors)\n", length(isolated)))
}

W <- nb2mat(adj_list, style = "B", zero.policy = TRUE)
rownames(W) <- sf_counties$STATE_PUMA
colnames(W) <- sf_counties$STATE_PUMA

cat("✓\n")

# Validate subsample size
min_size <- min(sapply(data_by_puma, nrow))

if (args$subsample_size != -1 && args$subsample_size > min_size) {
    stop(sprintf(
        "Subsample size (%d) exceeds minimum PUMA size (%d)",
        args$subsample_size, min_size
    ))
}

# ========== Generate datasets ==========

if (args$subsample_size == -1) {
    cat("\nNote: Subsample size set to -1. Generating 1 dataset with ALL observations.\n")
    args$num_datasets <- 1
}

cat(sprintf("\n=== Generating Datasets ===\n"))

# Save full dataset as data (original format for compatibility)
data <- lapply(data_by_puma, function(df) {
    df$LPINCP
})

# Save full dataset as CSV (long format)
cat("Saving full dataset as CSV... ")
full_data_long <- do.call(rbind, lapply(names(data_by_puma), function(puma_id) {
    data.frame(
        COD_PUMA = puma_id,
        NAME_PUMA = paste(sf_counties$State[sf_counties$STATE_PUMA == puma_id],
            sf_counties$Name[sf_counties$STATE_PUMA == puma_id],
            sep = " - "
        ),
        log_income = data_by_puma[[puma_id]]$LPINCP,
        stringsAsFactors = FALSE
    )
}))

write.csv(full_data_long, file.path(args$output_dir, "full_dataset.csv"), row.names = FALSE)
cat("✓\n")

# ========== CREATE ENHANCED full_dataset_covariates ==========

cat("\nBuilding and saving full_dataset_covariates with ALL PUMS variables...\n")

# Create list of dataframes, one per PUMA, with ALL PUMS variables
full_dataset_covariates <- list()

for (puma_id in names(data_by_puma)) {
    # Get all observations for this PUMA from census_data
    puma_subset <- census_data[census_data$STATE_PUMA == puma_id, ]

    # Select only PUMS variables (ensures consistency across PUMAs)
    puma_covariates <- puma_subset[, intersect(all_pums_vars, colnames(puma_subset)), drop = FALSE]

    # Add any missing PUMS variables as NA columns
    for (var in setdiff(all_pums_vars, colnames(puma_covariates))) {
        puma_covariates[[var]] <- NA
    }

    # Reorder to standard order
    puma_covariates <- puma_covariates[, intersect(all_pums_vars, colnames(puma_covariates)), drop = FALSE]

    full_dataset_covariates[[puma_id]] <- puma_covariates
}

names(full_dataset_covariates) <- names(data_by_puma)

# Verify structure
cat("\nVerification:\n")
cat(sprintf("  Number of PUMAs: %d\n", length(full_dataset_covariates)))
cat(sprintf("  Columns per PUMA: %d\n", ncol(full_dataset_covariates[[1]])))
cat(sprintf("  Observations per PUMA: %s\n", paste(
    sapply(full_dataset_covariates, nrow),
    collapse = " | "
)))

# Save as RDS (list format for LA.R)
saveRDS(
    full_dataset_covariates,
    file.path(args$output_dir, "full_dataset_covariates.rds")
)

cat(sprintf("✓ Saved full_dataset_covariates.rds (list format)\n"))

# Save as CSV (long format with COD_PUMA column)
cat("Saving full_dataset_covariates as CSV... ")
full_dataset_covariates_df <- do.call(rbind, lapply(names(full_dataset_covariates), function(puma_id) {
    df <- full_dataset_covariates[[puma_id]]
    df$COD_PUMA <- puma_id
    # Put COD_PUMA first, then all other columns
    df <- df[, c("COD_PUMA", setdiff(colnames(df), "COD_PUMA"))]
    return(df)
}))

write.csv(full_dataset_covariates_df, file.path(args$output_dir, "full_dataset_covariates.csv"), row.names = FALSE)
cat(sprintf("✓ (%d rows, %d columns)\n", nrow(full_dataset_covariates_df), ncol(full_dataset_covariates_df)))

# Generate subsampled datasets
set.seed(230196)

cat(sprintf(
    "Generating %d subsampled datasets (n=%s per PUMA)...\n",
    args$num_datasets, if (args$subsample_size == -1) "ALL" else args$subsample_size
))

for (i in 1:args$num_datasets) {
    subsample_list <- lapply(names(data_by_puma), function(puma_id) {
        if (args$subsample_size == -1) {
            sampled_values <- data_by_puma[[puma_id]]$LPINCP
        } else {
            sampled_values <- data_by_puma[[puma_id]]$LPINCP[
                sample(1:nrow(data_by_puma[[puma_id]]), args$subsample_size)
            ]
        }

        data.frame(
            COD_PUMA = puma_id,
            NAME_PUMA = paste(sf_counties$State[sf_counties$STATE_PUMA == puma_id],
                sf_counties$Name[sf_counties$STATE_PUMA == puma_id],
                sep = " - "
            ),
            log_income = sampled_values,
            stringsAsFactors = FALSE
        )
    })

    subsample_df <- do.call(rbind, subsample_list)
    filename <- file.path(args$output_dir, sprintf("data_%03d.csv", i))
    write.csv(subsample_df, filename, row.names = FALSE)

    if (i %% 5 == 0 || i == args$num_datasets) {
        cat(sprintf(" Progress: %d/%d datasets\n", i, args$num_datasets))
    }
}

# ========== Save spatial data ==========

cat("\nSaving spatial data...\n")

geometry_sf <- sf_counties %>%
    select(STATE_PUMA, ST, PUMA, State, Name, geometry) %>%
    rename(
        COD_PUMA = STATE_PUMA,
        COD_STATE = ST,
        COD_PUMA_ONLY = PUMA,
        NAME_STATE = State,
        NAME_PUMA = Name
    )

geometry_dir <- file.path(args$output_dir, "counties-pumas")
dir.create(geometry_dir, showWarnings = FALSE, recursive = TRUE)

st_write(geometry_sf, file.path(geometry_dir, "counties-pumas.shp"),
    append = FALSE, quiet = TRUE
)
cat(" ✓ Shapefile\n")

saveRDS(geometry_sf, file.path(args$output_dir, "geometry.rds"))
cat(" ✓ Geometry RDS\n")

saveRDS(W, file.path(args$output_dir, "adj_matrix.rds"))
cat(" ✓ Adjacency matrix RDS\n")

# ========== Cleanup ==========

if (!args$keep_raw) {
    cat("\nCleaning up raw files... ")
    unlink("raw", recursive = TRUE)
    cat("✓\n")
}

# ========== Summary ==========

cat("\n")
cat(strrep("=", 70), "\n", sep = "")
cat("✓ WORKFLOW COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat(sprintf(
    "States processed: %d (%s)\n",
    length(states_to_process), paste(states_to_process, collapse = ", ")
))

if (filter_la_bay_area) {
    cat("Region: LA Bay Area (Los Angeles, Orange, Ventura counties)\n")
}

cat(sprintf("PUMAs: %d\n", nrow(geometry_sf)))
cat(sprintf("Total observations: %s\n", format(nrow(census_data), big.mark = ",")))
cat(sprintf("Datasets generated: %d\n", args$num_datasets))

cat("\n")
cat("Output files:\n")
cat(sprintf(" %s/full_dataset.csv (all observations in long format)\n", args$output_dir))
cat(sprintf(" %s/full_dataset_covariates.rds (COMPLETE covariates as list, for LA.R)\n", args$output_dir))
cat(sprintf(" %s/full_dataset_covariates.csv (COMPLETE covariates as CSV)\n", args$output_dir))
cat(sprintf(" %s/data_001.csv ... data_%03d.csv (subsampled observations)\n", args$output_dir, args$num_datasets))
cat(sprintf(" %s/geometry.rds\n", args$output_dir))
cat(sprintf(" %s/adj_matrix.rds\n", args$output_dir))
cat(sprintf(" %s/counties-pumas/counties-pumas.shp\n", args$output_dir))

cat("\n")
cat("CSV Format for covariates:\n")
cat(" - COD_PUMA: Unique PUMA identifier (e.g., 06_03701)\n")
cat(" - All other columns: PUMS variables as defined in pums_variables_definition\n")
cat("\n")
cat("RDS Format for LA.R:\n")
cat(" - List object where names are COD_PUMA values\n")
cat(" - Each element is a dataframe containing all PUMS variables for that PUMA\n")