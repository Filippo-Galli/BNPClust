# Script to download and process Italian municipality data
# Simplified version of the original mechanism in spatialCMC by TeoGiane
# Author: Filippo Galli

# Load required libraries ====
library(dplyr)
library(sf)
library(readxl)
library(tidyr)

# Create necessary directories ====
dir.create("raw", showWarnings = FALSE, recursive = TRUE)
dir.create("input/municipalities", showWarnings = FALSE, recursive = TRUE)

cat("Starting data download and processing...\n\n")

# ========== Download function ==========
download_with_progress <- function(url, destfile, quiet = FALSE) {
    if (!quiet) cat(sprintf(" Downloading: %s\n", basename(url)))
    
    if (nzchar(Sys.which("aria2c"))) {
        dir_name <- dirname(destfile)
        file_name <- basename(destfile)
        
        cmd <- sprintf("aria2c -x 4 -s 4 --allow-overwrite=true --file-allocation=none -d '%s' -o '%s' '%s' %s", 
                       dir_name, file_name, url, if(quiet) "-q" else "")
        
        status <- system(cmd, ignore.stdout = quiet, ignore.stderr = quiet)
        if (status == 0 && file.exists(destfile)) {
            if (!quiet) cat(sprintf(" ✓ Downloaded (aria2c): %.2f MB\n", file.size(destfile) / 1024^2))
            return(TRUE)
        }
    }
    
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

# Step 1: Download raw data ====
cat("Step 1: Downloading raw data files...\n")

download_with_progress(
  url = "https://www.istat.it/storage/cartografia/confini_amministrativi/non_generalizzati/Limiti01012020.zip",
  destfile = "raw/Limiti01012020.zip"
)

download_with_progress(
  url = "https://www.istat.it/storage/dati_mortalita/dicembre-2025/decessi-comunali-giornalieri-provvisori_4-18122025.zip",
  destfile = "raw/Dataset-decessi-comunali.zip"
)

download_with_progress(
  url = "https://demo.istat.it/data/posas/POSAS_2021_it_Comuni.zip",
  destfile = "raw/POSAS_2021.zip"
)

download_with_progress(
  url = "https://www1.finanze.gov.it/finanze/analisi_stat/public/v_4_0_0/contenuti/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.zip",
  destfile = "raw/Redditi_2020.zip"
)

cat("Downloads complete!\n\n")

# Step 2: Unzip files ====
cat("Step 2: Extracting files...\n")

unzip("raw/Limiti01012020.zip", exdir = "raw")
unzip("raw/Dataset-decessi-comunali.zip", exdir = "raw")
unzip("raw/POSAS_2021.zip", exdir = "raw")
unzip("raw/Redditi_2020.zip", exdir = "raw")

cat("Extraction complete!\n\n")

# Step 3: Import shapefiles ====
cat("Step 3: Importing geographic data (shapefiles)...\n")

mun_sf <- read_sf("raw/Limiti01012020/Com01012020/Com01012020_WGS84.shp")

cat("Imported", nrow(mun_sf), "municipalities\n\n")

# Step 4: Merge Monteciccardo into Pesaro ====
cat("Step 4: Merging Monteciccardo into Pesaro municipality...\n")

merged_mun <- mun_sf %>%
  filter(COMUNE %in% c("Monteciccardo", "Pesaro")) %>%
  summarise(
    COD_RIP = 3, COD_REG = 11, COD_PROV = 41, COD_CM = 0, COD_UTS = 41,
    PRO_COM = 41044, PRO_COM_T = "041044", COMUNE = "Pesaro", COMUNE_A = NA, CC_UTS = 1,
    SHAPE_LENG = as.numeric(st_length(st_boundary(st_union(geometry)))),
    SHAPE_AREA = as.numeric(st_area(st_union(geometry))),
    geometry = st_union(geometry)
  )

mun_sf <- mun_sf %>%
  filter(!COMUNE %in% c("Monteciccardo", "Pesaro")) %>%
  bind_rows(merged_mun)

cat("Municipalities after merge:", nrow(mun_sf), "\n\n")

# Step 5: Process Income Data (for full_dataset.csv only) ====
cat("Step 5: Processing income data...\n")

income_cols <- names(read.csv(
  "raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.csv",
  sep = ";", header = TRUE, nrows = 0
))

income_df_raw <- read.csv(
  "raw/Redditi_e_principali_variabili_IRPEF_su_base_comunale_CSV_2020.csv",
  sep = ";", skip = 1, header = FALSE
) %>%
  select(-c(ncol(.) - 1, ncol(.))) %>%
  `colnames<-`(income_cols) %>%
  filter(Regione != "Mancante/errata")

income_distribution <- income_df_raw %>%
  select(
    Codice.Istat.Comune, Denominazione.Comune,
    matches("Reddito.complessivo.*Frequenza")
  ) %>%
  rename_with(
    ~ gsub("Reddito.complessivo.", "", .x) %>%
      gsub("\\.+Frequenza", "", .) %>%
      gsub("\\.", " ", .) %>%
      gsub("da ([0-9]+) a ([0-9]+)", "\\1-\\2", .) %>%
      gsub("oltre ([0-9]+)", "\\1+", .) %>%
      gsub("minore o uguale a zero", "0-", .) %>%
      gsub(" euro", "", .),
    matches("Frequenza")
  ) %>%
  mutate(across(matches("^[0-9]|^0-"), ~ replace(., is.na(.), 0))) %>%
  rename(COD_MUN = Codice.Istat.Comune, NAME_MUN = Denominazione.Comune) %>%
  mutate(COD_MUN = as.numeric(as.character(COD_MUN)))

cat("Processed income data for", nrow(income_distribution), "municipalities\n\n")

# Step 6: Process POSAS (Population/Demographic Structure) Data ====
cat("Step 6: Processing POSAS demographic data...\n")

posas_files <- list.files("raw", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
posas_file <- posas_files[grepl("POSAS|posas", posas_files, ignore.case = TRUE)][1]

if (!is.na(posas_file)) {
  cat(sprintf("Reading POSAS: %s\n", basename(posas_file)))
  
  # Skip first line (title), read header from second line
  posas_raw <- read.csv(posas_file, sep = ";", stringsAsFactors = FALSE, 
                        skip = 1, na.strings = "")
  
  # Standardize column names
  names(posas_raw) <- tolower(names(posas_raw))
  
  # Rename key columns
  posas_raw <- posas_raw %>%
    rename(
      cod_mun = codice.comune,
      name_mun = comune,
      age = età
    ) %>%
    mutate(cod_mun = as.numeric(as.character(cod_mun)))
  
  # Keep ONLY age 999 (summary row) - excludes individual age breakdowns
  posas_demo <- posas_raw %>%
    filter(age == 999) %>%
    select(-age, -name_mun) %>%  # Remove age column (always 999) and name
    rename_with(~toupper(gsub("\\.", "_", .x)), -cod_mun) %>%  # Uppercase all except cod_mun
    rename(COD_MUN = cod_mun)  # <-- ADD THIS LINE: Rename to uppercase for consistency
  cat(sprintf("POSAS processed: %d municipalities with summary demographics (age 999)\n", 
              nrow(posas_demo)))
  cat(sprintf("Variables included: %s\n", paste(names(posas_demo)[-1], collapse = ", ")))
} else {
  warning("POSAS file not found!")
  posas_demo <- data.frame(cod_mun = integer())
}

# Step 7: Process Decessi (Mortality) Data ====
cat("Step 7: Processing mortality (decessi) data...\n")

decessi_files <- list.files("raw", pattern = "\\.(csv|txt)$", recursive = TRUE, full.names = TRUE)
decessi_file <- decessi_files[grepl("comuni_giornaliero", decessi_files, ignore.case = TRUE)][1]

if (!is.na(decessi_file)) {
  cat(sprintf("Reading Decessi: %s\n", basename(decessi_file)))
  
  # Try to detect separator by reading first line
  first_line <- readLines(decessi_file, n = 1, encoding = "latin1")
  cat("First line preview:", substr(first_line, 1, 100), "...\n")
  
  # Count separators to determine which one to use
  n_semicolons <- length(gregexpr(";", first_line)[[1]])
  n_commas <- length(gregexpr(",", first_line)[[1]])
  
  cat(sprintf("Detected separators: %d semicolons, %d commas\n", n_semicolons, n_commas))
  
  # Choose separator based on detection
  if (n_semicolons > n_commas) {
    cat("Using semicolon separator\n")
    decessi_raw <- read.delim(decessi_file, sep = ";", stringsAsFactors = FALSE, 
                              fileEncoding = "latin1", 
                              na.strings = c("", "NA", "n.d.", "n/a"),
                              quote = "\"")
  } else {
    cat("Using comma separator\n")
    decessi_raw <- read.csv(decessi_file, stringsAsFactors = FALSE, 
                            fileEncoding = "latin1", 
                            na.strings = c("", "NA", "n.d.", "n/a"),
                            quote = "\"")
  }
  
  cat(sprintf("Raw data: %d rows, %d columns\n", nrow(decessi_raw), ncol(decessi_raw)))
  
  # Check for required columns
  if (!"COD_PROVCOM" %in% names(decessi_raw)) {
    stop("COD_PROVCOM column not found in mortality data")
  }
  
  # Convert municipality code to numeric
  decessi_raw <- decessi_raw %>%
    mutate(COD_MUN = as.numeric(as.character(COD_PROVCOM)))
  
  # Find death columns (T_11 to T_25)
  t_cols <- grep("^T_[0-9]+", names(decessi_raw), value = TRUE)
  cat(sprintf("Found %d death total columns: %s\n", length(t_cols), 
              paste(head(t_cols, 5), collapse = ", ")))
  
  # Check for GE (date) column for filtering
  if ("GE" %in% names(decessi_raw)) {
    cat("Processing dates from GE column...\n")
    decessi_processed <- decessi_raw %>%
      mutate(
        date = as.Date(sprintf("2025%04d", as.numeric(GE)), format = "%Y%m%d"),
        month = as.numeric(format(date, "%m"))  # Base R alternative to lubridate::month()
      ) %>%
      filter(month %in% c(1, 2, 3, 4, 5, 6))  # Filter months (adjust as needed)
  } else {
    decessi_processed <- decessi_raw
  }
  
  # Convert death columns to numeric
  decessi_processed <- decessi_processed %>%
    mutate(across(all_of(t_cols), ~as.numeric(as.character(.))))
  
  # Aggregate deaths by municipality across age classes (CL_ETA) if present
  if ("CL_ETA" %in% names(decessi_processed)) {
    cat("Aggregating across age classes (CL_ETA)...\n")
    decessi_df <- decessi_processed %>%
      group_by(COD_MUN) %>%
      summarise(
        tot_decessi = sum(across(all_of(t_cols)), na.rm = TRUE),
        n_age_classes = n_distinct(CL_ETA),
        .groups = "drop"
      )
  } else {
    decessi_df <- decessi_processed %>%
      group_by(COD_MUN) %>%
      summarise(
        tot_decessi = sum(across(all_of(t_cols)), na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Rename to uppercase for consistency
  decessi_df <- decessi_df %>%
    rename(TOT_DECESSI = tot_decessi)

  if("n_age_classes" %in% names(decessi_df)) {
    decessi_df <- decessi_df %>% rename(N_AGE_CLASSES = n_age_classes)
  }
  
  cat(sprintf("Mortality processed: %d municipalities, total deaths: %.0f\n", 
              nrow(decessi_df), sum(decessi_df$TOT_DECESSI, na.rm = TRUE)))
  print(head(decessi_df, 3))
  
} else {
  warning("Decessi file not found!")
  decessi_df <- data.frame(COD_MUN = integer(), TOT_DECESSI = integer())
}

# Step 8: Create full_dataset_covariates (NO GEOMETRY!) ====
cat("Step 8: Creating covariates dataset (POSAS + Decessi)...\n")

# Start with municipality list
base_df <- mun_sf %>%
  st_drop_geometry() %>%  # remove geometry for covariates
  select(PRO_COM, COMUNE, COD_REG, COD_PROV) %>%
  rename(COD_MUN = PRO_COM, NAME_MUN = COMUNE) %>%
  mutate(COD_MUN = as.numeric(as.character(COD_MUN))) %>%
  as.data.frame()  # Ensure it's a regular data.frame

# Merge POSAS covariates
if (nrow(posas_demo) > 0 && "COD_MUN" %in% names(posas_demo)) {
  # Select only covariate columns from POSAS (exclude duplicate NAME_MUN if present)
  posas_covariates <- posas_demo %>%
    select(-any_of("NAME_MUN"))  # We'll use the name from base_df
  
  full_dataset_covariates <- base_df %>%
    left_join(posas_covariates, by = "COD_MUN")
  
  cat(sprintf("Merged POSAS: %d columns added\n", ncol(posas_covariates) - 1))
} else {
  full_dataset_covariates <- base_df
  warning("POSAS data not merged - check file structure")
}

# Merge Decessi covariates
if (nrow(decessi_df) > 0 && "COD_MUN" %in% names(decessi_df)) {
  full_dataset_covariates <- full_dataset_covariates %>%
    left_join(decessi_df, by = "COD_MUN")
  
  cat(sprintf("Merged Decessi: %d columns added\n", ncol(decessi_df) - 1))
} else {
  warning("Decessi data not merged - check file structure")
}

# Replace NA with 0 for mortality counts (if applicable)
if ("TOT_DECESSI" %in% names(full_dataset_covariates)) {
  full_dataset_covariates <- full_dataset_covariates %>%
    mutate(TOT_DECESSI = replace_na(TOT_DECESSI, 0))
}

cat(sprintf("Final covariates dataset: %d municipalities, %d columns\n\n",
            nrow(full_dataset_covariates), ncol(full_dataset_covariates)))

# Check for missing data
missing_posas <- sum(is.na(full_dataset_covariates[[ncol(base_df) + 1]]))  # First POSAS column
if (missing_posas > 0) {
  cat(sprintf("Warning: %d municipalities missing POSAS data\n", missing_posas))
}

# Step 9: Prepare geometry for spatial files ====
cat("Step 9: Preparing geometry reference...\n")

geometry_sf <- mun_sf %>%
  select(PRO_COM, COD_REG, COD_PROV, COMUNE, geometry) %>%
  rename(COD_MUN = PRO_COM, NAME_MUN = COMUNE) %>%
  mutate(COD_MUN = as.numeric(as.character(COD_MUN)))

cat("Geometry prepared\n\n")

# Step 10: Compute adjacency matrix ====
cat("Step 10: Computing spatial adjacency matrix...\n")

neighbors <- st_touches(geometry_sf, sparse = FALSE)
adj_matrix <- matrix(as.integer(neighbors), nrow = nrow(geometry_sf))
rownames(adj_matrix) <- geometry_sf$COD_MUN
colnames(adj_matrix) <- geometry_sf$COD_MUN

if (!isSymmetric(adj_matrix)) {
  warning("Adjacency matrix is not symmetric!")
} else {
  cat("Adjacency matrix is symmetric ✓\n")
}

cat("Adjacency matrix:", nrow(adj_matrix), "x", ncol(adj_matrix), "\n\n")

# Step 11: Save outputs ====
cat("Step 11: Saving final datasets...\n")

# Save income distribution
write.csv(income_distribution, "input/municipalities/full_dataset.csv", row.names = FALSE)
cat("  ✓ Income distribution saved (full_dataset.csv)\n")

# Save covariates (verify no geometry!)
if (inherits(full_dataset_covariates, "sf")) {
  full_dataset_covariates <- st_drop_geometry(full_dataset_covariates)
}
write.csv(full_dataset_covariates, "input/municipalities/full_dataset_covariates.csv", row.names = FALSE)
cat(sprintf("  ✓ Covariates saved (full_dataset_covariates.csv): %d rows, %d cols\n", 
            nrow(full_dataset_covariates), ncol(full_dataset_covariates)))

# Save geometry as RDS
saveRDS(geometry_sf, "input/municipalities/geometry.rds")
cat("  ✓ Geometry RDS saved\n")

# Save shapefile
geometry_dir <- "input/municipalities/geometry"
dir.create(geometry_dir, showWarnings = FALSE, recursive = TRUE)
st_write(geometry_sf, file.path(geometry_dir, "municipalities.shp"), 
         append = FALSE, quiet = TRUE)
cat("  ✓ Shapefile saved\n")

# Save adjacency matrix
saveRDS(adj_matrix, "input/municipalities/adj_matrix.rds")
cat("  ✓ Adjacency matrix saved\n")

cat("\n✓ All processing complete!\n\n")

# Summary ====
cat("=== OUTPUT FILES ===\n")
cat("1. input/municipalities/full_dataset.csv - Income distribution brackets\n")
cat("2. input/municipalities/full_dataset_covariates.csv - POSAS demographics + Mortality covariates\n")
cat("3. input/municipalities/geometry.rds - Geographic boundaries\n")
cat("4. input/municipalities/geometry/municipalities.shp - Shapefile boundaries\n")
cat("5. input/municipalities/adj_matrix.rds - Spatial adjacency matrix\n\n")

cat("Covariates dataset preview:\n")
print(head(full_dataset_covariates[, 1:min(5, ncol(full_dataset_covariates))]))

# Optional cleanup
# unlink("raw", recursive = TRUE)