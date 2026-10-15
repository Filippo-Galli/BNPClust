# Data Processing Scripts - Quick Reference Guide

This document explains how to use the two data processing scripts for spatial analysis.

---

## 📊 Script Overview

| Script | Purpose | Data Source | Geography |
|--------|---------|-------------|-----------|
| **ACF.R** | US Census income data | American Community Survey 2020 | US States (PUMAs) |
| **municipalities.R** | Italian income data | ISTAT (Italian Statistics) | Italian Municipalities |

---

##  File Structure

Both scripts produce the same file structure:

```
input/[folder]/
├── full_dataset.csv       # Main data file (CSV format)
├── geometry.rds           # Geometry in RDS format
├── adj_matrix.rds         # Adjacency matrix (RDS format)
└── geometry/              # Folder with shapefiles
    ├── [name].shp
    ├── [name].shx
    ├── [name].dbf
    └── [name].prj
```

---

## 🇺🇸 ACF.R - US Census Data

### What It Does
Downloads and processes US Census income data at the PUMA (Public Use Microdata Area) level for spatial analysis.

### Quick Start

**Basic usage (California only):**
```r
source("ACF.R")
```

**Custom states:**
```r
Rscript ACF.R --states "California,Texas,New York"
```

**Entire region:**
```r
Rscript ACF.R --states "West"
```

**All parameters:**
```r
Rscript ACF.R \
  --states "California,Nevada,Oregon" \
  --num-datasets 20 \
  --subsample-size 150 \
  --output-dir "input/CA/" \
  --keep-raw
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--states` | `"California"` | States or regions to process |
| `--num-datasets` | `10` | Number of subsampled datasets to generate |
| `--subsample-size` | `100` | Observations per PUMA in each subsample |
| `--output-dir` | `"input/CA/"` | Output directory |
| `--keep-raw` | `FALSE` | Keep downloaded raw files |

### Available Regions

- **West** (13 states): California, Nevada, Oregon, Washington, Arizona, Idaho, Utah, Montana, Wyoming, Colorado, New Mexico, Alaska, Hawaii
- **Midwest** (12 states): North Dakota, South Dakota, Nebraska, Kansas, Minnesota, Iowa, Missouri, Wisconsin, Illinois, Michigan, Indiana, Ohio
- **South** (14 states): Texas, Oklahoma, Arkansas, Louisiana, Mississippi, Alabama, Tennessee, Kentucky, West Virginia, Virginia, North Carolina, South Carolina, Georgia, Florida
- **Northeast** (12 states): Maryland, Delaware, District of Columbia, Pennsylvania, New Jersey, New York, Connecticut, Rhode Island, Massachusetts, Vermont, New Hampshire, Maine
- **USA**: All 51 states/districts

### Output Files

```
input/CA/                      # Output directory
├── full_dataset.csv           # Complete dataset with summary statistics
├── data_001.csv               # Subsampled dataset 1
├── data_002.csv               # Subsampled dataset 2
├── ...
├── data_010.csv               # Subsampled dataset N
├── geometry.rds               # Geometry in RDS format
├── adj_matrix.rds             # Spatial adjacency matrix (W)
└── geometry/
    ├── pumas.shp              # Geographic boundaries
    ├── pumas.shx
    ├── pumas.dbf
    └── pumas.prj
```

### Loading the Data

```r
# Load packages
library(sf)
library(dplyr)

# Load full dataset
data <- read.csv("input/CA/full_dataset.csv")

# Load a subsampled dataset
subsample <- read.csv("input/CA/data_001.csv")

# Load geometry (RDS format)
geo <- readRDS("input/CA/geometry.rds")

# Or load shapefile
geo <- st_read("input/CA/geometry/pumas.shp")

# Load adjacency matrix
adj <- readRDS("input/CA/adj_matrix.rds")

# Join for spatial analysis
full <- left_join(geo, data, by = "COD_PUMA")

# Quick visualization
plot(full["mean_log_income"])
```

### Example Analyses

**1. Single state analysis:**
```r
Rscript ACF.R --states "California" --num-datasets 50
```

**2. Multi-state comparison:**
```r
Rscript ACF.R --states "California,Texas,Florida,New York"
```

**3. Regional study:**
```r
Rscript ACF.R --states "West" --num-datasets 100 --subsample-size 200
```

**4. Full US analysis:**
```r
Rscript ACF.R --states "USA" --num-datasets 10
```

### Requirements

```r
install.packages(c("sf", "dplyr", "spdep", "argparser"))
```

---

## 🇮🇹 municipalities.R - Italian Municipality Data

### What It Does
Downloads and processes Italian income distribution data and geographic boundaries for all Italian municipalities.

### Quick Start

**Just run it:**
```r
source("municipalities.R")
```

**Or from command line:**
```r
Rscript municipalities.R
```

**That's it!** No parameters needed - processes all Italian municipalities automatically.

### What It Downloads

1. **Geographic boundaries** (ISTAT Limiti amministrativi 2020)
2. **Mortality data** (Daily deaths by municipality - not currently used in output)
3. **Population data** (POSAS 2021 - not currently used in output)
4. **Income tax data** (2020 tax returns by income bracket)

### Output Files

```
input/municipalities/
├── full_dataset.csv       # Income distribution by municipality
├── geometry.rds           # Geographic boundaries (sf object)
├── adj_matrix.rds         # Spatial adjacency matrix
└── geometry/
    ├── municipalities.shp  # Shapefile boundaries
    ├── municipalities.shx
    ├── municipalities.dbf
    └── municipalities.prj
```

### Loading the Data

```r
# Load packages
library(sf)
library(dplyr)

# Load income distribution data
data <- read.csv("input/municipalities/full_dataset.csv")

# Load geometry (RDS format)
geo <- readRDS("input/municipalities/geometry.rds")

# Or load shapefile
geo <- st_read("input/municipalities/geometry/municipalities.shp")

# Load adjacency matrix
adj <- readRDS("input/municipalities/adj_matrix.rds")

# Join for spatial visualization
full <- left_join(geo, data, by = "COD_MUN")

# Quick map
plot(full["0-10000"])  # Plot income bracket 0-10k euro
```


### Special Processing Notes

**Monteciccardo-Pesaro Merger:**
The script automatically merges Monteciccardo into Pesaro municipality to reflect the 2017 administrative change. Final output has 7,903 municipalities (down from 7,904).

### Requirements

```r
install.packages(c("sf", "dplyr"))
```

---

## 🚑 Common Issues & Troubleshooting

### ACF.R Issues

**"Failed to download state data"**
- Check internet connection
- Census servers may be temporarily down
- Script will continue with successfully downloaded states

**"No matching PUMAs found"**
- This shouldn't happen with standard 2020 data
- Check that state names are spelled correctly
- Ensure shapefiles downloaded properly

**Memory issues**
- Process fewer states at once
- Reduce `--subsample-size`
- Close other applications

**CSV file issues**
- If you need the raw log income vectors, they're stored as list columns
- Use `data$log_income_data[[1]]` to access individual PUMA vectors

### municipalities.R Issues

**"Cannot download file"**
- ISTAT servers may be slow or down
- Try again later
- Check internet connection

**"Adjacency matrix is not symmetric"**
- This is a warning, not an error
- The script will still complete
- May indicate topological issues in the shapefile

**File encoding issues**
- Italian characters should display correctly
- If not, try: `Sys.setlocale("LC_ALL", "it_IT.UTF-8")`

---

## 💾 Cleanup

**ACF.R:**
```r
# Keep raw files
Rscript ACF.R --states "California" --keep-raw

# Or clean up manually
unlink("raw", recursive = TRUE)
```

**municipalities_modified.R:**
```r
# Uncomment the last lines in the script, or run:
unlink("raw", recursive = TRUE)
```

---

## 🔗 Common Workflow with Both Datasets

```r
library(sf)
library(dplyr)

# Load US data
us_data <- read.csv("input/CA/full_dataset.csv")
us_geo <- readRDS("input/CA/geometry.rds")
us_adj <- readRDS("input/CA/adj_matrix.rds")

# Load Italian data
it_data <- read.csv("input/municipalities/full_dataset.csv")
it_geo <- readRDS("input/municipalities/geometry.rds")
it_adj <- readRDS("input/municipalities/adj_matrix.rds")

# Both datasets now follow the same pattern!
# - CSV for data
# - RDS for geometry
# - RDS for adjacency
# - geometry/ folder for shapefiles

# Join and analyze
us_full <- left_join(us_geo, us_data, by = "COD_PUMA")
it_full <- left_join(it_geo, it_data, by = "COD_MUN")

# Plot both
plot(us_full["mean_log_income"], main = "US Income")
plot(it_full["0-10000"], main = "Italian Income (0-10k bracket)")
```

---

## 📚 Additional Resources

**ACF.R Data Sources:**
- [US Census ACS PUMS](https://www.census.gov/programs-surveys/acs/microdata.html)
- [IPUMS PUMA Shapefiles](https://usa.ipums.org/usa/)

**municipalities.R Data Sources:**
- [ISTAT Geographic Boundaries](https://www.istat.it/it/archivio/222527)
- [ISTAT Income Data](https://www1.finanze.gov.it/finanze/analisi_stat/)
