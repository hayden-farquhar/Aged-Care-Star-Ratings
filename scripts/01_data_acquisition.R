# =============================================================================
# 01_data_acquisition.R
# Downloads all raw data for the Aged Care Star Ratings project
#
# Data sources:
#   1. Star Ratings quarterly data extracts (health.gov.au) — 11 quarters
#   2. SEIFA 2021 IRSD at SA2 level (ABS) — for aggregation to SA3
#   3. SEIFA 2021 population distributions at SA3 level (ABS)
#   4. Remoteness Areas 2021 allocation file (ABS)
#   5. SA1 main structure allocation file (ABS) — SA1-SA2-SA3 hierarchy
#   6. SA3 boundary shapefiles, GDA2020 (ABS)
#   7. Census 2021 population by age at SA3 — MANUAL DOWNLOAD (instructions below)
#
# Notes:
#   - Star Ratings URLs on health.gov.au are occasionally moved/re-uploaded.
#     If a download fails, check the publication page for the current link:
#     https://www.health.gov.au/resources/collections/star-ratings-quarterly-data-extracts
#   - Census DataPacks require interactive selection on the ABS website and
#     cannot be downloaded programmatically. See instructions at end of script.
#   - All downloads go to data/raw/ or data/reference/ or data/spatial/
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(httr)
library(glue)

# -- Project paths ------------------------------------------------------------

project_root <- here::here()
dir_raw      <- file.path(project_root, "data", "raw")
dir_ref      <- file.path(project_root, "data", "reference")
dir_spatial  <- file.path(project_root, "data", "spatial")

# Create subdirectories
dir.create(file.path(dir_raw, "star_ratings"), recursive = TRUE, showWarnings = FALSE)
dir.create(dir_ref, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_spatial, recursive = TRUE, showWarnings = FALSE)

# -- Helper function ----------------------------------------------------------

download_file <- function(url, destfile, description = NULL) {
  if (file.exists(destfile)) {
    message(glue("  [SKIP] Already exists: {basename(destfile)}"))
    return(invisible(TRUE))
  }
  label <- if (!is.null(description)) description else basename(destfile)
  message(glue("  Downloading: {label} ..."))
  tryCatch({
    resp <- GET(url, write_disk(destfile, overwrite = TRUE),
                progress(), timeout(300),
                user_agent("Mozilla/5.0"),
                config(http_version = 2))
    if (status_code(resp) == 200) {
      size_mb <- round(file.size(destfile) / 1e6, 2)
      message(glue("  [OK] {basename(destfile)} ({size_mb} MB)"))
      return(invisible(TRUE))
    } else {
      message(glue("  [FAIL] HTTP {status_code(resp)}: {url}"))
      unlink(destfile)
      return(invisible(FALSE))
    }
  }, error = function(e) {
    message(glue("  [ERROR] {e$message}"))
    unlink(destfile)
    return(invisible(FALSE))
  })
}

# =============================================================================
# 1. STAR RATINGS QUARTERLY DATA EXTRACTS
# =============================================================================
#
# 11 quarterly extracts from May 2023 to October 2025
# Source: https://www.health.gov.au/resources/collections/star-ratings-quarterly-data-extracts
#
# Each file is an .xlsx with sheets: Notes, Star Ratings, Detailed data
# (May 2023 has a slightly different structure — only 2 sheets, 14 columns)
#
# URL WARNING: The subfolder paths on health.gov.au are not predictable — the
# Department occasionally re-uploads files to new subfolders. If a URL fails,
# visit the publication page to find the current download link:
#   https://www.health.gov.au/resources/publications/star-ratings-quarterly-data-extract-{month}-{year}
# =============================================================================

message("\n=== 1. Star Ratings Quarterly Data Extracts ===\n")

star_ratings_files <- list(
  list(
    label    = "May 2023",
    filename = "star-ratings-quarterly-data-extract-may-2023.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2024-02/star-ratings-quarterly-data-extract-may-2023.xlsx"
  ),
  list(
    label    = "August 2023",
    filename = "star-ratings-quarterly-data-extract-august-2023.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2024-02/star-ratings-quarterly-data-extract-august-2023.xlsx"
  ),
  list(
    label    = "December 2023",
    filename = "star-ratings-quarterly-data-extract-december-2023.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2023-12/star-ratings-quarterly-data-extract-december-2023.xlsx"
  ),
  list(
    label    = "February 2024",
    filename = "star-ratings-quarterly-data-extract-february-2024.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2024-03/star-ratings-quarterly-data-extract-february-2024_1.xlsx"
  ),
  list(
    label    = "May 2024",
    filename = "star-ratings-quarterly-data-extract-may-2024.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-08/star-ratings-quarterly-data-extract-may-2024_0.xlsx"
  ),
  list(
    label    = "July 2024",
    filename = "star-ratings-quarterly-data-extract-july-2024.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-08/star-ratings-quarterly-data-extract-july-2024_0.xlsx"
  ),
  list(
    label    = "November 2024",
    filename = "star-ratings-quarterly-data-extract-november-2024.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-08/star-ratings-quarterly-data-extract-november-2024_0.xlsx"
  ),
  list(
    label    = "February 2025",
    filename = "star-ratings-quarterly-data-extract-february-2025.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-07/star-ratings-quarterly-data-extract-february-2025_1.xlsx"
  ),
  list(
    label    = "May 2025",
    filename = "star-ratings-quarterly-data-extract-may-2025.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-09/star-ratings-quarterly-data-extract-may-2025.xlsx"
  ),
  list(
    label    = "August 2025",
    filename = "star-ratings-quarterly-data-extract-august-2025.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-09/star-ratings-quarterly-data-extract-august-2025_0.xlsx"
  ),
  list(
    label    = "October 2025",
    filename = "star-ratings-quarterly-data-extract-october-2025.xlsx",
    url      = "https://www.health.gov.au/sites/default/files/2025-10/star-ratings-quarterly-data-extract-october-2025_0.xlsx"
  )
)

results_sr <- vapply(star_ratings_files, function(f) {
  download_file(
    url         = f$url,
    destfile    = file.path(dir_raw, "star_ratings", f$filename),
    description = glue("Star Ratings — {f$label}")
  )
}, logical(1))

n_ok <- sum(results_sr)
n_total <- length(results_sr)
message(glue("\nStar Ratings: {n_ok}/{n_total} files downloaded successfully."))
if (n_ok < n_total) {
  message("For failed downloads, check current URLs at:")
  message("  https://www.health.gov.au/resources/collections/star-ratings-quarterly-data-extracts")
}


# =============================================================================
# 2. SEIFA 2021 — IRSD AT SA2 LEVEL
# =============================================================================
#
# SEIFA is not published as a single index at SA3 level. We download SA2-level
# IRSD scores and aggregate to SA3 using population-weighted averages in the
# cleaning script (02_data_cleaning.R).
#
# Also download the SA3 population distributions file as supplementary.
# =============================================================================

message("\n=== 2. SEIFA 2021 (IRSD) ===\n")

seifa_base <- "https://www.abs.gov.au/statistics/people/people-and-communities/socio-economic-indexes-areas-seifa-australia/2021"

download_file(
  url      = file.path(seifa_base, "Statistical%20Area%20Level%202%2C%20Indexes%2C%20SEIFA%202021.xlsx"),
  destfile = file.path(dir_ref, "seifa_2021_sa2_indexes.xlsx"),
  description = "SEIFA 2021 — SA2 Indexes"
)

download_file(
  url      = file.path(seifa_base, "Statistical%20Area%20Level%203%2C%20Population%20Distributions%2C%20SEIFA%202021.xlsx"),
  destfile = file.path(dir_ref, "seifa_2021_sa3_population_distributions.xlsx"),
  description = "SEIFA 2021 — SA3 Population Distributions"
)


# =============================================================================
# 3. REMOTENESS AREAS 2021 ALLOCATION FILE
# =============================================================================
#
# Maps every SA1 to its Remoteness Area category (Major Cities, Inner Regional,
# Outer Regional, Remote, Very Remote). We join this to the SA1 main structure
# file to derive SA3-level remoteness classification.
# =============================================================================

message("\n=== 3. Remoteness Areas 2021 ===\n")

asgs_base <- "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/allocation-files"

download_file(
  url      = file.path(asgs_base, "RA_2021_AUST.xlsx"),
  destfile = file.path(dir_ref, "remoteness_areas_2021.xlsx"),
  description = "Remoteness Areas 2021 — SA1 allocation"
)


# =============================================================================
# 4. SA1 MAIN STRUCTURE ALLOCATION FILE
# =============================================================================
#
# Maps every SA1 to its parent SA2, SA3, SA4, GCCSA, and State/Territory.
# Needed to link Remoteness Areas and SEIFA data up to SA3 level.
# =============================================================================

message("\n=== 4. SA1 Main Structure Allocation ===\n")

download_file(
  url      = file.path(asgs_base, "SA1_2021_AUST.xlsx"),
  destfile = file.path(dir_ref, "sa1_2021_allocation.xlsx"),
  description = "SA1 2021 Main Structure — SA1/SA2/SA3/SA4 hierarchy"
)


# =============================================================================
# 5. SA3 BOUNDARY SHAPEFILES (GDA2020)
# =============================================================================
#
# Digital boundary files for mapping and spatial analysis.
# Using GDA2020 (current Australian datum).
# =============================================================================

message("\n=== 5. SA3 Boundary Shapefiles ===\n")

shapefile_url <- "https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/SA3_2021_AUST_SHP_GDA2020.zip"
shapefile_zip <- file.path(dir_spatial, "SA3_2021_AUST_SHP_GDA2020.zip")

download_file(
  url         = shapefile_url,
  destfile    = shapefile_zip,
  description = "SA3 2021 Boundary Shapefile (GDA2020)"
)

# Unzip if download succeeded
if (file.exists(shapefile_zip)) {
  shapefile_dir <- file.path(dir_spatial, "SA3_2021_AUST_SHP_GDA2020")
  if (!dir.exists(shapefile_dir)) {
    message("  Extracting shapefile ...")
    unzip(shapefile_zip, exdir = shapefile_dir)
    message(glue("  [OK] Extracted to {basename(shapefile_dir)}/"))
  } else {
    message("  [SKIP] Shapefile already extracted.")
  }
}


# =============================================================================
# 6. CENSUS 2021 — POPULATION BY AGE AT SA3 (MANUAL DOWNLOAD)
# =============================================================================
#
# The ABS Census DataPacks require interactive selection and cannot be
# downloaded programmatically. Follow these steps:
#
# 1. Go to: https://www.abs.gov.au/census/find-census-data/datapacks
# 2. Select: 2021 Census > General Community Profile > SA3
# 3. Download the Australia-wide zip file (~12 MB)
# 4. Save the zip to: data/reference/
# 5. The relevant table is G04 (Age by Sex) — contains 5-year age groups
#    which we sum for 65+ and 85+ populations in the cleaning script.
#
# Alternatively, the "Regional population by age and sex" ERP data is at:
# https://www.abs.gov.au/statistics/people/population/regional-population-age-and-sex/latest-release
# =============================================================================

message("\n=== 6. Census 2021 Population by Age (SA3) ===\n")
message("  [MANUAL] Census DataPacks require interactive download from ABS website.")
message("  Instructions:")
message("    1. Go to: https://www.abs.gov.au/census/find-census-data/datapacks")
message("    2. Select: 2021 > General Community Profile > Statistical Area Level 3")
message("    3. Download Australia-wide zip")
message("    4. Save to: data/reference/")

census_files <- list.files(dir_ref, pattern = "census|datapack|GCP", ignore.case = TRUE)
if (length(census_files) > 0) {
  message(glue("  [OK] Found Census file(s): {paste(census_files, collapse = ', ')}"))
} else {
  message("  [PENDING] No Census DataPack file found yet.")
}


# =============================================================================
# DOWNLOAD SUMMARY
# =============================================================================

message("\n=== Download Summary ===\n")

count_files <- function(dir, pattern = NULL) {
  if (!dir.exists(dir)) return(0)
  length(list.files(dir, pattern = pattern, recursive = FALSE))
}

n_sr  <- count_files(file.path(dir_raw, "star_ratings"), "[.]xlsx$")
n_ref <- count_files(dir_ref, "[.]xlsx$")
n_sp  <- count_files(dir_spatial)
message(glue("  Star Ratings extracts: {n_sr}/11"))
message(glue("  Reference files:       {n_ref} .xlsx files"))
message(glue("  Spatial files:         {n_sp} items"))

message("\n=== Next Steps ===\n")
message("  1. If any Star Ratings downloads failed, manually download from:")
message("     https://www.health.gov.au/resources/collections/star-ratings-quarterly-data-extracts")
message("  2. Download Census 2021 DataPack (SA3) — see instructions above")
message("  3. Run 02_data_cleaning.R to build the longitudinal panel")
message("\nDone.")
