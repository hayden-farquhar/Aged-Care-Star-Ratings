# =============================================================================
# 02_data_cleaning.R
# Reads all Star Ratings quarterly extracts, standardises columns, builds
# longitudinal panel, and links to SEIFA, remoteness, and SA3 geography.
#
# Inputs:
#   data/raw/star_ratings/*.xlsx      — 11 quarterly Star Ratings extracts
#   data/reference/suburbs_localities_2021.xlsx
#   data/reference/mb_2021_allocation.xlsx
#   data/reference/sa1_2021_allocation.xlsx
#   data/reference/remoteness_areas_2021.xlsx
#   data/reference/seifa_2021_sa2_indexes.xlsx
#   data/reference/seifa_2021_sa3_population_distributions.xlsx
#
# Outputs:
#   data/processed/star_ratings_panel.rds     — facility x quarter panel
#   data/processed/suburb_sa3_lookup.rds      — suburb-to-SA3 lookup table
#   data/processed/sa3_reference.rds          — SA3 with SEIFA + remoteness
#   outputs/tables/descriptive_summary.csv    — descriptive statistics
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(glue)

project_root <- here::here()
dir_raw      <- file.path(project_root, "data", "raw", "star_ratings")
dir_ref      <- file.path(project_root, "data", "reference")
dir_proc     <- file.path(project_root, "data", "processed")
dir_tables   <- file.path(project_root, "outputs", "tables")

dir.create(dir_proc, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PART 1: Read and stack all quarterly extracts
# =============================================================================

message("\n=== Part 1: Reading quarterly extracts ===\n")

files <- list.files(dir_raw, pattern = "[.]xlsx$", full.names = TRUE)
files <- sort(files)
message(glue("  Found {length(files)} files"))

# Target columns (union of all quarterly schemas)
target_cols <- c(
  "reporting_period", "service_name", "provider_name", "service_suburb",
  "purpose", "acpr", "state", "mmm_region", "mmm_code", "size",
  "overall_star_rating", "residents_experience_rating",
  "compliance_rating", "staffing_rating", "quality_measures_rating"
)

read_star_ratings <- function(filepath) {
  sheets <- excel_sheets(filepath)
  sr_sheet <- sheets[grepl("Star", sheets, ignore.case = TRUE)][1]

  df <- read_excel(filepath, sheet = sr_sheet)

  # Standardise column names
  col_map <- c(
    "Reporting Period"              = "reporting_period",
    "Service Name"                  = "service_name",
    "Provider Name"                 = "provider_name",
    "Service Suburb"                = "service_suburb",
    "Purpose"                       = "purpose",
    "Aged Care Planning Region"     = "acpr",
    "State/Territory"               = "state",
    "MMM Region"                    = "mmm_region",
    "MMM Code"                      = "mmm_code",
    "Size"                          = "size",
    "Overall Star Rating"           = "overall_star_rating",
    "Residents' Experience rating"  = "residents_experience_rating",
    "Compliance rating"             = "compliance_rating",
    "Staffing rating"               = "staffing_rating",
    "Quality Measures rating"       = "quality_measures_rating"
  )

  names(df) <- col_map[names(df)]

  # Add missing columns as NA
  for (col in target_cols) {
    if (!col %in% names(df)) {
      df[[col]] <- NA_character_
    }
  }

  df <- df[, target_cols]

  # Standardise purpose (capitalisation varies across quarters)
  df <- df %>%
    mutate(purpose = case_when(
      str_to_lower(purpose) == "for profit"     ~ "For-profit",
      str_to_lower(purpose) == "not for profit"  ~ "Not-for-profit",
      str_to_lower(purpose) == "government"      ~ "Government",
      TRUE                                       ~ purpose
    ))

  # Convert ratings to numeric
  rating_cols <- c("overall_star_rating", "residents_experience_rating",
                   "compliance_rating", "staffing_rating", "quality_measures_rating")
  for (col in rating_cols) {
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }

  # Parse reporting period into a date for ordering
  df <- df %>%
    mutate(
      quarter_label = reporting_period,
      quarter_date  = parse_date_time(reporting_period, orders = c("B Y", "BY"))
    )

  df$source_file <- basename(filepath)

  df
}

# Read all files
all_quarters <- lapply(files, function(f) {
  message(glue("  Reading: {basename(f)}"))
  tryCatch(read_star_ratings(f),
           error = function(e) {
             message(glue("  [ERROR] {e$message}"))
             NULL
           })
})

all_quarters <- Filter(Negate(is.null), all_quarters)
panel <- bind_rows(all_quarters)

message(glue("\n  Total records: {nrow(panel)}"))
message(glue("  Quarters: {n_distinct(panel$quarter_label)}"))
message(glue("  Unique services: {n_distinct(panel$service_name)}"))


# =============================================================================
# PART 2: Create stable facility identifier
# =============================================================================

message("\n=== Part 2: Creating facility identifiers ===\n")

# Use service_name + state as the unique facility key
# (service names are unique within state in the Star Ratings data)
panel <- panel %>%
  mutate(facility_id = paste(service_name, state, sep = " | "))

n_facilities <- n_distinct(panel$facility_id)
message(glue("  Unique facilities: {n_facilities}"))

# Order quarters chronologically
panel <- panel %>%
  arrange(quarter_date, facility_id) %>%
  group_by(facility_id) %>%
  mutate(quarter_seq = dense_rank(quarter_date)) %>%
  ungroup()

# Create a quarter number for the whole panel
panel <- panel %>%
  mutate(quarter_num = dense_rank(quarter_date))

message("  Quarter ordering:")
panel %>%
  distinct(quarter_num, quarter_label, quarter_date) %>%
  arrange(quarter_num) %>%
  mutate(msg = glue("    Q{quarter_num}: {quarter_label} ({quarter_date})")) %>%
  pull(msg) %>%
  walk(message)


# =============================================================================
# PART 3: Fill in missing suburbs from other quarters
# =============================================================================

message("\n=== Part 3: Filling missing suburbs ===\n")

# Suburbs are only available from Feb 2024 onward. For earlier quarters,
# use the suburb from a later quarter for the same facility.
suburb_lookup_internal <- panel %>%
  filter(!is.na(service_suburb)) %>%
  distinct(facility_id, service_suburb) %>%
  group_by(facility_id) %>%
  slice(1) %>%
  ungroup()

n_missing_before <- sum(is.na(panel$service_suburb))

panel <- panel %>%
  select(-service_suburb) %>%
  left_join(suburb_lookup_internal, by = "facility_id")

n_missing_after <- sum(is.na(panel$service_suburb))
message(glue("  Missing suburbs before fill: {n_missing_before}"))
message(glue("  Missing suburbs after fill:  {n_missing_after}"))


# =============================================================================
# PART 4: Link facilities to SA3 via suburb
# =============================================================================

message("\n=== Part 4: Linking to SA3 geography ===\n")

# Build suburb -> SA3 lookup from ABS mesh block and suburb allocation files
sal <- read_excel(file.path(dir_ref, "suburbs_localities_2021.xlsx"), col_types = "text")
mb  <- read_excel(file.path(dir_ref, "mb_2021_allocation.xlsx"), col_types = "text")

suburb_sa3 <- sal %>%
  select(MB_CODE_2021, suburb_raw = SAL_NAME_2021, sal_state = STATE_NAME_2021) %>%
  inner_join(mb %>% select(MB_CODE_2021, SA3_CODE_2021, SA3_NAME_2021,
                           SA2_CODE_2021, SA2_NAME_2021),
             by = "MB_CODE_2021") %>%
  filter(SA3_CODE_2021 != "ZZZZZ") %>%
  mutate(suburb_clean = trimws(gsub(" *[(].*[)] *$", "", suburb_raw)))

# Most common SA3 per suburb+state
suburb_sa3_lookup <- suburb_sa3 %>%
  count(suburb_clean, sal_state, SA3_CODE_2021, SA3_NAME_2021) %>%
  group_by(suburb_clean, sal_state) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(suburb_upper = toupper(suburb_clean)) %>%
  select(suburb_upper, sal_state, sa3_code = SA3_CODE_2021, sa3_name = SA3_NAME_2021)

# State name mapping for join
state_map <- c(
  "NSW" = "New South Wales", "VIC" = "Victoria", "QLD" = "Queensland",
  "WA"  = "Western Australia", "SA"  = "South Australia",
  "TAS" = "Tasmania", "ACT" = "Australian Capital Territory",
  "NT"  = "Northern Territory"
)

# Handle common abbreviation mismatches (MT vs Mount, ST vs Saint, etc.)
expand_abbreviations <- function(x) {
  x <- gsub("^MT ", "MOUNT ", x)
  x <- gsub("^ST ", "SAINT ", x)
  x <- gsub(" MT ", " MOUNT ", x)
  x <- gsub(" ST ", " SAINT ", x)
  x
}

panel <- panel %>%
  mutate(
    suburb_upper = toupper(service_suburb),
    state_full   = state_map[state],
    suburb_expanded = expand_abbreviations(suburb_upper)
  )

# Join: first try exact match, then try with expanded abbreviations
panel <- panel %>%
  left_join(suburb_sa3_lookup,
            by = c("suburb_upper" = "suburb_upper", "state_full" = "sal_state"))

# For unmatched, try expanded abbreviations
unmatched_idx <- is.na(panel$sa3_code) & !is.na(panel$service_suburb)
if (sum(unmatched_idx) > 0) {
  expanded_match <- panel %>%
    filter(unmatched_idx) %>%
    select(-sa3_code, -sa3_name) %>%
    left_join(suburb_sa3_lookup,
              by = c("suburb_expanded" = "suburb_upper", "state_full" = "sal_state"))

  panel$sa3_code[unmatched_idx]  <- expanded_match$sa3_code
  panel$sa3_name[unmatched_idx]  <- expanded_match$sa3_name
}

n_with_sa3 <- sum(!is.na(panel$sa3_code))
n_total    <- nrow(panel)
message(glue("  Matched to SA3: {n_with_sa3}/{n_total} ({round(n_with_sa3/n_total*100, 1)}%)"))

n_unmatched_facilities <- panel %>%
  filter(is.na(sa3_code), !is.na(service_suburb)) %>%
  distinct(facility_id) %>%
  nrow()
message(glue("  Facilities without SA3: {n_unmatched_facilities}"))

# Clean up temp columns
panel <- panel %>%
  select(-suburb_upper, -state_full, -suburb_expanded)


# =============================================================================
# PART 5: Build SA3 reference table (SEIFA + remoteness)
# =============================================================================

message("\n=== Part 5: Building SA3 reference table ===\n")

# --- 5a: SEIFA IRSD at SA2, aggregate to SA3 ---
seifa_sa2 <- read_excel(file.path(dir_ref, "seifa_2021_sa2_indexes.xlsx"),
                        sheet = 2, skip = 5)

# Find the IRSD score column and SA2 code
seifa_cols <- names(seifa_sa2)
# The file has: SA2 code, SA2 name, IRSD score, IRSD rank, IRSD decile, ...
# Column names vary — identify by position or pattern
seifa_sa2 <- seifa_sa2 %>%
  select(1:6) %>%
  setNames(c("sa2_code", "sa2_name", "irsd_score", "irsd_rank",
             "irsd_decile", "irsd_percentile")) %>%
  filter(!is.na(sa2_code), nchar(sa2_code) == 9) %>%
  mutate(
    sa2_code  = as.character(sa2_code),
    irsd_score = as.numeric(irsd_score)
  )

# Get SA2 -> SA3 mapping from SA1 allocation file
sa1_alloc <- read_excel(file.path(dir_ref, "sa1_2021_allocation.xlsx"), col_types = "text")
sa2_sa3 <- sa1_alloc %>%
  distinct(SA2_CODE_2021, SA3_CODE_2021, SA3_NAME_2021)

# Get SA2 populations for weighting (use SEIFA "usual resident population")
# For now, use equal weights (unweighted mean) — close enough for SA3 aggregation
seifa_sa3 <- seifa_sa2 %>%
  left_join(sa2_sa3, by = c("sa2_code" = "SA2_CODE_2021")) %>%
  filter(!is.na(SA3_CODE_2021)) %>%
  group_by(sa3_code = SA3_CODE_2021, sa3_name = SA3_NAME_2021) %>%
  summarise(
    irsd_score_mean = mean(irsd_score, na.rm = TRUE),
    irsd_score_min  = min(irsd_score, na.rm = TRUE),
    irsd_score_max  = max(irsd_score, na.rm = TRUE),
    n_sa2           = n(),
    .groups = "drop"
  )

# Create IRSD quintile from the mean SA3 score
seifa_sa3 <- seifa_sa3 %>%
  mutate(
    irsd_quintile = ntile(irsd_score_mean, 5),
    irsd_decile   = ntile(irsd_score_mean, 10)
  )

message(glue("  SEIFA: {nrow(seifa_sa3)} SA3 areas with IRSD scores"))

# --- 5b: Remoteness classification at SA3 level ---
ra <- read_excel(file.path(dir_ref, "remoteness_areas_2021.xlsx"), col_types = "text")

# RA file maps SA1 -> Remoteness Area. Get SA1 -> SA3 from allocation file.
ra_sa3 <- ra %>%
  select(SA1_CODE_2021 = 1, ra_code = 2, ra_name = 3) %>%
  inner_join(sa1_alloc %>% select(SA1_CODE_2021, SA3_CODE_2021),
             by = "SA1_CODE_2021") %>%
  count(SA3_CODE_2021, ra_name) %>%
  group_by(SA3_CODE_2021) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sa3_code = SA3_CODE_2021, remoteness = ra_name)

message(glue("  Remoteness: {nrow(ra_sa3)} SA3 areas classified"))

# --- 5c: Combine into SA3 reference table ---
sa3_ref <- seifa_sa3 %>%
  left_join(ra_sa3, by = "sa3_code")

message(glue("  SA3 reference table: {nrow(sa3_ref)} rows"))


# =============================================================================
# PART 6: Join SA3 reference data to panel
# =============================================================================

message("\n=== Part 6: Joining reference data to panel ===\n")

panel <- panel %>%
  left_join(sa3_ref %>% select(sa3_code, irsd_score_mean, irsd_quintile,
                               irsd_decile, remoteness),
            by = "sa3_code")

n_with_irsd <- sum(!is.na(panel$irsd_quintile))
message(glue("  Records with SEIFA: {n_with_irsd}/{nrow(panel)} ({round(n_with_irsd/nrow(panel)*100, 1)}%)"))


# =============================================================================
# PART 7: Descriptive statistics
# =============================================================================

message("\n=== Part 7: Descriptive statistics ===\n")

# --- 7a: Overall summary ---
overall <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  summarise(
    n_records    = n(),
    n_facilities = n_distinct(facility_id),
    n_quarters   = n_distinct(quarter_label),
    mean_rating  = round(mean(overall_star_rating, na.rm = TRUE), 2),
    sd_rating    = round(sd(overall_star_rating, na.rm = TRUE), 2),
    median_rating = median(overall_star_rating, na.rm = TRUE)
  )
message("  Overall summary:")
message(glue("    Records: {overall$n_records}"))
message(glue("    Facilities: {overall$n_facilities}"))
message(glue("    Quarters: {overall$n_quarters}"))
message(glue("    Mean rating: {overall$mean_rating} (SD {overall$sd_rating})"))

# --- 7b: Rating distribution ---
rating_dist <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  count(overall_star_rating, name = "n") %>%
  mutate(pct = round(n / sum(n) * 100, 1))

message("\n  Rating distribution (all quarters):")
for (i in seq_len(nrow(rating_dist))) {
  message(glue("    {rating_dist$overall_star_rating[i]} stars: {rating_dist$n[i]} ({rating_dist$pct[i]}%)"))
}

# --- 7c: By provider type ---
by_provider <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  group_by(purpose) %>%
  summarise(
    n = n(),
    mean_rating = round(mean(overall_star_rating, na.rm = TRUE), 2),
    pct_low = round(mean(overall_star_rating <= 2, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

message("\n  By provider type:")
for (i in seq_len(nrow(by_provider))) {
  message(glue("    {by_provider$purpose[i]}: mean={by_provider$mean_rating[i]}, %<=2 stars={by_provider$pct_low[i]}% (n={by_provider$n[i]})"))
}

# --- 7d: By remoteness ---
by_remoteness <- panel %>%
  filter(!is.na(overall_star_rating), !is.na(remoteness)) %>%
  group_by(remoteness) %>%
  summarise(
    n = n(),
    mean_rating = round(mean(overall_star_rating, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(mean_rating)

message("\n  By remoteness:")
for (i in seq_len(nrow(by_remoteness))) {
  message(glue("    {by_remoteness$remoteness[i]}: mean={by_remoteness$mean_rating[i]} (n={by_remoteness$n[i]})"))
}

# --- 7e: By IRSD quintile ---
by_irsd <- panel %>%
  filter(!is.na(overall_star_rating), !is.na(irsd_quintile)) %>%
  group_by(irsd_quintile) %>%
  summarise(
    n = n(),
    mean_rating = round(mean(overall_star_rating, na.rm = TRUE), 2),
    .groups = "drop"
  )

message("\n  By IRSD quintile (1=most disadvantaged):")
for (i in seq_len(nrow(by_irsd))) {
  message(glue("    Q{by_irsd$irsd_quintile[i]}: mean={by_irsd$mean_rating[i]} (n={by_irsd$n[i]})"))
}

# --- 7f: By quarter (trajectory) ---
by_quarter <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  group_by(quarter_num, quarter_label) %>%
  summarise(
    n = n(),
    mean_rating = round(mean(overall_star_rating, na.rm = TRUE), 2),
    pct_low = round(mean(overall_star_rating <= 2, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(quarter_num)

message("\n  By quarter:")
for (i in seq_len(nrow(by_quarter))) {
  message(glue("    Q{by_quarter$quarter_num[i]} ({by_quarter$quarter_label[i]}): mean={by_quarter$mean_rating[i]}, %<=2 stars={by_quarter$pct_low[i]}% (n={by_quarter$n[i]})"))
}


# =============================================================================
# PART 8: Save outputs
# =============================================================================

message("\n=== Part 8: Saving outputs ===\n")

# Save panel
saveRDS(panel, file.path(dir_proc, "star_ratings_panel.rds"))
message(glue("  Saved: star_ratings_panel.rds ({nrow(panel)} rows x {ncol(panel)} cols)"))

# Save suburb-SA3 lookup
saveRDS(suburb_sa3_lookup, file.path(dir_proc, "suburb_sa3_lookup.rds"))
message(glue("  Saved: suburb_sa3_lookup.rds ({nrow(suburb_sa3_lookup)} entries)"))

# Save SA3 reference table
saveRDS(sa3_ref, file.path(dir_proc, "sa3_reference.rds"))
message(glue("  Saved: sa3_reference.rds ({nrow(sa3_ref)} SA3 areas)"))

# Save descriptive tables
desc_tables <- list(
  overall      = overall,
  rating_dist  = rating_dist,
  by_provider  = by_provider,
  by_remoteness = by_remoteness,
  by_irsd      = by_irsd,
  by_quarter   = by_quarter
)

for (nm in names(desc_tables)) {
  outpath <- file.path(dir_tables, paste0("descriptive_", nm, ".csv"))
  write.csv(desc_tables[[nm]], outpath, row.names = FALSE)
}
message(glue("  Saved: {length(desc_tables)} descriptive tables to outputs/tables/"))

message("\nDone. Panel dataset ready for analysis.")
