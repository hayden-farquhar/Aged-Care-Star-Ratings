# =============================================================================
# 06_international_comparison.R
# International comparison: Australian Star Ratings vs US CMS Five-Star
#
# Analyses:
#   1. Download US CMS Five-Star Provider Info dataset
#   2. Compute descriptive distribution of US ratings (overall + sub-categories)
#   3. Compare Australian and US rating distributions (descriptive only)
#   4. Compare ownership type distributions
#   5. Generate comparison figures and tables
#
# Inputs:
#   data/processed/star_ratings_panel.rds    (Australian data)
#   US CMS Provider Info CSV (downloaded)
#
# Outputs:
#   outputs/figures/  — comparison bar charts, distribution plots
#   outputs/tables/   — comparison summary tables
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(httr)
library(glue)
library(purrr)

project_root <- here::here()
dir_raw  <- file.path(project_root, "data", "raw")
dir_proc <- file.path(project_root, "data", "processed")
dir_fig  <- file.path(project_root, "outputs", "figures")
dir_tab  <- file.path(project_root, "outputs", "tables")

dir.create(file.path(dir_raw, "cms"), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Part 1: Download US CMS Five-Star data
# =============================================================================

message("\n=== Part 1: Download US CMS Five-Star data ===\n")

cms_file <- file.path(dir_raw, "cms", "NH_ProviderInfo.csv")

if (!file.exists(cms_file)) {
  # CMS Provider Info dataset — landing page: https://data.cms.gov/provider-data/dataset/4pq5-n9py
  # Try the January 2026 extract first, fall back to API endpoint
  cms_urls <- c(
    "https://data.cms.gov/provider-data/sites/default/files/resources/816c17cdfc690511f78287c6bb8267c0_1769652359/NH_ProviderInfo_Jan2026.csv",
    "https://data.cms.gov/provider-data/api/1/datastore/query/4pq5-n9py?limit=0&format=csv"
  )

  downloaded <- FALSE
  for (url in cms_urls) {
    message(glue("  Trying: {url}"))
    tryCatch({
      resp <- GET(url,
                  user_agent("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)"),
                  config(http_version = 2, followlocation = TRUE),
                  write_disk(cms_file, overwrite = TRUE),
                  timeout(120))
      if (status_code(resp) == 200 && file.size(cms_file) > 10000) {
        downloaded <- TRUE
        message(glue("  Downloaded: {round(file.size(cms_file)/1e6, 1)} MB"))
        break
      } else {
        message(glue("  HTTP {status_code(resp)}, trying next URL..."))
        file.remove(cms_file)
      }
    }, error = function(e) {
      message(glue("  Error: {e$message}"))
      if (file.exists(cms_file)) file.remove(cms_file)
    })
  }

  if (!downloaded) {
    stop("Could not download CMS data. Please download manually from:\n",
         "https://data.cms.gov/provider-data/dataset/4pq5-n9py\n",
         "Save as: ", cms_file)
  }
} else {
  message(glue("  CMS data already downloaded: {round(file.size(cms_file)/1e6, 1)} MB"))
}

# =============================================================================
# Part 2: Load and clean US data
# =============================================================================

message("\n=== Part 2: Load and clean US data ===\n")

cms_raw <- read_csv(cms_file, show_col_types = FALSE)
message(glue("  Raw CMS records: {nrow(cms_raw)}"))
message(glue("  Columns: {ncol(cms_raw)}"))

# Identify key columns (CMS column names vary slightly between extracts)
# Look for star rating columns
rating_cols <- names(cms_raw)[grepl("rating|star", names(cms_raw), ignore.case = TRUE)]
message(glue("  Rating-related columns: {paste(rating_cols, collapse = ', ')}"))

# Standardise column names — find the key ones
# Exclude "Chain Average" columns — we want facility-level ratings
find_col <- function(pattern, exclude_pattern = "chain") {
  matches <- names(cms_raw)[grepl(pattern, names(cms_raw), ignore.case = TRUE)]
  if (!is.null(exclude_pattern)) {
    matches <- matches[!grepl(exclude_pattern, matches, ignore.case = TRUE)]
  }
  if (length(matches) == 0) return(NA_character_)
  matches[1]
}

col_overall     <- find_col("^overall.*rating")
col_health      <- find_col("^health.inspection.*rating")
col_staffing    <- find_col("^staffing.*rating")
col_qm          <- find_col("^qm.*rating")
col_ownership   <- find_col("ownership", exclude_pattern = NULL)
col_state       <- find_col("^provider.state|^state$", exclude_pattern = NULL)
col_beds        <- find_col("certified.beds|number.*beds", exclude_pattern = NULL)
col_name        <- find_col("provider.name|facility.name", exclude_pattern = NULL)
col_id          <- find_col("federal.provider|provider.number|cms.certification", exclude_pattern = NULL)

message(glue("  Overall rating column: {col_overall}"))
message(glue("  Health inspection column: {col_health}"))
message(glue("  Staffing column: {col_staffing}"))
message(glue("  QM column: {col_qm}"))
message(glue("  Ownership column: {col_ownership}"))

# Build clean US dataset
us <- cms_raw %>%
  transmute(
    facility_id   = .data[[col_id]],
    facility_name = if (!is.na(col_name)) .data[[col_name]] else NA_character_,
    state         = if (!is.na(col_state)) .data[[col_state]] else NA_character_,
    overall_rating      = as.integer(.data[[col_overall]]),
    health_inspection   = as.integer(.data[[col_health]]),
    staffing_rating     = as.integer(.data[[col_staffing]]),
    quality_measures    = as.integer(.data[[col_qm]]),
    ownership   = if (!is.na(col_ownership)) .data[[col_ownership]] else NA_character_,
    beds        = if (!is.na(col_beds)) as.integer(.data[[col_beds]]) else NA_integer_
  ) %>%
  filter(!is.na(overall_rating))

# Standardise ownership
us <- us %>%
  mutate(
    ownership_std = case_when(
      grepl("for profit|for-profit|proprietary", ownership, ignore.case = TRUE) &
        !grepl("non|not", ownership, ignore.case = TRUE) ~ "For-profit",
      grepl("non.?profit|not.?for.?profit|voluntary", ownership, ignore.case = TRUE) ~ "Not-for-profit",
      grepl("government|state|county|city|federal|tribal", ownership, ignore.case = TRUE) ~ "Government",
      TRUE ~ "Other"
    )
  )

message(glue("  US facilities with overall rating: {nrow(us)}"))
message(glue("  US ownership breakdown:"))
us %>% count(ownership_std) %>% mutate(pct = round(n/sum(n)*100, 1)) %>%
  { walk(seq_len(nrow(.)), function(i) message(glue("    {.$ownership_std[i]}: {.$n[i]} ({.$pct[i]}%)"))) }

# =============================================================================
# Part 3: Load Australian data (latest quarter)
# =============================================================================

message("\n=== Part 3: Load Australian data ===\n")

panel <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
latest_q <- max(panel$quarter_num)
latest_label <- panel %>% filter(quarter_num == latest_q) %>% pull(quarter_label) %>% unique()

au <- panel %>%
  filter(quarter_num == latest_q, !is.na(overall_star_rating)) %>%
  transmute(
    facility_id,
    overall_rating            = overall_star_rating,
    residents_experience      = residents_experience_rating,
    compliance_rating         = compliance_rating,
    staffing_rating           = staffing_rating,
    quality_measures          = quality_measures_rating,
    ownership_std             = purpose
  )

message(glue("  AU facilities (latest quarter, {latest_label}): {nrow(au)}"))
message(glue("  AU ownership breakdown:"))
au %>% count(ownership_std) %>% mutate(pct = round(n/sum(n)*100, 1)) %>%
  { walk(seq_len(nrow(.)), function(i) message(glue("    {.$ownership_std[i]}: {.$n[i]} ({.$pct[i]}%)"))) }

# =============================================================================
# Part 4: Descriptive comparison — overall ratings
# =============================================================================

message("\n=== Part 4: Rating distribution comparison ===\n")

# Overall rating distributions
au_dist <- au %>%
  count(overall_rating) %>%
  mutate(pct = n / sum(n) * 100, country = "Australia")

us_dist <- us %>%
  count(overall_rating) %>%
  mutate(pct = n / sum(n) * 100, country = "United States")

comparison <- bind_rows(au_dist, us_dist) %>%
  select(country, overall_rating, n, pct)

message("  Overall rating distributions:")
for (r in 1:5) {
  au_pct <- au_dist %>% filter(overall_rating == r) %>% pull(pct) %>% { if(length(.) == 0) 0 else round(., 1) }
  us_pct <- us_dist %>% filter(overall_rating == r) %>% pull(pct) %>% { if(length(.) == 0) 0 else round(., 1) }
  message(glue("    {r} star: AU = {au_pct}%, US = {us_pct}%"))
}

au_mean <- round(mean(au$overall_rating, na.rm = TRUE), 2)
us_mean <- round(mean(us$overall_rating, na.rm = TRUE), 2)
au_sd   <- round(sd(au$overall_rating, na.rm = TRUE), 2)
us_sd   <- round(sd(us$overall_rating, na.rm = TRUE), 2)
message(glue("  AU mean = {au_mean} (SD = {au_sd}), US mean = {us_mean} (SD = {us_sd})"))

# Sub-category distributions
message("\n  Sub-category means:")

# US sub-categories
us_health_mean <- round(mean(us$health_inspection, na.rm = TRUE), 2)
us_staff_mean  <- round(mean(us$staffing_rating, na.rm = TRUE), 2)
us_qm_mean     <- round(mean(us$quality_measures, na.rm = TRUE), 2)
message(glue("    US: Health Inspection = {us_health_mean}, Staffing = {us_staff_mean}, QM = {us_qm_mean}"))

# AU sub-categories
au_re_mean   <- round(mean(au$residents_experience, na.rm = TRUE), 2)
au_comp_mean <- round(mean(au$compliance_rating, na.rm = TRUE), 2)
au_staff_mean <- round(mean(au$staffing_rating, na.rm = TRUE), 2)
au_qm_mean   <- round(mean(au$quality_measures, na.rm = TRUE), 2)
message(glue("    AU: Residents' Experience = {au_re_mean}, Compliance = {au_comp_mean}, Staffing = {au_staff_mean}, QM = {au_qm_mean}"))

# Ownership comparison
message("\n  Mean rating by ownership:")
au_own <- au %>% group_by(ownership_std) %>%
  summarise(mean_rating = round(mean(overall_rating, na.rm = TRUE), 2),
            n = n(), .groups = "drop") %>%
  mutate(country = "Australia")

us_own <- us %>%
  filter(ownership_std != "Other") %>%
  group_by(ownership_std) %>%
  summarise(mean_rating = round(mean(overall_rating, na.rm = TRUE), 2),
            n = n(), .groups = "drop") %>%
  mutate(country = "United States")

own_comparison <- bind_rows(au_own, us_own)
for (i in seq_len(nrow(own_comparison))) {
  message(glue("    {own_comparison$country[i]} — {own_comparison$ownership_std[i]}: {own_comparison$mean_rating[i]} (n={own_comparison$n[i]})"))
}

# =============================================================================
# Part 5: Save comparison tables
# =============================================================================

message("\n=== Part 5: Saving comparison tables ===\n")

# Table 1: Overall distribution comparison
dist_wide <- comparison %>%
  pivot_wider(names_from = country, values_from = c(n, pct), names_sep = "_") %>%
  arrange(overall_rating)
write_csv(dist_wide, file.path(dir_tab, "intl_rating_distribution.csv"))
message("  Saved: intl_rating_distribution.csv")

# Table 2: Summary comparison
summary_tab <- tibble(
  metric = c("Number of facilities", "Mean overall rating", "SD overall rating",
             "1-star (%)", "2-star (%)", "3-star (%)", "4-star (%)", "5-star (%)",
             "For-profit mean", "Not-for-profit mean", "Government mean",
             "For-profit share (%)",
             "Rating scale", "Year introduced", "Sub-categories"),
  australia = c(
    nrow(au), au_mean, au_sd,
    round(sum(au$overall_rating == 1) / nrow(au) * 100, 1),
    round(sum(au$overall_rating == 2) / nrow(au) * 100, 1),
    round(sum(au$overall_rating == 3) / nrow(au) * 100, 1),
    round(sum(au$overall_rating == 4) / nrow(au) * 100, 1),
    round(sum(au$overall_rating == 5) / nrow(au) * 100, 1),
    au_own %>% filter(ownership_std == "For-profit") %>% pull(mean_rating),
    au_own %>% filter(ownership_std == "Not-for-profit") %>% pull(mean_rating),
    au_own %>% filter(ownership_std == "Government") %>% pull(mean_rating),
    round(sum(au$ownership_std == "For-profit") / nrow(au) * 100, 1),
    "1-5 stars", "2022",
    "Residents' Experience, Compliance, Staffing, Quality Measures"
  ),
  united_states = c(
    nrow(us), us_mean, us_sd,
    round(sum(us$overall_rating == 1) / nrow(us) * 100, 1),
    round(sum(us$overall_rating == 2) / nrow(us) * 100, 1),
    round(sum(us$overall_rating == 3) / nrow(us) * 100, 1),
    round(sum(us$overall_rating == 4) / nrow(us) * 100, 1),
    round(sum(us$overall_rating == 5) / nrow(us) * 100, 1),
    us_own %>% filter(ownership_std == "For-profit") %>% pull(mean_rating),
    us_own %>% filter(ownership_std == "Not-for-profit") %>% pull(mean_rating),
    us_own %>% filter(ownership_std == "Government") %>% pull(mean_rating),
    round(sum(us$ownership_std == "For-profit") / nrow(us) * 100, 1),
    "1-5 stars", "2008",
    "Health Inspection, Staffing, Quality Measures"
  )
)
write_csv(summary_tab, file.path(dir_tab, "intl_summary_comparison.csv"))
message("  Saved: intl_summary_comparison.csv")

# Table 3: Ownership comparison
write_csv(own_comparison, file.path(dir_tab, "intl_ownership_comparison.csv"))
message("  Saved: intl_ownership_comparison.csv")

# =============================================================================
# Part 6: Comparison figures
# =============================================================================

message("\n=== Part 6: Generating comparison figures ===\n")

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Figure 1: Side-by-side bar chart of overall rating distributions
fig1_data <- bind_rows(
  au %>% count(overall_rating) %>% mutate(pct = n/sum(n)*100, country = "Australia"),
  us %>% count(overall_rating) %>% mutate(pct = n/sum(n)*100, country = "United States")
)

p1 <- ggplot(fig1_data, aes(x = factor(overall_rating), y = pct, fill = country)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Australia" = "#2171B5", "United States" = "#CB181D")) +
  labs(
    title = "Overall Star Rating Distribution: Australia vs United States",
    subtitle = glue("Australia: {nrow(au)} facilities ({latest_label}) | US: {nrow(us)} facilities (Jan 2026)"),
    x = "Overall Star Rating", y = "Percentage of Facilities (%)",
    fill = "Country"
  ) +
  theme_pub +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggsave(file.path(dir_fig, "intl_rating_distribution.png"), p1,
       width = 9, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: intl_rating_distribution.png")

# Figure 2: Ownership x rating comparison
fig2_data <- bind_rows(
  au %>% filter(ownership_std %in% c("For-profit", "Not-for-profit", "Government")) %>%
    mutate(country = "Australia"),
  us %>% filter(ownership_std %in% c("For-profit", "Not-for-profit", "Government")) %>%
    rename(overall_rating_val = overall_rating) %>%
    mutate(overall_rating = overall_rating_val, country = "United States") %>%
    select(overall_rating, ownership_std, country)
) %>%
  select(overall_rating, ownership_std, country)

p2 <- ggplot(fig2_data, aes(x = ownership_std, y = overall_rating, fill = country)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
               outlier.size = 0.5, outlier.alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3,
               position = position_dodge(width = 0.8), show.legend = FALSE) +
  scale_fill_manual(values = c("Australia" = "#2171B5", "United States" = "#CB181D")) +
  labs(
    title = "Star Rating by Ownership Type: Australia vs United States",
    subtitle = "Diamonds indicate mean rating",
    x = "Ownership Type", y = "Overall Star Rating",
    fill = "Country"
  ) +
  theme_pub +
  scale_y_continuous(breaks = 1:5, limits = c(0.5, 5.5))

ggsave(file.path(dir_fig, "intl_ownership_comparison.png"), p2,
       width = 9, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: intl_ownership_comparison.png")

# Figure 3: Sub-category comparison (grouped bar chart)
subcat_data <- tibble(
  country = rep(c("Australia", "United States"), each = 3),
  subcategory = rep(c("Staffing", "Quality Measures", "Compliance / Health Inspection"), 2),
  mean_rating = c(au_staff_mean, au_qm_mean, au_comp_mean,
                  us_staff_mean, us_qm_mean, us_health_mean)
) %>%
  mutate(subcategory = factor(subcategory, levels = c("Staffing", "Quality Measures",
                                                       "Compliance / Health Inspection")))

p3 <- ggplot(subcat_data, aes(x = subcategory, y = mean_rating, fill = country)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", mean_rating)),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Australia" = "#2171B5", "United States" = "#CB181D")) +
  labs(
    title = "Mean Sub-category Ratings: Australia vs United States",
    subtitle = "Comparable sub-categories only (Residents' Experience is AU-specific)",
    x = NULL, y = "Mean Star Rating",
    fill = "Country"
  ) +
  theme_pub +
  scale_y_continuous(limits = c(0, 5.5), breaks = 1:5, expand = expansion(mult = c(0, 0.05)))

ggsave(file.path(dir_fig, "intl_subcategory_comparison.png"), p3,
       width = 9, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: intl_subcategory_comparison.png")

# Figure 4: Stacked distribution by ownership
fig4_data <- bind_rows(
  au %>% filter(ownership_std %in% c("For-profit", "Not-for-profit", "Government")) %>%
    group_by(ownership_std) %>%
    count(overall_rating) %>%
    mutate(pct = n / sum(n) * 100, country = "Australia") %>% ungroup(),
  us %>% filter(ownership_std %in% c("For-profit", "Not-for-profit", "Government")) %>%
    group_by(ownership_std) %>%
    count(overall_rating) %>%
    mutate(pct = n / sum(n) * 100, country = "United States") %>% ungroup()
)

p4 <- ggplot(fig4_data, aes(x = factor(overall_rating), y = pct, fill = country)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ ownership_std) +
  scale_fill_manual(values = c("Australia" = "#2171B5", "United States" = "#CB181D")) +
  labs(
    title = "Rating Distribution by Ownership Type: Australia vs United States",
    x = "Overall Star Rating", y = "Percentage (%)",
    fill = "Country"
  ) +
  theme_pub +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggsave(file.path(dir_fig, "intl_distribution_by_ownership.png"), p4,
       width = 12, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: intl_distribution_by_ownership.png")

# =============================================================================
# Part 7: Methodological differences table
# =============================================================================

message("\n=== Part 7: Methodological differences ===\n")

method_diff <- tibble(
  feature = c(
    "Year introduced",
    "Facility type",
    "Number of facilities",
    "Rating scale",
    "Sub-categories",
    "Overall rating method",
    "Health inspection / Compliance",
    "Staffing basis",
    "Quality measures basis",
    "Consumer experience",
    "Inspection rating method",
    "For-profit share",
    "Data update frequency"
  ),
  australia = c(
    "Late 2022",
    "Residential aged care (all types)",
    as.character(nrow(au)),
    "1-5 stars",
    "Residents' Experience, Compliance, Staffing, Quality Measures",
    "Weighted composite of 4 sub-categories",
    "Based on regulatory audit outcomes",
    "Care minutes per resident day (mandated minimum)",
    "National Aged Care Mandatory Quality Indicator Programme",
    "Yes — Residents' Experience sub-category (consumer survey)",
    "Absolute standards (pass/fail compliance)",
    paste0(round(sum(au$ownership_std == "For-profit") / nrow(au) * 100, 0), "%"),
    "Quarterly"
  ),
  united_states = c(
    "2008",
    "Medicare/Medicaid-certified nursing homes (SNFs)",
    as.character(nrow(us)),
    "1-5 stars",
    "Health Inspection, Staffing, Quality Measures",
    "Starts from health inspection, adjusted by staffing and QM",
    "State survey agency inspections (unannounced)",
    "Payroll-Based Journal (PBJ) — verified payroll hours",
    "MDS assessments + Medicare claims data",
    "No — no direct consumer experience sub-category",
    "Within-state relative ranking (percentile cutoffs)",
    paste0(round(sum(us$ownership_std == "For-profit") / nrow(us) * 100, 0), "%"),
    "Quarterly (approx. Jan, Apr, Jul, Oct)"
  )
)
write_csv(method_diff, file.path(dir_tab, "intl_methodological_differences.csv"))
message("  Saved: intl_methodological_differences.csv")

# Key methodological distinctions for manuscript
message("\n  Key methodological differences:")
message("    1. US health inspection rating is RELATIVE (within-state percentiles);")
message("       AU compliance is ABSOLUTE (standards-based)")
message("    2. US overall = health inspection base + staffing/QM adjustments;")
message("       AU overall = weighted composite of 4 sub-categories")
message("    3. AU includes consumer experience survey; US does not")
message("    4. US for-profit share (~73%) is much higher than AU (~34%)")
message("    5. AU system is 3 years old; US system is 17+ years old")
message("    6. US uses relative ranking so distribution is ~fixed by design;")
message("       AU uses absolute standards so the distribution can shift over time")

# =============================================================================
# Summary
# =============================================================================

message("\n=== Phase 5 Complete ===\n")
message(glue("  Australian facilities: {nrow(au)} ({latest_label})"))
message(glue("  US facilities: {nrow(us)} (Jan 2026)"))
message(glue("  AU mean overall: {au_mean} (SD {au_sd})"))
message(glue("  US mean overall: {us_mean} (SD {us_sd})"))
message(glue("  AU for-profit share: {round(sum(au$ownership_std == 'For-profit') / nrow(au) * 100, 1)}%"))
message(glue("  US for-profit share: {round(sum(us$ownership_std == 'For-profit') / nrow(us) * 100, 1)}%"))
message("\n  Figures: outputs/figures/")
message("  Tables: outputs/tables/")
message("\nDone.")
