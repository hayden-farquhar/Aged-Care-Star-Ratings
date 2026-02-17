# =============================================================================
# 05_equity_synthesis.R
# Equity synthesis and quality desert identification
#
# Analyses:
#   1. Compute distance from each SA3 centroid to nearest 3+ star facility
#   2. Identify quality deserts (nearest 3+ star facility >50km)
#   3. Test overlap: quality deserts x disadvantage x aged population
#   4. Facility-level equity analysis (disadvantage x provider type x rating)
#   5. Quality desert maps with disadvantage overlay
#   6. Equity synthesis summary
#
# Inputs:
#   data/processed/star_ratings_panel.rds
#   data/processed/sa3_reference.rds
#   data/spatial/SA3_2021_AUST_SHP_GDA2020/
#
# Outputs:
#   outputs/maps/     — quality desert maps
#   outputs/figures/  — equity plots
#   outputs/tables/   — desert list, overlap analysis
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tmap)
library(glue)
library(purrr)

project_root <- here::here()
dir_proc   <- file.path(project_root, "data", "processed")
dir_figs   <- file.path(project_root, "outputs", "figures")
dir_tables <- file.path(project_root, "outputs", "tables")
dir_maps   <- file.path(project_root, "outputs", "maps")

save_map <- function(map, filename, ...) {
  tmap_save(map, filename = filename, device = ragg::agg_png, ...)
}

# =============================================================================
# PART 1: Load data and prepare facility locations
# =============================================================================

message("\n=== Part 1: Data preparation ===\n")

panel   <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
sa3_ref <- readRDS(file.path(dir_proc, "sa3_reference.rds"))

shp <- st_read(file.path(project_root, "data", "spatial",
                          "SA3_2021_AUST_SHP_GDA2020"), quiet = TRUE)

# Use most recent quarter
latest_q <- max(panel$quarter_num)
latest <- panel %>%
  filter(quarter_num == latest_q, !is.na(sa3_code), !is.na(overall_star_rating))

message(glue("  Latest quarter: {unique(latest$quarter_label)}"))
message(glue("  Facilities with ratings and SA3: {nrow(latest)}"))

# Get SA3 centroids (projected to Australian Albers for distance calculations)
shp_proj <- st_transform(shp, crs = 3577)
sa3_centroids <- st_centroid(shp_proj)

# Create facility-level spatial data using SA3 centroids as proxy locations
# (we don't have exact facility coordinates — use SA3 centroid)
facility_sf <- latest %>%
  left_join(
    sa3_centroids %>%
      select(SA3_CODE21) %>%
      mutate(centroid_geom = geometry),
    by = c("sa3_code" = "SA3_CODE21")
  ) %>%
  filter(!is.na(centroid_geom)) %>%
  st_as_sf(sf_column_name = "centroid_geom")

st_crs(facility_sf) <- 3577

message(glue("  Facilities with spatial data: {nrow(facility_sf)}"))
message(glue("  Facilities rated 3+ stars: {sum(facility_sf$overall_star_rating >= 3)}"))


# =============================================================================
# PART 2: Distance to nearest 3+ star facility from each SA3
# =============================================================================

message("\n=== Part 2: Distance to nearest adequate facility ===\n")

# Identify facilities rated 3+ stars (adequate)
adequate_facilities <- facility_sf %>% filter(overall_star_rating >= 3)
message(glue("  Adequate facilities (3+ stars): {nrow(adequate_facilities)}"))

# Compute distance from each SA3 centroid to nearest adequate facility
sa3_centroids_with_data <- sa3_centroids %>%
  left_join(sa3_ref, by = c("SA3_CODE21" = "sa3_code"))

# Distance matrix: each SA3 centroid to all adequate facilities
dist_matrix <- st_distance(sa3_centroids_with_data, adequate_facilities)

# Nearest adequate facility distance (in km)
sa3_centroids_with_data$dist_nearest_adequate_km <- apply(dist_matrix, 1, min) / 1000

# Also compute distance to nearest facility of any rating
all_dist <- st_distance(sa3_centroids_with_data, facility_sf)
sa3_centroids_with_data$dist_nearest_any_km <- apply(all_dist, 1, min) / 1000

# Count facilities within 50km
count_within <- function(dist_row, threshold_m) {
  sum(dist_row <= threshold_m)
}
sa3_centroids_with_data$n_adequate_within_50km <- apply(dist_matrix, 1, count_within, threshold_m = 50000)
sa3_centroids_with_data$n_any_within_50km <- apply(all_dist, 1, count_within, threshold_m = 50000)

message("  Distance to nearest adequate (3+) facility:")
dist_summary <- summary(sa3_centroids_with_data$dist_nearest_adequate_km)
message(glue("    Median: {round(dist_summary['Median'], 1)} km"))
message(glue("    Mean: {round(dist_summary['Mean'], 1)} km"))
message(glue("    Max: {round(dist_summary['Max.'], 1)} km"))
message(glue("    SA3s with nearest >50km: {sum(sa3_centroids_with_data$dist_nearest_adequate_km > 50)}"))
message(glue("    SA3s with nearest >100km: {sum(sa3_centroids_with_data$dist_nearest_adequate_km > 100)}"))


# =============================================================================
# PART 3: Identify quality deserts
# =============================================================================

message("\n=== Part 3: Quality deserts ===\n")

# Quality desert: SA3 where nearest 3+ star facility is >50km away
desert_threshold <- 50  # km

sa3_desert <- sa3_centroids_with_data %>%
  st_drop_geometry() %>%
  mutate(
    is_quality_desert = dist_nearest_adequate_km > desert_threshold,
    is_facility_desert = dist_nearest_any_km > desert_threshold
  )

n_quality_deserts   <- sum(sa3_desert$is_quality_desert, na.rm = TRUE)
n_facility_deserts  <- sum(sa3_desert$is_facility_desert, na.rm = TRUE)

message(glue("  Quality deserts (nearest 3+ star >50km): {n_quality_deserts} SA3 areas"))
message(glue("  Facility deserts (nearest any facility >50km): {n_facility_deserts} SA3 areas"))

# Population in quality deserts
pop_in_deserts <- sa3_desert %>%
  filter(is_quality_desert) %>%
  summarise(
    total_pop    = sum(pop_total, na.rm = TRUE),
    pop_65plus   = sum(pop_65plus, na.rm = TRUE),
    pop_85plus   = sum(pop_85plus, na.rm = TRUE)
  )

pop_total_all <- sum(sa3_desert$pop_total, na.rm = TRUE)
pop_65_all    <- sum(sa3_desert$pop_65plus, na.rm = TRUE)

message(glue("  Population in quality deserts: {format(pop_in_deserts$total_pop, big.mark = ',')} ({round(pop_in_deserts$total_pop/pop_total_all*100, 1)}% of total)"))
message(glue("  65+ population in quality deserts: {format(pop_in_deserts$pop_65plus, big.mark = ',')} ({round(pop_in_deserts$pop_65plus/pop_65_all*100, 1)}% of 65+ total)"))


# =============================================================================
# PART 4: Overlap with disadvantage and remoteness
# =============================================================================

message("\n=== Part 4: Quality deserts x disadvantage x remoteness ===\n")

# --- 4a: Quality deserts by IRSD quintile ---
desert_by_irsd <- sa3_desert %>%
  filter(!is.na(irsd_quintile)) %>%
  group_by(irsd_quintile) %>%
  summarise(
    n_sa3 = n(),
    n_deserts = sum(is_quality_desert, na.rm = TRUE),
    pct_desert = round(mean(is_quality_desert, na.rm = TRUE) * 100, 1),
    mean_dist = round(mean(dist_nearest_adequate_km, na.rm = TRUE), 1),
    .groups = "drop"
  )

message("  Quality deserts by IRSD quintile (1=most disadvantaged):")
for (i in seq_len(nrow(desert_by_irsd))) {
  message(glue("    Q{desert_by_irsd$irsd_quintile[i]}: {desert_by_irsd$n_deserts[i]}/{desert_by_irsd$n_sa3[i]} SA3s are deserts ({desert_by_irsd$pct_desert[i]}%), mean dist = {desert_by_irsd$mean_dist[i]} km"))
}

write.csv(desert_by_irsd, file.path(dir_tables, "quality_deserts_by_irsd.csv"),
          row.names = FALSE)

# Chi-squared test: is desert status associated with IRSD quintile?
desert_irsd_table <- sa3_desert %>%
  filter(!is.na(irsd_quintile)) %>%
  mutate(irsd_q = factor(irsd_quintile)) %>%
  with(table(irsd_q, is_quality_desert))

chi_irsd <- chisq.test(desert_irsd_table)
message(glue("\n  Chi-squared (desert x IRSD quintile): X2 = {round(chi_irsd$statistic, 2)}, p = {format.pval(chi_irsd$p.value, digits = 3)}"))

# --- 4b: Quality deserts by remoteness ---
desert_by_remote <- sa3_desert %>%
  filter(!is.na(remoteness)) %>%
  group_by(remoteness) %>%
  summarise(
    n_sa3 = n(),
    n_deserts = sum(is_quality_desert, na.rm = TRUE),
    pct_desert = round(mean(is_quality_desert, na.rm = TRUE) * 100, 1),
    mean_dist = round(mean(dist_nearest_adequate_km, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_desert))

message("\n  Quality deserts by remoteness:")
for (i in seq_len(nrow(desert_by_remote))) {
  message(glue("    {desert_by_remote$remoteness[i]}: {desert_by_remote$n_deserts[i]}/{desert_by_remote$n_sa3[i]} SA3s ({desert_by_remote$pct_desert[i]}%), mean dist = {desert_by_remote$mean_dist[i]} km"))
}

write.csv(desert_by_remote, file.path(dir_tables, "quality_deserts_by_remoteness.csv"),
          row.names = FALSE)

# --- 4c: Quality deserts by aged population proportion ---
desert_by_aged <- sa3_desert %>%
  filter(!is.na(pct_65plus)) %>%
  mutate(aged_tertile = ntile(pct_65plus, 3)) %>%
  group_by(aged_tertile) %>%
  summarise(
    n_sa3 = n(),
    n_deserts = sum(is_quality_desert, na.rm = TRUE),
    pct_desert = round(mean(is_quality_desert, na.rm = TRUE) * 100, 1),
    mean_pct_65plus = round(mean(pct_65plus, na.rm = TRUE), 1),
    .groups = "drop"
  )

message("\n  Quality deserts by aged population tertile:")
for (i in seq_len(nrow(desert_by_aged))) {
  message(glue("    T{desert_by_aged$aged_tertile[i]} (mean {desert_by_aged$mean_pct_65plus[i]}% aged 65+): {desert_by_aged$n_deserts[i]}/{desert_by_aged$n_sa3[i]} ({desert_by_aged$pct_desert[i]}%)"))
}

# --- 4d: "Double disadvantage" — desert AND high disadvantage ---
double_disadv <- sa3_desert %>%
  filter(!is.na(irsd_quintile)) %>%
  mutate(
    double_disadvantage = is_quality_desert & irsd_quintile <= 2,
    triple_disadvantage = is_quality_desert & irsd_quintile <= 2 & !is.na(pct_65plus) & pct_65plus > median(pct_65plus, na.rm = TRUE)
  )

n_double <- sum(double_disadv$double_disadvantage, na.rm = TRUE)
n_triple <- sum(double_disadv$triple_disadvantage, na.rm = TRUE)

message(glue("\n  Double disadvantage (desert + IRSD Q1-Q2): {n_double} SA3 areas"))
message(glue("  Triple disadvantage (desert + IRSD Q1-Q2 + above-median aged pop): {n_triple} SA3 areas"))

# Export quality desert list
desert_list <- sa3_desert %>%
  filter(is_quality_desert) %>%
  select(SA3_CODE21, SA3_NAME21, STE_NAME21,
         dist_nearest_adequate_km, dist_nearest_any_km,
         n_adequate_within_50km, n_any_within_50km,
         irsd_score_mean, irsd_quintile, remoteness,
         pop_total, pop_65plus, pct_65plus) %>%
  arrange(desc(dist_nearest_adequate_km))

write.csv(desert_list, file.path(dir_tables, "quality_desert_sa3_list.csv"),
          row.names = FALSE)
message(glue("\n  Saved: quality_desert_sa3_list.csv ({nrow(desert_list)} SA3 areas)"))


# =============================================================================
# PART 5: Quality desert maps
# =============================================================================

message("\n=== Part 5: Quality desert maps ===\n")

tmap_mode("plot")

# Prepare spatial data for mapping
shp_desert <- shp_proj %>%
  left_join(
    sa3_desert %>% select(SA3_CODE21, dist_nearest_adequate_km, is_quality_desert,
                          n_adequate_within_50km, irsd_quintile, remoteness,
                          pop_65plus, pct_65plus),
    by = "SA3_CODE21"
  )

# --- 5a: Distance to nearest adequate facility ---
map_dist <- tm_shape(shp_desert) +
  tm_fill(
    fill = "dist_nearest_adequate_km",
    fill.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
    fill.legend = tm_legend(title = "Distance to\nnearest 3+ star\nfacility (km)")
  ) +
  tm_borders(col = "grey70", lwd = 0.2) +
  tm_title("Distance to Nearest Adequate (3+ Star) Aged Care Facility") +
  tm_layout(frame = FALSE)

save_map(map_dist, file.path(dir_maps, "map_distance_to_adequate.png"),
         width = 10, height = 8, dpi = 300)
message("  Saved: map_distance_to_adequate.png")

# --- 5b: Quality deserts highlighted ---
shp_desert <- shp_desert %>%
  mutate(desert_label = case_when(
    is_quality_desert ~ "Quality desert (>50km)",
    !is.na(dist_nearest_adequate_km) ~ "Within 50km of 3+ star facility",
    TRUE ~ NA_character_
  ))

map_deserts <- tm_shape(shp_desert) +
  tm_fill(
    fill = "desert_label",
    fill.scale = tm_scale_categorical(
      values = c("Quality desert (>50km)" = "#d73027",
                 "Within 50km of 3+ star facility" = "#91cf60")
    ),
    fill.legend = tm_legend(title = "Quality Desert Status")
  ) +
  tm_borders(col = "grey70", lwd = 0.2) +
  tm_title("Aged Care Quality Deserts") +
  tm_layout(frame = FALSE)

save_map(map_deserts, file.path(dir_maps, "map_quality_deserts.png"),
         width = 10, height = 8, dpi = 300)
message("  Saved: map_quality_deserts.png")

# --- 5c: Quality deserts with IRSD disadvantage overlay ---
shp_desert_overlay <- shp_desert %>%
  mutate(
    desert_x_disadvantage = case_when(
      is_quality_desert & !is.na(irsd_quintile) & irsd_quintile <= 2 ~ "Desert + High Disadvantage",
      is_quality_desert ~ "Desert only",
      !is.na(irsd_quintile) & irsd_quintile <= 2 ~ "High Disadvantage only",
      TRUE ~ "Neither"
    )
  )

map_overlay <- tm_shape(shp_desert_overlay) +
  tm_fill(
    fill = "desert_x_disadvantage",
    fill.scale = tm_scale_categorical(
      values = c(
        "Desert + High Disadvantage" = "#d73027",
        "Desert only" = "#fc8d59",
        "High Disadvantage only" = "#4575b4",
        "Neither" = "#e0e0e0"
      ),
      levels = c("Desert + High Disadvantage", "Desert only",
                 "High Disadvantage only", "Neither")
    ),
    fill.legend = tm_legend(title = "Equity Classification")
  ) +
  tm_borders(col = "grey70", lwd = 0.2) +
  tm_title("Quality Deserts and Socioeconomic Disadvantage") +
  tm_layout(frame = FALSE)

save_map(map_overlay, file.path(dir_maps, "map_quality_desert_disadvantage.png"),
         width = 10, height = 8, dpi = 300)
message("  Saved: map_quality_desert_disadvantage.png")

# --- 5d: Number of adequate facilities within 50km ---
shp_access <- shp_desert %>%
  filter(!is.na(n_adequate_within_50km))

map_access <- tm_shape(shp_proj) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) +
  tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_access) +
  tm_fill(
    fill = "n_adequate_within_50km",
    fill.scale = tm_scale_continuous(values = "brewer.yl_gn_bu"),
    fill.legend = tm_legend(title = "3+ star facilities\nwithin 50km")
  ) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title("Access to Adequate Aged Care Facilities") +
  tm_layout(frame = FALSE)

save_map(map_access, file.path(dir_maps, "map_access_adequate_facilities.png"),
         width = 10, height = 8, dpi = 300)
message("  Saved: map_access_adequate_facilities.png")


# =============================================================================
# PART 6: Equity figures
# =============================================================================

message("\n=== Part 6: Equity figures ===\n")

# --- 6a: Distance to adequate facility by IRSD quintile ---
plot_data <- sa3_desert %>%
  filter(!is.na(irsd_quintile), !is.na(dist_nearest_adequate_km))

p_dist_irsd <- ggplot(plot_data, aes(x = factor(irsd_quintile),
                                      y = dist_nearest_adequate_km)) +
  geom_boxplot(fill = "#a6cee3", outlier.alpha = 0.3) +
  scale_y_log10() +
  labs(
    x = "IRSD Quintile (1 = Most Disadvantaged)",
    y = "Distance to nearest 3+ star facility (km, log scale)",
    title = "Access to Adequate Aged Care by Disadvantage Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(dir_figs, "distance_by_irsd_quintile.png"), p_dist_irsd,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: distance_by_irsd_quintile.png")

# --- 6b: Distance by remoteness ---
plot_remote <- sa3_desert %>%
  filter(!is.na(remoteness), !is.na(dist_nearest_adequate_km))

p_dist_remote <- ggplot(plot_remote, aes(x = reorder(remoteness, dist_nearest_adequate_km),
                                          y = dist_nearest_adequate_km)) +
  geom_boxplot(fill = "#b2df8a", outlier.alpha = 0.3) +
  scale_y_log10() +
  labs(
    x = "Remoteness Category",
    y = "Distance to nearest 3+ star facility (km, log scale)",
    title = "Access to Adequate Aged Care by Remoteness"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(dir_figs, "distance_by_remoteness.png"), p_dist_remote,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: distance_by_remoteness.png")

# --- 6c: Scatter: % 65+ population vs distance to adequate facility ---
p_aged_dist <- ggplot(plot_data, aes(x = pct_65plus, y = dist_nearest_adequate_km)) +
  geom_point(aes(colour = factor(irsd_quintile)), alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, colour = "black", linewidth = 0.8) +
  scale_y_log10() +
  scale_colour_brewer(palette = "RdYlGn", direction = -1, name = "IRSD Quintile\n(1=Most Disadv.)") +
  labs(
    x = "% Population Aged 65+",
    y = "Distance to nearest 3+ star facility (km, log scale)",
    title = "Aged Population and Access to Quality Aged Care"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "right")

ggsave(file.path(dir_figs, "aged_population_vs_distance.png"), p_aged_dist,
       width = 9, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: aged_population_vs_distance.png")


# =============================================================================
# PART 7: Equity synthesis — combining geographic and longitudinal findings
# =============================================================================

message("\n=== Part 7: Equity synthesis ===\n")

# Combine key equity indicators into a summary table
equity_summary <- tibble(
  indicator = c(
    "Total SA3 areas analysed",
    "SA3 areas that are quality deserts (>50km to 3+ star)",
    "SA3 areas that are facility deserts (>50km to any facility)",
    "Population in quality deserts",
    "65+ population in quality deserts",
    "Quality deserts in IRSD Q1-Q2 (most disadvantaged)",
    "Quality deserts in Remote/Very Remote areas",
    "Double disadvantage (desert + high disadvantage)",
    "Triple disadvantage (desert + high disadvantage + high aged pop)",
    "Persistent low performers (<=2 stars, 3+ quarters)",
    "For-profit share of persistent low performers"
  ),
  value = c(
    nrow(sa3_desert),
    n_quality_deserts,
    n_facility_deserts,
    format(pop_in_deserts$total_pop, big.mark = ","),
    format(pop_in_deserts$pop_65plus, big.mark = ","),
    sum(sa3_desert$is_quality_desert & !is.na(sa3_desert$irsd_quintile) & sa3_desert$irsd_quintile <= 2, na.rm = TRUE),
    sum(sa3_desert$is_quality_desert & sa3_desert$remoteness %in% c("Remote Australia", "Very Remote Australia"), na.rm = TRUE),
    n_double,
    n_triple,
    28,
    "50% (14/28)"
  )
)

write.csv(equity_summary, file.path(dir_tables, "equity_synthesis_summary.csv"),
          row.names = FALSE)
message("  Saved: equity_synthesis_summary.csv")

# Print full summary
message("\n  Equity Synthesis Summary:")
for (i in seq_len(nrow(equity_summary))) {
  message(glue("    {equity_summary$indicator[i]}: {equity_summary$value[i]}"))
}


# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Phase 4 Complete ===\n")
message("  Maps: outputs/maps/")
message("  Figures: outputs/figures/")
message("  Tables: outputs/tables/")
message("\nDone.")
