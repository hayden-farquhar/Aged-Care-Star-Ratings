# =============================================================================
# 03_spatial_analysis.R
# Cross-sectional geographic analysis of aged care Star Ratings
#
# Analyses (using most recent quarter):
#   1. SA3-level aggregation of facility ratings
#   2. Choropleth maps of overall and sub-category ratings
#   3. Global Moran's I for spatial autocorrelation
#   4. Bivariate analysis: SA3 mean rating vs SEIFA IRSD
#   5. Geographically Weighted Regression (GWR)
#   6. Sub-domain decomposition (repeat spatial tests for each sub-category)
#
# Inputs:
#   data/processed/star_ratings_panel.rds
#   data/processed/sa3_reference.rds
#   data/spatial/SA3_2021_AUST_SHP_GDA2020/
#
# Outputs:
#   outputs/figures/    — choropleth maps, correlation plots
#   outputs/tables/     — Moran's I results, GWR summary, bivariate tables
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tmap)
library(spdep)
library(GWmodel)
library(glue)
library(purrr)

project_root <- here::here()
dir_proc     <- file.path(project_root, "data", "processed")
dir_figs     <- file.path(project_root, "outputs", "figures")
dir_tables   <- file.path(project_root, "outputs", "tables")
dir_maps     <- file.path(project_root, "outputs", "maps")

dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_maps, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PART 1: Load data and aggregate to SA3
# =============================================================================

message("\n=== Part 1: SA3-level aggregation ===\n")

panel   <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
sa3_ref <- readRDS(file.path(dir_proc, "sa3_reference.rds"))

# Use most recent quarter for cross-sectional analysis
latest_q     <- max(panel$quarter_num)
latest_label <- panel %>% filter(quarter_num == latest_q) %>% pull(quarter_label) %>% unique()
message(glue("  Using latest quarter: {latest_label} (Q{latest_q})"))

latest <- panel %>%
  filter(quarter_num == latest_q, !is.na(sa3_code))

# Aggregate to SA3
rating_vars <- c("overall_star_rating", "residents_experience_rating",
                 "compliance_rating", "staffing_rating", "quality_measures_rating")

sa3_ratings <- latest %>%
  filter(!is.na(overall_star_rating)) %>%
  group_by(sa3_code, sa3_name) %>%
  summarise(
    n_facilities = n(),
    across(all_of(rating_vars), list(
      mean = ~mean(.x, na.rm = TRUE),
      sd   = ~sd(.x, na.rm = TRUE),
      min  = ~min(.x, na.rm = TRUE),
      max  = ~max(.x, na.rm = TRUE)
    ), .names = "{.col}_{.fn}"),
    pct_low     = mean(overall_star_rating <= 2, na.rm = TRUE) * 100,
    pct_high    = mean(overall_star_rating >= 4, na.rm = TRUE) * 100,
    n_for_profit     = sum(purpose == "For-profit", na.rm = TRUE),
    n_not_for_profit = sum(purpose == "Not-for-profit", na.rm = TRUE),
    n_government     = sum(purpose == "Government", na.rm = TRUE),
    .groups = "drop"
  )

# Join with SA3 reference data (SEIFA, remoteness, population)
sa3_ratings <- sa3_ratings %>%
  left_join(sa3_ref, by = c("sa3_code", "sa3_name"))

message(glue("  SA3 areas with facilities: {nrow(sa3_ratings)}"))
message(glue("  Total facilities: {sum(sa3_ratings$n_facilities)}"))
message(glue("  SA3 areas with >=3 facilities: {sum(sa3_ratings$n_facilities >= 3)}"))

# Save SA3-level aggregated data
write.csv(sa3_ratings, file.path(dir_tables, "sa3_ratings_latest_quarter.csv"),
          row.names = FALSE)


# =============================================================================
# PART 2: Load shapefile and join
# =============================================================================

message("\n=== Part 2: Preparing spatial data ===\n")

shp <- st_read(file.path(project_root, "data", "spatial",
                          "SA3_2021_AUST_SHP_GDA2020"), quiet = TRUE)

# Join ratings to shapefile
shp_ratings <- shp %>%
  left_join(sa3_ratings, by = c("SA3_CODE21" = "sa3_code"))

# Create a version with only SA3s that have facilities
shp_with_data <- shp_ratings %>% filter(!is.na(n_facilities))

message(glue("  SA3 areas in shapefile: {nrow(shp)}"))
message(glue("  SA3 areas with rating data: {nrow(shp_with_data)}"))

# Project to Australian Albers (EPSG:3577) for distance-based analyses
shp_proj <- st_transform(shp_with_data, crs = 3577)


# =============================================================================
# PART 3: Choropleth maps
# =============================================================================

message("\n=== Part 3: Generating choropleth maps ===\n")

tmap_mode("plot")

# --- 3a: Overall Star Rating map ---
map_overall <- tm_shape(shp, projection = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) +
  tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_with_data, projection = 3577) +
  tm_fill(
    fill = "overall_star_rating_mean",
    fill.scale = tm_scale_continuous(
      values = "viridis",
      midpoint = 3.5
    ),
    fill.legend = tm_legend(title = "Mean Star Rating")
  ) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title(glue("Overall Star Rating by SA3 — {latest_label}")) +
  tm_layout(frame = FALSE)

tmap_save(map_overall,
          filename = file.path(dir_maps, "map_overall_star_rating.png"),
          width = 10, height = 8, dpi = 300)
tmap_save(map_overall,
          filename = file.path(dir_maps, "map_overall_star_rating.pdf"),
          width = 10, height = 8)
message("  Saved: map_overall_star_rating.png/pdf")

# --- 3b: Sub-category rating maps ---
subcategory_labels <- c(
  "residents_experience_rating_mean" = "Residents' Experience",
  "compliance_rating_mean"           = "Compliance",
  "staffing_rating_mean"             = "Staffing",
  "quality_measures_rating_mean"     = "Quality Measures"
)

for (var in names(subcategory_labels)) {
  label <- subcategory_labels[var]
  fname <- gsub("_mean$", "", var) %>% gsub("_", "-", .)

  map_sub <- tm_shape(shp, projection = 3577) +
    tm_fill(fill = "grey90", fill_alpha = 0.5) +
    tm_borders(col = "grey80", lwd = 0.2) +
    tm_shape(shp_with_data, projection = 3577) +
    tm_fill(
      fill = var,
      fill.scale = tm_scale_continuous(
        values = "viridis",
        midpoint = 3.5
      ),
      fill.legend = tm_legend(title = paste("Mean", label, "Rating"))
    ) +
    tm_borders(col = "grey60", lwd = 0.3) +
    tm_title(glue("{label} Rating by SA3 — {latest_label}")) +
    tm_layout(frame = FALSE)

  tmap_save(map_sub,
            filename = file.path(dir_maps, glue("map_{fname}.png")),
            width = 10, height = 8, dpi = 300)
  message(glue("  Saved: map_{fname}.png"))
}

# --- 3c: IRSD quintile map with rating overlay ---
map_irsd <- tm_shape(shp, projection = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) +
  tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_with_data, projection = 3577) +
  tm_fill(
    fill = "irsd_quintile",
    fill.scale = tm_scale_continuous(
      values = "-RdYlGn",
      midpoint = 3
    ),
    fill.legend = tm_legend(title = "IRSD Quintile\n(1=Most Disadvantaged)")
  ) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title("SEIFA IRSD Quintile by SA3") +
  tm_layout(frame = FALSE)

tmap_save(map_irsd,
          filename = file.path(dir_maps, "map_irsd_quintile.png"),
          width = 10, height = 8, dpi = 300)
message("  Saved: map_irsd_quintile.png")


# =============================================================================
# PART 4: Global Moran's I — spatial autocorrelation
# =============================================================================

message("\n=== Part 4: Spatial autocorrelation (Moran's I) ===\n")

# Use SA3 areas with at least 3 facilities for robust mean estimates
sa3_robust <- shp_proj %>%
  filter(n_facilities >= 3, !is.na(overall_star_rating_mean))

message(glue("  SA3 areas with >=3 facilities: {nrow(sa3_robust)}"))

# Create spatial weights (queen contiguity)
nb_queen <- poly2nb(sa3_robust, queen = TRUE)

# Check for islands (disconnected areas)
n_islands <- sum(card(nb_queen) == 0)
message(glue("  Islands (no neighbours): {n_islands}"))

# If islands exist, use k-nearest neighbours instead
if (n_islands > 0) {
  message("  Using k-nearest neighbours (k=5) due to islands")
  coords <- st_centroid(sa3_robust) %>% st_coordinates()
  nb_knn <- knn2nb(knearneigh(coords, k = 5))
  listw <- nb2listw(nb_knn, style = "W")
} else {
  listw <- nb2listw(nb_queen, style = "W")
}

# Global Moran's I for overall rating
moran_overall <- moran.test(sa3_robust$overall_star_rating_mean, listw)
message(glue("  Overall Star Rating — Moran's I: {round(moran_overall$estimate['Moran I statistic'], 4)}, p = {format.pval(moran_overall$p.value, digits = 3)}"))

# Moran's I for each sub-category
moran_results <- tibble(
  variable = c("Overall Star Rating", "Residents' Experience",
               "Compliance", "Staffing", "Quality Measures"),
  column   = c("overall_star_rating_mean", "residents_experience_rating_mean",
               "compliance_rating_mean", "staffing_rating_mean",
               "quality_measures_rating_mean")
) %>%
  rowwise() %>%
  mutate(
    morans_i = {
      vals <- sa3_robust[[column]]
      if (sum(!is.na(vals)) > 10) {
        test <- moran.test(vals, listw, na.action = na.exclude)
        test$estimate["Moran I statistic"]
      } else NA_real_
    },
    p_value = {
      vals <- sa3_robust[[column]]
      if (sum(!is.na(vals)) > 10) {
        test <- moran.test(vals, listw, na.action = na.exclude)
        test$p.value
      } else NA_real_
    },
    interpretation = case_when(
      is.na(p_value) ~ "Insufficient data",
      p_value > 0.05 ~ "No significant spatial autocorrelation",
      morans_i > 0   ~ "Significant positive spatial autocorrelation (clustering)",
      morans_i < 0   ~ "Significant negative spatial autocorrelation (dispersion)",
      TRUE           ~ "No spatial autocorrelation"
    )
  ) %>%
  ungroup()

message("\n  Moran's I results:")
for (i in seq_len(nrow(moran_results))) {
  message(glue("    {moran_results$variable[i]}: I = {round(moran_results$morans_i[i], 4)}, p = {format.pval(moran_results$p_value[i], digits = 3)} — {moran_results$interpretation[i]}"))
}

write.csv(moran_results %>% select(-column),
          file.path(dir_tables, "morans_i_results.csv"), row.names = FALSE)
message("  Saved: morans_i_results.csv")

# Moran scatterplot for overall rating
png(file.path(dir_figs, "moran_scatterplot_overall.png"),
    width = 7, height = 6, units = "in", res = 300)
moran.plot(sa3_robust$overall_star_rating_mean, listw,
           xlab = "Mean Overall Star Rating (SA3)",
           ylab = "Spatially Lagged Mean Rating",
           main = "Moran Scatterplot — Overall Star Rating")
dev.off()
message("  Saved: moran_scatterplot_overall.png")


# =============================================================================
# PART 5: Bivariate analysis — Star Rating vs SEIFA IRSD
# =============================================================================

message("\n=== Part 5: Bivariate analysis (Rating vs IRSD) ===\n")

# Use SA3 areas with >=3 facilities and valid IRSD
sa3_bivar <- sa3_ratings %>%
  filter(n_facilities >= 3, !is.na(irsd_score_mean), !is.na(overall_star_rating_mean))

message(glue("  SA3 areas for bivariate analysis: {nrow(sa3_bivar)}"))

# --- 5a: Pearson correlation ---
cor_test <- cor.test(sa3_bivar$irsd_score_mean, sa3_bivar$overall_star_rating_mean)
message(glue("  Pearson r = {round(cor_test$estimate, 4)}, p = {format.pval(cor_test$p.value, digits = 3)}"))

# --- 5b: Spearman correlation ---
cor_spearman <- cor.test(sa3_bivar$irsd_score_mean, sa3_bivar$overall_star_rating_mean,
                         method = "spearman")
message(glue("  Spearman rho = {round(cor_spearman$estimate, 4)}, p = {format.pval(cor_spearman$p.value, digits = 3)}"))

# --- 5c: Mean rating by IRSD quintile (ANOVA) ---
sa3_bivar$irsd_q_factor <- factor(sa3_bivar$irsd_quintile,
                                  labels = paste0("Q", 1:5, " (", c("Most", "2nd", "3rd", "4th", "Least"), " disadvantaged)"))

aov_result <- aov(overall_star_rating_mean ~ irsd_q_factor, data = sa3_bivar)
aov_summary <- summary(aov_result)
f_stat <- aov_summary[[1]]["irsd_q_factor", "F value"]
p_val  <- aov_summary[[1]]["irsd_q_factor", "Pr(>F)"]
message(glue("  ANOVA F = {round(f_stat, 2)}, p = {format.pval(p_val, digits = 3)}"))

by_quintile_sa3 <- sa3_bivar %>%
  group_by(irsd_quintile) %>%
  summarise(
    n_sa3          = n(),
    n_facilities   = sum(n_facilities),
    mean_rating    = round(mean(overall_star_rating_mean, na.rm = TRUE), 3),
    sd_rating      = round(sd(overall_star_rating_mean, na.rm = TRUE), 3),
    .groups = "drop"
  )

message("\n  Mean SA3-level rating by IRSD quintile:")
for (i in seq_len(nrow(by_quintile_sa3))) {
  message(glue("    Q{by_quintile_sa3$irsd_quintile[i]}: mean = {by_quintile_sa3$mean_rating[i]} (SD {by_quintile_sa3$sd_rating[i]}), n_SA3 = {by_quintile_sa3$n_sa3[i]}, n_facilities = {by_quintile_sa3$n_facilities[i]}"))
}

write.csv(by_quintile_sa3, file.path(dir_tables, "rating_by_irsd_quintile_sa3.csv"),
          row.names = FALSE)

# --- 5d: Correlation for each sub-category ---
subcats <- c("residents_experience_rating_mean", "compliance_rating_mean",
             "staffing_rating_mean", "quality_measures_rating_mean")
subcat_labels <- c("Residents' Experience", "Compliance", "Staffing", "Quality Measures")

subcat_cors <- tibble(
  subcategory = c("Overall", subcat_labels),
  pearson_r   = NA_real_,
  pearson_p   = NA_real_,
  spearman_rho = NA_real_,
  spearman_p  = NA_real_
)

# Overall
subcat_cors$pearson_r[1]   <- cor_test$estimate
subcat_cors$pearson_p[1]   <- cor_test$p.value
subcat_cors$spearman_rho[1] <- cor_spearman$estimate
subcat_cors$spearman_p[1]  <- cor_spearman$p.value

for (i in seq_along(subcats)) {
  vals <- sa3_bivar[[subcats[i]]]
  valid <- !is.na(vals) & !is.na(sa3_bivar$irsd_score_mean)
  if (sum(valid) > 10) {
    ct <- cor.test(sa3_bivar$irsd_score_mean[valid], vals[valid])
    cs <- cor.test(sa3_bivar$irsd_score_mean[valid], vals[valid], method = "spearman")
    subcat_cors$pearson_r[i + 1]    <- ct$estimate
    subcat_cors$pearson_p[i + 1]    <- ct$p.value
    subcat_cors$spearman_rho[i + 1] <- cs$estimate
    subcat_cors$spearman_p[i + 1]   <- cs$p.value
  }
}

message("\n  Correlation of SA3 mean rating with IRSD (all sub-categories):")
for (i in seq_len(nrow(subcat_cors))) {
  message(glue("    {subcat_cors$subcategory[i]}: r = {round(subcat_cors$pearson_r[i], 4)}, p = {format.pval(subcat_cors$pearson_p[i], digits = 3)}"))
}

write.csv(subcat_cors, file.path(dir_tables, "irsd_rating_correlations.csv"),
          row.names = FALSE)
message("  Saved: irsd_rating_correlations.csv")

# --- 5e: Scatterplot ---
p_scatter <- ggplot(sa3_bivar, aes(x = irsd_score_mean, y = overall_star_rating_mean)) +
  geom_point(aes(size = n_facilities), alpha = 0.5, colour = "#2c7fb8") +
  geom_smooth(method = "lm", se = TRUE, colour = "red", linewidth = 0.8) +
  scale_size_continuous(range = c(1, 6), name = "Number of\nfacilities") +
  labs(
    x = "SEIFA IRSD Score (higher = less disadvantaged)",
    y = "Mean Overall Star Rating",
    title = glue("SA3 Mean Star Rating vs Socioeconomic Disadvantage — {latest_label}"),
    subtitle = glue("Pearson r = {round(cor_test$estimate, 3)}, p = {format.pval(cor_test$p.value, digits = 3)}")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(dir_figs, "scatter_rating_vs_irsd.png"), p_scatter,
       width = 8, height = 6, dpi = 300)
message("  Saved: scatter_rating_vs_irsd.png")

# --- 5f: Box plot by IRSD quintile ---
p_box <- ggplot(sa3_bivar, aes(x = factor(irsd_quintile), y = overall_star_rating_mean)) +
  geom_boxplot(fill = "#a6cee3", outlier.alpha = 0.4) +
  geom_jitter(aes(size = n_facilities), alpha = 0.3, width = 0.15) +
  scale_size_continuous(range = c(0.5, 4), name = "Facilities") +
  labs(
    x = "IRSD Quintile (1 = Most Disadvantaged)",
    y = "Mean Overall Star Rating (SA3)",
    title = glue("Star Rating by SEIFA Disadvantage Quintile — {latest_label}"),
    subtitle = glue("ANOVA F = {round(f_stat, 2)}, p = {format.pval(p_val, digits = 3)}")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(dir_figs, "boxplot_rating_by_irsd.png"), p_box,
       width = 8, height = 6, dpi = 300)
message("  Saved: boxplot_rating_by_irsd.png")

# --- 5g: By remoteness ---
sa3_remote <- sa3_ratings %>%
  filter(n_facilities >= 3, !is.na(remoteness), !is.na(overall_star_rating_mean))

by_remoteness_sa3 <- sa3_remote %>%
  group_by(remoteness) %>%
  summarise(
    n_sa3        = n(),
    n_facilities = sum(n_facilities),
    mean_rating  = round(mean(overall_star_rating_mean, na.rm = TRUE), 3),
    sd_rating    = round(sd(overall_star_rating_mean, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(mean_rating)

message("\n  Mean SA3-level rating by remoteness:")
for (i in seq_len(nrow(by_remoteness_sa3))) {
  message(glue("    {by_remoteness_sa3$remoteness[i]}: mean = {by_remoteness_sa3$mean_rating[i]} (n_SA3 = {by_remoteness_sa3$n_sa3[i]})"))
}

write.csv(by_remoteness_sa3, file.path(dir_tables, "rating_by_remoteness_sa3.csv"),
          row.names = FALSE)


# =============================================================================
# PART 6: Geographically Weighted Regression (GWR)
# =============================================================================

message("\n=== Part 6: Geographically Weighted Regression ===\n")

# Prepare spatial data for GWR
# Need: SpatialPointsDataFrame with centroids
sa3_gwr <- sa3_robust %>%
  filter(!is.na(irsd_score_mean))

message(glue("  SA3 areas for GWR: {nrow(sa3_gwr)}"))

# Convert to sp format (GWmodel requires sp objects)
sa3_gwr_sp <- as(sa3_gwr, "Spatial")
centroids_sf <- st_centroid(sa3_gwr)
coords <- st_coordinates(centroids_sf)
sa3_gwr_sp_pts <- SpatialPointsDataFrame(
  coords = coords,
  data = as.data.frame(sa3_gwr) %>% select(-geometry),
  proj4string = CRS(proj4string(sa3_gwr_sp))
)

# Select optimal bandwidth using cross-validation
message("  Selecting optimal bandwidth (adaptive) ...")
bw_adapt <- tryCatch({
  bw.gwr(overall_star_rating_mean ~ irsd_score_mean,
         data = sa3_gwr_sp_pts,
         approach = "AICc",
         kernel = "bisquare",
         adaptive = TRUE)
}, error = function(e) {
  message(glue("  [WARNING] Bandwidth selection failed: {e$message}"))
  message("  Using default bandwidth of 30 neighbours")
  30
})

message(glue("  Optimal bandwidth: {bw_adapt} neighbours"))

# Fit GWR model
message("  Fitting GWR model ...")
gwr_model <- tryCatch({
  gwr.basic(overall_star_rating_mean ~ irsd_score_mean,
            data = sa3_gwr_sp_pts,
            bw = bw_adapt,
            kernel = "bisquare",
            adaptive = TRUE)
}, error = function(e) {
  message(glue("  [ERROR] GWR failed: {e$message}"))
  NULL
})

if (!is.null(gwr_model)) {
  message("  GWR model fitted successfully")

  # Extract results
  gwr_results <- as.data.frame(gwr_model$SDF)
  gwr_results$sa3_code <- sa3_gwr$SA3_CODE21

  # Summary statistics of local coefficients
  message("\n  GWR local coefficient summary (IRSD effect):")
  irsd_coefs <- gwr_results$irsd_score_mean
  message(glue("    Mean:   {round(mean(irsd_coefs), 6)}"))
  message(glue("    Min:    {round(min(irsd_coefs), 6)}"))
  message(glue("    Max:    {round(max(irsd_coefs), 6)}"))
  message(glue("    SD:     {round(sd(irsd_coefs), 6)}"))
  message(glue("    % positive: {round(mean(irsd_coefs > 0) * 100, 1)}%"))
  message(glue("    % negative: {round(mean(irsd_coefs < 0) * 100, 1)}%"))

  # Compare GWR vs OLS
  ols_model <- lm(overall_star_rating_mean ~ irsd_score_mean, data = sa3_gwr)
  message(glue("\n  OLS R-squared: {round(summary(ols_model)$r.squared, 4)}"))
  message(glue("  GWR R-squared (mean local): {round(mean(gwr_results$Local_R2), 4)}"))
  message(glue("  GWR AICc: {round(gwr_model$GW.diagnostic$AICc, 2)}"))
  message(glue("  OLS AICc: {round(AIC(ols_model), 2)}"))

  # Save GWR coefficient summary
  gwr_summary <- tibble(
    statistic = c("Mean", "Median", "Min", "Max", "SD",
                  "% Positive", "% Negative",
                  "OLS R-squared", "GWR mean local R-squared",
                  "GWR AICc", "OLS AIC", "Bandwidth (neighbours)"),
    value = c(
      round(mean(irsd_coefs), 6),
      round(median(irsd_coefs), 6),
      round(min(irsd_coefs), 6),
      round(max(irsd_coefs), 6),
      round(sd(irsd_coefs), 6),
      round(mean(irsd_coefs > 0) * 100, 1),
      round(mean(irsd_coefs < 0) * 100, 1),
      round(summary(ols_model)$r.squared, 4),
      round(mean(gwr_results$Local_R2), 4),
      round(gwr_model$GW.diagnostic$AICc, 2),
      round(AIC(ols_model), 2),
      bw_adapt
    )
  )

  write.csv(gwr_summary, file.path(dir_tables, "gwr_summary.csv"), row.names = FALSE)
  message("  Saved: gwr_summary.csv")

  # Map of local IRSD coefficients
  gwr_sf <- st_sf(
    sa3_code = sa3_gwr$SA3_CODE21,
    irsd_coef = gwr_results$irsd_score_mean,
    local_r2 = gwr_results$Local_R2,
    geometry = st_geometry(sa3_gwr)
  )

  map_gwr <- tm_shape(shp, projection = 3577) +
    tm_fill(fill = "grey90", fill_alpha = 0.5) +
    tm_borders(col = "grey80", lwd = 0.2) +
    tm_shape(gwr_sf) +
    tm_fill(
      fill = "irsd_coef",
      fill.scale = tm_scale_continuous(
        values = "RdBu",
        midpoint = 0
      ),
      fill.legend = tm_legend(title = "Local IRSD\nCoefficient")
    ) +
    tm_borders(col = "grey60", lwd = 0.3) +
    tm_title("GWR Local Coefficients: IRSD Effect on Star Rating") +
    tm_layout(frame = FALSE)

  tmap_save(map_gwr,
            filename = file.path(dir_maps, "map_gwr_irsd_coefficients.png"),
            width = 10, height = 8, dpi = 300)
  message("  Saved: map_gwr_irsd_coefficients.png")

  # Map of local R-squared
  map_r2 <- tm_shape(shp, projection = 3577) +
    tm_fill(fill = "grey90", fill_alpha = 0.5) +
    tm_borders(col = "grey80", lwd = 0.2) +
    tm_shape(gwr_sf) +
    tm_fill(
      fill = "local_r2",
      fill.scale = tm_scale_continuous(values = "YlOrRd"),
      fill.legend = tm_legend(title = "Local R-squared")
    ) +
    tm_borders(col = "grey60", lwd = 0.3) +
    tm_title("GWR Local R-squared: Model Fit Across Australia") +
    tm_layout(frame = FALSE)

  tmap_save(map_r2,
            filename = file.path(dir_maps, "map_gwr_local_r2.png"),
            width = 10, height = 8, dpi = 300)
  message("  Saved: map_gwr_local_r2.png")
}


# =============================================================================
# PART 7: Facility density and provider type maps
# =============================================================================

message("\n=== Part 7: Supplementary maps ===\n")

# --- 7a: Facility density (facilities per 1000 aged population) ---
shp_density <- shp_with_data %>%
  filter(!is.na(pop_65plus), pop_65plus > 0) %>%
  mutate(facilities_per_1000_65plus = n_facilities / pop_65plus * 1000)

map_density <- tm_shape(shp, projection = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) +
  tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_density) +
  tm_fill(
    fill = "facilities_per_1000_65plus",
    fill.scale = tm_scale_continuous(values = "YlGnBu"),
    fill.legend = tm_legend(title = "Facilities per\n1,000 pop 65+")
  ) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title(glue("Aged Care Facility Density — {latest_label}")) +
  tm_layout(frame = FALSE)

tmap_save(map_density,
          filename = file.path(dir_maps, "map_facility_density.png"),
          width = 10, height = 8, dpi = 300)
message("  Saved: map_facility_density.png")

# --- 7b: Proportion for-profit map ---
shp_fp <- shp_with_data %>%
  mutate(pct_for_profit = n_for_profit / n_facilities * 100)

map_fp <- tm_shape(shp, projection = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) +
  tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_fp) +
  tm_fill(
    fill = "pct_for_profit",
    fill.scale = tm_scale_continuous(values = "OrRd"),
    fill.legend = tm_legend(title = "% For-profit\nfacilities")
  ) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title(glue("For-profit Facility Share by SA3 — {latest_label}")) +
  tm_layout(frame = FALSE)

tmap_save(map_fp,
          filename = file.path(dir_maps, "map_pct_for_profit.png"),
          width = 10, height = 8, dpi = 300)
message("  Saved: map_pct_for_profit.png")


# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Phase 2 Complete ===\n")
message("  Maps saved to: outputs/maps/")
message("  Figures saved to: outputs/figures/")
message("  Tables saved to: outputs/tables/")

message("\n  Key results:")
message(glue("    Moran's I (overall): {round(moran_results$morans_i[1], 4)} (p = {format.pval(moran_results$p_value[1], digits = 3)})"))
message(glue("    IRSD-rating correlation: r = {round(cor_test$estimate, 4)} (p = {format.pval(cor_test$p.value, digits = 3)})"))
if (!is.null(gwr_model)) {
  message(glue("    GWR shows {'spatially varying' } IRSD effect (range: {round(min(irsd_coefs), 5)} to {round(max(irsd_coefs), 5)})"))
}
message("\nDone.")
