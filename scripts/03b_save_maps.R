# 03b_save_maps.R — Re-save all maps using ragg device (cairo unavailable)

library(sf)
library(dplyr)
library(tmap)
library(glue)
library(GWmodel)
library(sp)

project_root <- here::here()
dir_proc <- file.path(project_root, "data", "processed")
dir_maps <- file.path(project_root, "outputs", "maps")

panel   <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
sa3_ref <- readRDS(file.path(dir_proc, "sa3_reference.rds"))

latest_q     <- max(panel$quarter_num)
latest_label <- panel %>% filter(quarter_num == latest_q) %>% pull(quarter_label) %>% unique()
latest       <- panel %>% filter(quarter_num == latest_q, !is.na(sa3_code))

rating_vars <- c("overall_star_rating", "residents_experience_rating",
                 "compliance_rating", "staffing_rating", "quality_measures_rating")

sa3_ratings <- latest %>%
  filter(!is.na(overall_star_rating)) %>%
  group_by(sa3_code, sa3_name) %>%
  summarise(
    n_facilities = n(),
    across(all_of(rating_vars), list(mean = ~mean(.x, na.rm = TRUE)), .names = "{.col}_mean"),
    n_for_profit = sum(purpose == "For-profit", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(sa3_ref, by = c("sa3_code", "sa3_name"))

shp <- st_read(file.path(project_root, "data", "spatial", "SA3_2021_AUST_SHP_GDA2020"), quiet = TRUE)
shp_with_data <- shp %>%
  left_join(sa3_ratings, by = c("SA3_CODE21" = "sa3_code")) %>%
  filter(!is.na(n_facilities))

tmap_mode("plot")

save_map <- function(map, filename, ...) {
  tmap_save(map, filename = filename, device = ragg::agg_png, ...)
}

# 1. Overall Star Rating
message("Saving overall star rating map...")
m1 <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_with_data, crs = 3577) +
  tm_fill(fill = "overall_star_rating_mean",
          fill.scale = tm_scale_continuous(values = "viridis", midpoint = 3.5),
          fill.legend = tm_legend(title = "Mean Star Rating")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title(paste("Overall Star Rating by SA3 --", latest_label)) +
  tm_layout(frame = FALSE)
save_map(m1, file.path(dir_maps, "map_overall_star_rating.png"), width = 10, height = 8, dpi = 300)

# 2. Sub-category maps
subcats <- list(
  list(var = "residents_experience_rating_mean", label = "Residents Experience", fname = "residents-experience"),
  list(var = "compliance_rating_mean",           label = "Compliance",           fname = "compliance"),
  list(var = "staffing_rating_mean",             label = "Staffing",             fname = "staffing"),
  list(var = "quality_measures_rating_mean",     label = "Quality Measures",     fname = "quality-measures")
)
for (sc in subcats) {
  message(glue("Saving {sc$label} map..."))
  m <- tm_shape(shp, crs = 3577) +
    tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
    tm_shape(shp_with_data, crs = 3577) +
    tm_fill(fill = sc$var,
            fill.scale = tm_scale_continuous(values = "viridis", midpoint = 3.5),
            fill.legend = tm_legend(title = paste("Mean", sc$label))) +
    tm_borders(col = "grey60", lwd = 0.3) +
    tm_title(paste(sc$label, "Rating by SA3 --", latest_label)) +
    tm_layout(frame = FALSE)
  save_map(m, file.path(dir_maps, paste0("map_", sc$fname, "-rating.png")), width = 10, height = 8, dpi = 300)
}

# 3. IRSD quintile
message("Saving IRSD quintile map...")
m_irsd <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_with_data, crs = 3577) +
  tm_fill(fill = "irsd_quintile",
          fill.scale = tm_scale_continuous(values = "brewer.rd_yl_gn", midpoint = 3),
          fill.legend = tm_legend(title = "IRSD Quintile\n(1=Most Disadvantaged)")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title("SEIFA IRSD Quintile by SA3") +
  tm_layout(frame = FALSE)
save_map(m_irsd, file.path(dir_maps, "map_irsd_quintile.png"), width = 10, height = 8, dpi = 300)

# 4. Facility density
message("Saving facility density map...")
shp_density <- shp_with_data %>%
  filter(!is.na(pop_65plus), pop_65plus > 0) %>%
  mutate(facilities_per_1000_65plus = n_facilities / pop_65plus * 1000)
m_den <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_density) +
  tm_fill(fill = "facilities_per_1000_65plus",
          fill.scale = tm_scale_continuous(values = "brewer.yl_gn_bu"),
          fill.legend = tm_legend(title = "Facilities per\n1,000 pop 65+")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title(paste("Aged Care Facility Density --", latest_label)) +
  tm_layout(frame = FALSE)
save_map(m_den, file.path(dir_maps, "map_facility_density.png"), width = 10, height = 8, dpi = 300)

# 5. For-profit share
message("Saving for-profit share map...")
shp_fp <- shp_with_data %>% mutate(pct_for_profit = n_for_profit / n_facilities * 100)
m_fp <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_fp) +
  tm_fill(fill = "pct_for_profit",
          fill.scale = tm_scale_continuous(values = "brewer.or_rd"),
          fill.legend = tm_legend(title = "% For-profit")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title(paste("For-profit Facility Share by SA3 --", latest_label)) +
  tm_layout(frame = FALSE)
save_map(m_fp, file.path(dir_maps, "map_pct_for_profit.png"), width = 10, height = 8, dpi = 300)

# 6. GWR maps
message("Saving GWR maps...")
sa3_gwr_data <- shp_with_data %>% filter(n_facilities >= 3, !is.na(irsd_score_mean))
sa3_gwr_proj <- st_transform(sa3_gwr_data, crs = 3577)
centroids_sf <- st_centroid(sa3_gwr_proj)
coords <- st_coordinates(centroids_sf)
sa3_sp <- SpatialPointsDataFrame(
  coords = coords,
  data = as.data.frame(sa3_gwr_proj) %>% select(-geometry),
  proj4string = CRS("+init=epsg:3577")
)

bw <- bw.gwr(overall_star_rating_mean ~ irsd_score_mean,
             data = sa3_sp, approach = "AICc", kernel = "bisquare", adaptive = TRUE)
gwr_fit <- gwr.basic(overall_star_rating_mean ~ irsd_score_mean,
                     data = sa3_sp, bw = bw, kernel = "bisquare", adaptive = TRUE)
gwr_sf <- st_sf(
  irsd_coef = as.data.frame(gwr_fit$SDF)$irsd_score_mean,
  local_r2  = as.data.frame(gwr_fit$SDF)$Local_R2,
  geometry  = st_geometry(sa3_gwr_proj)
)

m_gwr <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(gwr_sf) +
  tm_fill(fill = "irsd_coef",
          fill.scale = tm_scale_continuous(values = "brewer.rd_bu", midpoint = 0),
          fill.legend = tm_legend(title = "Local IRSD\nCoefficient")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title("GWR: Local IRSD Effect on Star Rating") +
  tm_layout(frame = FALSE)
save_map(m_gwr, file.path(dir_maps, "map_gwr_irsd_coefficients.png"), width = 10, height = 8, dpi = 300)

m_r2 <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey90", fill_alpha = 0.5) + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(gwr_sf) +
  tm_fill(fill = "local_r2",
          fill.scale = tm_scale_continuous(values = "brewer.yl_or_rd"),
          fill.legend = tm_legend(title = "Local R-squared")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title("GWR: Local R-squared Across Australia") +
  tm_layout(frame = FALSE)
save_map(m_r2, file.path(dir_maps, "map_gwr_local_r2.png"), width = 10, height = 8, dpi = 300)

message("\nAll maps saved.")
