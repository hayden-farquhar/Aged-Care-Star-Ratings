# 04b_supplement.R — Run ordinal model and fix desert threshold sensitivity
# (ordinal package was installed during main run; desert dist had NA issue)

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(lmerTest)
library(glue)
library(sf)
library(spdep)
library(purrr)

project_root <- here::here()
dir_proc <- file.path(project_root, "data", "processed")
dir_tab  <- file.path(project_root, "outputs", "tables")

panel   <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
sa3_ref <- readRDS(file.path(dir_proc, "sa3_reference.rds"))

dat <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  mutate(
    quarter_c = quarter_num - median(unique(quarter_num)),
    purpose_f = relevel(factor(purpose), ref = "Not-for-profit")
  )

# =============================================================================
# D. Ordinal Sensitivity Check (re-run with package installed)
# =============================================================================

message("\n=== D. Ordinal sensitivity check ===\n")

library(ordinal)

dat_ord <- dat %>%
  mutate(rating_ordered = factor(overall_star_rating, ordered = TRUE))

message("  Fitting cumulative link mixed model...")
m_clmm <- tryCatch({
  clmm(rating_ordered ~ quarter_c * purpose_f + (1 | facility_id),
       data = dat_ord, link = "logit",
       control = clmm.control(maxIter = 200, maxLineIter = 100))
}, error = function(e) {
  message(glue("  Error with interaction: {e$message}"))
  message("  Trying without interaction...")
  tryCatch(
    clmm(rating_ordered ~ quarter_c + purpose_f + (1 | facility_id),
         data = dat_ord, link = "logit",
         control = clmm.control(maxIter = 200, maxLineIter = 100)),
    error = function(e2) { message(glue("  Also failed: {e2$message}")); NULL }
  )
})

if (!is.null(m_clmm)) {
  clmm_coefs <- summary(m_clmm)$coefficients
  message("  CLMM coefficients:")
  for (nm in rownames(clmm_coefs)) {
    est <- round(clmm_coefs[nm, "Estimate"], 3)
    pval <- if ("Pr(>|z|)" %in% colnames(clmm_coefs)) {
      format(clmm_coefs[nm, "Pr(>|z|)"], digits = 3)
    } else "NA"
    message(glue("    {nm}: B = {est}, p = {pval}"))
  }

  clmm_df <- as.data.frame(clmm_coefs) %>%
    tibble::rownames_to_column("term")
  write_csv(clmm_df, file.path(dir_tab, "val_ordinal_model_coefficients.csv"))
  message("  Saved: val_ordinal_model_coefficients.csv")

  # Compare with linear model
  message("\n  Direction comparison (ordinal vs linear):")
  m3 <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 | facility_id), data = dat)
  m3_coefs <- broom.mixed::tidy(m3, effects = "fixed")

  compare_terms <- c("quarter_c", "purpose_fFor-profit", "purpose_fGovernment")
  for (t in compare_terms) {
    lin_est <- m3_coefs %>% filter(term == t) %>% pull(estimate)
    ord_est <- if (t %in% rownames(clmm_coefs)) clmm_coefs[t, "Estimate"] else NA
    if (!is.na(ord_est)) {
      same_dir <- sign(lin_est) == sign(ord_est)
      message(glue("    {t}: linear = {round(lin_est, 3)}, ordinal = {round(ord_est, 3)}, same direction = {same_dir}"))
    }
  }
} else {
  message("  Ordinal model could not be fitted.")
}

# =============================================================================
# L. Quality Desert Threshold Sensitivity (fixed)
# =============================================================================

message("\n=== L. Quality desert threshold sensitivity (fixed) ===\n")

shp <- st_read(file.path(project_root, "data", "spatial", "SA3_2021_AUST_SHP_GDA2020"), quiet = TRUE)

latest <- panel %>%
  filter(quarter_num == max(quarter_num), !is.na(sa3_code), !is.na(overall_star_rating))

# All SA3 centroids
shp_proj <- st_transform(shp, crs = 3577)
sa3_centroids <- suppressWarnings(st_centroid(shp_proj))
all_coords <- st_coordinates(sa3_centroids)

# SA3s with at least one 3+ star facility
adequate_sa3_codes <- latest %>%
  filter(overall_star_rating >= 3) %>%
  distinct(sa3_code) %>%
  pull(sa3_code)

adequate_idx <- which(sa3_centroids$SA3_CODE21 %in% adequate_sa3_codes)
adequate_coords <- all_coords[adequate_idx, , drop = FALSE]

message(glue("  Total SA3 areas: {nrow(all_coords)}"))
message(glue("  SA3 areas with 3+ star facility: {nrow(adequate_coords)}"))

# Compute distances from every SA3 to nearest adequate SA3
dist_to_adequate <- rep(NA_real_, nrow(all_coords))
for (i in seq_len(nrow(all_coords))) {
  dists <- sqrt((all_coords[i, 1] - adequate_coords[, 1])^2 +
                (all_coords[i, 2] - adequate_coords[, 2])^2)
  dist_to_adequate[i] <- min(dists) / 1000  # metres to km
}

sa3_centroids$dist_adequate_km <- dist_to_adequate

# Now test different thresholds
thresholds <- c(30, 50, 75, 100, 150)
sensitivity <- tibble(
  threshold_km = thresholds,
  n_deserts = sapply(thresholds, function(t) sum(dist_to_adequate > t, na.rm = TRUE)),
  pct_sa3 = sapply(thresholds, function(t) round(sum(dist_to_adequate > t, na.rm = TRUE) / sum(!is.na(dist_to_adequate)) * 100, 1)),
  pop_65plus_affected = sapply(thresholds, function(t) {
    desert_codes <- sa3_centroids$SA3_CODE21[dist_to_adequate > t]
    sum(sa3_ref$pop_65plus[sa3_ref$sa3_code %in% desert_codes], na.rm = TRUE)
  })
)

message("  Quality desert sensitivity to distance threshold:")
for (i in seq_len(nrow(sensitivity))) {
  message(glue("    >{sensitivity$threshold_km[i]}km: {sensitivity$n_deserts[i]} SA3s ({sensitivity$pct_sa3[i]}%), 65+ pop = {format(sensitivity$pop_65plus_affected[i], big.mark = ',')}"))
}

write_csv(sensitivity, file.path(dir_tab, "val_desert_threshold_sensitivity.csv"))
message("  Saved: val_desert_threshold_sensitivity.csv")

message("\nDone.")
