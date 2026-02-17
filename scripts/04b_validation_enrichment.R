# =============================================================================
# 04b_validation_enrichment.R
# Statistical validation, sensitivity analyses, and enrichment
#
# Sections:
#   A. Mixed model diagnostics and residual checks
#   B. Random slopes model
#   C. Extended mixed models (additional covariates: size, state, remoteness, IRSD)
#   D. Ordinal sensitivity check (cumulative link mixed model)
#   E. Marginal and conditional R-squared
#   F. Non-linear time trend (quadratic)
#   G. Balanced panel sensitivity analysis
#   H. LISA cluster analysis (local Moran's I)
#   I. Spatial weights sensitivity (Moran's I with varying k)
#   J. Facility-count weighted IRSD regression
#   K. ANOVA post-hoc tests
#   L. Quality desert sensitivity (alternative thresholds)
#   M. Quality desert regression (logistic and linear)
#   N. Transition probability modelling
#   O. Multiple testing corrections
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(ggplot2)
library(readr)
library(lmerTest)
library(broom.mixed)
library(sf)
library(spdep)
library(glue)
library(purrr)

project_root <- here::here()
dir_proc <- file.path(project_root, "data", "processed")
dir_fig  <- file.path(project_root, "outputs", "figures")
dir_tab  <- file.path(project_root, "outputs", "tables")
dir_maps <- file.path(project_root, "outputs", "maps")

panel   <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
sa3_ref <- readRDS(file.path(dir_proc, "sa3_reference.rds"))

theme_pub <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

# Prepare modelling data
dat <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  mutate(
    quarter_c = quarter_num - median(unique(quarter_num)),
    purpose_f = relevel(factor(purpose), ref = "Not-for-profit"),
    state_f   = factor(state),
    remoteness_f = factor(remoteness),
    size_f    = factor(size),
    irsd_q_f  = factor(irsd_quintile)
  )

message(glue("Modelling dataset: {nrow(dat)} records, {n_distinct(dat$facility_id)} facilities"))

# =============================================================================
# A. Mixed Model Diagnostics
# =============================================================================

message("\n=== A. Mixed model diagnostics ===\n")

m3 <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 | facility_id),
           data = dat)

# Extract residuals and fitted values
dat_diag <- dat %>%
  mutate(
    fitted_vals  = fitted(m3),
    residuals    = resid(m3),
    std_resid    = resid(m3, scaled = TRUE)
  )

# A1: Residuals vs fitted
p_rvf <- ggplot(dat_diag, aes(x = fitted_vals, y = residuals)) +
  geom_point(alpha = 0.05, size = 0.5) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, colour = "blue", linewidth = 0.8) +
  labs(title = "Residuals vs Fitted Values",
       subtitle = "Overall Star Rating mixed model (M3)",
       x = "Fitted Values", y = "Residuals") +
  theme_pub
ggsave(file.path(dir_fig, "diag_residuals_vs_fitted.png"), p_rvf,
       width = 8, height = 6, dpi = 300, device = ragg::agg_png)

# A2: QQ plot of residuals
p_qq <- ggplot(dat_diag, aes(sample = std_resid)) +
  stat_qq(alpha = 0.1, size = 0.5) +
  stat_qq_line(colour = "red") +
  labs(title = "QQ Plot of Standardised Residuals",
       subtitle = "Overall Star Rating mixed model (M3)",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_pub
ggsave(file.path(dir_fig, "diag_qq_residuals.png"), p_qq,
       width = 7, height = 6, dpi = 300, device = ragg::agg_png)

# A3: QQ plot of random effects
re <- ranef(m3)$facility_id[, 1]
p_re_qq <- ggplot(data.frame(re = re), aes(sample = re)) +
  stat_qq(alpha = 0.3, size = 0.5) +
  stat_qq_line(colour = "red") +
  labs(title = "QQ Plot of Facility Random Intercepts",
       subtitle = glue("{length(re)} facility-level intercepts"),
       x = "Theoretical Quantiles", y = "Random Intercepts") +
  theme_pub
ggsave(file.path(dir_fig, "diag_qq_random_effects.png"), p_re_qq,
       width = 7, height = 6, dpi = 300, device = ragg::agg_png)

# A4: Residuals over time
p_time <- dat_diag %>%
  group_by(quarter_num) %>%
  summarise(mean_resid = mean(residuals), sd_resid = sd(residuals), .groups = "drop") %>%
  ggplot(aes(x = quarter_num, y = mean_resid)) +
  geom_point(size = 3) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_resid - sd_resid, ymax = mean_resid + sd_resid), alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  labs(title = "Mean Residual by Quarter",
       subtitle = "Checking for systematic time patterns in residuals",
       x = "Quarter", y = "Mean Residual (± 1 SD)") +
  theme_pub
ggsave(file.path(dir_fig, "diag_residuals_over_time.png"), p_time,
       width = 8, height = 5, dpi = 300, device = ragg::agg_png)

message("  Saved 4 diagnostic plots")

# Shapiro-Wilk on random effects (sample if >5000)
sw_re <- shapiro.test(if (length(re) > 5000) sample(re, 5000) else re)
message(glue("  Shapiro-Wilk on random effects: W = {round(sw_re$statistic, 4)}, p = {format(sw_re$p.value, digits = 3)}"))

# Shapiro-Wilk on residuals (sample)
resid_sample <- sample(resid(m3), min(5000, length(resid(m3))))
sw_resid <- shapiro.test(resid_sample)
message(glue("  Shapiro-Wilk on residuals (n=5000 sample): W = {round(sw_resid$statistic, 4)}, p = {format(sw_resid$p.value, digits = 3)}"))

# =============================================================================
# B. Random Slopes Model
# =============================================================================

message("\n=== B. Random slopes model ===\n")

message("  Fitting random slopes model (may take a minute)...")
m_rs <- tryCatch({
  lmer(overall_star_rating ~ quarter_c * purpose_f + (1 + quarter_c | facility_id),
       data = dat, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
}, warning = function(w) {
  message(glue("  Warning: {w$message}"))
  suppressWarnings(
    lmer(overall_star_rating ~ quarter_c * purpose_f + (1 + quarter_c | facility_id),
         data = dat, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))
  )
}, error = function(e) {
  message(glue("  Error fitting random slopes: {e$message}"))
  NULL
})

if (!is.null(m_rs)) {
  # Compare random intercept vs random slopes
  anova_ri_rs <- anova(m3, m_rs)
  message(glue("  Random intercept AIC: {round(AIC(m3), 1)}"))
  message(glue("  Random slopes AIC: {round(AIC(m_rs), 1)}"))
  message(glue("  LRT chi2: {round(anova_ri_rs$Chisq[2], 1)}, p = {format(anova_ri_rs$`Pr(>Chisq)`[2], digits = 3)}"))

  # Variance components
  vc <- as.data.frame(VarCorr(m_rs))
  message("  Variance components (random slopes model):")
  for (i in seq_len(nrow(vc))) {
    message(glue("    {vc$grp[i]} / {vc$var1[i]}{ifelse(is.na(vc$var2[i]), '', paste0(' x ', vc$var2[i]))}: {round(vc$vcov[i], 5)}"))
  }

  # Save coefficients
  coefs_rs <- tidy(m_rs, effects = "fixed", conf.int = TRUE)
  write_csv(coefs_rs, file.path(dir_tab, "val_random_slopes_coefficients.csv"))
  message("  Saved: val_random_slopes_coefficients.csv")

  # Histogram of facility-level slopes
  re_rs <- ranef(m_rs)$facility_id
  p_slopes <- ggplot(re_rs, aes(x = quarter_c)) +
    geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
    geom_vline(xintercept = 0, colour = "red", linetype = "dashed") +
    labs(title = "Distribution of Facility-Level Improvement Rates",
         subtitle = glue("Random slopes for time (n = {nrow(re_rs)} facilities)"),
         x = "Facility-Level Slope (stars per quarter)", y = "Count") +
    theme_pub
  ggsave(file.path(dir_fig, "val_random_slopes_distribution.png"), p_slopes,
         width = 8, height = 5, dpi = 300, device = ragg::agg_png)
  message("  Saved: val_random_slopes_distribution.png")
}

# =============================================================================
# C. Extended Mixed Models (Additional Covariates)
# =============================================================================

message("\n=== C. Extended mixed models ===\n")

# C1: Model with size
dat_size <- dat %>% filter(!is.na(size_f))
m_size <- lmer(overall_star_rating ~ quarter_c * purpose_f + size_f + (1 | facility_id),
               data = dat_size)
message("  Model with size:")
coefs_size <- tidy(m_size, effects = "fixed", conf.int = TRUE)
coefs_size %>%
  filter(grepl("size", term)) %>%
  { walk(seq_len(nrow(.)), function(i) {
    message(glue("    {.$term[i]}: B = {round(.$estimate[i], 3)}, p = {format(.$p.value[i], digits = 3)}"))
  })}

# C2: Model with state
m_state <- lmer(overall_star_rating ~ quarter_c * purpose_f + state_f + (1 | facility_id),
                data = dat)
message("  Model with state:")
coefs_state <- tidy(m_state, effects = "fixed", conf.int = TRUE)
coefs_state %>%
  filter(grepl("state", term)) %>%
  { walk(seq_len(nrow(.)), function(i) {
    message(glue("    {.$term[i]}: B = {round(.$estimate[i], 3)}, p = {format(.$p.value[i], digits = 3)}"))
  })}

# C3: Model with remoteness
dat_remote <- dat %>% filter(!is.na(remoteness_f))
m_remote <- lmer(overall_star_rating ~ quarter_c * purpose_f + remoteness_f + (1 | facility_id),
                 data = dat_remote)
message("  Model with remoteness:")
coefs_remote <- tidy(m_remote, effects = "fixed", conf.int = TRUE)
coefs_remote %>%
  filter(grepl("remote", term, ignore.case = TRUE)) %>%
  { walk(seq_len(nrow(.)), function(i) {
    message(glue("    {.$term[i]}: B = {round(.$estimate[i], 3)}, p = {format(.$p.value[i], digits = 3)}"))
  })}

# C4: Model with IRSD quintile
dat_irsd <- dat %>% filter(!is.na(irsd_q_f))
m_irsd <- lmer(overall_star_rating ~ quarter_c * purpose_f + irsd_q_f + (1 | facility_id),
               data = dat_irsd)
message("  Model with IRSD quintile:")
coefs_irsd <- tidy(m_irsd, effects = "fixed", conf.int = TRUE)
coefs_irsd %>%
  filter(grepl("irsd", term, ignore.case = TRUE)) %>%
  { walk(seq_len(nrow(.)), function(i) {
    message(glue("    {.$term[i]}: B = {round(.$estimate[i], 3)}, p = {format(.$p.value[i], digits = 3)}"))
  })}

# C5: Full model — all covariates
dat_full <- dat %>% filter(!is.na(size_f), !is.na(remoteness_f), !is.na(irsd_q_f))
m_full <- lmer(overall_star_rating ~ quarter_c * purpose_f + size_f + state_f +
                 remoteness_f + irsd_q_f + (1 | facility_id),
               data = dat_full)
message("\n  Full model (all covariates):")
coefs_full <- tidy(m_full, effects = "fixed", conf.int = TRUE)
write_csv(coefs_full, file.path(dir_tab, "val_full_model_coefficients.csv"))
message("  Saved: val_full_model_coefficients.csv")

# Compare provider type effect: base vs full
base_fp <- tidy(m3, effects = "fixed") %>% filter(term == "purpose_fFor-profit") %>% pull(estimate)
full_fp <- coefs_full %>% filter(term == "purpose_fFor-profit") %>% pull(estimate)
base_gov <- tidy(m3, effects = "fixed") %>% filter(term == "purpose_fGovernment") %>% pull(estimate)
full_gov <- coefs_full %>% filter(term == "purpose_fGovernment") %>% pull(estimate)

message(glue("  For-profit effect: base = {round(base_fp, 3)}, adjusted = {round(full_fp, 3)} (change = {round((full_fp - base_fp) / abs(base_fp) * 100, 1)}%)"))
message(glue("  Government effect: base = {round(base_gov, 3)}, adjusted = {round(full_gov, 3)} (change = {round((full_gov - base_gov) / abs(base_gov) * 100, 1)}%)"))

# =============================================================================
# D. Ordinal Sensitivity Check
# =============================================================================

message("\n=== D. Ordinal sensitivity check ===\n")

# Check if ordinal package is available
if (requireNamespace("ordinal", quietly = TRUE)) {
  library(ordinal)

  dat_ord <- dat %>%
    mutate(rating_ordered = factor(overall_star_rating, ordered = TRUE))

  message("  Fitting cumulative link mixed model (this may take a while)...")
  m_clmm <- tryCatch({
    clmm(rating_ordered ~ quarter_c * purpose_f + (1 | facility_id),
         data = dat_ord, link = "logit")
  }, error = function(e) {
    message(glue("  Error: {e$message}"))
    # Try simpler model
    message("  Trying simpler model without interaction...")
    tryCatch(
      clmm(rating_ordered ~ quarter_c + purpose_f + (1 | facility_id),
           data = dat_ord, link = "logit"),
      error = function(e2) { message(glue("  Also failed: {e2$message}")); NULL }
    )
  })

  if (!is.null(m_clmm)) {
    clmm_coefs <- summary(m_clmm)$coefficients
    message("  CLMM coefficients:")
    for (nm in rownames(clmm_coefs)) {
      est <- round(clmm_coefs[nm, "Estimate"], 3)
      pval <- if ("Pr(>|z|)" %in% colnames(clmm_coefs)) format(clmm_coefs[nm, "Pr(>|z|)"], digits = 3) else "NA"
      message(glue("    {nm}: B = {est}, p = {pval}"))
    }
    # Save
    clmm_df <- as.data.frame(clmm_coefs) %>%
      tibble::rownames_to_column("term")
    write_csv(clmm_df, file.path(dir_tab, "val_ordinal_model_coefficients.csv"))
    message("  Saved: val_ordinal_model_coefficients.csv")

    # Key comparison: direction and significance should match linear model
    message("\n  Ordinal vs Linear comparison (direction and significance):")
    message("    If signs and significance match, linear approximation is adequate.")
  }
} else {
  message("  ordinal package not available — installing...")
  install.packages("ordinal", repos = "https://cloud.r-project.org", quiet = TRUE)
  message("  Installed. Please re-run this section.")
}

# =============================================================================
# E. Marginal and Conditional R-squared
# =============================================================================

message("\n=== E. R-squared (marginal and conditional) ===\n")

if (requireNamespace("performance", quietly = TRUE)) {
  library(performance)
  r2_m3 <- r2_nakagawa(m3)
  message(glue("  Base model (M3): R2_marginal = {round(r2_m3$R2_marginal, 4)}, R2_conditional = {round(r2_m3$R2_conditional, 4)}"))

  if (!is.null(m_rs)) {
    r2_rs <- r2_nakagawa(m_rs)
    message(glue("  Random slopes: R2_marginal = {round(r2_rs$R2_marginal, 4)}, R2_conditional = {round(r2_rs$R2_conditional, 4)}"))
  }

  r2_full <- r2_nakagawa(m_full)
  message(glue("  Full model:    R2_marginal = {round(r2_full$R2_marginal, 4)}, R2_conditional = {round(r2_full$R2_conditional, 4)}"))

  r2_tab <- tibble(
    model = c("M3 (base)", "Random slopes", "Full (all covariates)"),
    R2_marginal = c(r2_m3$R2_marginal,
                    if (!is.null(m_rs)) r2_rs$R2_marginal else NA,
                    r2_full$R2_marginal),
    R2_conditional = c(r2_m3$R2_conditional,
                       if (!is.null(m_rs)) r2_rs$R2_conditional else NA,
                       r2_full$R2_conditional)
  )
  write_csv(r2_tab, file.path(dir_tab, "val_r_squared.csv"))
  message("  Saved: val_r_squared.csv")
} else {
  message("  Installing performance package...")
  install.packages("performance", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(performance)
  r2_m3 <- r2_nakagawa(m3)
  message(glue("  Base model (M3): R2_marginal = {round(r2_m3$R2_marginal, 4)}, R2_conditional = {round(r2_m3$R2_conditional, 4)}"))
}

# =============================================================================
# F. Non-linear Time Trend
# =============================================================================

message("\n=== F. Non-linear time trend ===\n")

dat <- dat %>% mutate(quarter_c2 = quarter_c^2)

m_quad <- lmer(overall_star_rating ~ quarter_c + quarter_c2 + purpose_f +
                 quarter_c:purpose_f + quarter_c2:purpose_f + (1 | facility_id),
               data = dat)

anova_lin_quad <- anova(m3, m_quad)
message(glue("  Linear model AIC: {round(AIC(m3), 1)}"))
message(glue("  Quadratic model AIC: {round(AIC(m_quad), 1)}"))
message(glue("  LRT chi2: {round(anova_lin_quad$Chisq[2], 1)}, p = {format(anova_lin_quad$`Pr(>Chisq)`[2], digits = 3)}"))

quad_coef <- tidy(m_quad, effects = "fixed", conf.int = TRUE)
quad_term <- quad_coef %>% filter(term == "quarter_c2")
message(glue("  Quadratic term: B = {round(quad_term$estimate, 4)}, p = {format(quad_term$p.value, digits = 3)}"))
if (quad_term$estimate < 0) {
  message("  Interpretation: Negative quadratic = decelerating improvement (plateau)")
} else {
  message("  Interpretation: Positive quadratic = accelerating improvement")
}

write_csv(quad_coef, file.path(dir_tab, "val_quadratic_model_coefficients.csv"))
message("  Saved: val_quadratic_model_coefficients.csv")

# =============================================================================
# G. Balanced Panel Sensitivity
# =============================================================================

message("\n=== G. Balanced panel sensitivity ===\n")

# Restrict to facilities present in all 11 quarters
balanced_ids <- dat %>%
  group_by(facility_id) %>%
  summarise(n_q = n_distinct(quarter_num), .groups = "drop") %>%
  filter(n_q == 11) %>%
  pull(facility_id)

dat_balanced <- dat %>% filter(facility_id %in% balanced_ids)
n_balanced <- n_distinct(dat_balanced$facility_id)
message(glue("  Balanced panel: {n_balanced} facilities (of {n_distinct(dat$facility_id)})"))

m_balanced <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 | facility_id),
                   data = dat_balanced)

coefs_balanced <- tidy(m_balanced, effects = "fixed", conf.int = TRUE)
write_csv(coefs_balanced, file.path(dir_tab, "val_balanced_panel_coefficients.csv"))
message("  Saved: val_balanced_panel_coefficients.csv")

# Compare key coefficients
coefs_base <- tidy(m3, effects = "fixed")
for (t in c("quarter_c", "purpose_fFor-profit", "purpose_fGovernment")) {
  b_base <- coefs_base %>% filter(term == t) %>% pull(estimate)
  b_bal  <- coefs_balanced %>% filter(term == t) %>% pull(estimate)
  message(glue("  {t}: base = {round(b_base, 4)}, balanced = {round(b_bal, 4)}"))
}

# =============================================================================
# H. LISA Cluster Analysis
# =============================================================================

message("\n=== H. LISA cluster analysis ===\n")

shp <- st_read(file.path(project_root, "data", "spatial", "SA3_2021_AUST_SHP_GDA2020"), quiet = TRUE)

latest_q <- max(panel$quarter_num)
latest <- panel %>% filter(quarter_num == latest_q, !is.na(sa3_code), !is.na(overall_star_rating))

sa3_ratings <- latest %>%
  group_by(sa3_code, sa3_name) %>%
  summarise(n_facilities = n(),
            mean_rating = mean(overall_star_rating, na.rm = TRUE),
            .groups = "drop") %>%
  filter(n_facilities >= 3)

shp_data <- shp %>%
  left_join(sa3_ratings, by = c("SA3_CODE21" = "sa3_code")) %>%
  filter(!is.na(n_facilities))

shp_proj <- st_transform(shp_data, crs = 3577)
coords <- st_coordinates(st_centroid(suppressWarnings(shp_proj)))
knn <- knearneigh(coords, k = 5)
nb <- knn2nb(knn)
lw <- nb2listw(nb, style = "W")

# Local Moran's I
local_moran <- localmoran(shp_proj$mean_rating, lw)
shp_proj$li <- local_moran[, "Ii"]
shp_proj$li_p <- local_moran[, "Pr(z != E(Ii))"]

# LISA cluster classification
z_rating <- scale(shp_proj$mean_rating)[, 1]
lag_rating <- lag.listw(lw, z_rating)

shp_proj$lisa_cluster <- case_when(
  shp_proj$li_p > 0.05 ~ "Not significant",
  z_rating > 0 & lag_rating > 0 ~ "High-High",
  z_rating < 0 & lag_rating < 0 ~ "Low-Low",
  z_rating > 0 & lag_rating < 0 ~ "High-Low",
  z_rating < 0 & lag_rating > 0 ~ "Low-High"
)

lisa_counts <- table(shp_proj$lisa_cluster)
message("  LISA cluster classification:")
for (cl in names(lisa_counts)) {
  message(glue("    {cl}: {lisa_counts[cl]} SA3 areas"))
}

# Save LISA results
lisa_df <- shp_proj %>%
  st_drop_geometry() %>%
  select(SA3_CODE21, SA3_NAME21, STE_NAME21, mean_rating, n_facilities,
         li, li_p, lisa_cluster) %>%
  arrange(li_p)
write_csv(lisa_df, file.path(dir_tab, "val_lisa_clusters.csv"))
message("  Saved: val_lisa_clusters.csv")

# LISA map
library(tmap)
tmap_mode("plot")

shp_proj$lisa_cluster <- factor(shp_proj$lisa_cluster,
  levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant"))

lisa_colours <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
                  "High-Low" = "#fdae61", "Low-High" = "#abd9e9",
                  "Not significant" = "grey90")

m_lisa <- tm_shape(shp, crs = 3577) +
  tm_fill(fill = "grey95") + tm_borders(col = "grey80", lwd = 0.2) +
  tm_shape(shp_proj) +
  tm_fill(fill = "lisa_cluster",
          fill.scale = tm_scale_categorical(values = lisa_colours),
          fill.legend = tm_legend(title = "LISA Cluster")) +
  tm_borders(col = "grey60", lwd = 0.3) +
  tm_title("Local Indicators of Spatial Association: Star Rating Clusters") +
  tm_layout(frame = FALSE)

tmap_save(m_lisa, file.path(dir_maps, "map_lisa_clusters.png"),
          device = ragg::agg_png, width = 10, height = 8, dpi = 300)
message("  Saved: map_lisa_clusters.png")

# Identify specific Low-Low areas for the manuscript
low_low <- lisa_df %>% filter(lisa_cluster == "Low-Low") %>%
  select(SA3_CODE21, SA3_NAME21, STE_NAME21, mean_rating, n_facilities, li_p)
if (nrow(low_low) > 0) {
  message(glue("\n  Low-Low clusters (low quality surrounded by low quality): {nrow(low_low)} SA3s"))
  head(low_low, 10) %>%
    { walk(seq_len(nrow(.)), function(i) {
      message(glue("    {.$SA3_NAME21[i]} ({.$STE_NAME21[i]}): mean = {round(.$mean_rating[i], 2)}, p = {format(.$li_p[i], digits = 3)}"))
    })}
}

# =============================================================================
# I. Spatial Weights Sensitivity
# =============================================================================

message("\n=== I. Spatial weights sensitivity ===\n")

k_values <- c(3, 5, 8, 10, 15)
moran_sensitivity <- tibble(k = integer(), morans_i = numeric(), p_value = numeric())

for (k in k_values) {
  knn_k <- knearneigh(coords, k = k)
  nb_k <- knn2nb(knn_k)
  lw_k <- nb2listw(nb_k, style = "W")
  mt <- moran.test(shp_proj$mean_rating, lw_k)
  moran_sensitivity <- bind_rows(moran_sensitivity,
    tibble(k = k, morans_i = mt$estimate["Moran I statistic"],
           p_value = mt$p.value))
}

message("  Moran's I sensitivity to k:")
for (i in seq_len(nrow(moran_sensitivity))) {
  message(glue("    k = {moran_sensitivity$k[i]}: I = {round(moran_sensitivity$morans_i[i], 3)}, p = {format(moran_sensitivity$p_value[i], digits = 3)}"))
}

write_csv(moran_sensitivity, file.path(dir_tab, "val_morans_i_sensitivity.csv"))
message("  Saved: val_morans_i_sensitivity.csv")

# =============================================================================
# J. Facility-Count Weighted IRSD Regression
# =============================================================================

message("\n=== J. Weighted IRSD regression ===\n")

sa3_with_irsd <- shp_proj %>%
  st_drop_geometry() %>%
  left_join(sa3_ref %>% select(sa3_code, irsd_score_mean),
            by = c("SA3_CODE21" = "sa3_code")) %>%
  filter(!is.na(irsd_score_mean))

# Unweighted
lm_uw <- lm(mean_rating ~ irsd_score_mean, data = sa3_with_irsd)
# Weighted by facility count
lm_wt <- lm(mean_rating ~ irsd_score_mean, data = sa3_with_irsd, weights = n_facilities)

message(glue("  Unweighted: B = {round(coef(lm_uw)[2], 5)}, R2 = {round(summary(lm_uw)$r.squared, 4)}, p = {format(summary(lm_uw)$coefficients[2,4], digits = 3)}"))
message(glue("  Weighted:   B = {round(coef(lm_wt)[2], 5)}, R2 = {round(summary(lm_wt)$r.squared, 4)}, p = {format(summary(lm_wt)$coefficients[2,4], digits = 3)}"))

# =============================================================================
# K. ANOVA Post-hoc Tests
# =============================================================================

message("\n=== K. ANOVA post-hoc (Tukey HSD) ===\n")

sa3_quintile <- sa3_with_irsd %>%
  left_join(sa3_ref %>% select(sa3_code, irsd_quintile),
            by = c("SA3_CODE21" = "sa3_code")) %>%
  filter(!is.na(irsd_quintile))

aov_fit <- aov(mean_rating ~ factor(irsd_quintile), data = sa3_quintile)
tukey <- TukeyHSD(aov_fit)
tukey_df <- as.data.frame(tukey$`factor(irsd_quintile)`) %>%
  tibble::rownames_to_column("comparison") %>%
  mutate(significant = `p adj` < 0.05)

message("  Tukey HSD pairwise comparisons:")
for (i in seq_len(nrow(tukey_df))) {
  sig <- if (tukey_df$significant[i]) "*" else ""
  message(glue("    {tukey_df$comparison[i]}: diff = {round(tukey_df$diff[i], 3)}, p.adj = {format(tukey_df$`p adj`[i], digits = 3)}{sig}"))
}

write_csv(tukey_df, file.path(dir_tab, "val_tukey_hsd_irsd.csv"))
message("  Saved: val_tukey_hsd_irsd.csv")

# =============================================================================
# L. Quality Desert Sensitivity
# =============================================================================

message("\n=== L. Quality desert threshold sensitivity ===\n")

# Reload equity data
latest_spatial <- panel %>%
  filter(quarter_num == max(quarter_num), !is.na(sa3_code), !is.na(overall_star_rating))

facility_locations <- latest_spatial %>%
  group_by(sa3_code) %>%
  summarise(max_rating = max(overall_star_rating), .groups = "drop")

sa3_centroids <- shp %>%
  st_transform(3577) %>%
  left_join(sa3_ref %>% select(sa3_code, irsd_quintile, remoteness, pop_65plus),
            by = c("SA3_CODE21" = "sa3_code"))

sa3_centroids <- suppressWarnings(st_centroid(sa3_centroids))

# Facilities rated 3+ stars
adequate_sa3 <- latest_spatial %>%
  filter(overall_star_rating >= 3) %>%
  distinct(sa3_code) %>%
  pull(sa3_code)

adequate_pts <- sa3_centroids %>% filter(SA3_CODE21 %in% adequate_sa3)

# Compute distance for all SA3s
all_centroids_coords <- st_coordinates(sa3_centroids)
adequate_coords <- st_coordinates(adequate_pts)

if (nrow(adequate_coords) > 0) {
  dist_matrix <- as.matrix(dist(rbind(all_centroids_coords, adequate_coords), method = "euclidean"))
  n_all <- nrow(all_centroids_coords)
  n_adeq <- nrow(adequate_coords)
  dist_to_adequate <- apply(dist_matrix[1:n_all, (n_all+1):(n_all+n_adeq), drop = FALSE], 1, min) / 1000

  sa3_centroids$dist_adequate_km <- dist_to_adequate

  thresholds <- c(30, 50, 75, 100, 150)
  sensitivity <- tibble(
    threshold_km = thresholds,
    n_deserts = sapply(thresholds, function(t) sum(dist_to_adequate > t)),
    pct_sa3 = sapply(thresholds, function(t) round(sum(dist_to_adequate > t) / length(dist_to_adequate) * 100, 1)),
    pop_65plus_affected = sapply(thresholds, function(t) {
      desert_sa3 <- sa3_centroids$SA3_CODE21[dist_to_adequate > t]
      sum(sa3_ref$pop_65plus[sa3_ref$sa3_code %in% desert_sa3], na.rm = TRUE)
    })
  )

  message("  Quality desert sensitivity to distance threshold:")
  for (i in seq_len(nrow(sensitivity))) {
    message(glue("    >{sensitivity$threshold_km[i]}km: {sensitivity$n_deserts[i]} SA3s ({sensitivity$pct_sa3[i]}%), 65+ pop = {format(sensitivity$pop_65plus_affected[i], big.mark = ',')}"))
  }

  write_csv(sensitivity, file.path(dir_tab, "val_desert_threshold_sensitivity.csv"))
  message("  Saved: val_desert_threshold_sensitivity.csv")
}

# =============================================================================
# M. Quality Desert Regression
# =============================================================================

message("\n=== M. Quality desert regression ===\n")

desert_data <- sa3_centroids %>%
  st_drop_geometry() %>%
  filter(!is.na(irsd_quintile), !is.na(remoteness)) %>%
  mutate(
    is_desert = as.integer(dist_adequate_km > 50),
    log_dist = log1p(dist_adequate_km),
    remoteness_f = factor(remoteness),
    irsd_q_f = factor(irsd_quintile)
  )

# M1: Linear model for log-distance
lm_dist <- lm(log_dist ~ irsd_q_f + remoteness_f, data = desert_data)
message("  Linear model: log(1 + distance) ~ IRSD quintile + remoteness")
lm_dist_summary <- summary(lm_dist)
message(glue("  R-squared: {round(lm_dist_summary$r.squared, 3)}, adj R2: {round(lm_dist_summary$adj.r.squared, 3)}"))

# Report IRSD coefficients controlling for remoteness
irsd_coefs <- coef(summary(lm_dist))
irsd_rows <- grep("irsd", rownames(irsd_coefs))
message("  IRSD quintile effects (controlling for remoteness):")
for (r in irsd_rows) {
  message(glue("    {rownames(irsd_coefs)[r]}: B = {round(irsd_coefs[r, 1], 3)}, p = {format(irsd_coefs[r, 4], digits = 3)}"))
}

# M2: Logistic model for desert status
if (sum(desert_data$is_desert) >= 5) {
  glm_desert <- glm(is_desert ~ irsd_q_f + remoteness_f, data = desert_data, family = binomial)
  message("\n  Logistic model: is_desert ~ IRSD quintile + remoteness")
  glm_summary <- summary(glm_desert)
  message(glue("  AIC: {round(glm_summary$aic, 1)}"))

  # Check if IRSD matters after controlling for remoteness
  glm_remote_only <- glm(is_desert ~ remoteness_f, data = desert_data, family = binomial)
  lr_test <- anova(glm_remote_only, glm_desert, test = "Chisq")
  message(glue("  LRT (adding IRSD to remoteness-only model): chi2 = {round(lr_test$Deviance[2], 2)}, p = {format(lr_test$`Pr(>Chi)`[2], digits = 3)}"))

  if (lr_test$`Pr(>Chi)`[2] > 0.05) {
    message("  -> IRSD does NOT significantly predict desert status after controlling for remoteness")
  } else {
    message("  -> IRSD DOES significantly predict desert status even after controlling for remoteness")
  }
} else {
  message("  Too few deserts for logistic regression (n < 5)")
}

# M3: Fisher's exact test (instead of chi-squared)
fisher_tab <- table(
  desert = desert_data$is_desert,
  irsd = desert_data$irsd_quintile
)
fisher_result <- fisher.test(fisher_tab, simulate.p.value = TRUE, B = 10000)
message(glue("\n  Fisher's exact test (desert x IRSD): p = {format(fisher_result$p.value, digits = 3)}"))

# =============================================================================
# N. Transition Probability Model
# =============================================================================

message("\n=== N. Transition probability model ===\n")

transitions <- dat %>%
  arrange(facility_id, quarter_num) %>%
  group_by(facility_id) %>%
  mutate(
    prev_rating = lag(overall_star_rating),
    improved = as.integer(overall_star_rating > prev_rating),
    declined = as.integer(overall_star_rating < prev_rating)
  ) %>%
  filter(!is.na(prev_rating)) %>%
  ungroup()

# Multinomial: outcome = improved/stable/declined
transitions <- transitions %>%
  mutate(
    transition = case_when(
      improved == 1 ~ "Improved",
      declined == 1 ~ "Declined",
      TRUE ~ "Stable"
    ),
    transition_f = relevel(factor(transition), ref = "Stable")
  )

# Use binomial logistic for improvement vs not (simpler, more interpretable)
glm_improve <- glm(improved ~ purpose_f + prev_rating + quarter_c,
                   data = transitions, family = binomial)

message("  Logistic model: P(improve) ~ provider type + previous rating + time")
improve_coefs <- coef(summary(glm_improve))
for (r in seq_len(nrow(improve_coefs))) {
  or <- round(exp(improve_coefs[r, 1]), 3)
  message(glue("    {rownames(improve_coefs)[r]}: OR = {or}, p = {format(improve_coefs[r, 4], digits = 3)}"))
}

# Same for decline
glm_decline <- glm(declined ~ purpose_f + prev_rating + quarter_c,
                   data = transitions, family = binomial)

message("\n  Logistic model: P(decline) ~ provider type + previous rating + time")
decline_coefs <- coef(summary(glm_decline))
for (r in seq_len(nrow(decline_coefs))) {
  or <- round(exp(decline_coefs[r, 1]), 3)
  message(glue("    {rownames(decline_coefs)[r]}: OR = {or}, p = {format(decline_coefs[r, 4], digits = 3)}"))
}

# Save transition model results
transition_models <- bind_rows(
  as.data.frame(improve_coefs) %>% tibble::rownames_to_column("term") %>%
    mutate(model = "P(improve)", OR = exp(Estimate)),
  as.data.frame(decline_coefs) %>% tibble::rownames_to_column("term") %>%
    mutate(model = "P(decline)", OR = exp(Estimate))
)
write_csv(transition_models, file.path(dir_tab, "val_transition_logistic_models.csv"))
message("  Saved: val_transition_logistic_models.csv")

# =============================================================================
# O. Multiple Testing Corrections
# =============================================================================

message("\n=== O. Multiple testing corrections ===\n")

# Collect all p-values from sub-category tests
p_values <- tibble(
  test = c(
    "Moran's I: Overall", "Moran's I: Residents' Experience",
    "Moran's I: Compliance", "Moran's I: Staffing", "Moran's I: Quality Measures",
    "IRSD cor: Overall", "IRSD cor: Residents' Experience",
    "IRSD cor: Compliance", "IRSD cor: Staffing", "IRSD cor: Quality Measures"
  ),
  p_original = c(
    2.57e-13, 5.81e-48, 0.069, 1.61e-07, 8.44e-07,
    0.022, 0.014, 0.618, 0.047, 0.950
  )
)

p_values <- p_values %>%
  mutate(
    p_bonferroni = pmin(p_original * n(), 1),
    p_holm = p.adjust(p_original, method = "holm"),
    p_fdr  = p.adjust(p_original, method = "fdr"),
    sig_original = p_original < 0.05,
    sig_bonferroni = p_bonferroni < 0.05,
    sig_fdr = p_fdr < 0.05
  )

message("  Multiple testing corrections (10 tests):")
for (i in seq_len(nrow(p_values))) {
  orig <- if (p_values$sig_original[i]) "*" else ""
  bonf <- if (p_values$sig_bonferroni[i]) "*" else ""
  fdr  <- if (p_values$sig_fdr[i]) "*" else ""
  message(glue("    {p_values$test[i]}: p = {format(p_values$p_original[i], digits = 3)}{orig} | Bonf = {format(p_values$p_bonferroni[i], digits = 3)}{bonf} | FDR = {format(p_values$p_fdr[i], digits = 3)}{fdr}"))
}

# Count how many survive correction
n_orig <- sum(p_values$sig_original)
n_bonf <- sum(p_values$sig_bonferroni)
n_fdr  <- sum(p_values$sig_fdr)
message(glue("\n  Significant at 0.05: Original = {n_orig}/10, Bonferroni = {n_bonf}/10, FDR = {n_fdr}/10"))

write_csv(p_values, file.path(dir_tab, "val_multiple_testing_corrections.csv"))
message("  Saved: val_multiple_testing_corrections.csv")

# =============================================================================
# Summary
# =============================================================================

message("\n=== Validation & Enrichment Complete ===\n")
message("  Diagnostic plots: 4 (residuals, QQ, random effects, time)")
message("  Validation figures: 1 (random slopes distribution)")
message("  Validation tables: ~12")
message("  LISA cluster map: 1")
message("\n  All outputs saved to outputs/figures/, outputs/tables/, outputs/maps/")
message("\nDone.")
