# =============================================================================
# 04_longitudinal_models.R
# Longitudinal trajectory analysis of aged care Star Ratings
#
# Analyses:
#   1. Mixed-effects models: Rating ~ quarter * provider_type + (1|facility)
#   2. Sub-category trajectory models
#   3. Persistently low-rated facilities (<=2 stars for 3+ consecutive quarters)
#   4. Transition probabilities (improve, maintain, decline) between quarters
#   5. Trajectory plots by provider type and remoteness
#
# Inputs:
#   data/processed/star_ratings_panel.rds
#
# Outputs:
#   outputs/figures/  — trajectory plots, transition heatmaps
#   outputs/tables/   — model coefficients, transition matrices, persistent low
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom.mixed)
library(glue)

project_root <- here::here()
dir_proc   <- file.path(project_root, "data", "processed")
dir_figs   <- file.path(project_root, "outputs", "figures")
dir_tables <- file.path(project_root, "outputs", "tables")


# =============================================================================
# PART 1: Load and prepare data
# =============================================================================

message("\n=== Part 1: Data preparation ===\n")

panel <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))

# Centre quarter_num for interpretability (intercept = midpoint)
panel <- panel %>%
  mutate(
    quarter_c  = quarter_num - median(unique(quarter_num)),
    purpose_f  = factor(purpose, levels = c("Not-for-profit", "For-profit", "Government"))
  )

# Summary
n_obs <- sum(!is.na(panel$overall_star_rating))
n_fac <- n_distinct(panel$facility_id[!is.na(panel$overall_star_rating)])
n_q   <- n_distinct(panel$quarter_num)
message(glue("  Observations with ratings: {n_obs}"))
message(glue("  Facilities: {n_fac}"))
message(glue("  Quarters: {n_q}"))
message(glue("  Observations per facility: mean = {round(n_obs/n_fac, 1)}"))

by_purpose <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  count(purpose_f, name = "n")
message("\n  By provider type:")
for (i in seq_len(nrow(by_purpose))) {
  message(glue("    {by_purpose$purpose_f[i]}: {by_purpose$n[i]}"))
}


# =============================================================================
# PART 2: Mixed-effects models — Overall Star Rating
# =============================================================================

message("\n=== Part 2: Mixed-effects models (overall rating) ===\n")

dat <- panel %>% filter(!is.na(overall_star_rating), !is.na(purpose_f))

# Model 1: Unconditional growth (time only)
message("  Fitting Model 1: unconditional growth ...")
m1 <- lmer(overall_star_rating ~ quarter_c + (1 | facility_id), data = dat)

# Model 2: Add provider type
message("  Fitting Model 2: + provider type ...")
m2 <- lmer(overall_star_rating ~ quarter_c + purpose_f + (1 | facility_id), data = dat)

# Model 3: Add interaction (differential trajectories)
message("  Fitting Model 3: + quarter x provider type ...")
m3 <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 | facility_id), data = dat)

# Model comparison
anova_12 <- anova(m1, m2)
anova_23 <- anova(m2, m3)

message("\n  Model comparison:")
message(glue("    M1 (time only) AIC: {round(AIC(m1), 1)}"))
message(glue("    M2 (+ provider) AIC: {round(AIC(m2), 1)}  (vs M1: chi2 = {round(anova_12$Chisq[2], 2)}, p = {format.pval(anova_12$`Pr(>Chisq)`[2], digits = 3)})"))
message(glue("    M3 (+ interaction) AIC: {round(AIC(m3), 1)}  (vs M2: chi2 = {round(anova_23$Chisq[2], 2)}, p = {format.pval(anova_23$`Pr(>Chisq)`[2], digits = 3)})"))

# Extract Model 3 coefficients (the full model)
m3_tidy <- tidy(m3, effects = "fixed", conf.int = TRUE)
message("\n  Model 3 fixed effects:")
for (i in seq_len(nrow(m3_tidy))) {
  message(glue("    {m3_tidy$term[i]}: B = {round(m3_tidy$estimate[i], 4)} [{round(m3_tidy$conf.low[i], 4)}, {round(m3_tidy$conf.high[i], 4)}], p = {format.pval(m3_tidy$p.value[i], digits = 3)}"))
}

# Variance components
vc <- as.data.frame(VarCorr(m3))
icc <- vc$vcov[1] / sum(vc$vcov)
message(glue("\n  ICC (facility): {round(icc, 3)}"))
message(glue("  Facility variance: {round(vc$vcov[1], 4)}"))
message(glue("  Residual variance: {round(vc$vcov[2], 4)}"))

# Save model results
write.csv(m3_tidy, file.path(dir_tables, "mixed_model_overall_coefficients.csv"),
          row.names = FALSE)

model_summary <- tibble(
  model    = c("M1: Time only", "M2: + Provider type", "M3: + Interaction"),
  AIC      = round(c(AIC(m1), AIC(m2), AIC(m3)), 1),
  BIC      = round(c(BIC(m1), BIC(m2), BIC(m3)), 1),
  logLik   = round(c(logLik(m1), logLik(m2), logLik(m3)), 1),
  chi2_vs_prev = c(NA, round(anova_12$Chisq[2], 2), round(anova_23$Chisq[2], 2)),
  p_vs_prev    = c(NA, anova_12$`Pr(>Chisq)`[2], anova_23$`Pr(>Chisq)`[2])
)
write.csv(model_summary, file.path(dir_tables, "mixed_model_comparison.csv"),
          row.names = FALSE)
message("  Saved: mixed_model_overall_coefficients.csv, mixed_model_comparison.csv")


# =============================================================================
# PART 3: Sub-category trajectory models
# =============================================================================

message("\n=== Part 3: Sub-category trajectory models ===\n")

subcat_vars <- c(
  "residents_experience_rating",
  "compliance_rating",
  "staffing_rating",
  "quality_measures_rating"
)
subcat_labels <- c("Residents' Experience", "Compliance", "Staffing", "Quality Measures")

subcat_results <- list()

for (i in seq_along(subcat_vars)) {
  var <- subcat_vars[i]
  label <- subcat_labels[i]

  dat_sub <- panel %>%
    filter(!is.na(.data[[var]]), !is.na(purpose_f))

  if (nrow(dat_sub) < 100) {
    message(glue("  {label}: insufficient data ({nrow(dat_sub)} obs), skipping"))
    next
  }

  message(glue("  Fitting: {label} ({nrow(dat_sub)} obs) ..."))

  m <- tryCatch(
    lmer(as.formula(paste(var, "~ quarter_c * purpose_f + (1 | facility_id)")),
         data = dat_sub),
    error = function(e) {
      message(glue("    [ERROR] {e$message}"))
      NULL
    }
  )

  if (!is.null(m)) {
    coefs <- tidy(m, effects = "fixed", conf.int = TRUE) %>%
      mutate(outcome = label)
    subcat_results[[label]] <- coefs

    # Time slope
    time_coef <- coefs %>% filter(term == "quarter_c")
    message(glue("    Time effect: B = {round(time_coef$estimate, 4)}, p = {format.pval(time_coef$p.value, digits = 3)}"))

    # Interaction terms
    int_terms <- coefs %>% filter(grepl(":", term))
    for (j in seq_len(nrow(int_terms))) {
      message(glue("    {int_terms$term[j]}: B = {round(int_terms$estimate[j], 4)}, p = {format.pval(int_terms$p.value[j], digits = 3)}"))
    }
  }
}

all_subcat_coefs <- bind_rows(subcat_results)
write.csv(all_subcat_coefs, file.path(dir_tables, "mixed_model_subcategory_coefficients.csv"),
          row.names = FALSE)
message("  Saved: mixed_model_subcategory_coefficients.csv")


# =============================================================================
# PART 4: Trajectory plots
# =============================================================================

message("\n=== Part 4: Trajectory plots ===\n")

# --- 4a: Mean trajectory by provider type ---
traj_data <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  group_by(quarter_num, quarter_label, quarter_date, purpose) %>%
  summarise(
    n = n(),
    mean_rating = mean(overall_star_rating, na.rm = TRUE),
    se = sd(overall_star_rating, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_traj <- ggplot(traj_data, aes(x = quarter_date, y = mean_rating,
                                 colour = purpose, group = purpose)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = mean_rating - 1.96 * se, ymax = mean_rating + 1.96 * se,
                  fill = purpose), alpha = 0.15, colour = NA) +
  scale_colour_manual(values = c("For-profit" = "#e41a1c",
                                  "Not-for-profit" = "#377eb8",
                                  "Government" = "#4daf4a")) +
  scale_fill_manual(values = c("For-profit" = "#e41a1c",
                                "Not-for-profit" = "#377eb8",
                                "Government" = "#4daf4a")) +
  labs(
    x = "Quarter",
    y = "Mean Overall Star Rating",
    colour = "Provider type",
    fill = "Provider type",
    title = "Star Rating Trajectories by Provider Type",
    subtitle = "Mean overall rating with 95% CI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(3, 4.5))

ggsave(file.path(dir_figs, "trajectory_by_provider_type.png"), p_traj,
       width = 10, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: trajectory_by_provider_type.png")

# --- 4b: Sub-category trajectories (all providers combined) ---
subcat_traj <- panel %>%
  select(quarter_num, quarter_date, facility_id,
         residents_experience_rating, compliance_rating,
         staffing_rating, quality_measures_rating) %>%
  pivot_longer(
    cols = c(residents_experience_rating, compliance_rating,
             staffing_rating, quality_measures_rating),
    names_to = "subcategory",
    values_to = "rating"
  ) %>%
  filter(!is.na(rating)) %>%
  mutate(subcategory = case_when(
    subcategory == "residents_experience_rating" ~ "Residents' Experience",
    subcategory == "compliance_rating"           ~ "Compliance",
    subcategory == "staffing_rating"             ~ "Staffing",
    subcategory == "quality_measures_rating"     ~ "Quality Measures"
  )) %>%
  group_by(quarter_num, quarter_date, subcategory) %>%
  summarise(
    n = n(),
    mean_rating = mean(rating, na.rm = TRUE),
    se = sd(rating, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_subcat <- ggplot(subcat_traj, aes(x = quarter_date, y = mean_rating,
                                     colour = subcategory, group = subcategory)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Quarter",
    y = "Mean Rating",
    colour = "Sub-category",
    title = "Star Rating Sub-category Trajectories",
    subtitle = "Mean ratings across all providers"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(2.0, 5.0))

ggsave(file.path(dir_figs, "trajectory_subcategories.png"), p_subcat,
       width = 10, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: trajectory_subcategories.png")

# --- 4c: Provider type x sub-category panel plot ---
subcat_prov_traj <- panel %>%
  select(quarter_num, quarter_date, purpose,
         residents_experience_rating, compliance_rating,
         staffing_rating, quality_measures_rating) %>%
  pivot_longer(
    cols = c(residents_experience_rating, compliance_rating,
             staffing_rating, quality_measures_rating),
    names_to = "subcategory",
    values_to = "rating"
  ) %>%
  filter(!is.na(rating), !is.na(purpose)) %>%
  mutate(subcategory = case_when(
    subcategory == "residents_experience_rating" ~ "Residents' Experience",
    subcategory == "compliance_rating"           ~ "Compliance",
    subcategory == "staffing_rating"             ~ "Staffing",
    subcategory == "quality_measures_rating"     ~ "Quality Measures"
  )) %>%
  group_by(quarter_num, quarter_date, purpose, subcategory) %>%
  summarise(mean_rating = mean(rating, na.rm = TRUE), .groups = "drop")

p_panel <- ggplot(subcat_prov_traj, aes(x = quarter_date, y = mean_rating,
                                         colour = purpose)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~subcategory, scales = "free_y") +
  scale_colour_manual(values = c("For-profit" = "#e41a1c",
                                  "Not-for-profit" = "#377eb8",
                                  "Government" = "#4daf4a")) +
  labs(
    x = "Quarter", y = "Mean Rating", colour = "Provider type",
    title = "Sub-category Trajectories by Provider Type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(dir_figs, "trajectory_subcategory_by_provider.png"), p_panel,
       width = 12, height = 8, dpi = 300, device = ragg::agg_png)
message("  Saved: trajectory_subcategory_by_provider.png")


# =============================================================================
# PART 5: Persistently low-rated facilities
# =============================================================================

message("\n=== Part 5: Persistently low-rated facilities ===\n")

# Identify facilities with <=2 stars for 3+ consecutive quarters
compute_max_consec_low <- function(ratings) {
  is_low <- ratings <= 2
  if (!any(is_low)) return(0L)
  runs <- rle(is_low)
  low_runs <- runs$lengths[runs$values]
  if (length(low_runs) == 0) return(0L)
  max(low_runs)
}

facility_consec <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  arrange(facility_id, quarter_num) %>%
  group_by(facility_id) %>%
  summarise(
    n_quarters     = n(),
    n_low          = sum(overall_star_rating <= 2),
    max_consec_low = compute_max_consec_low(overall_star_rating),
    mean_rating    = mean(overall_star_rating),
    first_rating   = first(overall_star_rating),
    last_rating    = last(overall_star_rating),
    purpose        = first(purpose),
    state          = first(state),
    remoteness     = first(remoteness),
    sa3_name       = first(sa3_name),
    .groups = "drop"
  )

persistent_low <- facility_consec %>% filter(max_consec_low >= 3)

message(glue("  Total facilities with rating data: {nrow(facility_consec)}"))
message(glue("  Ever rated <=2 stars: {sum(facility_consec$n_low > 0)}"))
message(glue("  Persistently low (>=3 consecutive quarters <=2 stars): {nrow(persistent_low)}"))

if (nrow(persistent_low) > 0) {
  message("\n  Persistent low performers by provider type:")
  pl_by_purpose <- persistent_low %>% count(purpose)
  for (i in seq_len(nrow(pl_by_purpose))) {
    message(glue("    {pl_by_purpose$purpose[i]}: {pl_by_purpose$n[i]}"))
  }

  message("\n  By state:")
  pl_by_state <- persistent_low %>% count(state) %>% arrange(desc(n))
  for (i in seq_len(nrow(pl_by_state))) {
    message(glue("    {pl_by_state$state[i]}: {pl_by_state$n[i]}"))
  }

  write.csv(persistent_low, file.path(dir_tables, "persistent_low_performers.csv"),
            row.names = FALSE)
  message("  Saved: persistent_low_performers.csv")
}

# Also identify improvers and decliners
facility_consec <- facility_consec %>%
  mutate(
    trajectory = case_when(
      last_rating > first_rating  ~ "Improved",
      last_rating < first_rating  ~ "Declined",
      last_rating == first_rating ~ "Stable",
      TRUE                        ~ "Unknown"
    ),
    change = last_rating - first_rating
  )

traj_summary <- facility_consec %>%
  count(trajectory) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

message("\n  Facility trajectory (first to last quarter):")
for (i in seq_len(nrow(traj_summary))) {
  message(glue("    {traj_summary$trajectory[i]}: {traj_summary$n[i]} ({traj_summary$pct[i]}%)"))
}

# By provider type
traj_by_purpose <- facility_consec %>%
  group_by(purpose) %>%
  summarise(
    n = n(),
    pct_improved = round(mean(trajectory == "Improved") * 100, 1),
    pct_stable   = round(mean(trajectory == "Stable") * 100, 1),
    pct_declined = round(mean(trajectory == "Declined") * 100, 1),
    mean_change  = round(mean(change, na.rm = TRUE), 3),
    .groups = "drop"
  )

message("\n  Trajectory by provider type:")
for (i in seq_len(nrow(traj_by_purpose))) {
  message(glue("    {traj_by_purpose$purpose[i]}: improved={traj_by_purpose$pct_improved[i]}%, stable={traj_by_purpose$pct_stable[i]}%, declined={traj_by_purpose$pct_declined[i]}%, mean change={traj_by_purpose$mean_change[i]}"))
}

write.csv(traj_by_purpose, file.path(dir_tables, "facility_trajectory_by_provider.csv"),
          row.names = FALSE)


# =============================================================================
# PART 6: Transition probabilities
# =============================================================================

message("\n=== Part 6: Transition probabilities ===\n")

# Compute transitions between consecutive quarters
transitions <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  arrange(facility_id, quarter_num) %>%
  group_by(facility_id) %>%
  mutate(
    prev_rating = lag(overall_star_rating),
    prev_quarter = lag(quarter_num)
  ) %>%
  ungroup() %>%
  filter(!is.na(prev_rating), quarter_num == prev_quarter + 1)

transitions <- transitions %>%
  mutate(
    direction = case_when(
      overall_star_rating > prev_rating  ~ "Improved",
      overall_star_rating == prev_rating ~ "Maintained",
      overall_star_rating < prev_rating  ~ "Declined"
    )
  )

# Overall transition probabilities
overall_trans <- transitions %>%
  count(direction) %>%
  mutate(prob = round(n / sum(n), 3))

message("  Overall transition probabilities:")
for (i in seq_len(nrow(overall_trans))) {
  message(glue("    P({overall_trans$direction[i]}) = {overall_trans$prob[i]} (n = {overall_trans$n[i]})"))
}

# Transition matrix (from rating X to rating Y)
trans_matrix <- transitions %>%
  count(prev_rating, overall_star_rating) %>%
  group_by(prev_rating) %>%
  mutate(prob = round(n / sum(n), 3)) %>%
  ungroup()

# Pivot to matrix form
trans_wide <- trans_matrix %>%
  select(prev_rating, overall_star_rating, prob) %>%
  pivot_wider(names_from = overall_star_rating, values_from = prob,
              names_prefix = "to_", values_fill = 0) %>%
  arrange(prev_rating)

message("\n  Transition matrix (probabilities):")
print(as.data.frame(trans_wide))

write.csv(trans_wide, file.path(dir_tables, "transition_matrix.csv"), row.names = FALSE)

# By provider type
trans_by_purpose <- transitions %>%
  group_by(purpose, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(purpose) %>%
  mutate(prob = round(n / sum(n), 3)) %>%
  ungroup()

trans_purpose_wide <- trans_by_purpose %>%
  select(purpose, direction, prob) %>%
  pivot_wider(names_from = direction, values_from = prob, values_fill = 0)

message("\n  Transition probabilities by provider type:")
print(as.data.frame(trans_purpose_wide))

write.csv(trans_purpose_wide, file.path(dir_tables, "transition_probs_by_provider.csv"),
          row.names = FALSE)

# --- 6b: Transition probability over time ---
trans_by_quarter <- transitions %>%
  group_by(quarter_num, quarter_label, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(quarter_num, quarter_label) %>%
  mutate(prob = n / sum(n)) %>%
  ungroup()

p_trans_time <- ggplot(trans_by_quarter, aes(x = quarter_num, y = prob,
                                              colour = direction, group = direction)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("Improved" = "#2ca02c", "Maintained" = "#1f77b4",
                                  "Declined" = "#d62728")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Quarter", y = "Transition Probability",
    colour = "Direction",
    title = "Rating Transition Probabilities Over Time",
    subtitle = "Probability of improving, maintaining, or declining between consecutive quarters"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom")

ggsave(file.path(dir_figs, "transition_probs_over_time.png"), p_trans_time,
       width = 9, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: transition_probs_over_time.png")

# --- 6c: Transition heatmap ---
p_heat <- ggplot(trans_matrix, aes(x = factor(overall_star_rating),
                                    y = factor(prev_rating),
                                    fill = prob)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", prob * 100)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#2c7fb8", name = "Probability") +
  labs(
    x = "Rating (next quarter)",
    y = "Rating (current quarter)",
    title = "Star Rating Transition Probabilities",
    subtitle = "All consecutive quarter-pairs pooled"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold")) +
  coord_equal()

ggsave(file.path(dir_figs, "transition_heatmap.png"), p_heat,
       width = 7, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: transition_heatmap.png")


# =============================================================================
# PART 7: Rating distribution shifts over time
# =============================================================================

message("\n=== Part 7: Rating distribution over time ===\n")

dist_time <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  count(quarter_num, quarter_label, overall_star_rating) %>%
  group_by(quarter_num) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_dist <- ggplot(dist_time, aes(x = factor(quarter_num), y = pct,
                                 fill = factor(overall_star_rating))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = c("1" = "#d73027", "2" = "#fc8d59", "3" = "#fee08b",
               "4" = "#91cf60", "5" = "#1a9850"),
    name = "Star Rating"
  ) +
  labs(
    x = "Quarter", y = "Percentage of facilities (%)",
    title = "Star Rating Distribution Over Time"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold"))

ggsave(file.path(dir_figs, "rating_distribution_over_time.png"), p_dist,
       width = 10, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: rating_distribution_over_time.png")


# =============================================================================
# SUMMARY
# =============================================================================

message("\n=== Phase 3 Complete ===\n")

# Print key model results
message("  Key mixed model results (Model 3: quarter_c * purpose_f):")
key_terms <- m3_tidy %>% filter(term != "(Intercept)")
for (i in seq_len(nrow(key_terms))) {
  sig <- if (key_terms$p.value[i] < 0.001) "***" else if (key_terms$p.value[i] < 0.01) "**" else if (key_terms$p.value[i] < 0.05) "*" else ""
  message(glue("    {key_terms$term[i]}: B = {round(key_terms$estimate[i], 4)}{sig}"))
}
message(glue("\n  ICC: {round(icc, 3)}"))
message(glue("  Persistent low performers: {nrow(persistent_low)}"))
message(glue("  Transition P(maintain): {overall_trans$prob[overall_trans$direction == 'Maintained']}"))

message("\nDone.")
