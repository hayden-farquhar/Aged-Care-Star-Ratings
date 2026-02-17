# =============================================================================
# 04c_final_analyses.R
# Final enrichment analyses before manuscript
#
#   A. IRSD × provider type interaction model
#   B. Predictors of persistent low performance
#   C. Three-level model (facility nested in SA3)
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(lmerTest)
library(broom.mixed)
library(glue)
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
    purpose_f = relevel(factor(purpose), ref = "Not-for-profit"),
    state_f   = factor(state),
    size_f    = factor(size),
    irsd_q_f  = factor(irsd_quintile),
    remoteness_f = factor(remoteness)
  )

# =============================================================================
# A. IRSD × Provider Type Interaction
# =============================================================================

message("\n=== A. IRSD x provider type interaction ===\n")

dat_irsd <- dat %>% filter(!is.na(irsd_q_f))

# Model 1: Main effects only (already done in validation)
m_main <- lmer(overall_star_rating ~ quarter_c * purpose_f + irsd_q_f + (1 | facility_id),
               data = dat_irsd)

# Model 2: Add IRSD x provider type interaction
m_interact <- lmer(overall_star_rating ~ quarter_c * purpose_f + irsd_q_f * purpose_f + (1 | facility_id),
                   data = dat_irsd)

# LRT comparison
anova_interact <- anova(m_main, m_interact)
message(glue("  Main effects AIC: {round(AIC(m_main), 1)}"))
message(glue("  With IRSD x provider interaction AIC: {round(AIC(m_interact), 1)}"))
message(glue("  LRT chi2: {round(anova_interact$Chisq[2], 1)}, p = {format(anova_interact$`Pr(>Chisq)`[2], digits = 3)}"))

coefs_interact <- tidy(m_interact, effects = "fixed", conf.int = TRUE)
write_csv(coefs_interact, file.path(dir_tab, "val_irsd_provider_interaction.csv"))
message("  Saved: val_irsd_provider_interaction.csv")

# Report interaction terms
interact_terms <- coefs_interact %>%
  filter(grepl("irsd_q_f.*purpose|purpose.*irsd_q_f", term))

if (nrow(interact_terms) > 0) {
  message("\n  IRSD x Provider Type interaction terms:")
  walk(seq_len(nrow(interact_terms)), function(i) {
    message(glue("    {interact_terms$term[i]}: B = {round(interact_terms$estimate[i], 3)}, p = {format(interact_terms$p.value[i], digits = 3)}"))
  })

  # Marginal means: for-profit penalty at each IRSD quintile
  message("\n  For-profit penalty by IRSD quintile:")
  fp_main <- coefs_interact %>% filter(term == "purpose_fFor-profit") %>% pull(estimate)
  for (q in 2:5) {
    term_name <- paste0("irsd_q_f", q, ":purpose_fFor-profit")
    interact_est <- interact_terms %>% filter(term == term_name) %>% pull(estimate)
    if (length(interact_est) == 0) interact_est <- 0
    total_fp_effect <- fp_main + interact_est
    message(glue("    Q{q}: for-profit penalty = {round(total_fp_effect, 3)} stars (base {round(fp_main, 3)} + interaction {round(interact_est, 3)})"))
  }
  # Q1 is reference for IRSD, so for-profit effect in Q1 = main effect
  message(glue("    Q1 (ref): for-profit penalty = {round(fp_main, 3)} stars"))
} else {
  message("  No interaction terms found (model may have used different parameterisation)")
}

# Compute marginal means by provider type and IRSD quintile for plotting
message("\n  Marginal mean ratings by provider type x IRSD quintile:")
marginal <- dat_irsd %>%
  group_by(purpose_f, irsd_q_f) %>%
  summarise(mean_rating = mean(overall_star_rating, na.rm = TRUE),
            n = n(), .groups = "drop")

walk(seq_len(nrow(marginal)), function(i) {
  message(glue("    {marginal$purpose_f[i]} x IRSD Q{marginal$irsd_q_f[i]}: {round(marginal$mean_rating[i], 2)} (n={marginal$n[i]})"))
})

# Interaction plot
library(ggplot2)
theme_pub <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, colour = "grey40"),
        panel.grid.minor = element_blank(), legend.position = "bottom")

p_interact <- marginal %>%
  ggplot(aes(x = irsd_q_f, y = mean_rating, colour = purpose_f, group = purpose_f)) +
  geom_point(aes(size = n), alpha = 0.7) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("For-profit" = "#e41a1c", "Not-for-profit" = "#377eb8", "Government" = "#4daf4a")) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(title = "Mean Star Rating by Provider Type and IRSD Quintile",
       subtitle = "Q1 = most disadvantaged, Q5 = least disadvantaged",
       x = "IRSD Quintile", y = "Mean Overall Star Rating",
       colour = "Provider Type") +
  theme_pub +
  scale_y_continuous(limits = c(3, 4.5))

ggsave(file.path(project_root, "outputs", "figures", "val_irsd_provider_interaction.png"), p_interact,
       width = 9, height = 6, dpi = 300, device = ragg::agg_png)
message("  Saved: val_irsd_provider_interaction.png")

# =============================================================================
# B. Predictors of Persistent Low Performance
# =============================================================================

message("\n=== B. Predictors of persistent low performance ===\n")

# Compute max consecutive quarters at <=2 stars for each facility
compute_max_consec_low <- function(ratings, threshold = 2) {
  low <- ratings <= threshold
  if (!any(low, na.rm = TRUE)) return(0L)
  low[is.na(low)] <- FALSE
  runs <- rle(low)
  max(runs$lengths[runs$values], 0L)
}

facility_summary <- dat %>%
  arrange(facility_id, quarter_num) %>%
  group_by(facility_id) %>%
  summarise(
    max_consec_low = compute_max_consec_low(overall_star_rating),
    n_quarters = n(),
    mean_rating = mean(overall_star_rating, na.rm = TRUE),
    purpose = first(purpose),
    state = first(state),
    size = first(size),
    remoteness = first(remoteness),
    irsd_quintile = first(irsd_quintile),
    sa3_code = first(sa3_code),
    .groups = "drop"
  ) %>%
  mutate(
    persistent_low = as.integer(max_consec_low >= 3),
    purpose_f = relevel(factor(purpose), ref = "Not-for-profit"),
    size_f = factor(size),
    state_f = factor(state),
    remoteness_f = factor(remoteness),
    irsd_q_f = factor(irsd_quintile)
  )

n_persistent <- sum(facility_summary$persistent_low, na.rm = TRUE)
message(glue("  Persistent low performers (>=3 consecutive quarters <=2 stars): {n_persistent}"))
message(glue("  Total facilities: {nrow(facility_summary)}"))

# Descriptive breakdown
message("\n  Breakdown by provider type:")
facility_summary %>%
  group_by(purpose) %>%
  summarise(n = n(), n_low = sum(persistent_low), pct = round(n_low/n*100, 1), .groups = "drop") %>%
  { walk(seq_len(nrow(.)), function(i) {
    message(glue("    {.$purpose[i]}: {.$n_low[i]}/{.$n[i]} ({.$pct[i]}%)"))
  })}

message("\n  Breakdown by state:")
facility_summary %>%
  group_by(state) %>%
  summarise(n = n(), n_low = sum(persistent_low), pct = round(n_low/n*100, 1), .groups = "drop") %>%
  filter(n_low > 0) %>%
  arrange(desc(n_low)) %>%
  { walk(seq_len(nrow(.)), function(i) {
    message(glue("    {.$state[i]}: {.$n_low[i]}/{.$n[i]} ({.$pct[i]}%)"))
  })}

# Logistic regression
dat_logistic <- facility_summary %>%
  filter(!is.na(size_f), !is.na(remoteness_f), !is.na(irsd_q_f))

message(glue("\n  Logistic regression (n = {nrow(dat_logistic)}, {sum(dat_logistic$persistent_low)} events)"))

# Check if we have enough events for all predictors
if (sum(dat_logistic$persistent_low) >= 10) {
  # Full model
  glm_persist <- glm(persistent_low ~ purpose_f + size_f + state_f + remoteness_f + irsd_q_f,
                     data = dat_logistic, family = binomial)

  persist_coefs <- coef(summary(glm_persist))
  message("\n  Logistic model: persistent_low ~ provider + size + state + remoteness + IRSD")
  for (r in seq_len(nrow(persist_coefs))) {
    or <- round(exp(persist_coefs[r, 1]), 2)
    ci_lo <- round(exp(persist_coefs[r, 1] - 1.96 * persist_coefs[r, 2]), 2)
    ci_hi <- round(exp(persist_coefs[r, 1] + 1.96 * persist_coefs[r, 2]), 2)
    pval <- format(persist_coefs[r, 4], digits = 3)
    message(glue("    {rownames(persist_coefs)[r]}: OR = {or} ({ci_lo}-{ci_hi}), p = {pval}"))
  }

  # Simpler model (provider type only) for clean reporting
  glm_simple <- glm(persistent_low ~ purpose_f, data = dat_logistic, family = binomial)

  # LRT: does adding covariates change the provider type effect?
  lr_test <- anova(glm_simple, glm_persist, test = "Chisq")
  message(glue("\n  LRT (simple vs full): chi2 = {round(lr_test$Deviance[2], 1)}, p = {format(lr_test$`Pr(>Chi)`[2], digits = 3)}"))

  # Compare for-profit OR
  or_simple <- round(exp(coef(glm_simple)["purpose_fFor-profit"]), 2)
  or_full   <- round(exp(coef(glm_persist)["purpose_fFor-profit"]), 2)
  message(glue("  For-profit OR: unadjusted = {or_simple}, adjusted = {or_full}"))

  # Save
  persist_df <- as.data.frame(persist_coefs) %>%
    tibble::rownames_to_column("term") %>%
    mutate(OR = exp(Estimate),
           OR_lower = exp(Estimate - 1.96 * `Std. Error`),
           OR_upper = exp(Estimate + 1.96 * `Std. Error`))
  write_csv(persist_df, file.path(dir_tab, "val_persistent_low_predictors.csv"))
  message("  Saved: val_persistent_low_predictors.csv")
} else {
  message("  Too few events for multivariable logistic regression")

  # At least do Fisher's exact test for provider type
  fisher_tab <- table(dat_logistic$persistent_low, dat_logistic$purpose_f)
  fisher_result <- fisher.test(fisher_tab)
  message(glue("  Fisher's exact (persistent low x provider type): p = {format(fisher_result$p.value, digits = 3)}"))
}

# =============================================================================
# C. Three-Level Model (Facility Nested in SA3)
# =============================================================================

message("\n=== C. Three-level model (facility nested in SA3) ===\n")

dat_sa3 <- dat %>% filter(!is.na(sa3_code))
n_sa3 <- n_distinct(dat_sa3$sa3_code)
n_fac <- n_distinct(dat_sa3$facility_id)
message(glue("  Data: {nrow(dat_sa3)} records, {n_fac} facilities, {n_sa3} SA3 areas"))

# Two-level (baseline)
m_2level <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 | facility_id),
                 data = dat_sa3)

# Three-level: facility nested in SA3
m_3level <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 | sa3_code/facility_id),
                 data = dat_sa3)

# Compare
anova_23 <- anova(m_2level, m_3level)
message(glue("  Two-level AIC: {round(AIC(m_2level), 1)}"))
message(glue("  Three-level AIC: {round(AIC(m_3level), 1)}"))
message(glue("  LRT chi2: {round(anova_23$Chisq[2], 1)}, p = {format(anova_23$`Pr(>Chisq)`[2], digits = 3)}"))

# Variance decomposition
vc_3 <- as.data.frame(VarCorr(m_3level))
total_var <- sum(vc_3$vcov)
message("\n  Variance decomposition (three-level model):")
for (i in seq_len(nrow(vc_3))) {
  pct <- round(vc_3$vcov[i] / total_var * 100, 1)
  grp_label <- if (vc_3$grp[i] == "facility_id:sa3_code") {
    "Between facilities (within SA3)"
  } else if (vc_3$grp[i] == "sa3_code") {
    "Between SA3 areas"
  } else {
    "Within facilities (residual)"
  }
  message(glue("    {grp_label}: {round(vc_3$vcov[i], 4)} ({pct}%)"))
}

# Compare fixed effects
coefs_2 <- tidy(m_2level, effects = "fixed")
coefs_3 <- tidy(m_3level, effects = "fixed", conf.int = TRUE)

message("\n  Fixed effect comparison (two-level vs three-level):")
for (t in c("quarter_c", "purpose_fFor-profit", "purpose_fGovernment",
            "quarter_c:purpose_fFor-profit", "quarter_c:purpose_fGovernment")) {
  est_2 <- coefs_2 %>% filter(term == t) %>% pull(estimate)
  est_3 <- coefs_3 %>% filter(term == t) %>% pull(estimate)
  se_2  <- coefs_2 %>% filter(term == t) %>% pull(std.error)
  se_3  <- coefs_3 %>% filter(term == t) %>% pull(std.error)
  p_3   <- coefs_3 %>% filter(term == t) %>% pull(p.value)
  message(glue("    {t}: 2L = {round(est_2, 4)} (SE {round(se_2, 4)}), 3L = {round(est_3, 4)} (SE {round(se_3, 4)}), p = {format(p_3, digits = 3)}"))
}

write_csv(coefs_3, file.path(dir_tab, "val_three_level_model_coefficients.csv"))
message("  Saved: val_three_level_model_coefficients.csv")

# Save variance components
vc_df <- vc_3 %>%
  mutate(pct_total = round(vcov / total_var * 100, 1),
         level = case_when(
           grp == "sa3_code" ~ "SA3 area",
           grp == "facility_id:sa3_code" ~ "Facility (within SA3)",
           TRUE ~ "Residual (within facility)"
         )) %>%
  select(level, variance = vcov, pct_total)
write_csv(vc_df, file.path(dir_tab, "val_three_level_variance.csv"))
message("  Saved: val_three_level_variance.csv")

# =============================================================================
# Summary
# =============================================================================

message("\n=== Final Analyses Complete ===\n")
message("  Tables saved: 4")
message("  Figures saved: 1")
message("\nDone.")
