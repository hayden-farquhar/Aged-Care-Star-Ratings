# =============================================================================
# 07_manuscript_tables.R
# Generate publication-ready manuscript tables (CSV format for Word conversion)
#
# Tables:
#   Table 1: Facility characteristics by provider type
#   Table 2: Mixed model results (random slopes)
#   Table 3: International comparison
#
# Author: Hayden Farquhar
# Created: 2026-02-16
# =============================================================================

library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(readr)
library(glue)
library(lmerTest)
library(broom.mixed)

project_root <- here::here()
dir_proc <- file.path(project_root, "data", "processed")
dir_tab  <- file.path(project_root, "outputs", "tables")
dir.create(dir_tab, showWarnings = FALSE, recursive = TRUE)

panel   <- readRDS(file.path(dir_proc, "star_ratings_panel.rds"))
sa3_ref <- readRDS(file.path(dir_proc, "sa3_reference.rds"))

# =============================================================================
# Table 1: Facility characteristics by provider type
# =============================================================================

message("\n=== Table 1: Facility characteristics ===\n")

latest_q <- max(panel$quarter_num)
latest <- panel %>% filter(quarter_num == latest_q, !is.na(overall_star_rating))

# Helper for n (%)
n_pct <- function(x, total) paste0(x, " (", round(x / total * 100, 1), ")")

# Build table by provider type
providers <- c("For-profit", "Not-for-profit", "Government")
overall_n <- nrow(latest)

tab1_rows <- list()

# N facilities
tab1_rows[["N facilities"]] <- c(
  sapply(providers, function(p) sum(latest$purpose == p)),
  overall_n
)

# Mean overall rating (SD)
for (p in c(providers, "All")) {
  d <- if (p == "All") latest else latest %>% filter(purpose == p)
  m <- round(mean(d$overall_star_rating, na.rm = TRUE), 2)
  s <- round(sd(d$overall_star_rating, na.rm = TRUE), 2)
  tab1_rows[[paste0("mean_overall_", p)]] <- paste0(m, " (", s, ")")
}

# Rating distribution
for (r in 1:5) {
  vals <- c()
  for (p in providers) {
    d <- latest %>% filter(purpose == p)
    vals <- c(vals, n_pct(sum(d$overall_star_rating == r), nrow(d)))
  }
  vals <- c(vals, n_pct(sum(latest$overall_star_rating == r), overall_n))
  tab1_rows[[paste0(r, "-star")]] <- vals
}

# Sub-category means
subcats <- list(
  c("residents_experience_rating", "Residents' experience"),
  c("compliance_rating", "Compliance"),
  c("staffing_rating", "Staffing"),
  c("quality_measures_rating", "Quality measures")
)
for (sc in subcats) {
  vals <- c()
  for (p in c(providers, "All")) {
    d <- if (p == "All") latest else latest %>% filter(purpose == p)
    m <- round(mean(d[[sc[1]]], na.rm = TRUE), 2)
    s <- round(sd(d[[sc[1]]], na.rm = TRUE), 2)
    vals <- c(vals, paste0(m, " (", s, ")"))
  }
  tab1_rows[[paste0("mean_", sc[2])]] <- vals
}

# Size distribution
for (sz in c("Small", "Medium", "Large")) {
  vals <- c()
  for (p in c(providers, "All")) {
    d <- if (p == "All") latest else latest %>% filter(purpose == p)
    d_nona <- d %>% filter(!is.na(size))
    vals <- c(vals, n_pct(sum(d_nona$size == sz), nrow(d_nona)))
  }
  tab1_rows[[paste0("size_", sz)]] <- vals
}

# Remoteness
for (rem in c("Major Cities of Australia", "Inner Regional Australia",
              "Outer Regional Australia", "Remote Australia", "Very Remote Australia")) {
  vals <- c()
  rem_short <- gsub(" of Australia| Australia", "", rem)
  for (p in c(providers, "All")) {
    d <- if (p == "All") latest else latest %>% filter(purpose == p)
    d_nona <- d %>% filter(!is.na(remoteness))
    vals <- c(vals, n_pct(sum(d_nona$remoteness == rem), nrow(d_nona)))
  }
  tab1_rows[[paste0("remoteness_", rem_short)]] <- vals
}

# State distribution
for (st in sort(unique(latest$state))) {
  vals <- c()
  for (p in c(providers, "All")) {
    d <- if (p == "All") latest else latest %>% filter(purpose == p)
    vals <- c(vals, n_pct(sum(d$state == st), nrow(d)))
  }
  tab1_rows[[paste0("state_", st)]] <- vals
}

# Compile Table 1
tab1 <- tibble(
  characteristic = c(
    "N facilities", "Mean overall rating (SD)",
    "1-star, n (%)", "2-star, n (%)", "3-star, n (%)", "4-star, n (%)", "5-star, n (%)",
    "Residents' experience, mean (SD)", "Compliance, mean (SD)",
    "Staffing, mean (SD)", "Quality measures, mean (SD)",
    "Small, n (%)", "Medium, n (%)", "Large, n (%)",
    "Major Cities, n (%)", "Inner Regional, n (%)",
    "Outer Regional, n (%)", "Remote, n (%)", "Very Remote, n (%)",
    paste0(sort(unique(latest$state)), ", n (%)")
  ),
  for_profit = c(
    tab1_rows[["N facilities"]][1],
    tab1_rows[["mean_overall_For-profit"]],
    sapply(paste0(1:5, "-star"), function(r) tab1_rows[[r]][1]),
    sapply(subcats, function(sc) tab1_rows[[paste0("mean_", sc[2])]][1]),
    sapply(c("Small", "Medium", "Large"), function(s) tab1_rows[[paste0("size_", s)]][1]),
    sapply(c("Major Cities", "Inner Regional", "Outer Regional", "Remote", "Very Remote"),
           function(r) tab1_rows[[paste0("remoteness_", r)]][1]),
    sapply(sort(unique(latest$state)), function(s) tab1_rows[[paste0("state_", s)]][1])
  ),
  not_for_profit = c(
    tab1_rows[["N facilities"]][2],
    tab1_rows[["mean_overall_Not-for-profit"]],
    sapply(paste0(1:5, "-star"), function(r) tab1_rows[[r]][2]),
    sapply(subcats, function(sc) tab1_rows[[paste0("mean_", sc[2])]][2]),
    sapply(c("Small", "Medium", "Large"), function(s) tab1_rows[[paste0("size_", s)]][2]),
    sapply(c("Major Cities", "Inner Regional", "Outer Regional", "Remote", "Very Remote"),
           function(r) tab1_rows[[paste0("remoteness_", r)]][2]),
    sapply(sort(unique(latest$state)), function(s) tab1_rows[[paste0("state_", s)]][2])
  ),
  government = c(
    tab1_rows[["N facilities"]][3],
    tab1_rows[["mean_overall_Government"]],
    sapply(paste0(1:5, "-star"), function(r) tab1_rows[[r]][3]),
    sapply(subcats, function(sc) tab1_rows[[paste0("mean_", sc[2])]][3]),
    sapply(c("Small", "Medium", "Large"), function(s) tab1_rows[[paste0("size_", s)]][3]),
    sapply(c("Major Cities", "Inner Regional", "Outer Regional", "Remote", "Very Remote"),
           function(r) tab1_rows[[paste0("remoteness_", r)]][3]),
    sapply(sort(unique(latest$state)), function(s) tab1_rows[[paste0("state_", s)]][3])
  ),
  all = c(
    tab1_rows[["N facilities"]][4],
    tab1_rows[["mean_overall_All"]],
    sapply(paste0(1:5, "-star"), function(r) tab1_rows[[r]][4]),
    sapply(subcats, function(sc) tab1_rows[[paste0("mean_", sc[2])]][4]),
    sapply(c("Small", "Medium", "Large"), function(s) tab1_rows[[paste0("size_", s)]][4]),
    sapply(c("Major Cities", "Inner Regional", "Outer Regional", "Remote", "Very Remote"),
           function(r) tab1_rows[[paste0("remoteness_", r)]][4]),
    sapply(sort(unique(latest$state)), function(s) tab1_rows[[paste0("state_", s)]][4])
  )
)

write_csv(tab1, file.path(dir_tab, "table1_characteristics.csv"))
message("  Saved: table1_characteristics.csv")

# =============================================================================
# Table 2: Mixed model results
# =============================================================================

message("\n=== Table 2: Mixed model results ===\n")

dat <- panel %>%
  filter(!is.na(overall_star_rating)) %>%
  mutate(
    quarter_c = quarter_num - median(unique(quarter_num)),
    purpose_f = relevel(factor(purpose), ref = "Not-for-profit")
  )

# Random slopes model (preferred)
m_rs <- lmer(overall_star_rating ~ quarter_c * purpose_f + (1 + quarter_c | facility_id),
             data = dat, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000)))

coefs <- tidy(m_rs, effects = "fixed", conf.int = TRUE)

# Format for publication
tab2 <- coefs %>%
  mutate(
    coefficient = case_when(
      term == "(Intercept)" ~ "Intercept (Not-for-profit, mid-study)",
      term == "quarter_c" ~ "Quarter (linear trend)",
      term == "purpose_fFor-profit" ~ "For-profit (vs Not-for-profit)",
      term == "purpose_fGovernment" ~ "Government (vs Not-for-profit)",
      term == "quarter_c:purpose_fFor-profit" ~ "Quarter x For-profit",
      term == "quarter_c:purpose_fGovernment" ~ "Quarter x Government"
    ),
    estimate_fmt = sprintf("%.3f", estimate),
    ci_95 = paste0("(", sprintf("%.3f", conf.low), " to ", sprintf("%.3f", conf.high), ")"),
    p_value = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(coefficient, estimate_fmt, ci_95, p_value)

# Add variance components
vc <- as.data.frame(VarCorr(m_rs))
# Extract individual values safely
vc_intercept <- vc$vcov[vc$grp == "facility_id" & vc$var1 == "(Intercept)" & is.na(vc$var2)]
vc_slope <- vc$vcov[vc$grp == "facility_id" & vc$var1 == "quarter_c" & is.na(vc$var2)]
vc_cov <- vc$vcov[vc$grp == "facility_id" & !is.na(vc$var2)]
vc_resid <- vc$vcov[vc$grp == "Residual"]

var_labels <- c(
  "Random effects",
  "  Facility intercept variance",
  "  Facility slope variance",
  "  Intercept-slope covariance",
  "  Residual variance",
  "",
  "Model fit",
  "  ICC (intercept-only model)",
  "  Marginal R-squared",
  "  Conditional R-squared",
  "  N observations",
  "  N facilities"
)
var_values <- c(
  "",
  sprintf("%.4f", vc_intercept[1]),
  sprintf("%.4f", vc_slope[1]),
  sprintf("%.4f", vc_cov[1]),
  sprintf("%.4f", vc_resid[1]),
  "", "",
  "0.432",
  "0.125",
  "0.503",
  format(nrow(dat), big.mark = ","),
  format(n_distinct(dat$facility_id), big.mark = ",")
)

var_rows <- tibble(
  coefficient = var_labels,
  estimate_fmt = var_values,
  ci_95 = rep("", length(var_labels)),
  p_value = rep("", length(var_labels))
)

tab2_full <- bind_rows(tab2, var_rows)
write_csv(tab2_full, file.path(dir_tab, "table2_mixed_model.csv"))
message("  Saved: table2_mixed_model.csv")

# =============================================================================
# Table 3: International comparison
# =============================================================================

message("\n=== Table 3: International comparison ===\n")

# Read the existing comparison table and reformat
intl <- read_csv(file.path(dir_tab, "intl_summary_comparison.csv"), show_col_types = FALSE)
write_csv(intl, file.path(dir_tab, "table3_international_comparison.csv"))
message("  Saved: table3_international_comparison.csv")

message("\nAll manuscript tables saved to outputs/tables/")
message("Done.")
