# 07b_table1_subcategories.R
# Compute sub-category rating means (SD) by provider type for the Oct 2025 snapshot.
# Adds the four sub-category domains (residents' experience, compliance, staffing,
# quality measures) as rows for Table 1.

suppressMessages({library(dplyr)})

project_root <- here::here()
panel <- readRDS(file.path(project_root, "data", "processed", "star_ratings_panel.rds"))

latest_q <- max(panel$quarter_num)
latest <- panel %>% filter(quarter_num == latest_q, !is.na(overall_star_rating))

cat("Latest quarter_num:", latest_q, " n facilities:", nrow(latest), "\n")
cat("Reporting period label(s):", paste(unique(as.character(latest$reporting_period)), collapse="; "), "\n\n")

subcats <- c(
  "Residents' Experience" = "residents_experience_rating",
  "Compliance"            = "compliance_rating",
  "Staffing"              = "staffing_rating",
  "Quality Measures"      = "quality_measures_rating"
)

providers <- c("For-profit", "Not-for-profit", "Government", "All")

fmt <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

cat(sprintf("%-22s | %-16s | %-18s | %-14s | %-14s\n",
            "Sub-category", "For-profit", "Not-for-profit", "Government", "All"))
for (nm in names(subcats)) {
  col <- subcats[[nm]]
  vals <- sapply(providers, function(p) {
    d <- if (p == "All") latest else latest %>% filter(purpose == p)
    fmt(d[[col]])
  })
  cat(sprintf("%-22s | %-16s | %-18s | %-14s | %-14s\n", nm, vals[1], vals[2], vals[3], vals[4]))
}

# n per provider for header sanity
cat("\nn by provider: ",
    paste(sprintf("%s=%d", providers[1:3],
                  sapply(providers[1:3], function(p) sum(latest$purpose==p))), collapse=", "),
    ", All=", nrow(latest), "\n", sep="")
