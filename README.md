# Replication Code: Geographic Equity in Aged Care Quality

Analysis code for: **"Geographic Equity in Aged Care Quality: Spatial and Longitudinal Analysis of Australia's Star Ratings System"**

**Author:** Hayden Farquhar
**Contact:** hayden.farquhar@icloud.com
**ORCID:** [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X)

## Overview

This repository contains the R code to replicate all analyses in the manuscript. The study examines geographic variation, temporal trajectories, and equity dimensions of quality in Australian residential aged care using the Star Ratings system (May 2023 to October 2025).

Key analyses:
- **Spatial analysis:** Global and Local Moran's I, Geographically Weighted Regression (GWR), LISA cluster identification
- **Longitudinal modelling:** Linear mixed-effects models with facility random intercepts and slopes, cumulative link mixed models (CLMM)
- **Equity analysis:** Quality desert identification, distance-based access measures, remoteness-disadvantage decomposition
- **International comparison:** Descriptive comparison with the US CMS Five-Star system

## Requirements

**R version:** 4.3 or later (developed on R 4.5)

**Required packages:**

```r
install.packages(c(
  # Data wrangling
  "tidyverse", "readxl", "janitor", "glue",
  # Spatial analysis
  "sf", "spdep", "GWmodel", "tmap",
  # Longitudinal models
  "lme4", "lmerTest", "ordinal", "MuMIn",
  # Visualisation
  "ggplot2", "patchwork", "viridis", "RColorBrewer",
  # Data download
  "httr"
))
```

## Data Sources

All data are publicly available. The scripts download most files automatically; one file requires manual download.

| Source | Description | Access |
|--------|-------------|--------|
| [GEN Aged Care Data](https://www.gen-agedcaredata.gov.au/) | Star Ratings quarterly extracts (11 quarters) | Automatic download via `01_data_acquisition.R` |
| [ABS](https://www.abs.gov.au/) | SEIFA 2021 (IRSD at SA2 level) | Automatic download |
| [ABS](https://www.abs.gov.au/) | Remoteness Areas 2021 | Automatic download |
| [ABS](https://www.abs.gov.au/) | SA3 boundary shapefiles (GDA2020) | Automatic download |
| [ABS](https://www.abs.gov.au/) | SA1/SA2/SA3 allocation files | Automatic download |
| [ABS](https://www.abs.gov.au/) | Census 2021 population by age (G04B) | **Manual download** (see below) |
| [data.cms.gov](https://data.cms.gov/) | US CMS Five-Star Provider Information | Automatic download |

### Manual download: Census 2021 DataPack

The ABS Census DataPacks cannot be downloaded programmatically. To obtain the population-by-age data:

1. Go to https://www.abs.gov.au/census/find-census-data/datapacks
2. Select: **2021 Census DataPacks** > **General Community Profile**
3. Select geography: **Statistical Area Level 3 (SA3)**
4. Download the DataPack zip file
5. Extract the file `2021Census_G04B_AUST_SA3.csv` from the zip
6. Place it in `data/reference/2021Census_G04B_AUST_SA3.csv`

## Directory Structure

```
repository/
├── README.md                  # This file
├── LICENSE                    # MIT License
├── METADATA.txt               # Manuscript metadata
├── scripts/
│   ├── 01_data_acquisition.R  # Download all raw data
│   ├── 02_data_cleaning.R     # Build longitudinal panel dataset
│   ├── 03_spatial_analysis.R  # Moran's I, LISA, GWR
│   ├── 03b_save_maps.R        # Generate choropleth maps
│   ├── 04_longitudinal_models.R       # Mixed-effects models
│   ├── 04b_validation_enrichment.R    # Sensitivity analyses
│   ├── 04b_supplement.R       # Supplementary table generation
│   ├── 04c_final_analyses.R   # Interaction tests, persistent low performers
│   ├── 05_equity_synthesis.R  # Quality deserts, distance analysis
│   └── 06_international_comparison.R  # US CMS Five-Star comparison
├── data/                      # Created by scripts (not tracked in git)
│   ├── raw/                   # Downloaded source files
│   ├── reference/             # ABS reference files
│   ├── spatial/               # SA3 shapefiles
│   └── processed/             # Cleaned panel dataset
└── outputs/                   # Created by scripts (not tracked in git)
    ├── tables/                # Results tables (CSV)
    ├── figures/               # Plots (PNG)
    └── maps/                  # Choropleth maps (PNG)
```

## Running the Analysis

Scripts are numbered and should be run sequentially. All scripts assume the working directory is the repository root.

```r
setwd("/path/to/this/repository")
```

### Step 1: Data acquisition
```r
source("scripts/01_data_acquisition.R")
```
Downloads Star Ratings quarterly extracts, ABS reference files, shapefiles, and US CMS data. **Note:** You must manually download the Census DataPack file first (see above).

### Step 2: Data cleaning
```r
source("scripts/02_data_cleaning.R")
```
Standardises facility identifiers across quarters, links facilities to SA3 regions via suburb matching, and builds the longitudinal panel dataset (`data/processed/star_ratings_panel.rds`).

### Step 3: Spatial analysis
```r
source("scripts/03_spatial_analysis.R")
source("scripts/03b_save_maps.R")
```
Computes SA3-level mean ratings, tests spatial autocorrelation (Global Moran's I with KNN weights), identifies LISA clusters, runs Geographically Weighted Regression, and generates choropleth maps.

### Step 4: Longitudinal models
```r
source("scripts/04_longitudinal_models.R")
source("scripts/04b_validation_enrichment.R")
source("scripts/04b_supplement.R")
source("scripts/04c_final_analyses.R")
```
Fits linear mixed-effects models (random intercepts and slopes), ordinal sensitivity models (CLMM), three-level models, adjusted models, transition probability models, persistent low-performer analysis, and all sensitivity checks.

### Step 5: Equity synthesis
```r
source("scripts/05_equity_synthesis.R")
```
Computes distances to nearest adequately rated facility, identifies quality deserts, tests the independence of disadvantage and remoteness effects, and generates equity maps.

### Step 6: International comparison
```r
source("scripts/06_international_comparison.R")
```
Downloads and processes US CMS Five-Star data, computes descriptive comparison statistics, and generates comparison figures.

## Key Output Files

After running all scripts, the following key outputs are produced:

| File | Description |
|------|-------------|
| `data/processed/star_ratings_panel.rds` | Longitudinal panel (28,708 facility-quarter records) |
| `data/processed/sa3_reference.rds` | SA3 reference table with IRSD, remoteness, population |
| `outputs/tables/model_comparison_table.csv` | Mixed-effects model results |
| `outputs/tables/persistent_low_performers.csv` | 28 persistently low-performing facilities |
| `outputs/tables/lisa_clusters.csv` | LISA cluster assignments |
| `outputs/tables/quality_desert_sa3s.csv` | Quality desert SA3 areas |
| `outputs/maps/overall_star_rating_map.png` | Choropleth map of mean ratings |
| `outputs/maps/lisa_cluster_map.png` | LISA cluster map |
| `outputs/figures/quality_trajectories.png` | Rating trajectories by provider type |
| `outputs/maps/quality_desert_map.png` | Quality desert map |

## Notes

- **Star Ratings data updates:** The GEN website publishes new quarterly extracts periodically. URLs may change; `01_data_acquisition.R` includes the URLs valid as of February 2026. If downloads fail, check the [publication page](https://www.health.gov.au/resources/collections/star-ratings-quarterly-data-extracts) for current links.
- **Data period:** The analysis covers May 2023 to October 2025, the final quarter under the original Star Ratings methodology before the Staffing and Compliance rating redesigns of October--November 2025.
- **Computation time:** The full pipeline takes approximately 10--15 minutes on a modern laptop. The GWR step (in `03_spatial_analysis.R`) is the most computationally intensive.
- **Reproducibility:** Random number generation is not used in any analysis. All results are fully deterministic given the same input data.

## Citation

If you use this code, please cite the manuscript:

> Farquhar H. Geographic equity in aged care quality: spatial and longitudinal analysis of Australia's Star Ratings system. *Australian Health Review.* 2026. [DOI to be added]

## License

MIT License. See [LICENSE](LICENSE) for details.
