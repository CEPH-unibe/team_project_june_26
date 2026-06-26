# Heat–Health Warning Levels

Annual Team Project 2026

This repo turns district temperature forecasts into a daily heat **warning level
(0–3)**, based on the heat **attributable fraction (AF)** of mortality from
district-specific exposure–response curves (BLUP coefficients fitted upstream).

Two scripts, run in order:

- **`0_get_MMT_baselineAF.R`** — estimates each district's minimum-mortality
  temperature (MMT) and its baseline AF percentiles.
- **`1_calculate_warning_level.R`** — computes forecast AF and assigns a warning
  level to each district-day.

---

## Requirements

R (≥ 4.1) with four packages:

```r
install.packages(c("dlnm", "mixmeta", "dplyr", "lubridate"))
```

---

## Setup

Put the input files in one folder and set `data_dir` at the top of
`0_get_MMT_baselineAF.R` to point at it:

```r
data_dir <- "path/to/your/directory"
```

The file paths just below `data_dir` should all point to where your data actually
live — make them consistent (the same base folder for every input and output).

---

## Inputs & outputs

**Inputs** (place in your data folder):

| File | What it is |
|------|------------|
| `historical_temp_popw_mvngpop_2000_2024_district.csv` | Observed daily temperature by district. |
| `forcast_temp_popw_district.csv` | Forecast daily temperature by district. |
| `secondstage.rds` | Second-stage BLUP coefficients (one element per district). |

**Produced:**

| File | By | Contents |
|------|----|----------|
| `mmt.rds` | script 0 | Minimum-mortality temperature per district. |
| `baseline_summer.csv` | script 0 | Baseline AF percentiles per district. |
| `warning_level.csv` | script 1 | Final output: `BEZNAME`, `BEZNR`, `time`, `warning_level`. |

---

## Warning levels

Default ("fixed") cut-points on the attributable fraction:

| Level | AF range          |
|-------|-------------------|
| 0     | < 0.075           |
| 1     | 0.075 – 0.10      |
| 2     | 0.10 – 0.15       |
| 3     | ≥ 0.15            |

To use district-specific cut-points instead, set `threshold_method <- "baseline"`
in `1_calculate_warning_level.R` (uses the `p75 / p90 / p95` columns from
`baseline_summer.csv`).

---

## Notes

- `varfun` and `varper` (in script 0) must match the specification used to fit the
  BLUP coefficients.
- District `2604` is excluded (no exposure-response curve).
- The BLUP object has no names; district IDs are attached by position, assuming it
  is in the same order as the temperature files.
