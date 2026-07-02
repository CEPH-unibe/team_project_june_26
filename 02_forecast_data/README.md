# MeteoSwiss Forecast Processing Pipeline (ICON-CH2)

## Documentation of the Forecast Data

This repository contains a Python-based automation pipeline for retrieving, processing, and standardizing daily mean 2 m air temperature forecasts from MeteoSwiss using the **ICON-CH2 Numerical Weather Prediction (NWP) gridded model**.

To ensure strict temporal consistency, all operations are pinned to the **00:00 UTC model initialization run** of the execution day.

---

## Data Pipeline Architecture

```text
                           ┌──────────────────────────┐
                           │  MeteoSwiss Data Source  │
                           └─────────────┬────────────┘
                                         │
                                         ▼
                           ┌─────────────────────────┐
                           │   ICON-CH2 Grid (OGD)   │
                           └─────────────┬────────────┘
                                         │
                                         ▼
                           ┌─────────────────────────┐
                           │ Squeeze Dimensions      │
                           │ (Remove eps, z)         │
                           └─────────────┬────────────┘
                                         │
                                         ▼
                           ┌─────────────────────────┐
                           │ Spatial Grouping &      │
                           │ 24 h Block Mean         │
                           └─────────────┬────────────┘
                                         │
                               ┌─────────┴─────────┐
                               ▼                   ▼
                         ┌───────────┐       ┌───────────┐
                         │ Export:   │       │ Apply     │
                         │ Native    │       │ Conservative
                         │ NetCDF    │       │ Regridding │
                         └───────────┘       └─────┬─────┘
                                                   │
                                                   ▼
                                             ┌───────────┐
                                             │ Format to │
                                             │ 3D (T,Y,X)│
                                             └─────┬─────┘
                                                   │
                                                   ▼
                                             ┌───────────┐
                                             │ Export:   │
                                             │ R-Optimized
                                             │ NetCDF    │
                                             └───────────┘
```

---

## Detailed Technical Workflow

### 1. ICON-CH2 NWP Grid Processing (`T_2M`)

The gridded pipeline extracts surface temperature forecasts on an unstructured mesh, averages them temporally, and projects them into a standard coordinate system.

- **Earthkit Caching Control:** Configures `config.set("cache-policy", "temporary")` to manage high-volume GRIB/NetCDF downloads efficiently, enforcing automatic cleanup of temporary files to prevent local storage exhaustion.
- **Temporal Anchoring:** Explicitly builds a datetime object tracking the current day's **00:00 UTC** model run. This removes structural variance caused by executing the pipeline at different times of the day, ensuring a consistent forecast reference.
- **API Querying:** Generates ISO-8601 lead time strings (e.g., `P0DT0H`, `P1DT4H`) covering a 5-day forecast horizon (120 hourly steps) and submits the request to the `ogd-forecasting-icon-ch2` collection.
- **Dimensional Reduction:** Removes singleton dimensions (`eps`, representing ensemble members, and `z`, representing vertical levels) to isolate the native forecast arrays.
- **Native Spatial Aggregation:** Groups the 120 hourly forecasts into 24-hour blocks and computes the daily arithmetic mean directly on the native triangular mesh, minimizing interpolation bias.
- **Conservative Spatial Regridding:** Uses `regrid.iconremap` to project the unstructured grid onto a regular latitude/longitude grid (`EPSG:4326`) with dimensions **732 × 557** (approximately **2 km** spatial resolution). This conservative interpolation preserves atmospheric mass and reduces artifacts over complex Alpine terrain compared to standard bilinear interpolation.
- **3D Space-Time Reshaping for R:** Rearranges the data into a standard `(time, lat, lon)` layout. Latitude and longitude are stored as one-dimensional coordinate vectors, while time is encoded using CF-compliant metadata (`seconds since 1970-01-01 00:00:00`, standard calendar), ensuring compatibility with R packages such as `terra` and `raster`, as well as GDAL-based tools.

---

## Data Deliverables

Processed datasets are automatically written to the output directory using standardized filenames.

### Native Grid Output (currently not saved due to space issues)

#### `icon_ch2_native_daily_00utc_[YYYYMMDD].nc`

**Description**

Daily mean temperature fields aggregated directly on the native ICON-CH2 unstructured triangular grid.

**Dimensions**

```text
(ref_time: 1, lead_time: 5, cell: 283876)
```

---

### Regular Grid Output

#### `icon_ch2_2kmRegCons_daily_00utc_[YYYYMMDD].nc`

**Description**

Daily mean temperature fields conservatively remapped onto a regular 2 km latitude/longitude grid, optimized for analysis in R, GIS software, and other spatial processing workflows.

**Dimensions**

```text
(time: 5, lat: 557, lon: 732)
```