# Documentation of the forecast data

# MeteoSwiss Forecast Processing Pipeline (ICON-CH2 & Local Stations)

This repository contains a Python-based automation pipeline for retrieving, processing, and standardizing daily mean temperature forecasts ($2\text{m}$ air temperature) from MeteoSwiss. The script processes two complementary data streams: the **ICON-CH2 Numerical Weather Prediction (NWP) gridded model** and the **MOS/Local Stations point-based forecasts**. 

To ensure strict temporal consistency across both gridded and station data, all operations are pinned to the **00:00 UTC model initialization run** of the execution day.

---

## Data Pipeline Architecture

```
                               ┌──────────────────────────┐
                               │  MeteoSwiss Data Sources │
                               └─────────────┬────────────┘
                                             │
                    ┌────────────────────────┴────────────────────────┐
                    ▼                                                 ▼
       ┌─────────────────────────┐                       ┌─────────────────────────┐
       │   ICON-CH2 Grid (OGD)   │                       │   Local Stations (STAC) │
       └────────────┬────────────┘                       └────────────┬────────────┘
                    │                                                 │
                    ▼                                                 ▼
       ┌─────────────────────────┐                       ┌─────────────────────────┐
       │ Squeeze Dimensions      │                       │ Fetch 00 UTC Run Asset  │
       │ (Remove eps, z)         │                       └────────────┬────────────┘
       └────────────┬────────────┘                                    │
                    │                                                 ▼
                    ▼                                      ┌─────────────────────────┐
       ┌─────────────────────────┐                         │ Parse Semicolon CSV     │
       │ Spatial Grouping &      │                         │ (tre200h0 Parameter)    │
       │ 24h Block Mean          │                         └────────────┬────────────┘
       └────────────┬────────────┘                                    │
                    │                                                 ▼
         ┌──────────┴──────────┐                           ┌─────────────────────────┐
         ▼                     ▼                           │ Merge Station Metadata  │
   ┌───────────┐         ┌───────────┐                     │ (Lat, Lon, Elevation)   │
   │Export:    │         │ Apply     │                     └────────────┬────────────┘
   │Native Cell│         │Cons.      │                                │
   │NetCDF     │         │Regridding │                     ┌──────────┴──────────┐
   └───────────┘         └─────┬─────┘                     ▼                     ▼
                               │                     ┌───────────┐         ┌───────────┐
                               ▼                     │Export:    │         │ Filter QC │
                         ┌───────────┐               │Raw Hourly │         │ (24h count│
                         │Format to  │               │CSV        │         │ strictly) │
                         │3D (T,Y,X) │               └───────────┘         └─────┬─────┘
                         └─────┬─────┘                                           │
                               │                                                 ▼
                               ▼                                           ┌───────────┐
                         ┌───────────┐                                     │Export:    │
                         │Export: R- │                                     │Daily Mean │
                         │Optimized  │                                     │CSV        │
                         │NetCDF     │                                     └───────────┘
                         └───────────┘
```

---

## Detailed Technical Workflow

### 1. ICON-CH2 NWP Grid Processing (`T_2M`)
The gridded pipeline extracts surface temperature forecasts on an unstructured mesh, averages them temporally, and projects them into a standard coordinate system.

* **Earthkit Caching Control:** Configures `config.set("cache-policy", "temporary")` to manage high-volume GRIB/NetCDF downloads efficiently, enforcing automatic cleanup of temporary files to prevent local storage exhaustion.
* **Temporal Anchoring:** Explicitly builds a datetime object tracking today's `00:00 UTC` run. This removes structural variance caused by running the script at different times of the day, ensuring a consistent forecast base.
* **API Querying:** Generates systematic ISO-8601 lead time strings (`P0DT0H`, `P1DT4H`, etc.) extending across a 5-day horizon (120 hourly steps) and submits the data extraction request to the `ogd-forecasting-icon-ch2` collection.
* **Dimensional Reduction:** Squeezes out empty singleton dimensions (`eps` representing ensemble tracks, and `z` indicating vertical height layers) to isolate native data arrays.
* **Native Spatial Aggregation:** Groups the continuous 120-hour sequence into 24-hour steps via integer division. It computes the daily arithmetic mean directly on the native triangular mesh elements to minimize initial interpolation bias.
* **Conservative Spatial Regridding:** Utilizes `regrid.iconremap` to reproject the unstructured cells onto a regular lat/lon bounding box (`EPSG:4326`, $732 \times 557$ dimensions, $\sim 2\text{km}$ spatial resolution). Unlike generic bilinear methods, this conservative interpolation respects regional atmospheric boundaries, avoiding artificial thermal spikes across complex alpine topographies.
* **3D Spacetime Reshaping for R:** Rearranges the array axes into a pure 3D structure layout: `(time, lat, lon)`. Geographical dimensions are flattened to 1D coordinate vectors, and times are encoded using standard CF metadata specifications (`seconds since 1970-01-01 00:00:00` with a standard calendar). This ensures out-of-the-box compatibility with spatial packages in R (`terra`, `raster`) and GDAL tools.

### 2. Point-Based Local Station Forecast Processing (`tre200h0`)
The station-based pipeline fetches highly localized, statistically post-processed (MOS) meteorological station points through a web API.

* **Metadata Synchronization:** ingests the official global point identifier file (`ogd-local-forecasting_meta_point.csv`), parsing geographical definitions (Latitude, Longitude, Elevation) and mapping structural keys (`point_id`, `point_type_id`) for spatial accuracy.
* **STAC API Traversal:** Queries the MeteoSwiss SpatioTemporal Asset Catalog (STAC) endpoint to dynamically isolate today's specific catalog item (`YYYYMMDD-ch`). It loops across available asset keys to accurately locate and lock onto the `00:00 UTC` initialization run.
* **Parsing & Local Time Normalization:** Downloads the semicolon-delimited CSV for parameter `tre200h0` (Hourly mean of $2\text{m}$ air temperature). It converts raw UTC time strings into standard timestamp objects, then localizes them explicitly to Central European Time (`Europe/Zurich`) to maintain alignment with Swiss local civil calendars.
* **Completeness Quality Control (QC):** * **Hourly File:** Merges the time-series arrays directly with the geographical station attributes and saves a detailed hourly master table.
    * **Daily File:** Aggregates data by station and calendar date. To guarantee data integrity and eliminate edge-case distortions (such as partial observation blocks on the first or final forecast days), the script applies a strict look-back filter that drops any date group containing fewer than **24 valid hourly records**.

---

## Data Deliverables & Layout Specifications

Processed data products are continuously saved to the designated directory with standardized naming structures:

### Gridded Outputs (NetCDF Format)
1.  **`icon_ch2_native_daily_00utc_[YYYYMMDD].nc`**
    * **Description:** Daily aggregated mean fields preserved directly on the native unstructured triangular grid.
    * **Dimensions:** `(ref_time: 1, lead_time: 5, cell: 283876)`
2.  **`icon_ch2_2kmRegCons_daily_00utc_[YYYYMMDD].nc`**
    * **Description:** Reprojected regular grid using mass-conservative alpine remapping, explicitly structured for multi-temporal analysis in R or GIS environments.
    * **Dimensions:** `(time: 5, lat: 557, lon: 732)`

### Station Outputs (Semicolon CSV Format)
3.  **`MCH_local_forecast_[YYYYMMDD0000]_hourly.csv`**
    * **Description:** Semicolon-delimited file containing hourly raw localized forecasts for all Swiss meteorological stations merged with spatial attributes.
4.  **`MCH_local_forecast_[YYYYMMDD0000]_daily.csv`**
    * **Description:** Semicolon-delimited file detailing daily average station values, restricted strictly to complete 24-hour recording blocks.

