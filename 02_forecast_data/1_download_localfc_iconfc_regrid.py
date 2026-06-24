import datetime
from earthkit.data import config
from meteodatalab import ogd_api
from meteodatalab.operators import regrid
import numpy as np
import pandas as pd
from rasterio.crs import CRS
import xarray as xr

path_out = "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Forecast_team_data/Daily_MCH_forecasts/"
# 1. Configure Earthkit Caching
config.set("cache-policy", "temporary")

# 2. Hardcode the 00:00 UTC model run for today. 
# This is an important step because if you download "latest" then depending on the hour of download you will get forecasts initialized at different times.
today_00utc = datetime.datetime.now(datetime.timezone.utc).replace(
    hour=0, minute=0, second=0, microsecond=0
)

# The forecasts are hourly, but we need daily mean
# Generate lead times for 5 days (0 to 120 hours)
lead_times = [f"P0DT{h}H" if h < 24 else f"P{h//24}DT{h%24}H" for h in range(120)]


# 3. Request Data
req = ogd_api.Request(
    collection="ogd-forecasting-icon-ch2",
    variable="T_2M",
    ref_time=today_00utc,
    perturbed=False,
    lead_time=lead_times,
)
t2m = ogd_api.get_from_ogd(req)

# 4. Clean dimensions (remove eps, z)
t2m_cleaned = t2m.squeeze(dim=["eps", "z"], drop=True)
ref_time_val = t2m_cleaned["ref_time"].values

# 5. Group by 24-hour blocks and calculate the daily mean on native cells
hours = np.arange(120)
day_indices = hours // 24
day_indices[day_indices == 5] = 4


t2m_daily = t2m_cleaned.groupby(
    xr.DataArray(day_indices, dims="lead_time", name="day")
).mean(dim="lead_time")

# 6. Define coordinate tracking structure
target_lead_times = np.array([0,24,48,72,96], dtype="timedelta64[h]")
target_valid_times = ref_time_val + target_lead_times

# 7. Clean up attributes for NetCDF serialization compatibility
clean_attrs = {}
for k, v in t2m.attrs.items():
    if isinstance(v, (str, int, float, np.ndarray, list, tuple, bytes)):
        clean_attrs[k] = v
    else:
        clean_attrs[k] = str(v)

# 8. Reconstruct Final Dataset preserving native cells
# Reshape the data array to perfectly match the 3 specified dimensions (1, 5, 283876)
native_grid_array = t2m_daily.values.reshape(1, 5, -1)

ds_native = xr.Dataset(
    data_vars={"T_2M": (["ref_time", "lead_time", "cell"], native_grid_array)},
    coords={
        "ref_time": (["ref_time"], [ref_time_val.item()]),
        "lead_time": (["lead_time"], target_lead_times),
        "valid_time": (["ref_time", "lead_time"], target_valid_times[np.newaxis, :]),
        "lon": (["cell"], t2m_cleaned["lon"].values),
        "lat": (["cell"], t2m_cleaned["lat"].values),
    },
    attrs=clean_attrs,
)

# 9. Save the native to NetCDF on my own personal folder
path = "/Users/mp26i569/Documents/VS_CODE_projects/MeteoSwissData_annual_project/Data/"
output_filename = path+f"icon_ch2_native_daily_00utc_{today_00utc.strftime('%Y%m%d')}.nc"
ds_native.to_netcdf(output_filename)
print(f"Successfully processed and saved to {output_filename}")

# CREATES REGULAR FILE FROM THE INITIAL t2m hourly download

# 1. Define the regular lat/lon bounding box and resolution
# This domain fits the typical ICON-CH2 extent (approx. 2km resolution)
xmin, xmax, ymin, ymax = -0.817, 18.183, 41.183, 51.183
nx, ny = 732, 557  

destination = regrid.RegularGrid(CRS.from_epsg(4326), nx, ny, xmin, xmax, ymin, ymax)

#### 1. Apply regridding ####
"""regrid.iconremap: Dynamically scales or respects physical conservation laws between the triangular cells and your target square pixels. 
This prevents artificial cold or warm spikes in complex terrain like the Swiss Alps."""

output_file = path+f"icon_ch2_2kmRegCons_daily_00utc_{today_00utc.strftime('%Y%m%d')}.nc"

# 2. Slice the template down to match the dynamic lead_time steps of the current file
# Here I am basically taking the t2m template and replacing its values. I did this because I wanted 
# to match the metadata needed for the regrid function to work 
# regrid.iconremap(ds_native, destination) was not working
# but it was working when doing regrid.iconremap(t2m, destination)
da_matched = t2m.isel(lead_time=slice(0, ds_native.sizes["lead_time"])).copy(deep=False)

# 3. Assign coordinates explicitly mapping the 2D matrix layout for valid_time
da_matched = da_matched.assign_coords(
        ref_time=pd.to_datetime(ds_native["ref_time"].values),
        lead_time=ds_native["lead_time"].values,
        valid_time=(("ref_time", "lead_time"), ds_native["valid_time"].values)
    )

# 4. Inject the core data values
da_matched.values = ds_native["T_2M"].values[np.newaxis, :, :, np.newaxis, :]

# 5. Perform the spatial interpolation
t2m_regular = regrid.iconremap(da_matched, destination)

# Clean up ALL non-serializable dictionary or object attributes
def sanitize_attrs(obj):
    valid_types = (str, int, float, list, tuple, np.ndarray)
    bad_keys = [k for k, v in obj.attrs.items() if not isinstance(v, valid_types)]
    for k in bad_keys:
        obj.attrs[k] = str(obj.attrs[k])

# Sanitize main array attributes
sanitize_attrs(t2m_regular)

# Sanitize coordinate attributes (e.g., latitude, longitude, time)
for coord in t2m_regular.coords.values():
    sanitize_attrs(coord)

# ==================================================================
# NEW: Format to 3D Spacetime with (time, lat, lon) dimensions for R
# ==================================================================
# Extract 1D time values from the native dataset
ref_time_arr = pd.to_datetime(ds_native["ref_time"].values).to_numpy()
lead_time_vals = ds_native["lead_time"].values
valid_time_matrix = ref_time_arr[:, None] + lead_time_vals[None, :]
time_1d = valid_time_matrix.ravel()

# Squeeze out the singleton dimensions (eps, ref_time, z) from the numpy array
# Shape converts from (1, 1, 5, 1, 557, 732) -> (5, 557, 732)
raw_data = np.squeeze(t2m_regular.values, axis=(0, 1, 3))

# Extract the 1D geographic coordinate arrays generated by regrid
# (Since it's a regular grid, 2D lat/lon matrices are not needed as dimensions)
lat_1d = t2m_regular.coords["lat"].values[:, 0]  # Take a column slice across the grid rows
lon_1d = t2m_regular.coords["lon"].values[0, :]  # Take a row slice across the grid columns

# Build a clean, unencumbered DataArray using time, lat, and lon as explicit dimensions
t2m_r = xr.DataArray(
        data=raw_data,
        dims=["time", "lat", "lon"],
        coords={
            "time": time_1d,
            "lat": lat_1d,
            "lon": lon_1d,
        },
        attrs=t2m_regular.attrs
    )

# Explicitly enforce R spatial layout order
t2m_r = t2m_r.transpose("time", "lat", "lon")

# Time encoding properties required for continuous tracking in R/GDAL
time_encoding = {
        "time": {
            "units": "seconds since 1970-01-01 00:00:00",
            "calendar": "standard",
            "dtype": "int64",
        },
    }

# 7. Save this NetCDF to project folder and close file handles cleanly
t2m_r.to_netcdf(path_out+output_file, encoding=time_encoding)
ds_native.close()
    
print(f"Successfully processed and saved: {output_file}")


################################################################################

#Download the LOCAL stations MeteoSwiss forecast 

################################################################################

from datetime import datetime
from io import StringIO
from zoneinfo import ZoneInfo
import httpx
import pandas as pd

# ---------------------------------------------------------------------------
# 0. Setup Configuration & Metadata
# ---------------------------------------------------------------------------
STAC_BASE_URL = "https://data.geo.admin.ch/api/stac/v1"
COLLECTION_ID = "ch.meteoschweiz.ogd-local-forecasting"
POI_LIST_URL = f"https://data.geo.admin.ch/{COLLECTION_ID}/ogd-local-forecasting_meta_point.csv"
LOCAL_TZ = ZoneInfo("Europe/Zurich")

TARGET_PARAM = "tre200h0"  # Air temperature 2 m above ground; hourly mean
path = "/Users/mp26i569/Documents/VS_CODE_projects/MeteoSwissData_annual_project/Data/"

# Load and rename metadata for all stations globally
df_pois = pd.read_csv(POI_LIST_URL, sep=";", encoding="latin-1")
df_meta = df_pois[[
    "point_id", "point_type_id", "point_name", 
    "point_coordinates_wgs84_lat", "point_coordinates_wgs84_lon", "point_height_masl"
]].rename(columns={
    "point_name": "station_name",
    "point_coordinates_wgs84_lat": "latitude",
    "point_coordinates_wgs84_lon": "longitude",
    "point_height_masl": "elevation"
})

# ---------------------------------------------------------------------------
# Step 1 — Fetch today's STAC item and grab the 00 UTC run
# ---------------------------------------------------------------------------
today_id = f"{datetime.now(LOCAL_TZ).strftime('%Y%m%d')}-ch"
item_url = f"{STAC_BASE_URL}/collections/{COLLECTION_ID}/items/{today_id}"

with httpx.Client() as client:
    item = client.get(item_url)
    item.raise_for_status()
    stac_item = item.json()

assets = stac_item["assets"]

all_runs = sorted({key.split(".")[2] for key in assets if TARGET_PARAM in key})
if not all_runs:
    raise ValueError(f"No asset runs found for parameter: {TARGET_PARAM}")

utc00_run = all_runs[0]

# Verify run corresponds to 00 UTC
parsed_dt = datetime.strptime(utc00_run, "%Y%m%d%H%M")
if parsed_dt.hour != 0:
    raise ValueError(f"Invalid run selected: {utc00_run}. This is not a 00 UTC initialization.")

# Map the targeted parameter to its asset URL for the 00 UTC run
param_urls = {}
match = next((key for key in assets if TARGET_PARAM in key and utc00_run in key), None)
if match:
    param_urls[TARGET_PARAM] = assets[match]["href"]
else:
    raise ValueError(f"Could not locate asset URL matching {TARGET_PARAM} for run {utc00_run}")

# ---------------------------------------------------------------------------
# Step 2 — Download parameter CSV
# ---------------------------------------------------------------------------
raw_data = {}
with httpx.Client(timeout=30.0) as client:
    for param, url in param_urls.items():
        resp = client.get(url)
        if resp.status_code == 200:
            raw_data[param] = resp.content

# ---------------------------------------------------------------------------
# Step 3 — Parse CSV for all stations and merge metadata
# ---------------------------------------------------------------------------
def parse_parameter_csv(content: bytes, param: str, tz) -> pd.DataFrame:
    df = pd.read_csv(StringIO(content.decode("latin-1")), sep=";")
    if df.empty:
        return pd.DataFrame()

    time_col = next(c for c in df.columns if any(k in c.lower() for k in ("date", "time")))
    value_col = next(c for c in df.columns if c not in (time_col, "point_id", "point_type_id") and "type" not in c.lower())

    df[time_col] = pd.to_datetime(df[time_col].astype(int).astype(str), format="%Y%m%d%H%M", utc=True).dt.tz_convert(tz)
    df = df.rename(columns={time_col: "time", value_col: param})
    df[param] = pd.to_numeric(df[param], errors="coerce")
    
    return df[["point_id", "point_type_id", "time", param]]

if TARGET_PARAM in raw_data:
    df_parsed = parse_parameter_csv(raw_data[TARGET_PARAM], TARGET_PARAM, LOCAL_TZ)
    if not df_parsed.empty:
        # Merge the full forecast data with the structural station metadata
        df_all_stations = pd.merge(df_parsed, df_meta, on=["point_id", "point_type_id"], how="inner")
        df_all_stations = df_all_stations.sort_values(by=["station_name", "time"]).reset_index(drop=True)
        
        print(df_all_stations.head())
        
        # Save the hourly forecast dataset
        hourly_filename = f"MCH_local_forecast_{utc00_run}_hourly.csv"
        df_all_stations.to_csv(path+hourly_filename, index=False, sep=";")
        print(f"✓ Hourly file successfully saved as: {hourly_filename}")
        
        # ---------------------------------------------------------------------------
        # Step 4 — Calculate daily mean only for complete 24h days
        # ---------------------------------------------------------------------------
        df_daily = df_all_stations.copy()
        
        # Extract the calendar date directly from the local time stamp
        df_daily["date"] = df_daily["time"].dt.date
        
        # Group and aggregate both the mean value and the number of hourly observations
        group_cols = ["point_id", "point_type_id", "station_name", "latitude", "longitude", "elevation", "date"]
        df_daily_agg = df_daily.groupby(group_cols, as_index=False).agg(
            mean_value=(TARGET_PARAM, "mean"),
            hour_count=(TARGET_PARAM, "count")
        )
        
        # Filter groups to strictly retain days that have exactly 24 hours (00 to 23)
        df_daily_mean = df_daily_agg[df_daily_agg["hour_count"] == 24].copy()
        df_daily_mean = df_daily_mean.rename(columns={"mean_value": TARGET_PARAM}).drop(columns=["hour_count"])
        
        daily_filename = f"MCH_local_forecast_{utc00_run}_daily.csv"
        df_daily_mean.to_csv(path+daily_filename, index=False, sep=";")
        print(f"✓ Daily mean file successfully saved as: {daily_filename}")
