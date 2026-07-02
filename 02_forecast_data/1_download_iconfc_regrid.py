import datetime
from earthkit.data import config
from meteodatalab import ogd_api
from meteodatalab.operators import regrid
import numpy as np
import pandas as pd
from rasterio.crs import CRS
import xarray as xr

"""This script downloads and interpolates ICON MCH data using 
Meteodata-lab: a NumPy/Xarray-based Python library for processing and analyzing gridded meteorological data. 
See demos in: https://github.com/MeteoSwiss/opendata-nwp-demos"""

# ADAPT Path FOR YOUR DIRECTORY
path_out = "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Forecast_team_data/Daily_MCH_forecasts/"

# 1. Configure Earthkit Caching
config.set("cache-policy", "temporary")

# 2. Hardcode the 00:00 UTC model run for today. 
# This is an important step because if we download "latest" then depending on the 
# hour of download we will get forecasts initialized at different times.
today_00utc = datetime.datetime.now(datetime.timezone.utc).replace(
    hour=0, minute=0, second=0, microsecond=0
)

# This is the output filename
output_file = path_out+f"icon_ch2_2kmRegCons_daily_00utc_{today_00utc.strftime('%Y%m%d')}.nc"

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

# 9. Save the native to NetCDF on my own personal folder to have the native file as backup
# I omitted saving this file, as we do not need it for the webpage.

#path = "/Users/mp26i569/Documents/VS_CODE_projects/MeteoSwissData_annual_project/Data/"
#output_filename = path+f"icon_ch2_native_daily_00utc_{today_00utc.strftime('%Y%m%d')}.nc"
#ds_native.to_netcdf(output_filename)
#print(f"Successfully processed and saved to {output_filename}")


################################################################
# CREATE REGULAR FILE FROM THE INITIAL t2m hourly download

# 1. Define the regular lat/lon bounding box and resolution
# This domain fits the typical ICON-CH2 extent (approx. 2km resolution)
xmin, xmax, ymin, ymax = -0.817, 18.183, 41.183, 51.183
nx, ny = 732, 557  

destination = regrid.RegularGrid(CRS.from_epsg(4326), nx, ny, xmin, xmax, ymin, ymax)

#### 1. Apply regridding ####
"""regrid.iconremap: Dynamically scales or respects physical conservation laws between the triangular cells 
and your target square pixels. This prevents artificial cold or warm spikes in complex terrain like the Swiss Alps."""


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

# 7. Save this NetCDF to project folder 
t2m_r.to_netcdf(path_out+output_file, encoding=time_encoding)
ds_native.close()
print(f"Successfully processed and saved: {output_file}")