# ──────────────────────────────────────────────────────────────────────────────
# calculate_warning_level.R — forecast AF and assign heat-warning levels
# Annual Team Project 2026
#
# RUN AFTER get_MMT_baselineAF.R For each district and forecast day it computes the
# heat-attributable fraction (AF) from forecast temperatures, then classifies
# each district-day into a warning level (0–3).
#
# Inputs : forcast_temp_popw_district.csv
#          historical_temp_popw_mvngpop_2000_2024_district.csv  (defines the basis)
#          secondstage_blup.rds
#          mmt.rds                (from get_MMT_baselineAF.R)
#          baseline_summer.csv    (from get_MMT_baselineAF.R; baseline mode only)
# Output : warning_level.csv
# ──────────────────────────────────────────────────────────────────────────────

# ---- 1. Warning-level configuration ------------------------------------------
#   "fixed"    : the same AF cut-points for every district (fixed_cuts below).
#   "baseline" : district-specific percentiles from baseline_summer.csv
#                (columns dist_num, p75, p90, p95), produced by stage2_BLUP.R.
threshold_method <- "fixed"
fixed_cuts <- c(level1 = 0.075, level2 = 0.10, level3 = 0.15)


# ---- 2. Inputs ---------------------------------------------------------------
# Forecast temperatures (one row per district-day).
forecast_temp <- read.csv(f_forecast_temp) %>%
  select(-any_of("X")) %>%
  filter(BEZNR != exclude_districts)
forecast_temp_list <- split(forecast_temp, forecast_temp$BEZNR)

# Historical summer temperature — defines the spline support (knots/boundary)
# per district, so the basis matches the one used to fit the model.
temp_dist_sum <- read.csv(f_hist_temp) %>%
  select(-any_of("X")) %>%
  filter(month(time) %in% 5:9) %>%
  filter(BEZNR != exclude_districts)
temp_dist_list <- split(temp_dist_sum, temp_dist_sum$BEZNR)

# MMT per district, BLUP coefficients, and (optionally) baseline percentiles.
mmt_list <- readRDS(f_mmt)
res_blup <- readRDS(f_blup)
names(res_blup) <- names(forecast_temp_list)

baseline_summer <- if (file.exists(f_baseline)) read.csv(f_baseline) else NULL

# ---- 3. Forecast AF per district-day -----------------------------------------
af_fcst_list <- list()

# Calculate AF
for (k in seq(length(res_blup))) {
  
  dist_num <- names(res_blup)[k]
  
  # Extract temp data for selected district
  dist_temp <- temp_dist_list[[dist_num]]
  
  # Define argvar, centering and coef-vcov
  argvar <- list(fun = varfun,
                 knots = quantile(dist_temp$mean_value, varper / 100, na.rm = T), 
                 Bound = range(dist_temp$mean_value, na.rm = T))
  
  cen <- mmt_list$mintempcity[k]
  
  # EXTRACT PARAMETERS
  coef <- res_blup[[k]]$blup
  vcov <- res_blup[[k]]$vcov
  
  # Derive the centered basis
  bvar <- do.call(onebasis, c(list(x = forecast_temp_list[[dist_num]]$mean_value), argvar))
  cenvec <- do.call(onebasis, c(list(x = cen), argvar))
  bvarcen <- scale(bvar, center = cenvec, scale = F)
  
  
  af_dist_forecast <- (1 - exp(-bvarcen %*% coef))
  af_dsit_fcst_df <- cbind(forecast_temp_list[[dist_num]], af_dist_forecast)
  af_dsit_fcst_df$ind_heat_total <- af_dsit_fcst_df$mean_value > cen
  
  af_fcst_list[[dist_num]] <- af_dsit_fcst_df
  rm(af_dsit_fcst_df)
  
}

names(af_fcst_list) <- names(res_blup)


# ---- 4. Assign warning level -------------------------------------------------
assign_warning <- function(df, dist_id) {
  if (threshold_method == "baseline") {
    if (is.null(baseline_summer)) {
      stop("baseline_summer.csv not found — run stage2_BLUP.R or use threshold_method = 'fixed'.")
    }
    b  <- baseline_summer[as.character(baseline_summer$dist_num) == dist_id, ]
    c1 <- b$p75; c2 <- b$p90; c3 <- b$p95
  } else {
    c1 <- fixed_cuts[["level1"]]
    c2 <- fixed_cuts[["level2"]]
    c3 <- fixed_cuts[["level3"]]
  }

  # case_when evaluates top to bottom, so each line implies "and >= the previous
  # cut". NA AF values stay NA (there is no catch-all TRUE branch).
  df %>% mutate(warning_level = case_when(
    af_dist_forecast <  c1 ~ "0",
    af_dist_forecast <  c2 ~ "1",
    af_dist_forecast <  c3 ~ "2",
    af_dist_forecast >= c3 ~ "3"
  ))
}

# ---- 5. Combine & export -----------------------------------------------------
dist_fcst_warning <- Map(assign_warning, af_fcst_list, names(af_fcst_list)) %>%
  bind_rows() %>%
  select(BEZNAME, BEZNR, time, warning_level)

# row.names = FALSE avoids writing an extra index column (the "X" we drop on read).
write.csv(dist_fcst_warning, f_warning, row.names = FALSE)
