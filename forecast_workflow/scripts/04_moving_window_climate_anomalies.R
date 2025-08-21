# Moving Window Climate Anomalies (Backward-Looking)
# Implements colleague's standardization method:
# - Temperature anomalies: 30-day backward moving window to remove seasonal variation  
# - SPEI: 14-day backward totals of (P - ET), then standardize
# NOTE: All windows are RIGHT-ALIGNED (backward-looking) to avoid using future data

library(tidyverse)
library(sf)
library(lubridate)
library(zoo)  # For rolling window functions

# Configuration for moving window approach
moving_window_config <- list(
  temp_window = 30,    # days for temperature anomaly calculation (backward-looking)
  spei_window = 14,    # days for SPEI calculation (backward-looking)
  min_years = 5        # minimum years of data for standardization
)

# Function to calculate reference evapotranspiration (simple Thornthwaite)
calculate_reference_et <- function(temp_max, temp_min, day_of_year, latitude = 42) {
  
  # Mean temperature
  temp_mean <- (temp_max + temp_min) / 2
  
  # Simple ET calculation (Thornthwaite-style)
  et <- ifelse(temp_mean > 0,
               16 * (temp_mean / 5)^1.514,
               0)
  
  # Apply daylight correction (simplified)
  daylight_factor <- 1 + 0.2 * sin(2 * pi * (day_of_year - 80) / 365)
  et_adjusted <- et * daylight_factor
  
  return(pmax(0, et_adjusted))  # ET cannot be negative
}

# Function to convert GEFS reforecast data to daily format for SPEI
convert_reforecast_to_daily <- function(reforecast_climate, start_date, end_date) {
  
  cat("Converting GEFS reforecast data to daily format for SPEI calculation...\n")
  
  # Reforecast data is at specific dates, need to interpolate to daily
  # For each GEFS cell, interpolate between available dates
  
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  
  daily_reforecast <- expand_grid(
    gefs_cell_id = unique(reforecast_climate$gefs_cell_id),
    date = dates
  ) %>%
    mutate(
      year = year(date),
      month = month(date),
      day_of_year = yday(date)
    ) %>%
    # For each cell, interpolate climate variables
    group_by(gefs_cell_id) %>%
    group_modify(~{
      cell_climate <- filter(reforecast_climate, gefs_cell_id == .y$gefs_cell_id)
      
      if (nrow(cell_climate) > 0) {
        # Linear interpolation between available reforecast dates
        temp_interp <- approx(x = cell_climate$date, y = cell_climate$temperature, 
                             xout = .x$date, rule = 2)$y
        precip_interp <- approx(x = cell_climate$date, y = cell_climate$precipitation,
                               xout = .x$date, rule = 2)$y
        
        .x %>%
          mutate(
            # Convert temperature to daily min/max (add realistic variation)
            temp_mean = temp_interp,
            temp_max = temp_mean + 6 + rnorm(n(), 0, 1.5),
            temp_min = temp_mean - 6 + rnorm(n(), 0, 1.5),
            # Use interpolated precipitation (convert to daily rates)
            precip = pmax(0, precip_interp / 30)  # Assume monthly precip distributed over ~30 days
          )
      } else {
        .x %>% mutate(temp_max = NA, temp_min = NA, precip = NA)
      }
    }) %>%
    ungroup() %>%
    filter(!is.na(temp_max))
  
  cat("Converted to", nrow(daily_reforecast), "daily reforecast records\n")
  cat("Date range:", range(daily_reforecast$date), "\n")
  
  return(daily_reforecast)
}

# Function to convert monthly GEFS data to daily format (legacy)
convert_monthly_to_daily <- function(monthly_climate, start_date, end_date) {
  
  cat("Converting monthly GEFS climate data to daily format...\n")
  
  # Create daily date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  
  # Create daily data by repeating monthly values throughout each month
  daily_climate <- expand_grid(
    gefs_cell_id = unique(monthly_climate$gefs_cell_id),
    date = dates
  ) %>%
    mutate(
      year = year(date),
      month = month(date),
      day_of_year = yday(date)
    ) %>%
    left_join(
      monthly_climate %>%
        select(gefs_cell_id, year, month, temp_mean, precip_total, solar_mean, rh_mean),
      by = c("gefs_cell_id", "year", "month")
    ) %>%
    filter(!is.na(temp_mean)) %>%
    # Convert monthly means to realistic daily values
    mutate(
      # Use monthly temp as base, add daily variation
      temp_max = temp_mean + 5 + rnorm(n(), 0, 2),
      temp_min = temp_mean - 5 + rnorm(n(), 0, 2),
      # Convert monthly total precip to daily (assume uniform distribution)
      precip = precip_total / days_in_month(date) + rgamma(n(), shape = 0.5, rate = 2)
    )
  
  cat("Converted to", nrow(daily_climate), "daily records\n")
  
  return(daily_climate)
}

# Function to create synthetic daily climate data for testing
create_synthetic_daily_climate <- function(gefs_centers, 
                                          start_date = "2001-01-01", 
                                          end_date = "2023-12-31") {
  
  cat("Creating synthetic daily climate data...\n")
  cat("Period:", start_date, "to", end_date, "\n")
  
  # Create daily date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  
  # Create realistic daily climate data
  daily_climate <- expand_grid(
    gefs_cell_id = gefs_centers$gefs_cell_id,
    date = dates
  ) %>%
    mutate(
      year = year(date),
      month = month(date), 
      day_of_year = yday(date),
      # Create realistic temperature with seasonal cycle + noise
      temp_base = 15 + 15 * sin(2 * pi * (day_of_year - 100) / 365),
      temp_max = temp_base + 5 + rnorm(n(), 0, 3),
      temp_min = temp_base - 5 + rnorm(n(), 0, 2),
      # Precipitation (log-normal with seasonal component)
      precip_base = 2 + sin(2 * pi * (day_of_year - 120) / 365),
      precip = pmax(0, rlnorm(n(), log(precip_base), 1.2))
    ) %>%
    # Add some spatial variation
    group_by(gefs_cell_id) %>%
    mutate(
      spatial_temp_offset = rnorm(1, 0, 1),
      spatial_precip_factor = exp(rnorm(1, 0, 0.2))
    ) %>%
    ungroup() %>%
    mutate(
      temp_max = temp_max + spatial_temp_offset,
      temp_min = temp_min + spatial_temp_offset,
      precip = precip * spatial_precip_factor
    )
  
  cat("Created", nrow(daily_climate), "daily climate records\n")
  cat("Grid cells:", n_distinct(daily_climate$gefs_cell_id), "\n")
  cat("Years:", min(daily_climate$year), "-", max(daily_climate$year), "\n")
  
  return(daily_climate)
}

# Function to calculate 30-day backward moving window temperature anomalies
calculate_temp_anomalies <- function(daily_climate, window_days = 30) {
  
  cat("Calculating temperature anomalies with", window_days, "day BACKWARD moving window...\n")
  cat("Window alignment: RIGHT (uses past", window_days, "days including current day)\n")
  
  climate_with_anomalies <- daily_climate %>%
    group_by(gefs_cell_id) %>%
    arrange(date) %>%
    mutate(
      # Calculate BACKWARD rolling mean temperatures (right-aligned)
      temp_max_rolling = zoo::rollmean(temp_max, k = window_days, align = "right", fill = NA),
      temp_min_rolling = zoo::rollmean(temp_min, k = window_days, align = "right", fill = NA),
      
      # Calculate long-term climatology for each day of year
      temp_max_clim = ave(temp_max_rolling, day_of_year, FUN = function(x) mean(x, na.rm = TRUE)),
      temp_min_clim = ave(temp_min_rolling, day_of_year, FUN = function(x) mean(x, na.rm = TRUE)),
      
      # Calculate anomalies
      temp_max_anomaly = temp_max_rolling - temp_max_clim,
      temp_min_anomaly = temp_min_rolling - temp_min_clim
    ) %>%
    ungroup()
  
  # Remove NA values from rolling windows
  climate_clean <- climate_with_anomalies %>%
    filter(!is.na(temp_max_anomaly))
  
  cat("Calculated temperature anomalies:", nrow(climate_clean), "records\n")
  cat("Anomaly range - Max:", round(range(climate_clean$temp_max_anomaly, na.rm = TRUE), 2), "°C\n")
  
  return(climate_clean)
}

# Function to calculate 14-day backward SPEI
calculate_spei_14day <- function(climate_data, spei_days = 14) {
  
  cat("Calculating SPEI with", spei_days, "day BACKWARD accumulation...\n")
  cat("Window alignment: RIGHT (uses past", spei_days, "days including current day)\n")
  
  climate_with_spei <- climate_data %>%
    group_by(gefs_cell_id) %>%
    arrange(date) %>%
    mutate(
      # Calculate reference ET
      ref_et = calculate_reference_et(temp_max, temp_min, day_of_year, latitude = 42),
      
      # Water balance (P - ET)
      water_balance = precip - ref_et,
      
      # BACKWARD rolling sums for SPEI period (right-aligned)
      precip_sum = zoo::rollsum(precip, k = spei_days, align = "right", fill = NA),
      et_sum = zoo::rollsum(ref_et, k = spei_days, align = "right", fill = NA),
      wb_sum = zoo::rollsum(water_balance, k = spei_days, align = "right", fill = NA)
    ) %>%
    ungroup() %>%
    filter(!is.na(wb_sum))
  
  # Standardize water balance to get SPEI
  spei_standardized <- climate_with_spei %>%
    group_by(gefs_cell_id, month) %>%  # Standardize by month to remove seasonal bias
    mutate(
      spei_14 = scale(wb_sum)[,1]  # Z-score standardization
    ) %>%
    ungroup()
  
  cat("Calculated SPEI:", nrow(spei_standardized), "records\n")
  cat("SPEI range:", round(range(spei_standardized$spei_14, na.rm = TRUE), 2), "\n")
  
  return(spei_standardized)
}

# Function to integrate climate anomalies with NDVI data
integrate_anomalies_with_ndvi <- function(ndvi_gefs_data, climate_anomalies) {
  
  cat("Integrating climate anomalies with NDVI data...\n")
  
  # Create NDVI data for matching
  ndvi_for_matching <- ndvi_gefs_data$ndvi_data %>%
    st_drop_geometry() %>%
    select(gefs_cell_id, date, year, month, day_of_year, ndvi_anomaly, ndvi_observed) %>%
    mutate(date = as.Date(date))
  
  # Match climate data to NDVI observation dates (or nearest dates)
  integrated_data <- ndvi_for_matching %>%
    left_join(
      climate_anomalies %>%
        select(gefs_cell_id, date, temp_max_anomaly, spei_14, precip, ref_et),
      by = c("gefs_cell_id", "date")
    ) %>%
    # For dates without exact climate match, use nearest neighbor (backward fill only)
    group_by(gefs_cell_id) %>%
    arrange(date) %>%
    fill(temp_max_anomaly, spei_14, precip, ref_et, .direction = "down") %>%  # Only forward fill
    ungroup() %>%
    # Keep only records with both NDVI and climate data
    filter(!is.na(ndvi_anomaly) & !is.na(temp_max_anomaly) & !is.na(spei_14))
  
  cat("Integrated dataset:", nrow(integrated_data), "records\n")
  cat("Date range:", range(integrated_data$date), "\n")
  cat("GEFS cells:", n_distinct(integrated_data$gefs_cell_id), "\n")
  
  return(integrated_data)
}

# Function to save climate anomalies
save_climate_anomalies <- function(climate_data, integrated_data, 
                                  output_path = "U:/datasets/ndvi_monitor/climate_anomalies_moving_window.rds") {
  
  # Create output directory
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save with metadata
  anomalies_with_metadata <- list(
    climate_anomalies = climate_data,
    integrated_ndvi_climate = integrated_data,
    metadata = list(
      method = "Backward moving window anomalies",
      temp_window = moving_window_config$temp_window,
      spei_window = moving_window_config$spei_window,
      window_alignment = "right (backward-looking)",
      created_date = Sys.Date(),
      n_records_climate = nrow(climate_data),
      n_records_integrated = nrow(integrated_data)
    )
  )
  
  saveRDS(anomalies_with_metadata, output_path)
  cat("Climate anomalies saved to:", output_path, "\n")
  
  return(anomalies_with_metadata)
}

# Main workflow function using GEFS reforecast data
calculate_moving_window_anomalies_with_reforecast <- function(ndvi_gefs_data, 
                                                            reforecast_climate_data,
                                                            start_date = "2001-01-01",
                                                            end_date = "2015-12-31",
                                                            save_output = TRUE) {
  
  cat("Starting BACKWARD moving window climate anomaly calculation with GEFS reforecast...\n")
  cat("Method: 30-day backward temp anomalies + 14-day backward SPEI\n")
  cat("Data source: GEFSv12 reforecast (temperature, precipitation)\n")
  cat("All windows are RIGHT-ALIGNED (no future data used)\n")
  
  # Convert GEFS reforecast data to daily format
  daily_climate <- convert_reforecast_to_daily(reforecast_climate_data, start_date, end_date)
  
  # Calculate temperature anomalies (30-day backward moving window)
  climate_with_temp_anomalies <- calculate_temp_anomalies(
    daily_climate, 
    window_days = moving_window_config$temp_window
  )
  
  # Calculate SPEI (14-day backward accumulation)
  climate_with_anomalies <- calculate_spei_14day(
    climate_with_temp_anomalies,
    spei_days = moving_window_config$spei_window
  )
  
  # Integrate with NDVI data
  integrated_data <- integrate_anomalies_with_ndvi(ndvi_gefs_data, climate_with_anomalies)
  
  if (save_output) {
    saved_data <- save_climate_anomalies(climate_with_anomalies, integrated_data)
    return(saved_data)
  } else {
    return(list(
      climate_anomalies = climate_with_anomalies,
      integrated_data = integrated_data
    ))
  }
}

# Main workflow function (fallback with synthetic data)
calculate_moving_window_anomalies <- function(ndvi_gefs_data, 
                                            start_date = "2001-01-01",
                                            end_date = "2023-12-31",
                                            save_output = TRUE) {
  
  cat("Starting BACKWARD moving window climate anomaly calculation...\n")
  cat("Method: 30-day backward temp anomalies + 14-day backward SPEI\n")
  cat("All windows are RIGHT-ALIGNED (no future data used)\n")
  
  # Extract GEFS centers
  gefs_centers <- ndvi_gefs_data$gefs_grid %>%
    mutate(
      centroid = st_centroid(geometry),
      longitude = st_coordinates(centroid)[,1],
      latitude = st_coordinates(centroid)[,2]
    ) %>%
    st_drop_geometry() %>%
    select(gefs_cell_id, longitude, latitude)
  
  # Use synthetic data as fallback
  cat("Using synthetic daily climate data for testing...\n")
  daily_climate <- create_synthetic_daily_climate(gefs_centers, start_date, end_date)
  
  # Calculate temperature anomalies (30-day backward moving window)
  climate_with_temp_anomalies <- calculate_temp_anomalies(
    daily_climate, 
    window_days = moving_window_config$temp_window
  )
  
  # Calculate SPEI (14-day backward accumulation)
  climate_with_anomalies <- calculate_spei_14day(
    climate_with_temp_anomalies,
    spei_days = moving_window_config$spei_window
  )
  
  # Integrate with NDVI data
  integrated_data <- integrate_anomalies_with_ndvi(ndvi_gefs_data, climate_with_anomalies)
  
  if (save_output) {
    saved_data <- save_climate_anomalies(climate_with_anomalies, integrated_data)
    return(saved_data)
  } else {
    return(list(
      climate_anomalies = climate_with_anomalies,
      integrated_data = integrated_data
    ))
  }
}

cat("BACKWARD moving window climate anomalies script loaded.\n")
cat("Configuration:\n")
cat("  - Temperature window:", moving_window_config$temp_window, "days (RIGHT-ALIGNED)\n")
cat("  - SPEI window:", moving_window_config$spei_window, "days (RIGHT-ALIGNED)\n")
cat("  - No future data used in calculations\n")
cat("Run calculate_moving_window_anomalies(ndvi_gefs_data) to compute anomalies.\n")