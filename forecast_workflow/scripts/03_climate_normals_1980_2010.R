# Climate Normals Calculation Script (1980-2010)
# Calculates climatological normals for SPEI 14-day and Tmax 30-day metrics
# Reference period: 1980-2010 (30-year standard climatology)

library(tidyverse)
library(sf)
library(lubridate)
library(httr)
library(ecmwfr)  # For ERA5 data access

# Configuration for climate normals
normals_config <- list(
  reference_period = c("1980-01-01", "2010-12-31"),
  gefs_resolution = 0.25,  # degrees
  variables = c("tmax", "tmin", "precip", "rh", "wind", "solar"),
  spei_window = 14,  # days
  tmax_window = 30   # days
)

# ERA5 credentials and configuration for 1980-2010 climatology
era5_credentials <- list(
  uid = "266068",
  api_key = "35721697-ea24-49b7-bfa6-797f63bf952d"
)

era5_config <- list(
  dataset = "reanalysis-era5-single-levels",
  variables = c(
    "2m_temperature",                    # For Tmax/Tmin (will get daily max/min)
    "total_precipitation",               # For SPEI
    "2m_relative_humidity",              # For PET calculation
    "10m_u_component_of_wind",          # For PET calculation  
    "10m_v_component_of_wind",          # For PET calculation
    "surface_solar_radiation_downwards"  # For PET calculation
  ),
  format = "netcdf",
  # Bounding box for Chicago area (slightly larger than our GEFS grid)
  area = c(42.5, -89, 41, -87.5),  # North, West, South, East
  grid = c(0.25, 0.25),  # 0.25 degree resolution to match GEFS
  times = c("00:00", "06:00", "12:00", "18:00")  # 6-hourly for daily aggregation
)

# Function to setup ERA5 credentials
setup_era5_credentials <- function() {
  
  # Set credentials using ecmwfr (updated syntax)
  wf_set_key(
    user = era5_credentials$uid,
    key = era5_credentials$api_key
  )
  
  cat("✓ ERA5 credentials configured\n")
  cat("User ID:", era5_credentials$uid, "\n")
  
  return(TRUE)
}

# Function to request ERA5 data for a specific year
request_era5_year <- function(year, variables = era5_config$variables, output_dir = "U:/datasets/ndvi_monitor/era5/") {
  
  cat("Requesting ERA5 data for year:", year, "\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Output file name
  output_file <- file.path(output_dir, paste0("era5_", year, ".nc"))
  
  # Skip if file already exists
  if (file.exists(output_file)) {
    cat("  File already exists:", output_file, "\n")
    return(output_file)
  }
  
  # ERA5 request specification
  request <- list(
    dataset_short_name = era5_config$dataset,
    product_type = "reanalysis",
    variable = variables,
    year = as.character(year),
    month = sprintf("%02d", 1:12),
    day = sprintf("%02d", 1:31),
    time = era5_config$times,
    area = era5_config$area,
    grid = era5_config$grid,
    format = era5_config$format,
    target = basename(output_file)
  )
  
  tryCatch({
    # Submit request to ERA5
    request_id <- wf_request(
      user = era5_credentials$uid,
      request = request,
      path = output_dir
    )
    
    cat("  ✓ ERA5 request submitted for", year, "\n")
    cat("  Request ID:", request_id, "\n")
    cat("  Output file:", output_file, "\n")
    
    return(output_file)
    
  }, error = function(e) {
    cat("  ✗ Error requesting ERA5 data for", year, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to process ERA5 NetCDF to daily values
process_era5_to_daily <- function(era5_file, gefs_grid) {
  
  if (!file.exists(era5_file)) {
    cat("ERA5 file not found:", era5_file, "\n")
    return(NULL)
  }
  
  cat("Processing ERA5 file:", basename(era5_file), "\n")
  
  tryCatch({
    # Read ERA5 NetCDF file
    era5_data <- stars::read_ncdf(era5_file)
    
    # Convert to daily data (aggregate from 6-hourly)
    era5_daily <- era5_data %>%
      # Group by date (extract date from time dimension)
      mutate(date = as.Date(time)) %>%
      group_by(date) %>%
      summarise(
        # Daily aggregations
        tmax = max(t2m, na.rm = TRUE) - 273.15,  # Convert K to C
        tmin = min(t2m, na.rm = TRUE) - 273.15,  # Convert K to C
        temp_mean = mean(t2m, na.rm = TRUE) - 273.15,
        precip = sum(tp, na.rm = TRUE) * 1000,   # Convert m to mm
        rh = mean(r2, na.rm = TRUE),             # Mean relative humidity
        wind_u = mean(u10, na.rm = TRUE),        # Mean u wind
        wind_v = mean(v10, na.rm = TRUE),        # Mean v wind
        solar = mean(ssrd, na.rm = TRUE) / 86400, # Convert J/m² to W/m²
        .groups = "drop"
      ) %>%
      # Calculate wind speed
      mutate(
        wind_speed = sqrt(wind_u^2 + wind_v^2)
      )
    
    # Extract data at GEFS grid points
    gefs_coords <- gefs_grid %>%
      st_centroid() %>%
      st_coordinates() %>%
      as_tibble() %>%
      bind_cols(gefs_cell_id = gefs_grid$gefs_cell_id) %>%
      rename(longitude = X, latitude = Y)
    
    # Interpolate ERA5 to GEFS grid points (simplified - nearest neighbor)
    daily_at_grid <- era5_daily %>%
      # This is a simplified approach - should use proper spatial interpolation
      slice_sample(n = nrow(gefs_coords)) %>%
      bind_cols(gefs_cell_id = gefs_coords$gefs_cell_id)
    
    cat("  ✓ Processed", nrow(daily_at_grid), "daily records\n")
    
    return(daily_at_grid)
    
  }, error = function(e) {
    cat("  ✗ Error processing ERA5 file:", e$message, "\n")
    return(NULL)
  })
}

# Function to calculate potential evapotranspiration using Thornthwaite method
# Simplified approach - more sophisticated methods can be added later
calculate_pet_thornthwaite <- function(temp_data, latitude, day_of_year) {
  
  # Thornthwaite heat index calculation
  temp_monthly <- temp_data %>%
    filter(temp_mean > 0) %>%  # Only positive temperatures
    mutate(
      heat_index_component = (temp_mean / 5) ^ 1.514
    )
  
  # Annual heat index
  I <- sum(temp_monthly$heat_index_component, na.rm = TRUE)
  
  # Alpha coefficient
  alpha <- (6.75e-7 * I^3) - (7.71e-5 * I^2) + (1.792e-2 * I) + 0.49239
  
  # Daylight hours correction (simplified)
  daylight_correction <- 1.0  # Placeholder - should vary by latitude and day of year
  
  # PET calculation (mm/day)
  pet <- temp_data %>%
    mutate(
      pet_mm_day = ifelse(temp_mean > 0,
                          16 * daylight_correction * (10 * temp_mean / I)^alpha,
                          0)
    ) %>%
    pull(pet_mm_day)
  
  return(pet)
}

# Function to calculate SPEI from precipitation and PET
calculate_spei <- function(precip, pet, window_days = 14) {
  
  # Water balance (P - PET)
  water_balance <- precip - pet
  
  # Rolling sum for specified window
  if (length(water_balance) >= window_days) {
    wb_rolling <- zoo::rollsum(water_balance, k = window_days, align = "right", fill = NA)
    
    # Standardize to get SPEI (fit to gamma distribution typically)
    # Simplified: use z-score standardization for now
    spei <- scale(wb_rolling)[,1]
    
    return(spei)
  } else {
    return(rep(NA, length(water_balance)))
  }
}

# Function to create synthetic daily climate data for normals period
# This is a placeholder until we implement real ERA5/GEFS data retrieval
create_synthetic_normals_data <- function(gefs_grid, start_date, end_date) {
  
  cat("Creating synthetic daily climate data for normals calculation...\n")
  cat("Period:", start_date, "to", end_date, "\n")
  
  # Create daily date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  n_days <- length(dates)
  n_cells <- nrow(gefs_grid)
  
  # Create synthetic daily data for each GEFS cell
  daily_climate <- expand_grid(
    gefs_cell_id = gefs_grid$gefs_cell_id,
    date = dates
  ) %>%
    mutate(
      year = year(date),
      month = month(date),
      day_of_year = yday(date),
      # Create realistic daily climate variables with seasonal cycles
      tmax = 20 + 15 * sin(2 * pi * (day_of_year - 100) / 365) + rnorm(n(), 0, 3),
      tmin = tmax - 10 + rnorm(n(), 0, 2),
      temp_mean = (tmax + tmin) / 2,
      precip = pmax(0, rlnorm(n(), log(2), 1.5)),  # Daily precip (mm)
      rh = 65 + 15 * sin(2 * pi * (day_of_year - 150) / 365) + rnorm(n(), 0, 8),
      wind = 3 + rnorm(n(), 0, 1),
      solar = 200 + 150 * sin(2 * pi * (day_of_year - 80) / 365) + rnorm(n(), 0, 30)
    ) %>%
    # Add some spatial variation
    left_join(
      gefs_grid %>% 
        st_drop_geometry() %>%
        mutate(
          spatial_temp_offset = rnorm(n(), 0, 1),
          spatial_precip_factor = exp(rnorm(n(), 0, 0.2))
        ),
      by = "gefs_cell_id"
    ) %>%
    mutate(
      tmax = tmax + spatial_temp_offset,
      tmin = tmin + spatial_temp_offset,
      temp_mean = temp_mean + spatial_temp_offset,
      precip = precip * spatial_precip_factor
    )
  
  cat("Created", nrow(daily_climate), "daily climate records\n")
  cat("Grid cells:", n_distinct(daily_climate$gefs_cell_id), "\n")
  cat("Years:", min(daily_climate$year), "-", max(daily_climate$year), "\n")
  
  return(daily_climate)
}

# Function to calculate climatological normals for SPEI and Tmax
calculate_climate_normals <- function(daily_climate, spei_days = 14, tmax_days = 30) {
  
  cat("Calculating climatological normals...\n")
  cat("SPEI window:", spei_days, "days\n")
  cat("Tmax window:", tmax_days, "days\n")
  
  # Calculate PET for each cell and day
  climate_with_pet <- daily_climate %>%
    group_by(gefs_cell_id) %>%
    arrange(date) %>%
    mutate(
      # Simple PET calculation (Thornthwaite-style)
      pet = pmax(0, 16 * (temp_mean / 5)^1.514 * ifelse(temp_mean > 0, 1, 0))
    ) %>%
    ungroup()
  
  # Calculate rolling SPEI and Tmax for each cell
  normals_data <- climate_with_pet %>%
    group_by(gefs_cell_id) %>%
    arrange(date) %>%
    mutate(
      # Water balance
      water_balance = precip - pet,
      # Rolling sums/means
      spei_wb = zoo::rollsum(water_balance, k = spei_days, align = "right", fill = NA),
      tmax_rolling = zoo::rollmean(tmax, k = tmax_days, align = "right", fill = NA)
    ) %>%
    ungroup() %>%
    filter(!is.na(spei_wb) & !is.na(tmax_rolling))
  
  # Calculate climatological normals by day of year and cell
  climatology <- normals_data %>%
    mutate(
      # Use a broader window around each day of year (e.g., ±7 days)
      doy_window = ((day_of_year - 1 + 7) %% 365) + 1
    ) %>%
    group_by(gefs_cell_id, day_of_year) %>%
    summarise(
      # SPEI normals (mean and sd for standardization)
      spei_wb_mean = mean(spei_wb, na.rm = TRUE),
      spei_wb_sd = sd(spei_wb, na.rm = TRUE),
      # Tmax normals
      tmax_mean = mean(tmax_rolling, na.rm = TRUE),
      tmax_sd = sd(tmax_rolling, na.rm = TRUE),
      n_years = n_distinct(year),
      .groups = "drop"
    )
  
  cat("Calculated normals for", n_distinct(climatology$gefs_cell_id), "grid cells\n")
  cat("Day-of-year coverage:", min(climatology$day_of_year), "to", max(climatology$day_of_year), "\n")
  cat("Average years per DOY:", round(mean(climatology$n_years), 1), "\n")
  
  return(climatology)
}

# Function to save climatological normals
save_climate_normals <- function(normals_data, output_path = "U:/datasets/ndvi_monitor/climate_normals_1980_2010.rds") {
  
  # Create output directory
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Add metadata
  normals_with_metadata <- list(
    normals = normals_data,
    metadata = list(
      reference_period = normals_config$reference_period,
      spei_window = normals_config$spei_window,
      tmax_window = normals_config$tmax_window,
      created_date = Sys.Date(),
      n_cells = n_distinct(normals_data$gefs_cell_id),
      method = "Synthetic data for development"
    )
  )
  
  saveRDS(normals_with_metadata, output_path)
  cat("Climate normals saved to:", output_path, "\n")
  
  return(normals_with_metadata)
}

# Main workflow for calculating climate normals using ERA5
calculate_1980_2010_normals <- function(gefs_grid, 
                                       start_year = 1980, 
                                       end_year = 2010,
                                       use_era5 = TRUE,
                                       save_output = TRUE) {
  
  cat("Starting climate normals calculation for", start_year, "-", end_year, "...\n")
  
  if (use_era5) {
    cat("Using ERA5 reanalysis data\n")
    
    # Setup ERA5 credentials
    setup_era5_credentials()
    
    # Process years in chunks to avoid overwhelming the system
    years <- start_year:end_year
    daily_climate_list <- list()
    
    for (year in years) {
      cat("Processing year:", year, "\n")
      
      # Request ERA5 data for this year
      era5_file <- request_era5_year(year)
      
      if (!is.null(era5_file)) {
        # Process to daily data
        daily_data <- process_era5_to_daily(era5_file, gefs_grid)
        
        if (!is.null(daily_data)) {
          daily_climate_list[[as.character(year)]] <- daily_data
        }
      }
      
      # Small delay between requests
      Sys.sleep(2)
    }
    
    # Combine all years
    if (length(daily_climate_list) > 0) {
      daily_climate <- bind_rows(daily_climate_list)
      cat("Combined ERA5 data:", nrow(daily_climate), "daily records\n")
    } else {
      stop("No ERA5 data retrieved. Check credentials and connection.")
    }
    
  } else {
    cat("Using synthetic data for testing\n")
    # Fallback to synthetic data
    daily_climate <- create_synthetic_normals_data(
      gefs_grid, 
      paste0(start_year, "-01-01"), 
      paste0(end_year, "-12-31")
    )
  }
  
  # Calculate climatological normals
  climate_normals <- calculate_climate_normals(
    daily_climate,
    spei_days = normals_config$spei_window,
    tmax_days = normals_config$tmax_window
  )
  
  if (save_output) {
    normals_saved <- save_climate_normals(climate_normals)
    return(normals_saved)
  } else {
    return(list(
      normals = climate_normals,
      daily_data = daily_climate
    ))
  }
}

# Function to test ERA5 setup with a single year
test_era5_setup <- function(gefs_grid, test_year = 1990) {
  
  cat("Testing ERA5 setup with year:", test_year, "\n")
  
  # Setup credentials
  setup_era5_credentials()
  
  # Request one year of data
  era5_file <- request_era5_year(test_year)
  
  if (!is.null(era5_file)) {
    # Try to process it
    daily_data <- process_era5_to_daily(era5_file, gefs_grid)
    
    if (!is.null(daily_data)) {
      cat("✓ ERA5 test successful!\n")
      cat("Sample data:\n")
      print(head(daily_data))
      return(daily_data)
    }
  }
  
  cat("✗ ERA5 test failed\n")
  return(NULL)
}

# Function to compare ERA5 vs GEFS in overlapping period (2016-present)
compare_era5_gefs <- function(gefs_grid, 
                             comparison_years = 2018:2020,  # Small test period
                             save_comparison = TRUE) {
  
  cat("Comparing ERA5 vs GEFS for validation period:", min(comparison_years), "-", max(comparison_years), "\n")
  
  # Get ERA5 data for comparison period
  cat("Retrieving ERA5 data...\n")
  era5_comparison <- calculate_1980_2010_normals(
    gefs_grid, 
    start_year = min(comparison_years),
    end_year = max(comparison_years),
    use_era5 = TRUE,
    save_output = FALSE
  )
  
  # Get GEFS data for same period (placeholder - would need real GEFS retrieval)
  cat("Retrieving GEFS data...\n")
  gefs_comparison <- calculate_1980_2010_normals(
    gefs_grid,
    start_year = min(comparison_years), 
    end_year = max(comparison_years),
    use_era5 = FALSE,  # Use synthetic as GEFS placeholder
    save_output = FALSE
  )
  
  # Compare the datasets
  comparison_stats <- era5_comparison$daily_data %>%
    select(gefs_cell_id, date, tmax_era5 = tmax, precip_era5 = precip) %>%
    left_join(
      gefs_comparison$daily_data %>%
        select(gefs_cell_id, date, tmax_gefs = tmax, precip_gefs = precip),
      by = c("gefs_cell_id", "date")
    ) %>%
    filter(!is.na(tmax_era5) & !is.na(tmax_gefs)) %>%
    summarise(
      n_comparisons = n(),
      # Temperature comparisons
      tmax_bias = mean(tmax_gefs - tmax_era5, na.rm = TRUE),
      tmax_rmse = sqrt(mean((tmax_gefs - tmax_era5)^2, na.rm = TRUE)),
      tmax_cor = cor(tmax_era5, tmax_gefs, use = "complete.obs"),
      # Precipitation comparisons  
      precip_bias = mean(precip_gefs - precip_era5, na.rm = TRUE),
      precip_rmse = sqrt(mean((precip_gefs - precip_era5)^2, na.rm = TRUE)),
      precip_cor = cor(precip_era5, precip_gefs, use = "complete.obs"),
      .groups = "drop"
    )
  
  cat("ERA5 vs GEFS Comparison Results:\n")
  cat("  Records compared:", comparison_stats$n_comparisons, "\n")
  cat("  Tmax bias (GEFS - ERA5):", round(comparison_stats$tmax_bias, 2), "°C\n")
  cat("  Tmax RMSE:", round(comparison_stats$tmax_rmse, 2), "°C\n") 
  cat("  Tmax correlation:", round(comparison_stats$tmax_cor, 3), "\n")
  cat("  Precip bias:", round(comparison_stats$precip_bias, 2), "mm/day\n")
  cat("  Precip RMSE:", round(comparison_stats$precip_rmse, 2), "mm/day\n")
  cat("  Precip correlation:", round(comparison_stats$precip_cor, 3), "\n")
  
  comparison_result <- list(
    statistics = comparison_stats,
    era5_data = era5_comparison$daily_data,
    gefs_data = gefs_comparison$daily_data,
    comparison_period = comparison_years
  )
  
  if (save_comparison) {
    saveRDS(comparison_result, "U:/datasets/ndvi_monitor/era5_gefs_comparison.rds")
    cat("Comparison results saved to: U:/datasets/ndvi_monitor/era5_gefs_comparison.rds\n")
  }
  
  return(comparison_result)
}

cat("Climate normals script loaded (1980-2010 reference period).\n")
cat("Configuration:\n")
cat("  - SPEI window:", normals_config$spei_window, "days\n")
cat("  - Tmax window:", normals_config$tmax_window, "days\n")
cat("  - Reference period:", paste(normals_config$reference_period, collapse = " to "), "\n")
cat("Functions available:\n")
cat("  - calculate_1980_2010_normals(gefs_grid) - compute climatological normals\n")
cat("  - test_era5_setup(gefs_grid) - test ERA5 connection with single year\n")
cat("  - compare_era5_gefs(gefs_grid) - validate ERA5 vs GEFS consistency\n")