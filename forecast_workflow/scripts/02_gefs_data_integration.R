# GEFS Climate Data Integration Script
# Integrates NOAA GEFS 0.25° climate data with NDVI grid structure

library(tidyverse)
library(sf)
library(lubridate)
library(httr)

# Custom GEFS functions for drought forecasting

# GEFS variable mapping for 0.25 degree data
gefs_variables <- list(
  "tmp2m" = "2m_above_ground/tmp/",     # 2m temperature 
  "rh2m" = "2m_above_ground/rh/",       # 2m relative humidity
  "apcp" = "surface/apcp/",             # precipitation
  "dswrf" = "surface/dswrf/"            # downward shortwave radiation
)

# Configuration for GEFS data
gefs_config <- list(
  resolution = 0.25,  # degrees (0.25° = ~28km)
  variables = names(gefs_variables),
  ensemble = "geavg",  # Ensemble mean
  base_url = "https://noaa-gefs-pds.s3.amazonaws.com",
  cycle = "00"         # Daily 00Z cycle
)

# Function to build GEFS URLs for specific variables and dates
build_gefs_urls <- function(date, variables = gefs_config$variables, 
                           ensemble = "geavg", cycle = "00") {
  
  date_str <- format(as.Date(date), "%Y%m%d")
  base_path <- paste0(gefs_config$base_url, "/gefs.", date_str, "/", cycle, "/atmos/pgrb2sp25/")
  
  urls <- list()
  for (var in variables) {
    if (var %in% names(gefs_variables)) {
      # Build URL for this variable
      url <- paste0(base_path, ensemble, ".t", cycle, "z.pgrb2sp25.0p25.f006")  # 6-hour forecast
      urls[[var]] <- url
    }
  }
  
  return(urls)
}

# Function to check GEFS data availability for a date
check_gefs_availability <- function(date) {
  urls <- build_gefs_urls(date, variables = "tmp2m")
  
  tryCatch({
    # Try to access the URL to see if data exists
    response <- httr::HEAD(urls[["tmp2m"]])
    return(httr::status_code(response) == 200)
  }, error = function(e) {
    return(FALSE)
  })
}

# Function to get GEFS grid cell coordinates from our spatial structure
extract_gefs_coordinates <- function(ndvi_gefs_data) {
  
  # Get center coordinates of each GEFS cell
  gefs_centers <- ndvi_gefs_data$gefs_grid %>%
    mutate(
      centroid = st_centroid(geometry),
      longitude = st_coordinates(centroid)[,1],
      latitude = st_coordinates(centroid)[,2]
    ) %>%
    st_drop_geometry() %>%
    select(gefs_cell_id, longitude, latitude)
  
  cat("Extracted coordinates for", nrow(gefs_centers), "GEFS cells\n")
  cat("Longitude range:", range(gefs_centers$longitude), "\n")
  cat("Latitude range:", range(gefs_centers$latitude), "\n")
  
  return(gefs_centers)
}

# Function to create sf points for GEFS cell centers
create_gefs_points <- function(gefs_centers) {
  
  gefs_points <- gefs_centers %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326")
  
  return(gefs_points)
}

# Simplified function to create synthetic climate data for testing
# (Replace with real GEFS data retrieval once we get the URLs working)
create_synthetic_climate_data <- function(gefs_centers, 
                                         start_date = "2023-01-01", 
                                         end_date = "2023-03-31") {
  
  cat("Creating synthetic climate data for testing...\n")
  
  # Create date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "month")
  
  # Create synthetic climate data for each GEFS cell and date
  climate_data <- expand_grid(
    gefs_cell_id = gefs_centers$gefs_cell_id,
    date = dates
  ) %>%
    mutate(
      year = year(date),
      month = month(date),
      # Create realistic synthetic climate variables
      temp_mean = 15 + 10 * sin(2 * pi * (month - 1) / 12) + rnorm(n(), 0, 2),  # Seasonal temp
      rh_mean = 65 + 10 * sin(2 * pi * (month - 6) / 12) + rnorm(n(), 0, 5),    # Seasonal humidity  
      precip_total = pmax(0, rlnorm(n(), log(50), 0.5)),                        # Log-normal precip
      solar_mean = 200 + 100 * sin(2 * pi * (month - 3) / 12) + rnorm(n(), 0, 20) # Seasonal solar
    ) %>%
    left_join(gefs_centers, by = "gefs_cell_id")
  
  cat("Created", nrow(climate_data), "synthetic climate records\n")
  cat("GEFS cells:", n_distinct(climate_data$gefs_cell_id), "\n")
  cat("Months:", n_distinct(climate_data$date), "\n")
  
  return(climate_data)
}

# Function to pull real GEFS data (placeholder for now)
get_real_gefs_data <- function(gefs_centers, start_date, end_date) {
  
  cat("Real GEFS data retrieval not yet implemented.\n")
  cat("Using synthetic data for development.\n")
  
  # Check if recent dates are available
  recent_date <- Sys.Date() - 7  # 7 days ago
  
  if (check_gefs_availability(recent_date)) {
    cat("GEFS data appears to be available for recent dates.\n")
  } else {
    cat("GEFS data availability check failed.\n")
  }
  
  # For now, return synthetic data
  return(create_synthetic_climate_data(gefs_centers, start_date, end_date))
}

# Function to match climate data to GEFS cells
match_climate_to_gefs <- function(climate_data, gefs_centers) {
  
  # Join climate data with GEFS cell IDs
  climate_matched <- climate_data %>%
    # Assuming climate data has FID or similar identifier matching row order
    mutate(gefs_cell_id = rep(gefs_centers$gefs_cell_id, 
                              length.out = nrow(climate_data))) %>%
    left_join(gefs_centers, by = "gefs_cell_id")
  
  return(climate_matched)
}

# Climate data is already monthly, so just pass through with validation
validate_climate_data <- function(climate_data) {
  
  # Check that we have the expected columns
  expected_cols <- c("gefs_cell_id", "date", "year", "month", "temp_mean", "rh_mean", "precip_total", "solar_mean")
  missing_cols <- setdiff(expected_cols, colnames(climate_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing climate data columns: ", paste(missing_cols, collapse = ", "))
  }
  
  cat("Validated", nrow(climate_data), "monthly climate records\n")
  cat("Variables: temp_mean, rh_mean, precip_total, solar_mean\n")
  
  return(climate_data)
}

# Function to integrate climate with NDVI data
integrate_climate_ndvi <- function(ndvi_gefs_data, climate_monthly) {
  
  # Create monthly NDVI aggregations
  ndvi_monthly <- ndvi_gefs_data$ndvi_data %>%
    st_drop_geometry() %>%
    mutate(month = month(date)) %>%
    group_by(gefs_cell_id, year, month) %>%
    summarise(
      ndvi_anomaly_mean = mean(ndvi_anomaly, na.rm = TRUE),
      ndvi_anomaly_sd = sd(ndvi_anomaly, na.rm = TRUE),
      n_pixels = n(),
      .groups = "drop"
    ) %>%
    mutate(date = as.Date(paste(year, month, "01", sep = "-")))
  
  # Join NDVI and climate data
  integrated_data <- ndvi_monthly %>%
    left_join(climate_monthly, by = c("gefs_cell_id", "year", "month", "date")) %>%
    filter(!is.na(temp_mean))  # Keep only records with climate data
  
  cat("Integrated dataset:", nrow(integrated_data), "records\n")
  cat("Date range:", range(integrated_data$date), "\n")
  cat("GEFS cells:", n_distinct(integrated_data$gefs_cell_id), "\n")
  
  return(integrated_data)
}

# Main workflow function
run_gefs_integration <- function(ndvi_gefs_data, 
                                 start_date = "2023-01-01",
                                 end_date = "2023-03-31",  # Start small for testing
                                 save_output = TRUE) {
  
  cat("Starting GEFS climate data integration...\n")
  
  # Extract GEFS grid coordinates
  gefs_centers <- extract_gefs_coordinates(ndvi_gefs_data)
  gefs_points <- create_gefs_points(gefs_centers)
  
  # Pull climate data (using synthetic for now)
  climate_data <- get_real_gefs_data(gefs_centers, start_date, end_date)
  
  if (!is.null(climate_data)) {
    # Validate climate data (already monthly from synthetic function)
    climate_validated <- validate_climate_data(climate_data)
    
    # Integrate with NDVI
    integrated_data <- integrate_climate_ndvi(ndvi_gefs_data, climate_validated)
    
    if (save_output) {
      saveRDS(integrated_data, "forecast_workflow/data/integrated_ndvi_climate.rds")
      cat("Integrated data saved to: forecast_workflow/data/integrated_ndvi_climate.rds\n")
    }
    
    return(list(
      climate_data = climate_validated,
      integrated_data = integrated_data,
      gefs_grid = ndvi_gefs_data$gefs_grid
    ))
  } else {
    stop("Failed to retrieve climate data")
  }
}

cat("GEFS integration script loaded.\n")
cat("GEFS resolution: 0.25° (~28km cells)\n") 
cat("Variables:", paste(gefs_config$variables, collapse = ", "), "\n")
cat("Run run_gefs_integration(ndvi_gefs_data) to integrate climate data.\n")