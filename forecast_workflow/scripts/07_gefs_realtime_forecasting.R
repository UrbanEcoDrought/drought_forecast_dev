# GEFS Real-time Data Extraction for Forecasting
# Extracts current GEFS forecasts from noaa-gefs-pds S3 bucket for operational forecasting
# Uses existing NOAA GEFS functions with modifications for single ensemble member

library(tidyverse)
library(sf)
library(lubridate)
library(gdalcubes)
library(stars)

# Load existing GEFS functions
source("../../00_NOAA_GEFS_functions.R")

# Configuration for real-time GEFS forecasting
realtime_config <- list(
  base_url = "https://noaa-gefs-pds.s3.amazonaws.com",
  ensemble = "gec00",  # Control run for consistency
  resolution = 0.25,   # degrees
  cycle = "00",        # 00Z cycle (most reliable)
  variables = c("TMP", "APCP", "DSWRF", "RH"),  # Core variables
  forecast_horizon = "003",  # 3-hour forecast (more common)
  retention_days = 30  # Typical retention period on AWS
)

# Function to get available GEFS forecast dates
get_available_forecast_dates <- function(max_days_back = 5) {
  
  cat("Checking available GEFS forecast dates...\n")
  
  available_dates <- c()
  
  # Try today first (GEFS runs 4x daily)
  test_dates <- c(Sys.Date(), Sys.Date() - seq_len(max_days_back))
  
  for (test_date in test_dates) {
    cat("  Checking date:", as.character(test_date), "\n")
    
    # Check if data exists for this date
    if (check_realtime_gefs_availability(test_date)) {
      available_dates <- c(available_dates, as.character(test_date))
    }
  }
  
  if (length(available_dates) > 0) {
    cat("Found", length(available_dates), "dates with available data\n")
    cat("Most recent:", min(available_dates), "\n")
    cat("Oldest:", max(available_dates), "\n")
  } else {
    cat("No available GEFS data found in last", max_days_back, "days\n")
    cat("NOTE: Real-time GEFS data may have limited availability or very short retention\n")
    cat("This is normal - operational forecasts are prioritized for real-time users\n")
  }
  
  return(as.Date(available_dates))
}

# Function to check real-time GEFS data availability
check_realtime_gefs_availability <- function(date) {
  
  date_str <- format(as.Date(date), "%Y%m%d")
  
  # Try multiple forecast horizons - real-time may have different availability
  horizons_to_try <- c("003", "006", "009", "012")
  
  for (horizon in horizons_to_try) {
    # Test URL structure - use official GEFS structure from documentation
    test_url <- paste0(
      realtime_config$base_url, "/gefs.", date_str, "/", 
      realtime_config$cycle, "/pgrb2a/", realtime_config$ensemble, 
      ".t", realtime_config$cycle, "z.pgrb2af", horizon
    )
    
    cat("    Testing URL:", test_url, "\n")
    
    tryCatch({
      response <- httr::HEAD(test_url)
      status <- httr::status_code(response)
      available <- status == 200
      
      cat("    Status:", status, "Available:", available, "\n")
      
      if (available) {
        cat("    Found data at horizon f", horizon, "\n")
        return(TRUE)
      }
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
    })
  }
  
  cat("    No data found for any forecast horizon\n")
  return(FALSE)
}

# Function to extract real-time GEFS data using existing functions
extract_realtime_gefs_data <- function(date, gefs_centers, variables = realtime_config$variables) {
  
  cat("Extracting real-time GEFS data for", as.character(date), "\n")
  
  # Convert GEFS centers to spatial points for extraction
  gefs_points <- gefs_centers %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326") %>%
    rowid_to_column("FID")  # Add FID column expected by GEFS functions
  
  tryCatch({
    # Use existing GEFS extraction functions
    gefs_data <- grib_extract(
      ens = realtime_config$ensemble,
      date = date,
      sites = gefs_points,
      bands = variables,
      horizon = realtime_config$forecast_horizon
    )
    
    if (!is.null(gefs_data) && nrow(gefs_data) > 0) {
      # Process extracted data
      climate_df <- gefs_data %>%
        as_tibble() %>%
        mutate(
          date = date,
          year = year(date),
          month = month(date),
          day_of_year = yday(date),
          # Convert units and rename variables
          temperature = if("TMP" %in% names(.)) TMP - 273.15 else NA,  # K to C
          precipitation = if("APCP" %in% names(.)) APCP else NA,       # mm
          solar = if("DSWRF" %in% names(.)) DSWRF else NA,             # W/m2
          humidity = if("RH" %in% names(.)) RH else NA                 # %
        ) %>%
        # Add GEFS cell IDs
        mutate(gefs_cell_id = gefs_centers$gefs_cell_id[FID]) %>%
        select(gefs_cell_id, date, year, month, day_of_year, 
               temperature, precipitation, solar, humidity)
      
      cat("Successfully extracted data for", nrow(climate_df), "GEFS cells\n")
      return(climate_df)
      
    } else {
      cat("No data returned from GEFS extraction\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("Error extracting real-time GEFS data:", e$message, "\n")
    return(NULL)
  })
}

# Function to get most recent GEFS forecast
get_latest_gefs_forecast <- function(gefs_centers, max_days_back = 7) {
  
  cat("Getting latest GEFS forecast data...\n")
  
  # Find most recent available date
  for (i in 1:max_days_back) {
    test_date <- Sys.Date() - i
    
    if (check_realtime_gefs_availability(test_date)) {
      cat("Using GEFS data from", as.character(test_date), "\n")
      
      forecast_data <- extract_realtime_gefs_data(test_date, gefs_centers)
      
      if (!is.null(forecast_data)) {
        return(forecast_data)
      }
    }
  }
  
  cat("No recent GEFS data available\n")
  return(NULL)
}

# Function to build forecast time series over multiple days
build_forecast_time_series <- function(gefs_centers, 
                                      start_date = Sys.Date() - 10,
                                      end_date = Sys.Date() - 1) {
  
  cat("Building GEFS forecast time series...\n")
  cat("Period:", as.character(start_date), "to", as.character(end_date), "\n")
  
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  forecast_data <- list()
  
  for (i in seq_along(dates)) {
    date <- dates[i]
    
    cat(paste0("[", i, "/", length(dates), "] "))
    
    if (check_realtime_gefs_availability(date)) {
      daily_data <- extract_realtime_gefs_data(date, gefs_centers)
      
      if (!is.null(daily_data)) {
        forecast_data[[i]] <- daily_data
      }
    } else {
      cat("No data available for", as.character(date), "\n")
    }
    
    # Small delay
    Sys.sleep(0.2)
  }
  
  if (length(forecast_data) > 0) {
    combined_data <- bind_rows(forecast_data)
    
    cat("\nForecast time series completed:\n")
    cat("  Records:", nrow(combined_data), "\n")
    cat("  Date range:", range(combined_data$date), "\n")
    cat("  GEFS cells:", n_distinct(combined_data$gefs_cell_id), "\n")
    
    return(combined_data)
  } else {
    cat("No forecast data was successfully extracted\n")
    return(NULL)
  }
}

# Function to save forecast data
save_forecast_data <- function(forecast_data, 
                              output_path = "U:/datasets/ndvi_monitor/gefs_realtime_forecast_data.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Add metadata
  forecast_with_metadata <- list(
    climate_data = forecast_data,
    metadata = list(
      source = "GEFS Real-time",
      ensemble = realtime_config$ensemble,
      variables = c("temperature", "precipitation", "solar", "humidity"),
      period = range(forecast_data$date),
      n_records = nrow(forecast_data),
      n_cells = n_distinct(forecast_data$gefs_cell_id),
      forecast_horizon = realtime_config$forecast_horizon,
      created_date = Sys.Date(),
      url_base = realtime_config$base_url
    )
  )
  
  saveRDS(forecast_with_metadata, output_path)
  cat("Forecast data saved to:", output_path, "\n")
  
  return(forecast_with_metadata)
}

# Main workflow function for forecasting pipeline
extract_gefs_realtime_forecast <- function(gefs_centers,
                                          mode = "latest",  # "latest" or "time_series"
                                          start_date = Sys.Date() - 10,
                                          end_date = Sys.Date() - 1,
                                          save_output = TRUE) {
  
  cat("Starting GEFS real-time forecast extraction...\n")
  cat("Data source:", realtime_config$base_url, "\n")
  cat("Ensemble:", realtime_config$ensemble, "\n")
  cat("Mode:", mode, "\n")
  
  if (mode == "latest") {
    # Get most recent forecast
    forecast_data <- get_latest_gefs_forecast(gefs_centers)
  } else {
    # Get time series of forecasts
    forecast_data <- build_forecast_time_series(gefs_centers, start_date, end_date)
  }
  
  if (!is.null(forecast_data) && nrow(forecast_data) > 0) {
    if (save_output) {
      result <- save_forecast_data(forecast_data)
      return(result)
    } else {
      return(forecast_data)
    }
  } else {
    stop("No forecast data was successfully extracted")
  }
}

cat("GEFS real-time forecasting script loaded.\n")
cat("Data source:", realtime_config$base_url, "\n")
cat("Ensemble:", realtime_config$ensemble, "\n")
cat("Variables:", paste(realtime_config$variables, collapse = ", "), "\n")
cat("Run extract_gefs_realtime_forecast(gefs_centers) for latest forecast.\n")
cat("Use get_available_forecast_dates() to check data availability.\n")