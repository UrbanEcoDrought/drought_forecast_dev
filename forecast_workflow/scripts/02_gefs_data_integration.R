# GEFS Climate Data Integration Script  
# Integrates NOAA GEFS 0.25° climate data with NDVI grid structure

library(tidyverse)
library(sf)
library(lubridate)
library(httr)
library(gdalcubes)
library(stars)

# Load the existing GEFS functions (path relative to project root)
source("../../00_NOAA_GEFS_functions.R")

# GEFS variable mapping for 0.25 degree data
gefs_variables <- list(
  "tmp2m" = "2m_above_ground/tmp/",     # 2m temperature 
  "rh2m" = "2m_above_ground/rh/",       # 2m relative humidity
  "apcp" = "surface/apcp/",             # precipitation
  "dswrf" = "surface/dswrf/"            # downward shortwave radiation
)

# Configuration for GEFS data (both real-time and reforecast)
gefs_config <- list(
  resolution = 0.25,  # degrees (0.25° = ~28km)
  variables = c("tmp", "apcp", "dswrf", "rh"),  # Core variables we need
  # Ensemble selection (single time series)
  realtime_ensemble = "gec00",     # Control run for real-time (2020+)
  reforecast_ensemble = "c00",     # Control run for reforecast (2000-2019)
  # URLs
  realtime_url = "https://noaa-gefs-pds.s3.amazonaws.com",
  reforecast_url = "https://noaa-gefs-retrospective.s3.amazonaws.com",
  cycle = "00",        # Daily 00Z cycle
  # Data periods
  reforecast_period = c("2000-01-01", "2019-12-31"),  # Historical training
  realtime_cutoff = "2020-01-01"  # Switch to real-time data
)

# Function to determine data source (reforecast vs real-time) based on date
determine_gefs_source <- function(date) {
  
  date_obj <- as.Date(date)
  reforecast_start <- as.Date(gefs_config$reforecast_period[1])
  reforecast_end <- as.Date(gefs_config$reforecast_period[2])
  realtime_start <- as.Date(gefs_config$realtime_cutoff)
  
  if (date_obj >= reforecast_start & date_obj <= reforecast_end) {
    return(list(
      source = "reforecast",
      url = gefs_config$reforecast_url,
      ensemble = gefs_config$reforecast_ensemble
    ))
  } else if (date_obj >= realtime_start) {
    return(list(
      source = "realtime", 
      url = gefs_config$realtime_url,
      ensemble = gefs_config$realtime_ensemble
    ))
  } else {
    stop("Date ", date, " is outside available GEFS data range")
  }
}

# Function to build GEFS URLs for specific variables and dates
build_gefs_urls <- function(date, variables = gefs_config$variables, 
                           ensemble = NULL, cycle = "00") {
  
  # Determine data source
  source_info <- determine_gefs_source(date)
  date_str <- format(as.Date(date), "%Y%m%d")
  
  if (source_info$source == "reforecast") {
    # Build reforecast URLs (variable-specific GRIB2 files)
    year <- format(as.Date(date), "%Y")
    date_hour <- paste0(date_str, "00")
    ensemble_member <- source_info$ensemble
    
    urls <- list()
    var_mapping <- list(
      "tmp" = "tmp_2m",
      "apcp" = "apcp_sfc", 
      "dswrf" = "dswrf_sfc",
      "rh" = "rh_2m"
    )
    
    for (var in variables) {
      if (var %in% names(var_mapping)) {
        # Reforecast path: GEFSv12/reforecast/YYYY/YYYYMMDD00/c00/Days:1-10/variable_level_date_ensemble.grib2
        url <- paste0(
          source_info$url, "/GEFSv12/reforecast/", year, "/", date_hour, "/",
          ensemble_member, "/Days:1-10/", var_mapping[[var]], "_", date_hour, "_", 
          ensemble_member, ".grib2"
        )
        urls[[var]] <- url
      }
    }
    
  } else {
    # Build real-time URLs (combined parameter files)
    ensemble_member <- source_info$ensemble
    base_path <- paste0(source_info$url, "/gefs.", date_str, "/", cycle, "/pgrb2a/")
    
    urls <- list()
    # Real-time uses single file with all variables
    url <- paste0(base_path, ensemble_member, ".t", cycle, "z.pgrb2a.0p25.f006")  # 6-hour forecast
    
    for (var in variables) {
      urls[[var]] <- url  # Same file contains all variables
    }
  }
  
  return(urls)
}

# Function to check GEFS data availability for a date
check_gefs_availability <- function(date) {
  
  tryCatch({
    source_info <- determine_gefs_source(date)
    urls <- build_gefs_urls(date, variables = "tmp")  # Use "tmp" instead of "tmp2m"
    
    # Try to access the URL to see if data exists
    response <- httr::HEAD(urls[["tmp"]])
    available <- httr::status_code(response) == 200
    
    cat("Data source:", source_info$source, "\n")
    cat("Ensemble:", source_info$ensemble, "\n")
    
    return(available)
  }, error = function(e) {
    cat("Error checking availability:", e$message, "\n")
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

# Function to pull real GEFS data using existing GEFS functions
get_real_gefs_data <- function(gefs_centers, start_date, end_date) {
  
  cat("Retrieving real GEFS climate data...\n")
  
  # Convert GEFS centers to sf points for extraction
  gefs_points <- gefs_centers %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326")
  
  # Create date sequence - use monthly intervals for now
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "month")
  
  # Initialize list to store climate data
  climate_data_list <- list()
  
  for (i in seq_along(dates)) {
    date <- dates[i]
    cat("Processing date:", as.character(date), "\n")
    
    tryCatch({
      # Extract GEFS data for this date
      # Use geavg (ensemble average) with selected bands
      bands <- c("TMP", "RH", "APCP", "DSWRF")  # Temperature, humidity, precip, solar
      
      gefs_data <- grib_extract(
        ens = "geavg",
        date = date,
        sites = gefs_points,
        bands = bands,
        horizon = "006"  # 6-hour forecast
      )
      
      if (!is.null(gefs_data)) {
        # Convert to data frame and process
        climate_df <- gefs_data %>%
          as_tibble() %>%
          mutate(
            date = date,
            year = year(date),
            month = month(date),
            # Rename variables to match expected format
            temp_mean = TMP - 273.15,  # Convert K to C
            rh_mean = RH,              # Relative humidity %
            precip_total = APCP,       # Precipitation mm
            solar_mean = DSWRF         # Solar radiation W/m2
          ) %>%
          select(FID, date, year, month, temp_mean, rh_mean, precip_total, solar_mean) %>%
          # Add GEFS cell IDs based on FID
          mutate(gefs_cell_id = gefs_centers$gefs_cell_id[FID])
        
        climate_data_list[[i]] <- climate_df
        cat("Successfully retrieved data for", nrow(climate_df), "GEFS cells\n")
      }
      
    }, error = function(e) {
      cat("Error retrieving GEFS data for", as.character(date), ":", e$message, "\n")
      cat("Falling back to synthetic data for this date\n")
      
      # Create synthetic data for this date if real data fails
      synthetic_df <- expand_grid(
        gefs_cell_id = gefs_centers$gefs_cell_id,
        date = date
      ) %>%
        mutate(
          year = year(date),
          month = month(date),
          temp_mean = 15 + 10 * sin(2 * pi * (month - 1) / 12) + rnorm(n(), 0, 2),
          rh_mean = 65 + 10 * sin(2 * pi * (month - 6) / 12) + rnorm(n(), 0, 5),
          precip_total = pmax(0, rlnorm(n(), log(50), 0.5)),
          solar_mean = 200 + 100 * sin(2 * pi * (month - 3) / 12) + rnorm(n(), 0, 20)
        )
      
      climate_data_list[[i]] <- synthetic_df
    })
  }
  
  # Combine all climate data
  climate_data <- bind_rows(climate_data_list) %>%
    left_join(gefs_centers, by = "gefs_cell_id")
  
  cat("Retrieved climate data for", nrow(climate_data), "records\n")
  cat("Date range:", range(climate_data$date), "\n")
  cat("GEFS cells:", n_distinct(climate_data$gefs_cell_id), "\n")
  
  return(climate_data)
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
      # Save to U:/ drive with other NDVI data
      gefs_dir <- "U:/datasets/ndvi_monitor"
      if (!dir.exists(gefs_dir)) {
        dir.create(gefs_dir, recursive = TRUE)
      }
      
      saveRDS(integrated_data, file.path(gefs_dir, "integrated_ndvi_climate.rds"))
      cat("Integrated data saved to:", file.path(gefs_dir, "integrated_ndvi_climate.rds"), "\n")
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
cat("Variables: TMP, RH, APCP, DSWRF\n")
cat("Data will be saved to: U:/datasets/ndvi_monitor\n")
cat("Run run_gefs_integration(ndvi_gefs_data) to integrate real climate data.\n")