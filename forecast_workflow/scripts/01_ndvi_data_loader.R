# NDVI Data Housing Script
# Loads and processes NDVI anomaly data for forecast modeling

library(tidyverse)
library(sf)
library(stars)
library(lubridate)

# Configuration
config <- list(
  data_path = "forecast_workflow/data/ndvi/",
  crs = "EPSG:4326",  # WGS84 for compatibility with GEFS
  temporal_resolution = "monthly",
  spatial_bounds = NULL  # Will be set based on data
)

# Function to load spatial NDVI anomaly data from 16-day window dataset
load_spatial_ndvi_anomalies <- function(file_path = NULL, format = "csv") {
  
  if (is.null(file_path)) {
    cat("Please specify file path for spatial NDVI anomaly data\n")
    cat("Expected format: CSV with columns [xy, x, y, yday, year, norm, norm_lwr, norm_upr, mean, lwr, upr, anoms_mean, anoms_lwr, anoms_upr]\n")
    return(NULL)
  }
  
  if (!file.exists(file_path)) {
    stop("NDVI data file not found: ", file_path)
  }
  
  # Load data based on format
  if (format == "csv") {
    ndvi_data <- read_csv(file_path, show_col_types = FALSE)
  } else if (format == "rds") {
    ndvi_data <- readRDS(file_path)
  } else {
    stop("Unsupported format. Use 'csv' or 'rds'")
  }
  
  # Validate expected columns for spatial anomaly dataset
  expected_cols <- c("xy", "x", "y", "yday", "year", "norm", "norm_lwr", "norm_upr", 
                     "mean", "lwr", "upr", "anoms_mean", "anoms_lwr", "anoms_upr")
  missing_cols <- setdiff(expected_cols, colnames(ndvi_data))
  
  if (length(missing_cols) > 0) {
    warning("Missing expected columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process data focusing on anomalies
  ndvi_processed <- ndvi_data %>%
    mutate(
      longitude = as.numeric(x),
      latitude = as.numeric(y),
      day_of_year = as.numeric(yday),
      year = as.numeric(year),
      # Create approximate date from year and day of year
      date = as.Date(paste(year, day_of_year), format = "%Y %j"),
      month = month(date),
      # Focus on anomaly variables
      ndvi_anomaly = as.numeric(anoms_mean),
      ndvi_anomaly_lwr = as.numeric(anoms_lwr), 
      ndvi_anomaly_upr = as.numeric(anoms_upr),
      # Keep observed NDVI for reference
      ndvi_observed = as.numeric(mean),
      ndvi_observed_lwr = as.numeric(lwr),
      ndvi_observed_upr = as.numeric(upr),
      # Keep climatological normal for reference  
      ndvi_normal = as.numeric(norm),
      ndvi_normal_lwr = as.numeric(norm_lwr),
      ndvi_normal_upr = as.numeric(norm_upr)
    ) %>%
    filter(!is.na(ndvi_anomaly)) %>%
    arrange(year, day_of_year, longitude, latitude)
  
  cat("Loaded", nrow(ndvi_processed), "spatial NDVI anomaly observations\n")
  cat("Date range:", as.character(min(ndvi_processed$date, na.rm = TRUE)), "to", as.character(max(ndvi_processed$date, na.rm = TRUE)), "\n")
  cat("Year range:", min(ndvi_processed$year), "to", max(ndvi_processed$year), "\n")
  cat("Day of year range:", min(ndvi_processed$day_of_year), "to", max(ndvi_processed$day_of_year), "\n")
  cat("Spatial extent:\n")
  cat("  Longitude:", range(ndvi_processed$longitude), "\n")
  cat("  Latitude:", range(ndvi_processed$latitude), "\n")
  cat("Anomaly statistics:\n")
  cat("  Mean anomaly:", round(mean(ndvi_processed$ndvi_anomaly, na.rm = TRUE), 4), "\n")
  cat("  Anomaly range:", range(ndvi_processed$ndvi_anomaly, na.rm = TRUE), "\n")
  
  return(ndvi_processed)
}

# Function to convert to spatial object
create_spatial_ndvi <- function(ndvi_data) {
  
  # Remove any NA coordinates first
  clean_data <- ndvi_data %>%
    filter(!is.na(longitude) & !is.na(latitude))
  
  # Create simple features with explicit CRS
  ndvi_sf <- clean_data %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326")
  
  return(ndvi_sf)
}

# Function to create GEFS grid overlay (0.25 degree cells containing multiple NDVI pixels)
create_gefs_grid_overlay <- function(ndvi_data, gefs_resolution = 0.25) {
  
  # Create GEFS grid based on data extent
  bbox <- st_bbox(ndvi_data)
  
  # Snap to GEFS grid boundaries
  lon_min <- floor(bbox[["xmin"]] / gefs_resolution) * gefs_resolution
  lon_max <- ceiling(bbox[["xmax"]] / gefs_resolution) * gefs_resolution
  lat_min <- floor(bbox[["ymin"]] / gefs_resolution) * gefs_resolution
  lat_max <- ceiling(bbox[["ymax"]] / gefs_resolution) * gefs_resolution
  
  cat("Grid bounds: lon", lon_min, "to", lon_max, ", lat", lat_min, "to", lat_max, "\n")
  
  # Calculate grid dimensions
  n_lon <- round((lon_max - lon_min) / gefs_resolution)
  n_lat <- round((lat_max - lat_min) / gefs_resolution)
  
  cat("Grid dimensions:", n_lon, "x", n_lat, "cells\n")
  
  # Create GEFS grid manually to avoid st_make_grid issues
  lon_seq <- seq(lon_min, lon_max - gefs_resolution, by = gefs_resolution)
  lat_seq <- seq(lat_min, lat_max - gefs_resolution, by = gefs_resolution)
  
  # Create grid centers
  grid_centers <- expand.grid(lon = lon_seq + gefs_resolution/2, 
                              lat = lat_seq + gefs_resolution/2)
  
  # Create polygons for each grid cell
  grid_polys <- list()
  for (i in 1:nrow(grid_centers)) {
    lon_c <- grid_centers$lon[i]
    lat_c <- grid_centers$lat[i]
    
    poly <- st_polygon(list(matrix(c(
      lon_c - gefs_resolution/2, lat_c - gefs_resolution/2,
      lon_c + gefs_resolution/2, lat_c - gefs_resolution/2, 
      lon_c + gefs_resolution/2, lat_c + gefs_resolution/2,
      lon_c - gefs_resolution/2, lat_c + gefs_resolution/2,
      lon_c - gefs_resolution/2, lat_c - gefs_resolution/2
    ), ncol = 2, byrow = TRUE)))
    
    grid_polys[[i]] <- poly
  }
  
  # Create sf object
  gefs_grid <- st_sf(
    gefs_cell_id = 1:length(grid_polys),
    geometry = st_sfc(grid_polys, crs = st_crs(ndvi_data))
  )
  
  # Assign each NDVI pixel to its parent GEFS cell
  ndvi_with_gefs <- st_join(ndvi_data, gefs_grid) %>%
    filter(!is.na(gefs_cell_id))
  
  cat("Created", nrow(gefs_grid), "GEFS grid cells at", gefs_resolution, "degree resolution\n")
  cat("Assigned", nrow(ndvi_with_gefs), "NDVI pixels to GEFS cells\n")
  
  return(list(
    ndvi_data = ndvi_with_gefs,
    gefs_grid = gefs_grid
  ))
}

# Function to summarize NDVI within GEFS cells (optional aggregation)
summarize_ndvi_by_gefs <- function(ndvi_gefs_data) {
  
  ndvi_summary <- ndvi_gefs_data$ndvi_data %>%
    st_drop_geometry() %>%
    group_by(gefs_cell_id, date) %>%
    summarise(
      ndvi_anomaly_mean = mean(ndvi_anomaly, na.rm = TRUE),
      ndvi_anomaly_sd = sd(ndvi_anomaly, na.rm = TRUE),
      ndvi_anomaly_median = median(ndvi_anomaly, na.rm = TRUE),
      n_pixels = n(),
      .groups = "drop"
    ) %>%
    left_join(ndvi_gefs_data$gefs_grid, by = "gefs_cell_id") %>%
    st_as_sf()
  
  return(ndvi_summary)
}

# Function to save processed data
save_processed_ndvi <- function(ndvi_data, output_path = "U:/datasets/ndvi_monitor/processed_ndvi_anomalies.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ndvi_data, output_path)
  cat("Processed NDVI data saved to:", output_path, "\n")
  
}

# Data paths for research group access
data_paths <- list(
  # Local mounted Google Drive path
  local_path = "G:/Shared drives/Urban Ecological Drought/data/spatial_NDVI_monitoring/16_day_window_spatial_data_with_anomalies.csv",
  # Google Drive link for research group sharing
  google_drive_link = "https://drive.google.com/file/d/1vC86OZe-Vl2ZNbzOU6qZPwNOo_lwiSuE/view",
  # Google Drive direct download link
  direct_download = "https://drive.google.com/uc?export=download&id=1vC86OZe-Vl2ZNbzOU6qZPwNOo_lwiSuE"
)

# Function to get NDVI data path based on availability
get_ndvi_data_path <- function() {
  if (file.exists(data_paths$local_path)) {
    cat("Using local mounted Google Drive path\n")
    return(data_paths$local_path)
  } else {
    cat("Local path not found. Use Google Drive link:\n")
    cat("View: ", data_paths$google_drive_link, "\n") 
    cat("Direct download: ", data_paths$direct_download, "\n")
    return(NULL)
  }
}

# Ready-to-use workflow for spatial NDVI anomalies
run_spatial_ndvi_processing <- function(save_output = TRUE) {
  
  # Get data path
  data_file <- get_ndvi_data_path()
  
  if (is.null(data_file)) {
    stop("Please download data from Google Drive link or mount drive locally")
  }
  
  cat("Processing spatial NDVI anomaly data...\n")
  
  # Load and process data
  ndvi_raw <- load_spatial_ndvi_anomalies(data_file)
  ndvi_spatial <- create_spatial_ndvi(ndvi_raw)
  ndvi_with_gefs <- create_gefs_grid_overlay(ndvi_spatial)
  
  if (save_output) {
    save_processed_ndvi(ndvi_with_gefs, "U:/datasets/ndvi_monitor/processed_spatial_ndvi_anomalies.rds")
  }
  
  return(ndvi_with_gefs)
}

cat("NDVI data housing script loaded with spatial anomaly data paths.\n")
cat("Run run_spatial_ndvi_processing() to process the 16-day window spatial anomaly data.\n")
cat("Focusing on: anoms_mean, anoms_lwr, anoms_upr (NDVI anomalies)\n")
cat("Google Drive access: ", data_paths$google_drive_link, "\n")