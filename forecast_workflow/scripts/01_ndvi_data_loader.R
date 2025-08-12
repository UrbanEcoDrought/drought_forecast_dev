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

# Function to load NDVI anomaly data
load_ndvi_anomalies <- function(file_path = NULL, format = "csv") {
  
  if (is.null(file_path)) {
    cat("Please specify file path for NDVI anomaly data\n")
    cat("Expected format: CSV with columns [date, longitude, latitude, ndvi_anomaly]\n")
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
  
  # Validate required columns
  required_cols <- c("date", "longitude", "latitude", "ndvi_anomaly")
  missing_cols <- setdiff(required_cols, colnames(ndvi_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process data
  ndvi_processed <- ndvi_data %>%
    mutate(
      date = as.Date(date),
      year = year(date),
      month = month(date),
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude),
      ndvi_anomaly = as.numeric(ndvi_anomaly)
    ) %>%
    filter(!is.na(ndvi_anomaly)) %>%
    arrange(date, longitude, latitude)
  
  cat("Loaded", nrow(ndvi_processed), "NDVI anomaly observations\n")
  cat("Date range:", min(ndvi_processed$date), "to", max(ndvi_processed$date), "\n")
  cat("Spatial extent:\n")
  cat("  Longitude:", range(ndvi_processed$longitude), "\n")
  cat("  Latitude:", range(ndvi_processed$latitude), "\n")
  
  return(ndvi_processed)
}

# Function to convert to spatial object
create_spatial_ndvi <- function(ndvi_data, target_crs = "EPSG:4326") {
  
  ndvi_sf <- ndvi_data %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(target_crs)
  
  return(ndvi_sf)
}

# Function to create GEFS grid overlay (0.25 degree cells containing multiple NDVI pixels)
create_gefs_grid_overlay <- function(ndvi_data, gefs_resolution = 0.25) {
  
  # Create GEFS grid based on data extent
  bbox <- st_bbox(ndvi_data)
  
  # Snap to GEFS grid boundaries
  lon_min <- floor(bbox$xmin / gefs_resolution) * gefs_resolution
  lon_max <- ceiling(bbox$xmax / gefs_resolution) * gefs_resolution
  lat_min <- floor(bbox$ymin / gefs_resolution) * gefs_resolution
  lat_max <- ceiling(bbox$ymax / gefs_resolution) * gefs_resolution
  
  # Create GEFS grid
  gefs_grid <- st_make_grid(
    c(lon_min, lat_min, lon_max, lat_max),
    cellsize = gefs_resolution,
    square = TRUE
  ) %>%
    st_sf(gefs_cell_id = 1:length(.))
  
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
save_processed_ndvi <- function(ndvi_data, output_path = "forecast_workflow/data/processed_ndvi_anomalies.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ndvi_data, output_path)
  cat("Processed NDVI data saved to:", output_path, "\n")
  
}

# Example usage (commented out - will be activated when data is provided)
# ndvi_raw <- load_ndvi_anomalies("path/to/your/ndvi_data.csv")
# ndvi_spatial <- create_spatial_ndvi(ndvi_raw)
# ndvi_gridded <- aggregate_to_grid(ndvi_spatial, grid_resolution = 0.25)
# save_processed_ndvi(ndvi_gridded)

cat("NDVI data housing script loaded. Ready to process anomaly data.\n")
cat("Run load_ndvi_anomalies('your_file_path.csv') to get started.\n")