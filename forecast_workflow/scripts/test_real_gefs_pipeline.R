# Test Script for Real GEFS Climate Data Pipeline
# Tests the integration of real NOAA GEFS data with the forecast workflow

library(tidyverse)
library(sf)
library(lubridate)

# Load the NDVI data and spatial structure
cat("Loading NDVI data and spatial structure...\n")
source("01_ndvi_data_loader.R")

# Load NDVI data with GEFS grid overlay from Google Drive
# The data loader script is already configured for Google Drive access
ndvi_gefs_data <- run_spatial_ndvi_processing(save_output = FALSE)

cat("NDVI data loaded:", nrow(ndvi_gefs_data$ndvi_data), "observations\n")
cat("GEFS grid:", nrow(ndvi_gefs_data$gefs_grid), "cells\n")

# Test real GEFS data integration
cat("\nTesting GEFS climate data integration...\n")
source("02_gefs_data_integration.R")

# Try to get real GEFS data from the reforecast period (2000-2019)
test_start <- "2010-01-01"  # Within reforecast period
test_end <- "2010-03-31"    # Within reforecast period

cat("Attempting to retrieve GEFS data for:", test_start, "to", test_end, "\n")

# Extract GEFS coordinates
gefs_centers <- ndvi_gefs_data$gefs_grid %>%
  mutate(
    centroid = st_centroid(geometry),
    longitude = st_coordinates(centroid)[,1],
    latitude = st_coordinates(centroid)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(gefs_cell_id, longitude, latitude)

cat("Testing with", nrow(gefs_centers), "GEFS cells\n")

# Test the real GEFS data retrieval function
tryCatch({
  
  # Test data availability check first (use reforecast period date)
  test_date <- as.Date("2010-01-15")  # Date within reforecast period
  cat("Checking GEFS availability for:", as.character(test_date), "\n")
  
  # Test building GEFS URLs
  test_urls <- build_gefs_urls(test_date, variables = "tmp")
  cat("Built test URLs:", length(test_urls), "variables\n")
  cat("Example URL:", test_urls[[1]], "\n")
  
  # Check URL accessibility
  availability <- check_gefs_availability(test_date)
  cat("GEFS data availability:", availability, "\n")
  
  if (availability) {
    cat("GEFS data appears to be accessible. Proceeding with data extraction...\n")
    
    # Try extracting real GEFS data for just 2 months
    real_climate <- get_real_gefs_data(
      gefs_centers, 
      start_date = test_start,
      end_date = "2024-02-28"  # Just 2 months for testing
    )
    
    if (!is.null(real_climate) && nrow(real_climate) > 0) {
      cat("SUCCESS: Retrieved real GEFS data!\n")
      cat("Climate data shape:", nrow(real_climate), "x", ncol(real_climate), "\n")
      cat("Variables:", paste(colnames(real_climate), collapse = ", "), "\n")
      
      # Show summary
      cat("\nClimate data summary:\n")
      print(summary(real_climate[c("temp_mean", "rh_mean", "precip_total", "solar_mean")]))
      
    } else {
      cat("No real GEFS data retrieved, falling back to synthetic\n")
    }
    
  } else {
    cat("GEFS data not accessible, will use synthetic data\n")
  }
  
}, error = function(e) {
  cat("ERROR in GEFS data retrieval:", e$message, "\n")
  cat("This is expected during development - using synthetic data as fallback\n")
})

cat("\nTest completed. Check output above for results.\n")
cat("If GEFS data was successfully retrieved, the pipeline is ready for real data.\n")
cat("If not, the pipeline will use synthetic data as fallback.\n")