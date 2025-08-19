# Test script for NDVI data loading
source('forecast_workflow/scripts/01_ndvi_data_loader.R')

# Test data loading
cat("Testing NDVI data loading...\n")
result <- run_spatial_ndvi_processing(save_output = FALSE)

if (!is.null(result)) {
  cat("\nData loading successful!\n")
  cat("NDVI data dimensions:", nrow(result$ndvi_data), "rows\n")
  cat("GEFS grid cells:", nrow(result$gefs_grid), "cells\n")
  
  # Quick summary of data
  data_summary <- result$ndvi_data %>%
    st_drop_geometry() %>%
    summarise(
      years = paste(range(year), collapse = "-"),
      spatial_pixels = n_distinct(paste(longitude, latitude)),
      gefs_cells = n_distinct(gefs_cell_id),
      mean_anomaly = round(mean(ndvi_anomaly, na.rm = TRUE), 4),
      anomaly_range = paste(round(range(ndvi_anomaly, na.rm = TRUE), 4), collapse = " to ")
    )
  
  cat("Summary:\n")
  cat("  Years:", data_summary$years, "\n")
  cat("  NDVI pixels:", data_summary$spatial_pixels, "\n") 
  cat("  GEFS cells:", data_summary$gefs_cells, "\n")
  cat("  Mean anomaly:", data_summary$mean_anomaly, "\n")
  cat("  Anomaly range:", data_summary$anomaly_range, "\n")
}