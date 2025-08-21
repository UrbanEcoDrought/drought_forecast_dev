# Build Training Dataset Script
# Creates NDVI-SPEI training dataset using GEFS reforecast data
# Implements: GEFS reforecast → 30-day temp anomalies + 14-day SPEI → NDVI anomalies

library(tidyverse)
library(sf)
library(lubridate)
library(zoo)
library(stars)
library(httr)

# Set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

cat("=== TRAINING DATASET BUILDER ===\n")
cat("Working directory:", getwd(), "\n")
cat("Start time:", Sys.time(), "\n\n")

# Configuration
build_config <- list(
  # Dataset scope
  training_start = "2001-01-01",
  training_end = "2015-12-31",
  
  # Development settings (set ndvi_sample_size = NULL for full dataset)
  ndvi_sample_size = 100000,  # Use 100K for development, NULL for production
  
  # GEFS extraction settings
  temporal_sampling = "2 weeks",
  use_production_extraction = TRUE,  # Use production function with caching
  
  # Output
  save_dataset = TRUE
)

cat("Build Configuration:\n")
cat("  Period:", build_config$training_start, "to", build_config$training_end, "\n")
cat("  NDVI sample:", ifelse(is.null(build_config$ndvi_sample_size), "FULL DATASET", format(build_config$ndvi_sample_size, big.mark = ",")), "\n")
cat("  Sampling interval:", build_config$temporal_sampling, "\n")
cat("  Production extraction:", build_config$use_production_extraction, "\n")
cat("  Save output:", build_config$save_dataset, "\n\n")

# Load required scripts
cat("Loading pipeline scripts...\n")
required_scripts <- c(
  "01_ndvi_data_loader.R",
  "04_moving_window_climate_anomalies.R",
  "06_gefs_reforecast_training_data.R",
  "08_build_training_dataset.R"
)

for (script in required_scripts) {
  cat("  Loading:", script, "\n")
  if (!file.exists(script)) {
    stop("Required script not found: ", script)
  }
  tryCatch({
    source(script)
  }, error = function(e) {
    stop("Failed to load ", script, ": ", e$message)
  })
}
cat("Scripts loaded successfully!\n\n")

# Preliminary checks
cat("=== PRELIMINARY CHECKS ===\n")

# Check if NDVI data file exists (Google Shared Drive)
ndvi_file_path <- "G:/Shared drives/Urban Ecological Drought/data/spatial_NDVI_monitoring/16_day_window_spatial_data_with_anomalies.csv"
if (file.exists(ndvi_file_path)) {
  cat("✓ NDVI data file found on Google Drive\n")
} else {
  cat("✗ NDVI data file not found:", ndvi_file_path, "\n")
  cat("  Please check Google Drive is mounted and accessible\n")
}

# Check if Google Drive is accessible
if (dir.exists("G:/Shared drives/")) {
  cat("✓ Google Shared Drive accessible\n")
} else {
  cat("✗ Google Shared Drive not accessible\n")
  cat("  Please mount Google Drive or check G: drive connection\n")
}

# Check output directory (U: drive for caching/output)
if (dir.exists("U:/")) {
  cat("✓ U: drive accessible for output\n")
} else {
  cat("✗ U: drive not accessible for output\n")
  cat("  Outputs will use local directory instead\n")
}

# Check required packages
required_packages <- c("tidyverse", "sf", "lubridate", "zoo", "stars", "httr")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  cat("✗ Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("  Install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n")
} else {
  cat("✓ All required packages installed\n")
}

cat("Preliminary checks completed\n\n")

# Build training dataset
cat("=== BUILDING TRAINING DATASET ===\n")
build_start_time <- Sys.time()

tryCatch({
  # Test NDVI data access first
  cat("Testing NDVI data access...\n")
  test_ndvi_path <- get_ndvi_data_path()
  if (is.null(test_ndvi_path)) {
    stop("NDVI data file not accessible. Please ensure Google Drive is mounted at G:/Shared drives/")
  }
  cat("✓ NDVI data path confirmed:", test_ndvi_path, "\n")
  
  # Run the complete pipeline
  cat("Calling build_complete_training_dataset()...\n")
  training_result <- build_complete_training_dataset(
    ndvi_sample_size = build_config$ndvi_sample_size,
    start_date = build_config$training_start,
    end_date = build_config$training_end,
    save_output = build_config$save_dataset
  )
  
  build_end_time <- Sys.time()
  build_duration <- difftime(build_end_time, build_start_time, units = "mins")
  
  # Success summary
  cat("\n🎉 TRAINING DATASET BUILD COMPLETED! 🎉\n")
  cat("Build duration:", round(build_duration, 1), "minutes\n\n")
  
  # Dataset summary
  training_data <- training_result$training_data
  cat("Dataset Summary:\n")
  cat("  Total records:", format(nrow(training_data), big.mark = ","), "\n")
  cat("  Date range:", as.character(range(training_data$date)), "\n")
  cat("  GEFS cells:", n_distinct(training_data$gefs_cell_id), "\n")
  cat("  Spatial extent: 0.25° GEFS grid cells\n\n")
  
  # Variable summary
  cat("Variables:\n")
  cat("  Response: ndvi_anomaly\n")
  cat("  Predictors: temp_max_anomaly, spei_14\n")
  cat("  Auxiliary:", setdiff(colnames(training_data), 
                              c("gefs_cell_id", "date", "year", "month", "day_of_year", 
                                "ndvi_anomaly", "temp_max_anomaly", "spei_14")), "\n\n")
  
  # Data quality checks
  cat("Data Quality:\n")
  cat("  NDVI anomalies - Range:", round(range(training_data$ndvi_anomaly, na.rm = TRUE), 3), "\n")
  cat("  Temp anomalies - Range:", round(range(training_data$temp_max_anomaly, na.rm = TRUE), 2), "°C\n")
  cat("  SPEI values - Range:", round(range(training_data$spei_14, na.rm = TRUE), 2), "\n")
  cat("  Missing values:", sum(is.na(training_data$ndvi_anomaly)), "NDVI,", 
      sum(is.na(training_data$temp_max_anomaly)), "Temp,", 
      sum(is.na(training_data$spei_14)), "SPEI\n\n")
  
  # File location
  if (build_config$save_dataset) {
    cat("Output Location:\n")
    cat("  File: U:/datasets/ndvi_monitor/ndvi_spei_training_dataset.rds\n")
    cat("  Load with: readRDS('U:/datasets/ndvi_monitor/ndvi_spei_training_dataset.rds')\n\n")
  }
  
  cat("STATUS: ✅ READY FOR MODEL TRAINING\n")
  cat("Next step: Run model training script\n")
  
}, error = function(e) {
  cat("\n❌ TRAINING DATASET BUILD FAILED\n")
  cat("Full error details:\n")
  cat("  Message:", e$message, "\n")
  cat("  Call:", deparse(e$call), "\n")
  
  # Print the full traceback for debugging
  cat("\nFull traceback:\n")
  traceback()
  
  cat("\nTroubleshooting tips:\n")
  cat("1. Check if NDVI data file exists: U:/datasets/ndvi_monitor/16_day_window_spatial_data_with_anomalies.csv\n")
  cat("2. Check network connection for GEFS data access\n")
  cat("3. Try reducing ndvi_sample_size for testing\n")
  cat("4. Check available disk space for caching\n")
  cat("5. Check if all required R packages are installed\n")
  
  stop("Build failed - see error details above")
})

cat("\n=== BUILD COMPLETED ===\n")
cat("End time:", Sys.time(), "\n")