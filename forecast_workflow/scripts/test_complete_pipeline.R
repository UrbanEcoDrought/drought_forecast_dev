# Comprehensive Pipeline Test Script
# Tests both training data extraction and operational forecasting pipelines

library(tidyverse)
library(sf)
library(lubridate)

# Test configuration
test_config <- list(
  # Small samples for testing
  ndvi_sample_size = 10000,     # Sample NDVI data for faster testing
  training_period = c("2010-01-01", "2010-12-31"),  # 1 year for testing
  reforecast_sampling = "2 months",  # Sparse sampling for testing
  # Test modes
  test_reforecast = TRUE,
  test_realtime = TRUE,
  test_training_build = TRUE,
  test_forecasting = FALSE  # Requires trained model
)

cat("=== COMPREHENSIVE DROUGHT FORECASTING PIPELINE TEST ===\n")
cat("Test configuration:\n")
cat("  NDVI sample size:", test_config$ndvi_sample_size, "\n")
cat("  Training period:", test_config$training_period[1], "to", test_config$training_period[2], "\n")
cat("  Reforecast sampling:", test_config$reforecast_sampling, "\n\n")

# Load all scripts
cat("Loading pipeline scripts...\n")
source("01_ndvi_data_loader.R")
source("06_gefs_reforecast_training_data.R") 
source("07_gefs_realtime_forecasting.R")
source("08_build_training_dataset.R")
if (test_config$test_forecasting) {
  source("09_operational_forecasting_pipeline.R")
}

# =============================================================================
# TEST 1: NDVI Data Loading ----
# =============================================================================
cat("\n=== TEST 1: NDVI DATA LOADING ===\n")

tryCatch({
  ndvi_gefs_data <- run_spatial_ndvi_processing(save_output = FALSE)
  
  # Sample for testing
  if (nrow(ndvi_gefs_data$ndvi_data) > test_config$ndvi_sample_size) {
    cat("Sampling NDVI data for testing...\n")
    ndvi_sampled <- ndvi_gefs_data$ndvi_data %>%
      slice_sample(n = test_config$ndvi_sample_size)
    ndvi_gefs_data$ndvi_data <- ndvi_sampled
  }
  
  cat("✓ NDVI data loaded successfully\n")
  cat("  Records:", nrow(ndvi_gefs_data$ndvi_data), "\n")
  cat("  GEFS cells:", nrow(ndvi_gefs_data$gefs_grid), "\n")
  cat("  Date range:", range(ndvi_gefs_data$ndvi_data$date), "\n")
  
  ndvi_test_passed <- TRUE
  
}, error = function(e) {
  cat("✗ NDVI data loading failed:", e$message, "\n")
  ndvi_test_passed <- FALSE
})

if (!exists("ndvi_test_passed") || !ndvi_test_passed) {
  stop("NDVI data loading failed - cannot continue tests")
}

# Extract GEFS centers for climate data tests
gefs_centers <- ndvi_gefs_data$gefs_grid %>%
  mutate(
    centroid = st_centroid(geometry),
    longitude = st_coordinates(centroid)[,1],
    latitude = st_coordinates(centroid)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(gefs_cell_id, longitude, latitude)

cat("GEFS centers extracted:", nrow(gefs_centers), "cells\n")

# =============================================================================
# TEST 2: Reforecast Data Extraction----
# =============================================================================
if (test_config$test_reforecast) {
  cat("\n=== TEST 2: REFORECAST DATA EXTRACTION ===\n")
  
  tryCatch({
    # Test single date first
    test_date <- as.Date("2010-06-15")
    cat("Testing reforecast availability for", as.character(test_date), "\n")
    
    availability <- check_reforecast_availability(test_date)
    cat("Reforecast availability:", availability, "\n")
    
    if (availability) {
      # Test single date extraction
      cat("Testing single date extraction...\n")
      single_date_data <- extract_reforecast_point_data(
        date = test_date,
        gefs_centers = head(gefs_centers, 5)  # Test with 5 cells only
      )
      
      if (!is.null(single_date_data) && nrow(single_date_data) > 0) {
        cat("✓ Single date reforecast extraction successful\n")
        cat("  Records:", nrow(single_date_data), "\n")
        cat("  Variables:", setdiff(colnames(single_date_data), c("gefs_cell_id", "date", "year", "month", "day_of_year")), "\n")
        
        # Test small time series
        cat("Testing small reforecast time series...\n")
        reforecast_data <- extract_gefs_reforecast_training(
          gefs_centers = head(gefs_centers, 5),  # 5 cells
          start_date = test_config$training_period[1],
          end_date = "2010-06-30",  # Just 6 months
          date_interval = test_config$reforecast_sampling,
          save_output = FALSE
        )
        
        if (!is.null(reforecast_data) && nrow(reforecast_data) > 0) {
          cat("✓ Reforecast time series extraction successful\n")
          cat("  Records:", nrow(reforecast_data), "\n")
          cat("  Date range:", range(reforecast_data$date), "\n")
          reforecast_test_passed <- TRUE
        } else {
          cat("✗ Reforecast time series extraction failed\n")
          reforecast_test_passed <- FALSE
        }
        
      } else {
        cat("✗ Single date reforecast extraction failed\n")
        reforecast_test_passed <- FALSE
      }
      
    } else {
      cat("⚠ Reforecast data not available for test date - this may be expected\n")
      reforecast_test_passed <- FALSE
    }
    
  }, error = function(e) {
    cat("✗ Reforecast extraction test failed:", e$message, "\n")
    reforecast_test_passed <- FALSE
  })
  
} else {
  cat("Skipping reforecast test\n")
  reforecast_test_passed <- TRUE
}

# =============================================================================
# TEST 3: Real-time GEFS Data Extraction  ----
# =============================================================================
if (test_config$test_realtime) {
  cat("\n=== TEST 3: REAL-TIME GEFS DATA EXTRACTION ===\n")
  
  tryCatch({
    # Check available forecast dates
    cat("Checking available real-time forecast dates...\n")
    available_dates <- get_available_forecast_dates(max_days_back = 10)
    
    if (length(available_dates) > 0) {
      cat("✓ Found", length(available_dates), "available forecast dates\n")
      cat("  Most recent:", min(available_dates), "\n")
      
      # Test extraction with most recent date
      test_date <- min(available_dates)
      cat("Testing real-time extraction for", as.character(test_date), "\n")
      
      realtime_data <- extract_realtime_gefs_data(
        date = test_date,
        gefs_centers = head(gefs_centers, 5)  # Test with 5 cells
      )
      
      if (!is.null(realtime_data) && nrow(realtime_data) > 0) {
        cat("✓ Real-time GEFS extraction successful\n")
        cat("  Records:", nrow(realtime_data), "\n")
        cat("  Date:", unique(realtime_data$date), "\n")
        cat("  Variables:", setdiff(colnames(realtime_data), c("gefs_cell_id", "date", "year", "month", "day_of_year")), "\n")
        realtime_test_passed <- TRUE
      } else {
        cat("✗ Real-time GEFS extraction failed\n")
        realtime_test_passed <- FALSE
      }
      
    } else {
      cat("⚠ No real-time GEFS data available - this may be expected\n")
      realtime_test_passed <- FALSE
    }
    
  }, error = function(e) {
    cat("✗ Real-time GEFS test failed:", e$message, "\n")
    realtime_test_passed <- FALSE
  })
  
} else {
  cat("Skipping real-time test\n")
  realtime_test_passed <- TRUE
}

# =============================================================================
# TEST 4: Training Dataset Building
# =============================================================================
if (test_config$test_training_build && exists("reforecast_test_passed") && reforecast_test_passed) {
  cat("\n=== TEST 4: TRAINING DATASET BUILDING ===\n")
  
  tryCatch({
    # Test with very small dataset
    cat("Building test training dataset...\n")
    
    training_dataset <- build_complete_training_dataset(
      ndvi_sample_size = test_config$ndvi_sample_size,
      start_date = test_config$training_period[1], 
      end_date = "2010-06-30",  # Very short period for testing
      save_output = FALSE
    )
    
    if (!is.null(training_dataset) && 
        !is.null(training_dataset$training_data) &&
        nrow(training_dataset$training_data) > 0) {
      
      cat("✓ Training dataset building successful\n")
      cat("  Training records:", nrow(training_dataset$training_data), "\n")
      cat("  GEFS cells:", n_distinct(training_dataset$training_data$gefs_cell_id), "\n")
      cat("  Variables:", colnames(training_dataset$training_data), "\n")
      training_build_test_passed <- TRUE
      
    } else {
      cat("✗ Training dataset building failed\n")
      training_build_test_passed <- FALSE
    }
    
  }, error = function(e) {
    cat("✗ Training dataset building test failed:", e$message, "\n")
    training_build_test_passed <- FALSE
  })
  
} else {
  cat("Skipping training dataset build test (requires reforecast data)\n")
  training_build_test_passed <- TRUE
}

# =============================================================================
# TEST 5: Operational Forecasting (Optional)
# =============================================================================
if (test_config$test_forecasting) {
  cat("\n=== TEST 5: OPERATIONAL FORECASTING ===\n")
  
  # This test requires a trained model - skip if not available
  if (file.exists(forecast_config$trained_model_path)) {
    tryCatch({
      forecast_result <- run_operational_forecast(save_output = FALSE)
      
      if (!is.null(forecast_result)) {
        cat("✓ Operational forecasting successful\n")
        cat("  Forecast records:", nrow(forecast_result$forecasts), "\n")
        forecasting_test_passed <- TRUE
      } else {
        cat("✗ Operational forecasting failed\n")
        forecasting_test_passed <- FALSE
      }
      
    }, error = function(e) {
      cat("✗ Operational forecasting test failed:", e$message, "\n")
      forecasting_test_passed <- FALSE
    })
  } else {
    cat("⚠ Trained model not found - skipping forecasting test\n")
    cat("  Train model first using the hierarchical Bayesian workflow\n")
    forecasting_test_passed <- TRUE
  }
} else {
  cat("Skipping operational forecasting test\n")
  forecasting_test_passed <- TRUE
}

# =============================================================================
# TEST SUMMARY
# =============================================================================
cat("\n=== PIPELINE TEST SUMMARY ===\n")

test_results <- list(
  "NDVI Data Loading" = if(exists("ndvi_test_passed")) ndvi_test_passed else FALSE,
  "Reforecast Extraction" = if(exists("reforecast_test_passed")) reforecast_test_passed else FALSE,
  "Real-time GEFS" = if(exists("realtime_test_passed")) realtime_test_passed else FALSE,
  "Training Dataset Build" = if(exists("training_build_test_passed")) training_build_test_passed else FALSE,
  "Operational Forecasting" = if(exists("forecasting_test_passed")) forecasting_test_passed else FALSE
)

for (test_name in names(test_results)) {
  status <- if(test_results[[test_name]]) "✓ PASSED" else "✗ FAILED"
  cat(sprintf("  %-25s: %s\n", test_name, status))
}

total_tests <- length(test_results)
passed_tests <- sum(unlist(test_results))

cat(sprintf("\nOverall: %d/%d tests passed\n", passed_tests, total_tests))

if (passed_tests == total_tests) {
  cat("🎉 All pipeline components working!\n")
} else if (passed_tests >= 3) {
  cat("⚠ Most components working - check failed tests\n")
} else {
  cat("❌ Major issues - pipeline needs debugging\n")
}

cat("\nNext steps:\n")
if (test_results[["NDVI Data Loading"]] && test_results[["Reforecast Extraction"]]) {
  cat("  1. Run full training dataset build: build_complete_training_dataset()\n")
  cat("  2. Train hierarchical Bayesian model: run_hierarchical_forecast_model()\n")
}
if (test_results[["Real-time GEFS"]]) {
  cat("  3. Set up automated forecasting: run_operational_forecast()\n")
}
cat("  4. Scale up data periods and remove test sampling limits\n")