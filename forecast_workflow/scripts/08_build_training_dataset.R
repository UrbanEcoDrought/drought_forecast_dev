# Unified Training Dataset Builder
# Combines NDVI anomaly data with GEFSv12 reforecast climate data for model training

library(tidyverse)
library(sf)
library(lubridate)

# Load required scripts
source("01_ndvi_data_loader.R")
source("06_gefs_reforecast_training_data.R")
source("04_moving_window_climate_anomalies.R")

# Configuration for training dataset
training_config <- list(
  # Training period (overlap of NDVI and reforecast availability)
  start_date = "2001-01-01",
  end_date = "2015-12-31",  # Training data
  validation_start = "2016-01-01", 
  validation_end = "2019-12-31",   # Validation/test data
  # Temporal sampling for training (2-week intervals for better coverage)
  temporal_sampling = "2 weeks",
  # Sample size for development (NULL = use all data)
  ndvi_sample_size = NULL,
  # Climate variables to include (removed humidity due to availability)
  climate_variables = c("temperature", "precipitation", "solar")
)

# Function to load and prepare NDVI data for training
prepare_ndvi_training_data <- function(sample_size = training_config$ndvi_sample_size) {
  
  cat("Loading NDVI data for training...\n")
  
  # Load NDVI data with GEFS grid
  ndvi_gefs_data <- run_spatial_ndvi_processing(save_output = FALSE)
  
  # Sample if requested
  if (!is.null(sample_size) && sample_size < nrow(ndvi_gefs_data$ndvi_data)) {
    cat("Sampling", sample_size, "NDVI observations for development...\n")
    
    ndvi_sampled <- ndvi_gefs_data$ndvi_data %>%
      slice_sample(n = sample_size)
    
    ndvi_gefs_data$ndvi_data <- ndvi_sampled
  }
  
  # Filter to training period
  ndvi_data_processed <- ndvi_gefs_data$ndvi_data %>%
    st_drop_geometry() %>%
    filter(
      date >= as.Date(training_config$start_date),
      date <= as.Date(training_config$end_date)
    )
  
  # Debug: Check available columns
  cat("Available columns after processing:", paste(colnames(ndvi_data_processed), collapse = ", "), "\n")
  
  # Select only available columns
  required_cols <- c("gefs_cell_id", "date", "year", "month", "day_of_year", 
                     "ndvi_anomaly", "ndvi_observed", "longitude", "latitude")
  available_cols <- intersect(required_cols, colnames(ndvi_data_processed))
  missing_cols <- setdiff(required_cols, colnames(ndvi_data_processed))
  
  if (length(missing_cols) > 0) {
    cat("Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(available_cols, collapse = ", "), "\n")
  }
  
  ndvi_training <- ndvi_data_processed %>%
    select(all_of(available_cols))
  
  cat("Filtered NDVI training data:\n")
  cat("  Records:", nrow(ndvi_training), "\n")
  cat("  Period:", range(ndvi_training$date), "\n")
  cat("  GEFS cells:", n_distinct(ndvi_training$gefs_cell_id), "\n")
  
  return(list(
    ndvi_data = ndvi_training,
    gefs_grid = ndvi_gefs_data$gefs_grid
  ))
}

# Function to extract reforecast climate data for training period
extract_training_climate_data <- function(gefs_centers, 
                                         start_date = training_config$start_date,
                                         end_date = training_config$end_date,
                                         sampling = training_config$temporal_sampling) {
  
  cat("Extracting reforecast climate data for training period...\n")
  cat("Period:", start_date, "to", end_date, "\n")
  cat("Sampling:", sampling, "\n")
  
  # Use reforecast extraction function
  climate_data <- extract_gefs_reforecast_training(
    gefs_centers = gefs_centers,
    start_date = start_date,
    end_date = end_date,
    date_interval = sampling,
    save_output = FALSE
  )
  
  return(climate_data)
}

# Function to align NDVI and climate data temporally
align_ndvi_climate_data <- function(ndvi_data, climate_data, 
                                   temporal_tolerance = 7) {  # days
  
  cat("Aligning NDVI and climate data...\n")
  cat("Temporal tolerance:", temporal_tolerance, "days\n")
  
  # For each NDVI observation, find the closest climate data within tolerance
  aligned_data <- ndvi_data %>%
    rowwise() %>%
    mutate(
      # Find closest climate date within tolerance for this GEFS cell
      closest_climate_date = {
        cell_climate <- filter(climate_data, gefs_cell_id == !!gefs_cell_id)
        if (nrow(cell_climate) > 0) {
          date_diffs <- abs(as.numeric(cell_climate$date - !!date))
          min_diff_idx <- which.min(date_diffs)
          if (date_diffs[min_diff_idx] <= temporal_tolerance) {
            cell_climate$date[min_diff_idx]
          } else {
            as.Date(NA)
          }
        } else {
          as.Date(NA)
        }
      }
    ) %>%
    ungroup() %>%
    filter(!is.na(closest_climate_date))
  
  # Join with climate data
  integrated_data <- aligned_data %>%
    left_join(
      climate_data %>%
        select(gefs_cell_id, date, all_of(training_config$climate_variables)),
      by = c("gefs_cell_id", "closest_climate_date" = "date")
    ) %>%
    # Remove records without climate data
    filter(!is.na(temperature)) %>%
    # Clean up
    select(-closest_climate_date)
  
  cat("Aligned training dataset:\n")
  cat("  Records:", nrow(integrated_data), "\n")
  cat("  NDVI-Climate pairs:", sum(!is.na(integrated_data$temperature)), "\n")
  cat("  Period:", range(integrated_data$date), "\n")
  cat("  GEFS cells:", n_distinct(integrated_data$gefs_cell_id), "\n")
  
  return(integrated_data)
}

# Function to compute climate anomalies for training data
compute_training_climate_anomalies <- function(integrated_data) {
  
  cat("Computing climate anomalies for training data...\n")
  
  # Calculate monthly climatologies for each GEFS cell
  climate_normals <- integrated_data %>%
    group_by(gefs_cell_id, month) %>%
    summarise(
      temp_normal = mean(temperature, na.rm = TRUE),
      precip_normal = mean(precipitation, na.rm = TRUE),
      solar_normal = mean(solar, na.rm = TRUE),
      humidity_normal = mean(humidity, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate anomalies
  training_with_anomalies <- integrated_data %>%
    left_join(climate_normals, by = c("gefs_cell_id", "month")) %>%
    mutate(
      temp_anomaly = temperature - temp_normal,
      precip_anomaly = precipitation - precip_normal,
      solar_anomaly = solar - solar_normal,
      humidity_anomaly = humidity - humidity_normal
    )
  
  cat("Climate anomalies computed\n")
  cat("Temperature anomaly range:", round(range(training_with_anomalies$temp_anomaly, na.rm = TRUE), 2), "°C\n")
  cat("Precipitation anomaly range:", round(range(training_with_anomalies$precip_anomaly, na.rm = TRUE), 2), "mm\n")
  
  return(training_with_anomalies)
}

# Function to save complete training dataset with SPEI
save_training_dataset_with_spei <- function(training_data, gefs_grid,
                                           output_path = "U:/datasets/ndvi_monitor/ndvi_spei_training_dataset.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Package with metadata
  training_dataset <- list(
    training_data = training_data,
    gefs_grid = gefs_grid,
    metadata = list(
      description = "Complete NDVI-SPEI training dataset for drought forecasting",
      ndvi_source = "16-day HLS anomalies from Google Drive",
      climate_source = "GEFSv12 reforecast (c00 ensemble)",
      climate_method = "30-day backward temp anomalies + 14-day backward SPEI",
      period = range(training_data$date),
      n_records = nrow(training_data),
      n_cells = n_distinct(training_data$gefs_cell_id),
      variables = list(
        response = "ndvi_anomaly",
        predictors = c("temp_max_anomaly", "spei_14"),
        auxiliary = c("precip", "ref_et"),
        covariates = c("longitude", "latitude", "year", "month", "day_of_year")
      ),
      spatial_resolution = "~4km NDVI pixels nested in 0.25° GEFS cells",
      temporal_resolution = training_config$temporal_sampling,
      window_alignment = "RIGHT (backward-looking, no future data)",
      created_date = Sys.Date()
    )
  )
  
  saveRDS(training_dataset, output_path)
  cat("Complete NDVI-SPEI training dataset saved to:", output_path, "\n")
  
  return(training_dataset)
}

# Function to save complete training dataset (legacy)
save_training_dataset <- function(training_data, gefs_grid,
                                 output_path = "U:/datasets/ndvi_monitor/complete_training_dataset.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Package with metadata
  training_dataset <- list(
    training_data = training_data,
    gefs_grid = gefs_grid,
    metadata = list(
      description = "Complete NDVI-Climate training dataset",
      ndvi_source = "16-day HLS anomalies from Google Drive",
      climate_source = "GEFSv12 reforecast (c00 ensemble)",
      period = range(training_data$date),
      n_records = nrow(training_data),
      n_cells = n_distinct(training_data$gefs_cell_id),
      variables = list(
        response = "ndvi_anomaly",
        predictors = c("temp_anomaly", "precip_anomaly", "solar_anomaly", "humidity_anomaly"),
        covariates = c("longitude", "latitude", "year", "month", "day_of_year")
      ),
      spatial_resolution = "~4km NDVI pixels nested in 0.25° GEFS cells",
      temporal_resolution = training_config$temporal_sampling,
      created_date = Sys.Date()
    )
  )
  
  saveRDS(training_dataset, output_path)
  cat("Complete training dataset saved to:", output_path, "\n")
  
  return(training_dataset)
}

# Main workflow function with SPEI integration
build_complete_training_dataset <- function(ndvi_sample_size = NULL,
                                          start_date = training_config$start_date,
                                          end_date = training_config$end_date,
                                          save_output = TRUE) {
  
  cat("Building complete NDVI-Climate-SPEI training dataset for drought forecasting...\n")
  cat("Training period:", start_date, "to", end_date, "\n")
  cat("Method: GEFS reforecast → 30-day temp anomalies + 14-day SPEI → NDVI anomalies\n")
  
  # Step 1: Load and prepare NDVI data
  ndvi_training <- prepare_ndvi_training_data(sample_size = ndvi_sample_size)
  
  # Step 2: Extract GEFS centers for reforecast extraction
  gefs_centers <- ndvi_training$gefs_grid %>%
    mutate(
      centroid = st_centroid(geometry),
      longitude = st_coordinates(centroid)[,1],
      latitude = st_coordinates(centroid)[,2]
    ) %>%
    st_drop_geometry() %>%
    select(gefs_cell_id, longitude, latitude)
  
  # Step 3: Extract reforecast climate data (temperature, precipitation)
  cat("\n=== EXTRACTING GEFS REFORECAST DATA ===\n")
  reforecast_climate <- extract_gefs_reforecast_production(
    gefs_centers = gefs_centers,
    start_date = start_date,
    end_date = end_date,
    date_interval = training_config$temporal_sampling,
    save_output = FALSE
  )
  
  # Step 4: Calculate moving window climate anomalies and SPEI
  cat("\n=== CALCULATING CLIMATE ANOMALIES AND SPEI ===\n")
  climate_anomalies_result <- calculate_moving_window_anomalies_with_reforecast(
    ndvi_gefs_data = ndvi_training,
    reforecast_climate_data = reforecast_climate,
    start_date = start_date,
    end_date = end_date,
    save_output = FALSE
  )
  
  # Extract the integrated NDVI-climate data
  training_data <- climate_anomalies_result$integrated_data
  
  cat("\n=== FINAL TRAINING DATASET ===\n")
  cat("Records:", nrow(training_data), "\n")
  cat("Variables:", colnames(training_data), "\n")
  cat("Date range:", range(training_data$date), "\n")
  cat("GEFS cells:", n_distinct(training_data$gefs_cell_id), "\n")
  
  if (save_output) {
    # Step 5: Save complete dataset with SPEI
    result <- save_training_dataset_with_spei(training_data, ndvi_training$gefs_grid)
    return(result)
  } else {
    return(list(
      training_data = training_data,
      gefs_grid = ndvi_training$gefs_grid,
      climate_anomalies = climate_anomalies_result$climate_anomalies
    ))
  }
}

cat("Unified training dataset builder loaded.\n")
cat("Training period:", training_config$start_date, "to", training_config$end_date, "\n")
cat("Temporal sampling:", training_config$temporal_sampling, "\n")
cat("Climate variables:", paste(training_config$climate_variables, collapse = ", "), "\n")
cat("Run build_complete_training_dataset() to create integrated NDVI-Climate dataset.\n")