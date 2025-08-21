# Operational Forecasting Pipeline
# Uses trained model with real-time GEFS data to generate NDVI anomaly forecasts

library(tidyverse)
library(sf)
library(lubridate)
library(brms)

# Load required scripts
source("07_gefs_realtime_forecasting.R")
source("05_hierarchical_bayesian_model.R")

# Configuration for operational forecasting
forecast_config <- list(
  # Model and data paths
  trained_model_path = "U:/datasets/ndvi_monitor/hierarchical_forecast_results.rds",
  training_dataset_path = "U:/datasets/ndvi_monitor/complete_training_dataset.rds",
  # Forecast parameters
  forecast_horizon_days = 30,  # How far ahead to forecast
  update_frequency = "weekly", # How often to update forecasts
  # Climate anomaly calculation
  climate_normal_years = 10,   # Years of historical data for normals
  # Output
  forecast_output_path = "U:/datasets/ndvi_monitor/operational_forecasts/"
)

# Function to load trained model and training data
load_trained_model <- function(model_path = forecast_config$trained_model_path,
                              training_path = forecast_config$training_dataset_path) {
  
  cat("Loading trained model and training data...\n")
  
  # Load model results
  if (file.exists(model_path)) {
    model_results <- readRDS(model_path)
    cat("Loaded trained model from:", model_path, "\n")
  } else {
    stop("Trained model not found at: ", model_path)
  }
  
  # Load training dataset (for climate normals)
  if (file.exists(training_path)) {
    training_data <- readRDS(training_path)
    cat("Loaded training dataset from:", training_path, "\n")
  } else {
    stop("Training dataset not found at: ", training_path)
  }
  
  return(list(
    model = model_results,
    training_data = training_data
  ))
}

# Function to get current climate data for forecasting
get_current_climate_data <- function(gefs_centers, mode = "latest") {
  
  cat("Getting current climate data for forecasting...\n")
  
  # Extract real-time GEFS data
  current_climate <- extract_gefs_realtime_forecast(
    gefs_centers = gefs_centers,
    mode = mode,
    save_output = FALSE
  )
  
  if (is.null(current_climate) || nrow(current_climate) == 0) {
    stop("Failed to retrieve current climate data")
  }
  
  cat("Retrieved current climate data:\n")
  cat("  Date:", unique(current_climate$date), "\n")
  cat("  GEFS cells:", n_distinct(current_climate$gefs_cell_id), "\n")
  
  return(current_climate)
}

# Function to compute climate anomalies for forecasting
compute_forecast_climate_anomalies <- function(current_climate, training_data) {
  
  cat("Computing climate anomalies for current conditions...\n")
  
  # Extract climate normals from training data
  climate_normals <- training_data$training_data %>%
    # Use recent years for more relevant normals
    filter(year >= (max(year) - forecast_config$climate_normal_years)) %>%
    group_by(gefs_cell_id, month) %>%
    summarise(
      temp_normal = mean(temperature, na.rm = TRUE),
      precip_normal = mean(precipitation, na.rm = TRUE),
      solar_normal = mean(solar, na.rm = TRUE),
      humidity_normal = mean(humidity, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate anomalies for current climate
  current_with_anomalies <- current_climate %>%
    left_join(climate_normals, by = c("gefs_cell_id", "month")) %>%
    mutate(
      temp_anomaly = temperature - temp_normal,
      precip_anomaly = precipitation - precip_normal,
      solar_anomaly = solar - solar_normal,
      humidity_anomaly = humidity - humidity_normal
    ) %>%
    # Remove records without normals
    filter(!is.na(temp_normal))
  
  cat("Climate anomalies computed for", nrow(current_with_anomalies), "records\n")
  
  return(current_with_anomalies)
}

# Function to prepare prediction data for model
prepare_forecast_prediction_data <- function(current_climate_anomalies, training_metadata) {
  
  cat("Preparing prediction data for forecasting model...\n")
  
  # Format data to match training dataset structure
  prediction_data <- current_climate_anomalies %>%
    mutate(
      # Add spatial coordinates if missing
      longitude = if(!"longitude" %in% names(.)) NA else longitude,
      latitude = if(!"latitude" %in% names(.)) NA else latitude
    ) %>%
    select(
      gefs_cell_id, date, year, month, day_of_year,
      longitude, latitude,
      temp_anomaly, precip_anomaly, solar_anomaly, humidity_anomaly
    ) %>%
    # Remove any rows with missing predictors
    filter(
      !is.na(temp_anomaly),
      !is.na(precip_anomaly), 
      !is.na(solar_anomaly),
      !is.na(humidity_anomaly)
    )
  
  cat("Prediction dataset prepared:\n")
  cat("  Records:", nrow(prediction_data), "\n")
  cat("  GEFS cells:", n_distinct(prediction_data$gefs_cell_id), "\n")
  
  return(prediction_data)
}

# Function to generate NDVI forecasts
generate_ndvi_forecasts <- function(model_fit, prediction_data) {
  
  cat("Generating NDVI anomaly forecasts...\n")
  
  tryCatch({
    # Generate predictions using the trained model
    forecasts <- predict(
      model_fit$model_fit,
      newdata = prediction_data,
      probs = c(0.05, 0.25, 0.75, 0.95),
      robust = TRUE
    ) %>%
      as_tibble() %>%
      bind_cols(prediction_data) %>%
      rename(
        ndvi_forecast = Estimate,
        forecast_error = Est.Error,
        forecast_lower_90 = Q5,
        forecast_lower_50 = Q25, 
        forecast_upper_50 = Q75,
        forecast_upper_90 = Q95
      )
    
    cat("Forecasts generated for", nrow(forecasts), "locations\n")
    cat("Mean forecast:", round(mean(forecasts$ndvi_forecast), 3), "\n")
    cat("Forecast range:", round(range(forecasts$ndvi_forecast), 3), "\n")
    
    return(forecasts)
    
  }, error = function(e) {
    cat("Error generating forecasts:", e$message, "\n")
    return(NULL)
  })
}

# Function to save operational forecasts
save_operational_forecasts <- function(forecasts, output_dir = forecast_config$forecast_output_path) {
  
  # Create timestamped output directory
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  forecast_dir <- file.path(output_dir, paste0("forecast_", timestamp))
  dir.create(forecast_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save forecast data
  forecast_file <- file.path(forecast_dir, "ndvi_forecasts.rds")
  
  forecast_package <- list(
    forecasts = forecasts,
    metadata = list(
      forecast_date = Sys.Date(),
      forecast_timestamp = Sys.time(),
      climate_source = "GEFS real-time",
      model_source = "Hierarchical Bayesian brms",
      n_locations = nrow(forecasts),
      forecast_horizon = "Current conditions -> NDVI response",
      variables = c("ndvi_forecast", "forecast_error", "forecast_lower_90", 
                   "forecast_lower_50", "forecast_upper_50", "forecast_upper_90"),
      created_by = "operational_forecasting_pipeline.R"
    )
  )
  
  saveRDS(forecast_package, forecast_file)
  
  # Also save as CSV for easy viewing
  csv_file <- file.path(forecast_dir, "ndvi_forecasts.csv")
  write_csv(forecasts, csv_file)
  
  cat("Operational forecasts saved:\n")
  cat("  Directory:", forecast_dir, "\n") 
  cat("  RDS file:", forecast_file, "\n")
  cat("  CSV file:", csv_file, "\n")
  
  return(forecast_package)
}

# Function to create forecast summary
create_forecast_summary <- function(forecasts) {
  
  # Overall summary
  overall_summary <- forecasts %>%
    summarise(
      forecast_date = unique(date),
      n_locations = n(),
      mean_forecast = mean(ndvi_forecast, na.rm = TRUE),
      median_forecast = median(ndvi_forecast, na.rm = TRUE),
      min_forecast = min(ndvi_forecast, na.rm = TRUE),
      max_forecast = max(ndvi_forecast, na.rm = TRUE),
      n_positive_anomalies = sum(ndvi_forecast > 0, na.rm = TRUE),
      n_negative_anomalies = sum(ndvi_forecast < 0, na.rm = TRUE),
      mean_uncertainty = mean(forecast_error, na.rm = TRUE)
    )
  
  # Summary by GEFS cell
  by_cell_summary <- forecasts %>%
    group_by(gefs_cell_id) %>%
    summarise(
      longitude = mean(longitude, na.rm = TRUE),
      latitude = mean(latitude, na.rm = TRUE),
      ndvi_forecast = mean(ndvi_forecast, na.rm = TRUE),
      forecast_uncertainty = mean(forecast_error, na.rm = TRUE),
      drought_risk = case_when(
        ndvi_forecast < -0.5 ~ "High",
        ndvi_forecast < -0.2 ~ "Moderate", 
        ndvi_forecast < 0.2 ~ "Low",
        TRUE ~ "Very Low"
      ),
      .groups = "drop"
    )
  
  return(list(
    overall = overall_summary,
    by_cell = by_cell_summary
  ))
}

# Main operational forecasting workflow
run_operational_forecast <- function(save_output = TRUE, 
                                   climate_mode = "latest") {
  
  cat("Starting operational NDVI drought forecasting...\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # Step 1: Load trained model and training data
  model_data <- load_trained_model()
  
  # Step 2: Extract GEFS centers from training data
  gefs_centers <- model_data$training_data$gefs_grid %>%
    mutate(
      centroid = st_centroid(geometry),
      longitude = st_coordinates(centroid)[,1],
      latitude = st_coordinates(centroid)[,2]
    ) %>%
    st_drop_geometry() %>%
    select(gefs_cell_id, longitude, latitude)
  
  # Step 3: Get current climate data
  current_climate <- get_current_climate_data(gefs_centers, mode = climate_mode)
  
  # Step 4: Compute climate anomalies
  climate_anomalies <- compute_forecast_climate_anomalies(
    current_climate, model_data$training_data
  )
  
  # Step 5: Prepare prediction data
  prediction_data <- prepare_forecast_prediction_data(
    climate_anomalies, model_data$training_data$metadata
  )
  
  # Step 6: Generate forecasts
  forecasts <- generate_ndvi_forecasts(model_data$model, prediction_data)
  
  if (!is.null(forecasts) && nrow(forecasts) > 0) {
    # Step 7: Create summary
    forecast_summary <- create_forecast_summary(forecasts)
    
    # Step 8: Save results
    if (save_output) {
      result <- save_operational_forecasts(forecasts)
      result$summary <- forecast_summary
      return(result)
    } else {
      return(list(
        forecasts = forecasts,
        summary = forecast_summary
      ))
    }
  } else {
    stop("Failed to generate forecasts")
  }
}

cat("Operational forecasting pipeline loaded.\n")
cat("Forecast horizon:", forecast_config$forecast_horizon_days, "days\n")
cat("Update frequency:", forecast_config$update_frequency, "\n")
cat("Output directory:", forecast_config$forecast_output_path, "\n")
cat("Run run_operational_forecast() to generate current NDVI forecasts.\n")