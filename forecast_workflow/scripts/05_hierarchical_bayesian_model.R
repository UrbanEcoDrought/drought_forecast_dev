# Hierarchical Bayesian NDVI Anomaly Forecast Model
# Uses spatially nested structure with climate forcing (Tmax + SPEI)
# Modern Bayesian approach with brms/Stan backend

library(tidyverse)
library(brms)
library(sf)
library(lubridate)

# Model configuration
model_config <- list(
  response = "ndvi_anomaly",
  climate_predictors = c("temp_max_anomaly", "spei_14"),
  spatial_hierarchy = "gefs_cell_id",
  temporal_structure = "autoregressive",
  priors = "weakly_informative",
  forecast_horizon = 30,  # days ahead
  ensemble_size = 1000    # posterior samples for uncertainty
)

# Function to prepare data for Bayesian modeling
prepare_model_data <- function(integrated_data, 
                              train_end_date = "2022-12-31",
                              test_start_date = "2023-01-01") {
  
  cat("Preparing data for hierarchical Bayesian modeling...\n")
  
  # Convert to proper date format and sort
  model_data <- integrated_data %>%
    mutate(
      date = as.Date(date),
      # Create time index for temporal modeling
      time_index = as.numeric(date - min(date, na.rm = TRUE)) + 1,
      # Standardize predictors for better convergence
      temp_max_anomaly_std = scale(temp_max_anomaly)[,1],
      spei_14_std = scale(spei_14)[,1],
      # Create seasonal components
      cos_doy = cos(2 * pi * day_of_year / 365.25),
      sin_doy = sin(2 * pi * day_of_year / 365.25),
      # Lag variables for autoregressive component
      gefs_cell_id = factor(gefs_cell_id)
    ) %>%
    arrange(gefs_cell_id, date)
  
  # Create lagged NDVI for temporal dependence
  model_data <- model_data %>%
    group_by(gefs_cell_id) %>%
    mutate(
      ndvi_anomaly_lag1 = lag(ndvi_anomaly, 1),
      ndvi_anomaly_lag2 = lag(ndvi_anomaly, 2)
    ) %>%
    ungroup() %>%
    filter(!is.na(ndvi_anomaly_lag1))  # Remove first observation per cell
  
  # Split into training and testing
  train_data <- model_data %>%
    filter(date <= as.Date(train_end_date))
  
  test_data <- model_data %>%
    filter(date >= as.Date(test_start_date))
  
  cat("Training data:", nrow(train_data), "observations\n")
  cat("Testing data:", nrow(test_data), "observations\n")
  cat("GEFS cells:", n_distinct(model_data$gefs_cell_id), "\n")
  cat("Date range:", range(model_data$date), "\n")
  
  return(list(
    full_data = model_data,
    train_data = train_data,
    test_data = test_data,
    model_formula = create_model_formula()
  ))
}

# Function to create hierarchical model formula
create_model_formula <- function() {
  
  # Hierarchical model with climate forcing
  formula_str <- "
    ndvi_anomaly ~ 
      # Fixed effects: climate forcing + seasonal + autoregressive
      temp_max_anomaly_std + 
      spei_14_std + 
      cos_doy + sin_doy +
      ndvi_anomaly_lag1 +
      
      # Random effects: hierarchical spatial structure
      (1 + temp_max_anomaly_std + spei_14_std | gefs_cell_id) +
      
      # Temporal correlation within cells
      ar(time = time_index, gr = gefs_cell_id, p = 1)
  "
  
  model_formula <- as.formula(formula_str)
  
  cat("Model specification:\n")
  cat("Response: NDVI anomalies\n")
  cat("Fixed effects: Climate forcing (Tmax + SPEI) + seasonality + lag\n") 
  cat("Random effects: Hierarchical by GEFS cell\n")
  cat("Temporal: AR(1) within cells\n")
  
  return(model_formula)
}

# Function to set up Bayesian priors
setup_priors <- function(model_data) {
  
  cat("Setting up weakly informative priors...\n")
  
  priors <- c(
    # Fixed effect priors
    prior(normal(0, 0.1), class = Intercept),  # NDVI anomalies centered at 0
    prior(normal(0, 0.05), class = b, coef = temp_max_anomaly_std),  # Climate effects
    prior(normal(0, 0.05), class = b, coef = spei_14_std),
    prior(normal(0, 0.02), class = b, coef = cos_doy),  # Seasonal effects  
    prior(normal(0, 0.02), class = b, coef = sin_doy),
    prior(normal(0.5, 0.2), class = b, coef = ndvi_anomaly_lag1),  # Autoregressive
    
    # Random effect variances
    prior(exponential(10), class = sd),  # Random effect standard deviations
    
    # Residual variance
    prior(exponential(20), class = sigma),  # Model residual standard deviation
    
    # AR parameter
    prior(normal(0, 0.3), class = ar)  # Temporal autocorrelation
  )
  
  return(priors)
}

# Function to fit hierarchical Bayesian model
fit_hierarchical_model <- function(model_data, 
                                  chains = 4, 
                                  iter = 2000,
                                  cores = 4,
                                  save_model = TRUE) {
  
  cat("Fitting hierarchical Bayesian NDVI model...\n")
  cat("Chains:", chains, ", Iterations:", iter, "\n")
  
  # Set up model components
  model_formula <- model_data$model_formula
  train_data <- model_data$train_data
  priors <- setup_priors(train_data)
  
  # Fit model with brms
  model_fit <- brm(
    formula = model_formula,
    data = train_data,
    family = gaussian(),
    prior = priors,
    chains = chains,
    iter = iter,
    cores = cores,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE),
    file = if(save_model) "U:/datasets/ndvi_monitor/hierarchical_ndvi_model" else NULL
  )
  
  cat("Model fitting complete!\n")
  
  return(model_fit)
}

# Function to generate forecasts
generate_forecasts <- function(model_fit, 
                              test_data,
                              forecast_horizon = 30,
                              n_samples = 1000) {
  
  cat("Generating probabilistic forecasts...\n")
  cat("Forecast horizon:", forecast_horizon, "days\n")
  
  # Generate posterior predictions
  forecasts <- posterior_predict(
    model_fit,
    newdata = test_data,
    ndraws = n_samples,
    allow_new_levels = FALSE
  )
  
  # Summarize forecasts
  forecast_summary <- test_data %>%
    mutate(
      forecast_mean = apply(forecasts, 2, mean),
      forecast_sd = apply(forecasts, 2, sd),
      forecast_lower = apply(forecasts, 2, quantile, 0.025),
      forecast_upper = apply(forecasts, 2, quantile, 0.975),
      forecast_lower_50 = apply(forecasts, 2, quantile, 0.25),
      forecast_upper_50 = apply(forecasts, 2, quantile, 0.75)
    )
  
  cat("Generated forecasts for", nrow(forecast_summary), "observations\n")
  
  return(list(
    forecasts = forecasts,
    summary = forecast_summary,
    model = model_fit
  ))
}

# Function to evaluate forecast performance
evaluate_forecasts <- function(forecast_results) {
  
  cat("Evaluating forecast performance...\n")
  
  forecast_summary <- forecast_results$summary
  
  # Calculate performance metrics
  performance <- forecast_summary %>%
    summarise(
      # Point forecast accuracy
      mae = mean(abs(ndvi_anomaly - forecast_mean), na.rm = TRUE),
      rmse = sqrt(mean((ndvi_anomaly - forecast_mean)^2, na.rm = TRUE)),
      bias = mean(forecast_mean - ndvi_anomaly, na.rm = TRUE),
      correlation = cor(ndvi_anomaly, forecast_mean, use = "complete.obs"),
      
      # Probabilistic forecast skill
      coverage_95 = mean(ndvi_anomaly >= forecast_lower & 
                        ndvi_anomaly <= forecast_upper, na.rm = TRUE),
      coverage_50 = mean(ndvi_anomaly >= forecast_lower_50 & 
                        ndvi_anomaly <= forecast_upper_50, na.rm = TRUE),
      
      # Sharpness (forecast precision)
      sharpness_95 = mean(forecast_upper - forecast_lower, na.rm = TRUE),
      sharpness_50 = mean(forecast_upper_50 - forecast_lower_50, na.rm = TRUE)
    )
  
  # Performance by GEFS cell
  performance_by_cell <- forecast_summary %>%
    group_by(gefs_cell_id) %>%
    summarise(
      n_obs = n(),
      mae = mean(abs(ndvi_anomaly - forecast_mean), na.rm = TRUE),
      rmse = sqrt(mean((ndvi_anomaly - forecast_mean)^2, na.rm = TRUE)),
      correlation = cor(ndvi_anomaly, forecast_mean, use = "complete.obs"),
      coverage_95 = mean(ndvi_anomaly >= forecast_lower & 
                        ndvi_anomaly <= forecast_upper, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Overall Performance:\n")
  print(performance)
  cat("\nPerformance by GEFS cell:\n")
  print(performance_by_cell)
  
  return(list(
    overall = performance,
    by_cell = performance_by_cell
  ))
}

# Function to save model results
save_model_results <- function(results, output_path = "U:/datasets/ndvi_monitor/hierarchical_forecast_results.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, output_path)
  cat("Model results saved to:", output_path, "\n")
  
}

# Main workflow function
run_hierarchical_forecast_model <- function(integrated_data,
                                           train_end = "2022-12-31", 
                                           fit_model = TRUE,
                                           generate_forecast = TRUE,
                                           save_output = TRUE) {
  
  cat("Starting hierarchical Bayesian NDVI forecast workflow...\n")
  
  # Prepare data
  model_data <- prepare_model_data(integrated_data, train_end_date = train_end)
  
  if (fit_model) {
    # Fit model
    model_fit <- fit_hierarchical_model(model_data)
    
    if (generate_forecast) {
      # Generate forecasts
      forecasts <- generate_forecasts(model_fit, model_data$test_data)
      
      # Evaluate performance
      performance <- evaluate_forecasts(forecasts)
      
      results <- list(
        model_data = model_data,
        model_fit = model_fit,
        forecasts = forecasts,
        performance = performance
      )
      
      if (save_output) {
        save_model_results(results)
      }
      
      return(results)
    } else {
      return(list(
        model_data = model_data,
        model_fit = model_fit
      ))
    }
  } else {
    return(list(
      model_data = model_data
    ))
  }
}

cat("Hierarchical Bayesian NDVI forecast model loaded.\n")
cat("Configuration:\n")
cat("  - Response: NDVI anomalies\n")
cat("  - Climate forcing: 30-day Tmax anomalies + 14-day SPEI\n")
cat("  - Spatial hierarchy: NDVI pixels nested in GEFS cells\n")
cat("  - Temporal: AR(1) + seasonality + lag structure\n")
cat("  - Backend: brms/Stan for full Bayesian inference\n")
cat("Run run_hierarchical_forecast_model(integrated_data) to start modeling.\n")