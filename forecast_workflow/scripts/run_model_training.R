# Model Training Script
# Trains hierarchical Bayesian model for NDVI-SPEI drought forecasting
# Uses brms/Stan for modern Bayesian inference with spatial hierarchy

library(tidyverse)
library(sf)
library(lubridate)
library(brms)        # Bayesian regression models using Stan
library(bayesplot)   # Bayesian model visualization
library(loo)         # Leave-one-out cross-validation
library(posterior)   # Posterior analysis

# Set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

cat("=== HIERARCHICAL BAYESIAN MODEL TRAINING ===\n")
cat("Working directory:", getwd(), "\n")
cat("Start time:", Sys.time(), "\n\n")

# Configuration
model_config <- list(
  # Data source
  training_data_path = "U:/datasets/ndvi_monitor/ndvi_spei_training_dataset.rds",
  
  # Model specification
  response_var = "ndvi_anomaly",
  predictors = c("temp_max_anomaly", "spei_14"),
  spatial_grouping = "gefs_cell_id",
  temporal_vars = c("year", "month"),
  
  # Model complexity
  include_spatial_smooth = FALSE,  # Add spatial GP smoother
  include_temporal_ar = TRUE,      # Add AR(1) temporal structure
  include_interactions = FALSE,    # Include predictor interactions
  
  # Stan/brms settings
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  
  # Model selection
  fit_multiple_models = TRUE,      # Fit several model variants
  use_cross_validation = TRUE,     # Evaluate with LOO-CV
  
  # Output
  save_models = TRUE,
  save_predictions = TRUE
)

cat("Model Configuration:\n")
cat("  Response variable:", model_config$response_var, "\n")
cat("  Predictors:", paste(model_config$predictors, collapse = ", "), "\n")
cat("  Spatial grouping:", model_config$spatial_grouping, "\n")
cat("  Include spatial smoother:", model_config$include_spatial_smooth, "\n")
cat("  Include temporal AR(1):", model_config$include_temporal_ar, "\n")
cat("  Stan chains:", model_config$chains, "x", model_config$iter, "iterations\n")
cat("  Fit multiple models:", model_config$fit_multiple_models, "\n\n")

# Load training dataset
cat("=== LOADING TRAINING DATA ===\n")
if (!file.exists(model_config$training_data_path)) {
  stop("Training dataset not found: ", model_config$training_data_path, 
       "\nRun 'run_build_training_data.R' first!")
}

training_dataset <- readRDS(model_config$training_data_path)
training_data <- training_dataset$training_data
gefs_grid <- training_dataset$gefs_grid

cat("Loaded training dataset:\n")
cat("  Records:", format(nrow(training_data), big.mark = ","), "\n")
cat("  Period:", range(training_data$date), "\n")
cat("  GEFS cells:", n_distinct(training_data$gefs_cell_id), "\n")
cat("  Variables:", paste(colnames(training_data), collapse = ", "), "\n\n")

# Data preprocessing for modeling
cat("=== DATA PREPROCESSING ===\n")

# Check for required variables
required_vars <- c(model_config$response_var, model_config$predictors, 
                   model_config$spatial_grouping, model_config$temporal_vars)
missing_vars <- setdiff(required_vars, colnames(training_data))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
}

# Prepare modeling dataset
model_data <- training_data %>%
  # Filter complete cases
  filter(complete.cases(select(., all_of(required_vars)))) %>%
  # Add temporal structure
  arrange(gefs_cell_id, date) %>%
  group_by(gefs_cell_id) %>%
  mutate(
    # Lagged NDVI for AR structure
    ndvi_lag1 = lag(ndvi_anomaly, 1),
    # Time index within each cell
    time_idx = row_number(),
    # Note: temp_max_anomaly and spei_14 are already standardized/anomalies from data pipeline
    # No additional standardization needed - preserves interpretable drought thresholds
  ) %>%
  ungroup() %>%
  # Add spatial coordinates
  left_join(
    gefs_grid %>%
      st_drop_geometry() %>%
      mutate(
        centroid = st_centroid(gefs_grid$geometry),
        longitude = st_coordinates(centroid)[,1],
        latitude = st_coordinates(centroid)[,2]
      ) %>%
      select(gefs_cell_id, longitude, latitude),
    by = "gefs_cell_id"
  ) %>%
  # Remove incomplete cases after lagging
  filter(!is.na(ndvi_lag1)) %>%
  # Create cell factor for random effects
  mutate(
    gefs_cell_factor = factor(gefs_cell_id),
    year_factor = factor(year),
    month_factor = factor(month)
  )

cat("Preprocessed data:\n")
cat("  Records after preprocessing:", format(nrow(model_data), big.mark = ","), "\n")
cat("  Complete cases:", sum(complete.cases(model_data)), "\n")
cat("  Cells with data:", n_distinct(model_data$gefs_cell_id), "\n\n")

# Define model formulas
cat("=== DEFINING MODEL FORMULAS ===\n")

# Base formula - using original predictors (already standardized in pipeline)
base_formula <- paste(model_config$response_var, "~", 
                      paste(model_config$predictors, collapse = " + "))

# Model variants to fit
model_formulas <- list(
  # Model 1: Basic linear model
  linear = as.formula(paste(base_formula, "+ (1|gefs_cell_factor)")),
  
  # Model 2: With temporal AR(1)
  ar1 = as.formula(paste(base_formula, "+ ndvi_lag1 + (1|gefs_cell_factor)")),
  
  # Model 3: With random slopes
  random_slopes = as.formula(paste(base_formula, "+ (temp_max_anomaly + spei_14|gefs_cell_factor)")),
  
  # Model 4: AR(1) + random slopes  
  full = as.formula(paste(base_formula, "+ ndvi_lag1 + (temp_max_anomaly + spei_14|gefs_cell_factor)"))
)

# Add spatial smoother if requested
if (model_config$include_spatial_smooth) {
  model_formulas$spatial = update(model_formulas$full, . ~ . + gp(longitude, latitude))
}

# Select models to fit
if (model_config$fit_multiple_models) {
  models_to_fit <- model_formulas
} else {
  models_to_fit <- list(main = model_formulas$ar1)  # Default to AR(1) model
}

cat("Models to fit:\n")
for (i in seq_along(models_to_fit)) {
  cat("  ", names(models_to_fit)[i], ":", deparse(models_to_fit[[i]]), "\n")
}
cat("\n")

# Fit models
cat("=== FITTING BAYESIAN MODELS ===\n")
fitted_models <- list()
model_summaries <- list()
fit_start_time <- Sys.time()

for (model_name in names(models_to_fit)) {
  cat("Fitting model:", model_name, "\n")
  cat("Formula:", deparse(models_to_fit[[model_name]]), "\n")
  
  model_start_time <- Sys.time()
  
  tryCatch({
    # Fit model with brms
    fit <- brm(
      formula = models_to_fit[[model_name]],
      data = model_data,
      family = gaussian(),
      prior = c(
        prior(normal(0, 0.5), class = Intercept),
        prior(normal(0, 0.25), class = b),
        prior(exponential(2), class = sd),
        prior(exponential(1), class = sigma)
      ),
      chains = model_config$chains,
      iter = model_config$iter,
      warmup = model_config$warmup,
      cores = model_config$cores,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      refresh = 200,  # Print progress every 200 iterations
      seed = 12345
    )
    
    model_end_time <- Sys.time()
    model_duration <- difftime(model_end_time, model_start_time, units = "mins")
    
    # Store results
    fitted_models[[model_name]] <- fit
    
    # Create summary
    model_summary <- list(
      model_name = model_name,
      formula = models_to_fit[[model_name]],
      fit_time = model_duration,
      n_obs = nobs(fit),
      r_squared = bayes_R2(fit),
      loo_cv = loo(fit)
    )
    model_summaries[[model_name]] <- model_summary
    
    cat("✓ Model", model_name, "completed in", round(model_duration, 1), "minutes\n")
    cat("  Bayesian R²:", round(mean(model_summary$r_squared[,1]), 3), 
        "(", round(sd(model_summary$r_squared[,1]), 3), ")\n")
    cat("  LOO-CV ELPD:", round(model_summary$loo_cv$estimates["elpd_loo", "Estimate"], 1), "\n\n")
    
  }, error = function(e) {
    cat("✗ Model", model_name, "failed:", e$message, "\n\n")
    fitted_models[[model_name]] <- NULL
  })
}

fit_end_time <- Sys.time()
total_fit_duration <- difftime(fit_end_time, fit_start_time, units = "mins")

# Model comparison
cat("=== MODEL COMPARISON ===\n")
successful_models <- fitted_models[!sapply(fitted_models, is.null)]

if (length(successful_models) > 1) {
  # LOO model comparison
  loo_compare_results <- loo_compare(lapply(model_summaries[names(successful_models)], 
                                           function(x) x$loo_cv))
  
  cat("LOO Cross-validation comparison:\n")
  print(loo_compare_results)
  cat("\n")
  
  # Best model
  best_model_name <- rownames(loo_compare_results)[1]
  best_model <- fitted_models[[best_model_name]]
  
  cat("🏆 Best model:", best_model_name, "\n")
  
} else if (length(successful_models) == 1) {
  best_model_name <- names(successful_models)[1]
  best_model <- successful_models[[1]]
  cat("Single model fitted:", best_model_name, "\n")
} else {
  stop("No models fitted successfully!")
}

# Save results
if (model_config$save_models) {
  cat("=== SAVING RESULTS ===\n")
  
  # Create output directory
  output_dir <- "U:/datasets/ndvi_monitor/trained_models/"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save best model
  best_model_path <- paste0(output_dir, "best_ndvi_spei_model.rds")
  saveRDS(list(
    model = best_model,
    model_name = best_model_name,
    model_data = model_data,
    config = model_config,
    training_dataset_metadata = training_dataset$metadata,
    fit_summary = model_summaries[[best_model_name]]
  ), best_model_path)
  
  cat("Best model saved to:", best_model_path, "\n")
  
  # Save all models if multiple fitted
  if (length(successful_models) > 1) {
    all_models_path <- paste0(output_dir, "all_ndvi_spei_models.rds")
    saveRDS(list(
      models = fitted_models,
      summaries = model_summaries,
      comparison = if(exists("loo_compare_results")) loo_compare_results else NULL,
      config = model_config
    ), all_models_path)
    
    cat("All models saved to:", all_models_path, "\n")
  }
}

# Final summary
cat("\n🎉 MODEL TRAINING COMPLETED! 🎉\n")
cat("Total training time:", round(total_fit_duration, 1), "minutes\n")
cat("Models fitted:", length(successful_models), "/", length(models_to_fit), "\n")
if (exists("best_model_name")) {
  cat("Best model:", best_model_name, "\n")
}

cat("\nNext steps:\n")
cat("1. Run model evaluation: run_model_evaluation.R\n")
cat("2. Generate forecasts: run_operational_forecast.R\n")
cat("3. Create visualizations and reports\n\n")

cat("STATUS: ✅ READY FOR MODEL EVALUATION\n")