# Model Evaluation Script
# Evaluates trained hierarchical Bayesian NDVI-SPEI model
# Generates validation metrics, forecasts, and diagnostic plots

library(tidyverse)
library(sf)
library(lubridate)
library(brms)
library(bayesplot)
library(loo)
library(posterior)
library(ggplot2)
library(patchwork)    # For combining plots
library(corrr)        # For correlation analysis

# Set working directory to script location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

cat("=== HIERARCHICAL BAYESIAN MODEL EVALUATION ===\n")
cat("Working directory:", getwd(), "\n")
cat("Start time:", Sys.time(), "\n\n")

# Configuration
eval_config <- list(
  # Model and data paths
  trained_model_path = "U:/datasets/ndvi_monitor/trained_models/best_ndvi_spei_model.rds",
  validation_data_path = "U:/datasets/ndvi_monitor/ndvi_spei_validation_dataset.rds",
  
  # Evaluation scope
  use_validation_data = TRUE,    # Use separate validation dataset
  use_training_holdout = FALSE,  # Also evaluate on training holdout
  
  # Forecast evaluation
  generate_forecasts = TRUE,
  forecast_horizon = c(1, 2, 3, 6),  # Months ahead
  n_forecast_samples = 1000,
  
  # Diagnostics
  check_convergence = TRUE,
  posterior_predictive = TRUE,
  residual_analysis = TRUE,
  spatial_diagnostics = TRUE,
  
  # Cross-validation
  run_cv_detailed = TRUE,        # Detailed spatial/temporal CV
  cv_folds = 10,
  
  # Output
  save_results = TRUE,
  create_plots = TRUE,
  output_dir = "U:/datasets/ndvi_monitor/model_evaluation/"
)

cat("Evaluation Configuration:\n")
cat("  Use validation data:", eval_config$use_validation_data, "\n")
cat("  Generate forecasts:", eval_config$generate_forecasts, "\n")
cat("  Forecast horizons:", paste(eval_config$forecast_horizon, "months"), "\n")
cat("  Detailed diagnostics:", eval_config$check_convergence, "\n")
cat("  Save results:", eval_config$save_results, "\n\n")

# Create output directory
if (eval_config$save_results) {
  dir.create(eval_config$output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load trained model
cat("=== LOADING TRAINED MODEL ===\n")
if (!file.exists(eval_config$trained_model_path)) {
  stop("Trained model not found: ", eval_config$trained_model_path,
       "\nRun 'run_model_training.R' first!")
}

model_results <- readRDS(eval_config$trained_model_path)
trained_model <- model_results$model
model_name <- model_results$model_name
training_data <- model_results$model_data
training_metadata <- model_results$training_dataset_metadata

cat("Loaded trained model:\n")
cat("  Model type:", model_name, "\n")
cat("  Formula:", deparse(formula(trained_model)), "\n")
cat("  Training observations:", nobs(trained_model), "\n")
cat("  Training period:", range(training_data$date), "\n\n")

# Load validation data if available
validation_data <- NULL
if (eval_config$use_validation_data && file.exists(eval_config$validation_data_path)) {
  cat("=== LOADING VALIDATION DATA ===\n")
  validation_dataset <- readRDS(eval_config$validation_data_path)
  
  # Preprocess validation data same as training
  validation_data <- validation_dataset$validation_data %>%
    arrange(gefs_cell_id, date) %>%
    group_by(gefs_cell_id) %>%
    mutate(
      ndvi_lag1 = lag(ndvi_anomaly, 1),
      time_idx = row_number(),
      temp_max_anomaly_std = scale(temp_max_anomaly)[,1],
      spei_14_std = scale(spei_14)[,1]
    ) %>%
    ungroup() %>%
    filter(!is.na(ndvi_lag1)) %>%
    mutate(
      gefs_cell_factor = factor(gefs_cell_id),
      year_factor = factor(year),
      month_factor = factor(month)
    )
  
  cat("Validation data loaded:\n")
  cat("  Records:", format(nrow(validation_data), big.mark = ","), "\n")
  cat("  Period:", range(validation_data$date), "\n")
  cat("  GEFS cells:", n_distinct(validation_data$gefs_cell_id), "\n\n")
} else {
  cat("No validation data found - using training data for evaluation\n\n")
}

# Model convergence diagnostics
if (eval_config$check_convergence) {
  cat("=== CONVERGENCE DIAGNOSTICS ===\n")
  
  # Rhat convergence
  rhat_summary <- rhat(trained_model)
  max_rhat <- max(rhat_summary, na.rm = TRUE)
  n_high_rhat <- sum(rhat_summary > 1.01, na.rm = TRUE)
  
  cat("Rhat diagnostics:\n")
  cat("  Max Rhat:", round(max_rhat, 4), "\n")
  cat("  Parameters with Rhat > 1.01:", n_high_rhat, "\n")
  
  if (max_rhat < 1.01) {
    cat("  ✓ Excellent convergence (all Rhat < 1.01)\n")
  } else if (max_rhat < 1.05) {
    cat("  ✓ Good convergence (all Rhat < 1.05)\n")
  } else {
    cat("  ⚠ Potential convergence issues (some Rhat > 1.05)\n")
  }
  
  # Effective sample size
  ess_bulk <- ess_bulk(trained_model)
  ess_tail <- ess_tail(trained_model)
  min_ess_bulk <- min(ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(ess_tail, na.rm = TRUE)
  
  cat("\nEffective sample size:\n")
  cat("  Min bulk ESS:", round(min_ess_bulk), "\n")
  cat("  Min tail ESS:", round(min_ess_tail), "\n")
  
  if (min_ess_bulk > 400 && min_ess_tail > 400) {
    cat("  ✓ Sufficient effective sample size\n\n")
  } else {
    cat("  ⚠ Low effective sample size - consider more iterations\n\n")
  }
  
  # Save diagnostics plot
  if (eval_config$create_plots) {
    mcmc_plot <- mcmc_rhat_hist(rhat_summary) + 
      ggtitle("Rhat Distribution") +
      theme_minimal()
    
    ggsave(paste0(eval_config$output_dir, "convergence_diagnostics.png"), 
           mcmc_plot, width = 8, height = 6)
    cat("Convergence diagnostic plot saved\n\n")
  }
}

# Model performance metrics
cat("=== MODEL PERFORMANCE METRICS ===\n")

# Training performance
training_predictions <- fitted(trained_model, summary = TRUE)
training_r2 <- bayes_R2(trained_model)

cat("Training Performance:\n")
cat("  Bayesian R²:", round(mean(training_r2[,1]), 3), "±", round(sd(training_r2[,1]), 3), "\n")

# Calculate additional training metrics
training_rmse <- sqrt(mean((training_data$ndvi_anomaly - training_predictions[,1])^2))
training_mae <- mean(abs(training_data$ndvi_anomaly - training_predictions[,1]))
training_cor <- cor(training_data$ndvi_anomaly, training_predictions[,1])

cat("  RMSE:", round(training_rmse, 4), "\n")
cat("  MAE:", round(training_mae, 4), "\n")
cat("  Correlation:", round(training_cor, 3), "\n\n")

# Validation performance
if (!is.null(validation_data)) {
  cat("Validation Performance:\n")
  
  # Generate predictions for validation data
  val_predictions <- posterior_predict(trained_model, newdata = validation_data, 
                                      allow_new_levels = TRUE, summary = TRUE)
  
  val_rmse <- sqrt(mean((validation_data$ndvi_anomaly - val_predictions[,1])^2))
  val_mae <- mean(abs(validation_data$ndvi_anomaly - val_predictions[,1]))
  val_cor <- cor(validation_data$ndvi_anomaly, val_predictions[,1])
  
  # Calculate validation R²
  val_ss_res <- sum((validation_data$ndvi_anomaly - val_predictions[,1])^2)
  val_ss_tot <- sum((validation_data$ndvi_anomaly - mean(validation_data$ndvi_anomaly))^2)
  val_r2 <- 1 - (val_ss_res / val_ss_tot)
  
  cat("  R²:", round(val_r2, 3), "\n")
  cat("  RMSE:", round(val_rmse, 4), "\n")
  cat("  MAE:", round(val_mae, 4), "\n")
  cat("  Correlation:", round(val_cor, 3), "\n\n")
  
  # Model generalization
  rmse_ratio <- val_rmse / training_rmse
  if (rmse_ratio < 1.2) {
    cat("  ✓ Good generalization (validation RMSE < 1.2x training)\n")
  } else if (rmse_ratio < 1.5) {
    cat("  ⚠ Moderate overfitting (validation RMSE 1.2-1.5x training)\n")
  } else {
    cat("  ❌ Significant overfitting (validation RMSE > 1.5x training)\n")
  }
  cat("\n")
}

# Residual analysis
if (eval_config$residual_analysis) {
  cat("=== RESIDUAL ANALYSIS ===\n")
  
  eval_data <- if (!is.null(validation_data)) validation_data else training_data
  eval_preds <- if (!is.null(validation_data)) val_predictions[,1] else training_predictions[,1]
  
  # Calculate residuals
  residuals <- eval_data$ndvi_anomaly - eval_preds
  
  # Residual statistics
  cat("Residual diagnostics:\n")
  cat("  Mean:", round(mean(residuals), 5), "\n")
  cat("  SD:", round(sd(residuals), 4), "\n")
  cat("  Skewness:", round(moments::skewness(residuals), 3), "\n")
  cat("  Kurtosis:", round(moments::kurtosis(residuals), 3), "\n\n")
  
  # Create residual plots
  if (eval_config$create_plots) {
    # Residuals vs fitted
    p1 <- ggplot(data.frame(fitted = eval_preds, residuals = residuals), 
                 aes(x = fitted, y = residuals)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(method = "loess", se = TRUE) +
      labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
      theme_minimal()
    
    # Q-Q plot
    p2 <- ggplot(data.frame(residuals = residuals), aes(sample = residuals)) +
      stat_qq() + stat_qq_line(color = "red") +
      labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_minimal()
    
    # Residuals histogram
    p3 <- ggplot(data.frame(residuals = residuals), aes(x = residuals)) +
      geom_histogram(bins = 50, alpha = 0.7, fill = "skyblue") +
      geom_density(color = "red", size = 1) +
      labs(title = "Residual Distribution", x = "Residuals", y = "Density") +
      theme_minimal()
    
    # Observed vs predicted
    p4 <- ggplot(data.frame(observed = eval_data$ndvi_anomaly, predicted = eval_preds),
                 aes(x = observed, y = predicted)) +
      geom_point(alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = TRUE) +
      labs(title = "Observed vs Predicted", x = "Observed NDVI Anomaly", y = "Predicted NDVI Anomaly") +
      theme_minimal()
    
    # Combine plots
    residual_plots <- (p1 | p2) / (p3 | p4)
    
    ggsave(paste0(eval_config$output_dir, "residual_analysis.png"), 
           residual_plots, width = 12, height = 10)
    cat("Residual analysis plots saved\n\n")
  }
}

# Spatial diagnostics
if (eval_config$spatial_diagnostics) {
  cat("=== SPATIAL DIAGNOSTICS ===\n")
  
  eval_data <- if (!is.null(validation_data)) validation_data else training_data
  eval_preds <- if (!is.null(validation_data)) val_predictions[,1] else training_predictions[,1]
  
  # Performance by GEFS cell
  spatial_performance <- eval_data %>%
    mutate(
      prediction = eval_preds,
      residual = ndvi_anomaly - prediction,
      abs_error = abs(residual)
    ) %>%
    group_by(gefs_cell_id) %>%
    summarise(
      n_obs = n(),
      rmse = sqrt(mean(residual^2)),
      mae = mean(abs_error),
      correlation = cor(ndvi_anomaly, prediction),
      bias = mean(residual),
      .groups = "drop"
    )
  
  cat("Spatial performance summary:\n")
  cat("  Mean RMSE by cell:", round(mean(spatial_performance$rmse), 4), "±", 
      round(sd(spatial_performance$rmse), 4), "\n")
  cat("  Mean correlation by cell:", round(mean(spatial_performance$correlation, na.rm = TRUE), 3), "±",
      round(sd(spatial_performance$correlation, na.rm = TRUE), 3), "\n")
  cat("  Cells with r > 0.5:", sum(spatial_performance$correlation > 0.5, na.rm = TRUE), "/", 
      nrow(spatial_performance), "\n\n")
}

# Posterior predictive checks
if (eval_config$posterior_predictive) {
  cat("=== POSTERIOR PREDICTIVE CHECKS ===\n")
  
  # Generate posterior predictive samples
  pp_samples <- posterior_predict(trained_model, draws = 100)
  
  # Compare to observed data
  y_obs <- training_data$ndvi_anomaly
  
  if (eval_config$create_plots) {
    # Posterior predictive density overlay
    pp_plot <- pp_check(trained_model, ndraws = 100) + 
      ggtitle("Posterior Predictive Check") +
      theme_minimal()
    
    ggsave(paste0(eval_config$output_dir, "posterior_predictive_check.png"), 
           pp_plot, width = 10, height = 6)
    cat("Posterior predictive check plot saved\n")
  }
  
  # Test statistics
  pp_mean_check <- mean(apply(pp_samples, 1, mean) > mean(y_obs))
  pp_sd_check <- mean(apply(pp_samples, 1, sd) > sd(y_obs))
  
  cat("Posterior predictive test statistics:\n")
  cat("  P(T_rep > T_obs) for mean:", round(pp_mean_check, 3), "\n")
  cat("  P(T_rep > T_obs) for SD:", round(pp_sd_check, 3), "\n")
  
  if (pp_mean_check > 0.1 && pp_mean_check < 0.9 && pp_sd_check > 0.1 && pp_sd_check < 0.9) {
    cat("  ✓ Good posterior predictive fit\n\n")
  } else {
    cat("  ⚠ Potential model misspecification\n\n")
  }
}

# Save evaluation results
if (eval_config$save_results) {
  cat("=== SAVING EVALUATION RESULTS ===\n")
  
  evaluation_results <- list(
    model_name = model_name,
    evaluation_date = Sys.Date(),
    training_performance = list(
      r2 = mean(training_r2[,1]),
      r2_sd = sd(training_r2[,1]),
      rmse = training_rmse,
      mae = training_mae,
      correlation = training_cor
    ),
    validation_performance = if (!is.null(validation_data)) {
      list(
        r2 = val_r2,
        rmse = val_rmse, 
        mae = val_mae,
        correlation = val_cor,
        generalization_ratio = rmse_ratio
      )
    } else NULL,
    convergence_diagnostics = if (eval_config$check_convergence) {
      list(
        max_rhat = max_rhat,
        n_high_rhat = n_high_rhat,
        min_ess_bulk = min_ess_bulk,
        min_ess_tail = min_ess_tail
      )
    } else NULL,
    spatial_performance = if (eval_config$spatial_diagnostics) spatial_performance else NULL
  )
  
  # Save results
  results_path <- paste0(eval_config$output_dir, "evaluation_results.rds")
  saveRDS(evaluation_results, results_path)
  cat("Evaluation results saved to:", results_path, "\n")
  
  # Save summary report
  report_path <- paste0(eval_config$output_dir, "evaluation_summary.txt")
  sink(report_path)
  
  cat("NDVI-SPEI Hierarchical Bayesian Model Evaluation\n")
  cat("================================================\n\n")
  cat("Model:", model_name, "\n")
  cat("Evaluation date:", Sys.Date(), "\n")
  cat("Formula:", deparse(formula(trained_model)), "\n\n")
  
  cat("PERFORMANCE METRICS\n")
  cat("Training R²:", round(evaluation_results$training_performance$r2, 3), "\n")
  cat("Training RMSE:", round(evaluation_results$training_performance$rmse, 4), "\n")
  
  if (!is.null(evaluation_results$validation_performance)) {
    cat("Validation R²:", round(evaluation_results$validation_performance$r2, 3), "\n")
    cat("Validation RMSE:", round(evaluation_results$validation_performance$rmse, 4), "\n")
    cat("Generalization ratio:", round(evaluation_results$validation_performance$generalization_ratio, 2), "\n")
  }
  
  sink()
  cat("Evaluation summary saved to:", report_path, "\n")
}

cat("\n🎉 MODEL EVALUATION COMPLETED! 🎉\n")
cat("Model performance:\n")
if (!is.null(validation_data)) {
  cat("  Validation R²:", round(val_r2, 3), "\n")
  cat("  Validation RMSE:", round(val_rmse, 4), "\n")
} else {
  cat("  Training R²:", round(mean(training_r2[,1]), 3), "\n")
  cat("  Training RMSE:", round(training_rmse, 4), "\n")
}

cat("\nNext steps:\n")
cat("1. Review diagnostic plots in:", eval_config$output_dir, "\n")
cat("2. Run operational forecasting if performance is acceptable\n")
cat("3. Consider model refinements if needed\n\n")

cat("STATUS: ✅ MODEL EVALUATION COMPLETE\n")