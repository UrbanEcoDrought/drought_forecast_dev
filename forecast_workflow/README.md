# NDVI Drought Forecast Workflow

Modern Bayesian forecasting workflow for predicting NDVI anomalies using GEFS climate data.

## Structure

- `data/` - Raw and processed data storage
- `models/` - Bayesian model specifications  
- `scripts/` - Data processing and modeling scripts
- `output/` - Forecast results and visualizations

## Getting Started

1. Load NDVI anomaly data using `scripts/01_ndvi_data_loader.R`
2. Process GEFS climate data (coming soon)
3. Fit spatial-temporal Bayesian model (coming soon)
4. Generate forecasts with uncertainty (coming soon)

## Recommended Modeling Approach

- **INLA** for fast spatial-temporal modeling
- **brms** for flexible regression structures
- 0.25° grid resolution to match GEFS data