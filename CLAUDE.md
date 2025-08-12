# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based drought forecasting research project focused on urban ecological systems. The codebase implements Bayesian models using JAGS to analyze NDVI (vegetation health) patterns and their relationship with climate variables like VPD (Vapor Pressure Deficit).

## Key Architecture

### Core Components

1. **NOAA GEFS Data Interface** (`00_NOAA_GEFS_functions.R`)
   - Custom functions for accessing NOAA Global Ensemble Forecast System data
   - Modified to use 0.25-degree resolution data (changed from 0.5-degree)
   - Functions: `grib_extract()`, `gefs_grib_collection()`, `gefs_view()`, `gefs_urls()`
   - Uses gdalcubes for raster data processing and AWS S3 for data access

2. **Time Series Modeling** (`0_random_walk_NDVI.R`, `00_ndvi_vpd_test.R`)
   - Bayesian state-space models implemented in JAGS
   - Random walk models for NDVI forecasting
   - Dynamic linear models incorporating climate drivers (VPD)
   - Uses monthly aggregated NDVI data from Landsat imagery

3. **Visualization** (`xx_gif_creator.R`)
   - Creates animated GIFs from forecast output images using magick package

### Data Flow

1. NDVI data is read from CSV files (Landsat-derived vegetation indices)
2. Climate data (VPD indices) comes from Morton drought indices
3. Data is aggregated to monthly timesteps
4. JAGS models fit Bayesian time series models
5. Forecasts are generated with uncertainty quantification
6. Results can be animated for visualization

## Key Dependencies

- **Bayesian Modeling**: `rjags`, `coda` for MCMC analysis
- **Time Series**: `ecoforecastR` for ecological forecasting utilities  
- **Geospatial**: `gdalcubes`, `sf`, `stars` for raster processing
- **Data**: `tidyverse`, `lubridate`, `readxl` for data manipulation
- **Visualization**: `ggplot2`, `magick` for plots and animations

## Data Sources

- **NDVI**: Landsat 5-9 derived vegetation indices from Google Earth Engine
- **Climate**: NOAA GEFS ensemble weather forecasts (0.25° resolution)
- **VPD**: Morton Arboretum drought indices dataset
- Data paths reference Google Shared Drive: "G:/Shared drives/Urban Ecological Drought"

## Model Structure

The project implements hierarchical Bayesian models:
- **Process Model**: Random walk or autoregressive structure for NDVI evolution
- **Observation Model**: Links observed NDVI to latent state with measurement error
- **Climate Drivers**: VPD index as predictor in dynamic linear models
- **Priors**: Weakly informative priors on regression coefficients and precision parameters

## Development Notes

- Models use JAGS textConnection() syntax for inline model specification
- MCMC convergence assessed with Gelman-Rubin diagnostics
- Forecast uncertainty propagated through ensemble predictions
- File paths are currently hardcoded to specific Google Drive locations
- The codebase appears to be in active development/testing phase

## New Forecast Workflow (In Development)

### Modern Bayesian Approach
- Moving away from JAGS to modern alternatives: **INLA** (recommended for spatial-temporal), **brms**, **Stan**
- Target: Spatially explicit NDVI anomaly forecasting using GEFS climate data
- Architecture: Hierarchical spatial model with nested grid structure

### Spatial Resolution Strategy
- **NDVI data**: Fine resolution (~4km pixels) 
- **GEFS climate**: Coarse resolution (0.25° = ~28km cells)
- **Structure**: Multiple NDVI pixels nested within each GEFS climate cell
- **Modeling**: Hierarchical Bayesian approach handling both spatial scales

### Workflow Status
- ✅ Created `forecast_workflow/` directory structure
- ✅ Built NDVI data loader with nested spatial handling (`scripts/01_ndvi_data_loader.R`)
- ⏸️ Waiting for NDVI anomaly data from colleague
- 🔄 Next: GEFS data integration, model specification, implementation