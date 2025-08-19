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

### Workflow Status (Updated 2025-01-19)
- ✅ Created `forecast_workflow/` directory structure
- ✅ Built NDVI data loader with nested spatial handling (`scripts/01_ndvi_data_loader.R`)
- ✅ Successfully loaded actual NDVI anomaly data (4.5M observations, 2001-2025)
- ✅ Created GEFS grid overlay (25 active cells, 0.25° resolution)
- ✅ Implemented moving window climate anomalies (`scripts/04_moving_window_climate_anomalies.R`)
- ✅ Built hierarchical Bayesian model framework (`scripts/05_hierarchical_bayesian_model.R`)
- ✅ Successfully tested brms/Stan model with 10K subset (correlation = 0.97)

### Current Status: Ready for Full-Scale Modeling

**Data Pipeline Complete:**
1. **NDVI Anomaly Data**: 4,566,052 observations from `16_day_window_spatial_data_with_anomalies.csv`
2. **Spatial Structure**: 25 GEFS cells (0.25°) containing multiple NDVI pixels (~4km)
3. **Climate Forcing**: 30-day backward Tmax anomalies + 14-day backward SPEI
4. **Integration**: Climate and NDVI data successfully merged at high temporal resolution

**Model Architecture Validated:**
- **Response**: NDVI anomalies (maintaining 16-day HLS resolution)
- **Predictors**: Standardized 30-day Tmax anomalies + 14-day SPEI (colleague's proven method)
- **Spatial Hierarchy**: Random effects by GEFS cell 
- **Temporal Structure**: AR(1) + lagged NDVI + optional spatial smoothers
- **Backend**: brms/Stan (modern Bayesian inference)

**Test Results:**
- Model converges successfully with Rtools 4.4 installation
- Strong predictive performance (r = 0.97 on 10K subset)
- Climate effects detectable even with synthetic data
- Ready to scale to full 4.1M observation dataset

**Next Session Tasks:**
1. Fit full hierarchical model on complete dataset (4.1M observations)
2. Generate probabilistic forecasts with uncertainty quantification
3. Evaluate spatial and temporal forecast performance
4. Consider adding spatial smoothers if needed: `s(longitude, latitude, bs = "gp")`
5. Replace synthetic climate with real GEFS/ERA5 data for production use

**File Structure:**
```
forecast_workflow/
├── scripts/
│   ├── 01_ndvi_data_loader.R              # ✅ NDVI spatial data processing
│   ├── 02_gefs_data_integration.R         # ✅ Synthetic climate (needs real GEFS)
│   ├── 03_climate_normals_1980_2010.R     # ⏸️ ERA5 normals (auth issues)
│   ├── 04_moving_window_climate_anomalies.R # ✅ Backward-looking anomalies
│   └── 05_hierarchical_bayesian_model.R   # ✅ brms/Stan model framework
└── data/                                   # Uses network drive: U:/datasets/ndvi_monitor/
```

**Key Decisions Made:**
- Using colleague's standardization method (30-day Tmax + 14-day SPEI) instead of ERA5 normals
- All moving windows are RIGHT-ALIGNED (backward-looking) to avoid data leakage
- Preserving high temporal resolution (16-day) rather than monthly aggregation
- Modern Bayesian approach (brms/Stan) instead of JAGS
- Hierarchical spatial structure (pixels nested in climate cells)