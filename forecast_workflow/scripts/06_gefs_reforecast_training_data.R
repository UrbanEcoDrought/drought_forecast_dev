# GEFS Reforecast Data Extraction for Training
# Extracts GEFSv12 reforecast data (2000-2019) to build training datasets
# Uses direct GRIB2 file access from noaa-gefs-retrospective S3 bucket

library(tidyverse)
library(sf)
library(lubridate)
library(stars)
library(httr)

# Configuration for GEFSv12 reforecast
reforecast_config <- list(
  base_url = "https://noaa-gefs-retrospective.s3.amazonaws.com",
  period = c("2000-01-01", "2019-12-31"),
  ensemble = "c00",  # Control run only
  resolution = 0.25, # degrees
  variables = list(
    temperature = "tmp_2m",
    precipitation = "apcp_sfc", 
    solar = "dswrf_sfc"
    # Note: humidity (rh_2m) removed due to limited availability
  ),
  forecast_day = 1,  # Use day-1 forecast (closest to observations)
  # Production settings
  batch_size = 50,   # Process dates in batches
  retry_attempts = 3, # Retry failed extractions
  cache_extractions = TRUE  # Cache successful extractions
)

# Function to build reforecast URLs for a specific date and variable
build_reforecast_url <- function(date, variable, ensemble = "c00") {
  
  date_obj <- as.Date(date)
  year <- format(date_obj, "%Y")
  date_hour <- paste0(format(date_obj, "%Y%m%d"), "00")
  
  # GEFSv12 reforecast path structure
  url <- paste0(
    reforecast_config$base_url,
    "/GEFSv12/reforecast/", year, "/", date_hour, "/",
    ensemble, "/Days:1-10/", variable, "_", date_hour, "_", ensemble, ".grib2"
  )
  
  # Debug URL construction
  cat("    Built URL:", url, "\n")
  
  return(url)
}

# Function to check reforecast data availability
check_reforecast_availability <- function(date, variable = "tmp_2m") {
  
  url <- build_reforecast_url(date, variable)
  
  tryCatch({
    response <- httr::HEAD(url)
    available <- httr::status_code(response) == 200
    return(available)
  }, error = function(e) {
    return(FALSE)
  })
}

# Function to extract reforecast data at specific points
extract_reforecast_point_data <- function(date, gefs_centers, variables = c("temperature", "precipitation", "solar")) {
  
  cat("Extracting reforecast data for", as.character(date), "\n")
  
  # Convert GEFS centers to spatial points
  gefs_points <- gefs_centers %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = "EPSG:4326")
  
  climate_data <- list()
  
  for (var in variables) {
    var_name <- reforecast_config$variables[[var]]
    url <- build_reforecast_url(date, var_name)
    
    cat("  Processing", var, "from", basename(url), "\n")
    
    tryCatch({
      # Read GRIB2 file directly using stars
      grib_url <- paste0("/vsicurl/", url)
      raster_data <- read_stars(grib_url)
      
      # Ensure CRS compatibility
      gefs_points_transformed <- st_transform(gefs_points, st_crs(raster_data))
      
      # Extract values at GEFS cell locations - ensure point extraction only
      extracted_values <- st_extract(raster_data, gefs_points_transformed, exact = TRUE)
      
      # Convert to data frame and ensure proper alignment
      values_df <- data.frame(
        gefs_cell_id = gefs_centers$gefs_cell_id
      )
      
      # Extract the actual values - handle different return formats
      extracted_df <- as.data.frame(extracted_values)
      
      # Find the data column (typically first non-geometry column)
      data_columns <- extracted_df[, !sapply(extracted_df, function(x) inherits(x, "sfc")), drop = FALSE]
      
      if (ncol(data_columns) > 0) {
        data_values <- data_columns[, 1]  # Take first data column
        if (length(data_values) == nrow(gefs_centers)) {
          values_df[[var]] <- data_values
        } else {
          # If we got too many values, take first N matching our points
          if (length(data_values) > nrow(gefs_centers)) {
            values_df[[var]] <- data_values[1:nrow(gefs_centers)]
            cat("    Warning: extracted", length(data_values), "values, using first", nrow(gefs_centers), "\n")
          } else {
            # If too few, pad with NAs
            values_df[[var]] <- c(data_values, rep(NA, nrow(gefs_centers) - length(data_values)))
            cat("    Warning: extracted only", length(data_values), "values, padding to", nrow(gefs_centers), "\n")
          }
        }
      } else {
        values_df[[var]] <- rep(NA, nrow(gefs_centers))
        cat("    Warning: no data columns found in extraction result\n")
      }
      
      climate_data[[var]] <- values_df
      cat("    Successfully extracted", nrow(values_df), "values\n")
      
    }, error = function(e) {
      cat("    Error extracting", var, ":", e$message, "\n")
      # Create NA values as fallback
      na_df <- data.frame(
        gefs_cell_id = gefs_centers$gefs_cell_id
      )
      na_df[[var]] <- rep(NA, nrow(gefs_centers))
      climate_data[[var]] <- na_df
    })
  }
  
  # Combine all variables - handle case where no data was extracted
  if (length(climate_data) == 0) {
    cat("  No climate data extracted successfully\n")
    return(NULL)
  }
  
  # Start with first non-empty dataset
  valid_data <- climate_data[lengths(climate_data) > 0]
  if (length(valid_data) == 0) {
    cat("  All climate extractions failed\n")
    return(NULL)
  }
  
  result <- valid_data %>%
    reduce(full_join, by = "gefs_cell_id") %>%
    mutate(
      date = date,
      year = year(date),
      month = month(date),
      day_of_year = yday(date)
    ) %>%
    relocate(gefs_cell_id, date, year, month, day_of_year)
  
  return(result)
}

# Function to build training dataset over date range
build_reforecast_training_data <- function(gefs_centers, 
                                          start_date = "2010-01-01", 
                                          end_date = "2010-12-31",
                                          date_interval = "month") {
  
  cat("Building reforecast training dataset...\n")
  cat("Period:", start_date, "to", end_date, "\n")
  cat("Sampling interval:", date_interval, "\n")
  
  # Create date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = date_interval)
  cat("Processing", length(dates), "dates\n")
  
  # Extract data for each date
  training_data <- list()
  
  for (i in seq_along(dates)) {
    date <- dates[i]
    
    cat(paste0("[", i, "/", length(dates), "] "))
    
    # Check availability first
    if (check_reforecast_availability(date)) {
      climate_data <- extract_reforecast_point_data(date, gefs_centers)
      training_data[[i]] <- climate_data
    } else {
      cat("Reforecast data not available for", as.character(date), "\n")
    }
    
    # Small delay to avoid overwhelming the server
    Sys.sleep(0.5)
  }
  
  # Combine all dates
  if (length(training_data) > 0) {
    full_dataset <- bind_rows(training_data)
    
    cat("\nTraining dataset completed:\n")
    cat("  Records:", nrow(full_dataset), "\n")
    cat("  Date range:", as.character(range(full_dataset$date)), "\n")
    cat("  GEFS cells:", n_distinct(full_dataset$gefs_cell_id), "\n")
    cat("  Variables:", setdiff(colnames(full_dataset), c("gefs_cell_id", "date", "year", "month", "day_of_year")), "\n")
    
    return(full_dataset)
  } else {
    stop("No reforecast data was successfully extracted")
  }
}

# Function to save training data
save_reforecast_training_data <- function(training_data, 
                                         output_path = "U:/datasets/ndvi_monitor/gefs_reforecast_training_data.rds") {
  
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  
  # Add metadata
  training_with_metadata <- list(
    climate_data = training_data,
    metadata = list(
      source = "GEFSv12 Reforecast",
      ensemble = reforecast_config$ensemble,
      variables = reforecast_config$variables,
      period = range(training_data$date),
      n_records = nrow(training_data),
      n_cells = n_distinct(training_data$gefs_cell_id),
      created_date = Sys.Date(),
      url_base = reforecast_config$base_url
    )
  )
  
  saveRDS(training_with_metadata, output_path)
  cat("Training data saved to:", output_path, "\n")
  
  return(training_with_metadata)
}

# Main workflow function
extract_gefs_reforecast_training <- function(gefs_centers,
                                            start_date = "2005-01-01",
                                            end_date = "2015-12-31", 
                                            date_interval = "month",
                                            save_output = TRUE) {
  
  cat("Starting GEFSv12 reforecast extraction for training...\n")
  cat("Data source:", reforecast_config$base_url, "\n")
  cat("Period:", reforecast_config$period[1], "to", reforecast_config$period[2], "\n")
  
  # Build training dataset
  training_data <- build_reforecast_training_data(
    gefs_centers = gefs_centers,
    start_date = start_date,
    end_date = end_date,
    date_interval = date_interval
  )
  
  if (save_output) {
    result <- save_reforecast_training_data(training_data)
    return(result)
  } else {
    return(training_data)
  }
}

# Production-scale reforecast extraction with batching and caching
extract_gefs_reforecast_production <- function(gefs_centers,
                                              start_date = "2001-01-01",
                                              end_date = "2015-12-31", 
                                              date_interval = "2 weeks",
                                              save_output = TRUE,
                                              cache_dir = "U:/datasets/ndvi_monitor/reforecast_cache/") {
  
  cat("=== PRODUCTION-SCALE GEFSv12 REFORECAST EXTRACTION ===\n")
  cat("Data source:", reforecast_config$base_url, "\n")
  cat("Training period:", start_date, "to", end_date, "\n")
  cat("Variables:", paste(names(reforecast_config$variables), collapse = ", "), "\n")
  
  # Create cache directory
  if (reforecast_config$cache_extractions && !dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Created cache directory:", cache_dir, "\n")
  }
  
  # Create date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = date_interval)
  total_dates <- length(dates)
  cat("Total dates to process:", total_dates, "\n")
  
  # Process in batches
  batch_size <- reforecast_config$batch_size
  n_batches <- ceiling(total_dates / batch_size)
  cat("Processing in", n_batches, "batches of", batch_size, "dates each\n\n")
  
  all_data <- list()
  successful_extractions <- 0
  
  for (batch_i in 1:n_batches) {
    cat("=== BATCH", batch_i, "of", n_batches, "===\n")
    
    # Get dates for this batch
    start_idx <- (batch_i - 1) * batch_size + 1
    end_idx <- min(batch_i * batch_size, total_dates)
    batch_dates <- dates[start_idx:end_idx]
    
    cat("Processing", length(batch_dates), "dates:", 
        as.character(min(batch_dates)), "to", as.character(max(batch_dates)), "\n")
    
    # Process each date in batch
    batch_data <- list()
    
    for (i in seq_along(batch_dates)) {
      date <- batch_dates[i]
      global_idx <- start_idx + i - 1
      
      cat(sprintf("[%d/%d] ", global_idx, total_dates))
      
      # Check cache first
      cache_file <- NULL
      if (reforecast_config$cache_extractions) {
        cache_file <- paste0(cache_dir, "reforecast_", format(date, "%Y%m%d"), ".rds")
        if (file.exists(cache_file)) {
          cat("Loading from cache:", format(date, "%Y-%m-%d"), "\n")
          batch_data[[i]] <- readRDS(cache_file)
          successful_extractions <- successful_extractions + 1
          next
        }
      }
      
      # Extract data with retries
      extracted_data <- NULL
      for (attempt in 1:reforecast_config$retry_attempts) {
        if (check_reforecast_availability(date)) {
          extracted_data <- extract_reforecast_point_data(date, gefs_centers)
          if (!is.null(extracted_data) && nrow(extracted_data) > 0) {
            break  # Success
          }
        }
        
        if (attempt < reforecast_config$retry_attempts) {
          cat("    Retry", attempt + 1, "for", format(date, "%Y-%m-%d"), "\n")
          Sys.sleep(1)  # Brief pause before retry
        }
      }
      
      if (!is.null(extracted_data) && nrow(extracted_data) > 0) {
        batch_data[[i]] <- extracted_data
        successful_extractions <- successful_extractions + 1
        
        # Cache successful extraction
        if (reforecast_config$cache_extractions && !is.null(cache_file)) {
          saveRDS(extracted_data, cache_file)
        }
      } else {
        cat("Failed to extract data for", format(date, "%Y-%m-%d"), "\n")
      }
      
      # Brief pause between extractions
      Sys.sleep(0.2)
    }
    
    # Combine batch data
    valid_batch_data <- batch_data[lengths(batch_data) > 0]
    if (length(valid_batch_data) > 0) {
      all_data <- c(all_data, valid_batch_data)
    }
    
    cat("Batch", batch_i, "completed:", length(valid_batch_data), "successful extractions\n")
    cat("Cumulative progress:", successful_extractions, "/", global_idx, 
        sprintf("(%.1f%%)\n\n", 100 * successful_extractions / global_idx))
  }
  
  # Combine all successful extractions
  if (length(all_data) > 0) {
    final_dataset <- bind_rows(all_data)
    
    cat("=== EXTRACTION SUMMARY ===\n")
    cat("Total successful extractions:", successful_extractions, "/", total_dates, 
        sprintf("(%.1f%%)\n", 100 * successful_extractions / total_dates))
    cat("Final dataset:\n")
    cat("  Records:", nrow(final_dataset), "\n")
    cat("  Date range:", as.character(range(final_dataset$date)), "\n")
    cat("  GEFS cells:", n_distinct(final_dataset$gefs_cell_id), "\n")
    cat("  Variables:", setdiff(colnames(final_dataset), c("gefs_cell_id", "date", "year", "month", "day_of_year")), "\n")
    
    if (save_output) {
      result <- save_reforecast_training_data(final_dataset)
      return(result)
    } else {
      return(final_dataset)
    }
    
  } else {
    stop("No reforecast data was successfully extracted")
  }
}

cat("GEFSv12 reforecast training data extraction script loaded.\n")
cat("Period available:", reforecast_config$period[1], "to", reforecast_config$period[2], "\n")
cat("Variables:", paste(names(reforecast_config$variables), collapse = ", "), "\n")
cat("Run extract_gefs_reforecast_production(gefs_centers) for full-scale extraction.\n")
cat("Run extract_gefs_reforecast_training(gefs_centers) for smaller-scale testing.\n")