# This script provides several functions for the data preparation
#
# These functions require the package "netcdf.conversions".
# Installation:
# install.packages("devtools")
# devtools::install_github("elisabethvogel/netcdf.conversions", quiet = TRUE)


# load required libraries
require(ncdf4)
require(raster)
require(plyr)
require(dplyr)
require(reshape2)
require(assertthat)
require(lubridate)
require(ncdf.tools)
require(netcdf.conversions)

project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"
data_path = file.path(project_path, "data")

source(sprintf("%s/code/other/file_information.R", project_path), local = TRUE)

# function: calculate_standardised_timeseries ----------------------------------------
#
#   This function is a function to standardise all climate and yield time series
#   Reads in:
#     - crop: "wheat", "maize", "rice" or "soybeans"
#     - crop_calendar: "agmip"
#     - irrigation: "combined", "irrigated", "rainfed"
#     - masked: "masked", "_masked", ""
#     - data_type: "extreme_indicators", "climate_data", "soil_moisture", "drought_indicators"
#     - data_set: depending on data_type, several options
#     - variable: name of the variable to detrend
#     - resolution: 0.5, 2.5, 3.75
#     - skip_existing_files: should files that exist be overwritten or should the calculation be skipped (can speed up process)
#     - verbose
#   Returns:
#     - filename of output netcdf when successful
calculate_standardised_timeseries = function(crop, crop_calendar,
                                             irrigation, masked, 
                                             data_type, data_set, 
                                             variable, 
                                             res, 
                                             skip_existing_files = FALSE,
                                             verbose = FALSE) {
  
  # Legacy variables
  # (specify that the data preparation is done for the whole growing season, without a time lag)
  time_range = "gs"
  lag = 0
  
  if (verbose) print("function: calculate_standardised_timeseries")
  
  resolution = as.character(res)
  
  # some tests
  assert_that(is.character(crop), crop %in% c("wheat", "maize", "rice", "soybeans"))
  assert_that(is.character(crop_calendar), crop_calendar %in% data_sets$crop_calendars)
  assert_that(is.character(irrigation), irrigation %in% c("irrigated", "rainfed", "combined"))
  assert_that(is.character(masked), masked %in% c("_masked", "masked", ""))
  assert_that(is.character(data_type), data_type %in% data_types$growing_season)
  print(data_set)
  print(data_type)
  print(data_sets[[data_type]])
  assert_that(is.character(data_set), data_set %in% data_sets[[data_type]])
  assert_that(is.character(variable), variable %in% variables[[data_type]][[data_set]])
  assert_that(resolution %in% c("0.5", "2.5", "3.75", "1.5"))
  assert_that(is.number(lag))
  
  if (verbose) {
    print(sprintf("Calculation of standardised time series for: %s", crop))
    print(sprintf("Crop calendar: %s, %s, %s", crop_calendar, irrigation, sub("_", "", masked)))
    print(sprintf("Data set: %s, %s", data_set, variable))
    print(sprintf("Resolution: %s", resolution))
  }
  
  # prepare data paths
  if (data_type == "yield") {
    input_dir = sprintf("%s/growing_season_data_detrended/%s_yield_data", data_path, data_set)
    if (data_set == "deepak") {
      input_file = sprintf("%s/%s_%s_yield_1961_2008_%sdeg.nc", input_dir, data_set, crop, resolution)
      sd_file = sprintf("%s/%s_%s_yield_1961_2008_%sdeg_stdev.nc", input_dir, data_set, crop, resolution)
      output_file = sprintf("%s/%s_%s_yield_1961_2008_%sdeg_standardised.nc", input_dir, data_set, crop, resolution)
    } else if (data_set == "iizumi") {
      input_file = sprintf("%s/%s_%s_yield_1982_2006_%sdeg.nc", input_dir, data_set, crop, resolution)
      sd_file = sprintf("%s/%s_%s_yield_1982_2006_%sdeg_stdev.nc", input_dir, data_set, crop, resolution)
      output_file = sprintf("%s/%s_%s_yield_1982_2006_%sdeg_standardised.nc", input_dir, data_set, crop, resolution)
    }
    
  } else {
    if (masked == "masked") masked = "_masked" # replace, if there is a typo
    input_dir = sprintf("%s/growing_season_data_detrended/%s/%s/%s_%s_%s%s/%s",
                        data_path, data_type, data_set, crop_calendar, crop, irrigation, masked,
                        ifelse(lag == 0, "original", sprintf("lag_%s", lag)))
    input_file = sprintf("%s/%s_%s_%s_%s_%s_%sdeg.nc", input_dir, crop, irrigation, data_set, variable, time_range, resolution)
    sd_file = sprintf("%s/%s_%s_%s_%s_%s_%sdeg_stdev.nc", input_dir, crop, irrigation, data_set, variable, time_range, resolution)
    output_file = sprintf("%s/%s_%s_%s_%s_%s_%sdeg_standardised.nc", input_dir, crop, irrigation, data_set, variable, time_range, resolution)
    
  }
  
  # test
  assert_that(is.readable(input_file))
  
  if (skip_existing_files & file.exists(output_file)) {
    print(sprintf("Skipped, because file already exists: %s", output_file))
    return(output_file)
  }
  
  
  ################################################################################
  # Start the calculation
  if (verbose) print("Reading in data...")
  # read in file
  data_frame = netcdf2dataframe(netcdf_file = input_file, variables = "all", time_format = "%Y")
  
  # remove columns with "pct" because relative detrended values are already in some way standardised
  # also the percentage values sometimes make problems, when they are infinitive etc.
  if (any(sapply(names(data_frame), function(x) grepl(x = x, pattern = "pct", fixed = TRUE))))
    data_frame = select(data_frame, -ends_with("pct"))
  
  if (verbose) print("Calculating standard deviation and standardised values...")
  # calculate stdev and standardised values
  # (it seems faster to do it this way - to duplicate standard deviation for each time step than
  # to calculate stdev once for each grid cell and to then merge it with the time series)
  data_frame_long = data_frame %>%
    melt(id.vars = c("lon", "lat", "time")) %>%
    group_by(lon, lat, variable) %>%
    mutate(
      value = ifelse(is.infinite(value), NA, value), # remove infinite values (the first increment detrended values are wrongly infinitive)
      n = sum(!is.na(value)),
      stdev = sd(value, na.rm = TRUE),
      standardised = value/stdev) %>%
    ungroup()
  
  data_frame_long$standardised[data_frame_long$stdev == 0] = 0 # if the standard deviation is 0 (i.e. same value every where, set all standardised values to 0, instead of having NaNs)
  
  # remove all standard deviations / standardisations where less than 5 values where available
  data_frame_long$stdev[data_frame_long$n < 5] = NA
  data_frame_long$standardised[data_frame_long$n < 5] = NA
  
  # retrieve unique values for stdev (it's duplicated for all time steps)
  stdev = dcast(data = data_frame_long, lon + lat ~ variable, value.var = "stdev", fun.aggregate = first)
  standardised = dcast(data = data_frame_long, lon + lat + time ~ variable, value.var = "standardised")
  
  # save netcdf files
  if (verbose) print("Writing files...")
  dataframe2netcdf(stdev, netcdf_file = sd_file, overwrite_existing = TRUE)
  dataframe2netcdf(standardised, netcdf_file = output_file, overwrite_existing = TRUE)
  
  if (verbose) print(sprintf("Standardised timeseries saved at: %s", output_file))
  return(output_file)
  
}