# This script provides several functions for the data preparation
#
# The following functions are included:
# detrending_growing_season_climate         - wrapper function for detrending the yield and climate data
# seriesFill                                - filling missing values in a time series, required for SSA detrending
# ssatrend                                  - for SSA detrending
# detrending.ssa, detrending.ssa.pct        - SSA detrending (abs/rel residuals)
# detrending.linear, detrending.linear.pct  - linear detrending (abs/rel residuals)
# detrending.quad, detrending.quad.pct      - quadratic detrending (abs/rel residuals)
# detrending.best.fit, detrending.best.fit.pct - best fit detrending (abs/rel residuals) - for yields
# detrending.incr, detrending.incr.pct      - annual increments (abs/rel residuals)
# detrending                                - wrapper function for all the other detrending types

# These functions require the package "netcdf.conversions".
# Installation:
# install.packages("devtools")
# devtools::install_github("elisabethvogel/netcdf.conversions", quiet = TRUE)


# load required libraries
library(ncdf4)
library(raster)
library(plyr)
library(dplyr)
library(reshape2)
library(assertthat)
library(lubridate)
library(ncdf.tools)
library(netcdf.conversions)


project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"
data_path = file.path(project_path, "data")

source(file.path(project_path, "code/other/file_information.R"), local = TRUE)

# function: detrend_growing_season_data ----------------------------------------
#
#   This function is a wrapper for the detrend_netcdf function to detrend all growing season data
#   Reads in:
#     - crop: "wheat", "maize", "rice" or "soybeans"
#     - crop_calendar: "agmip"
#     - irrigation: "combined", "irrigated", "rainfed"
#     - masked: "masked", "_masked", ""
#     - data_type: "extreme_indicators", "climate_data", "soil_moisture", "drought_indicators"
#     - data_set: depending on data_type, several options
#     - variable: name of the variable to detrend
#     - statistic: for gs climate data: default = mean, min and max 
#     - resolution: 0.5, 2.5, 3.75
#     - skip_existing_files: should files that exist be overwritten or should the calculation be skipped (can speed up process)
#     - detrending_types: which detrending_method to use, if NA then this will be set automatically
#     - skip_existing_files: should files that exist be overwritten or should the calculation be skipped (can speed up process)
#     - verbose
#   Returns:
#     - filename of output netcdf when successful

detrend_growing_season_data = function(crop, crop_calendar,
                                       irrigation, masked, 
                                       data_type, data_set, variable, 
                                       res, 
                                       detrending_types = NA,
                                       skip_existing_files = FALSE,
                                       verbose = FALSE) {
  
  # Legacy variables
  # (specify that the data preparation is done for the whole growing season, without a time lag)
  time_range = "gs"
  lag = 0
  
  if (verbose) print("function: detrend_growing_season_data")
  
  resolution = as.character(res)
  
  # some tests
  assert_that(is.character(crop), crop %in% c("wheat", "maize", "rice", "soybeans"))
  assert_that(is.character(crop_calendar), crop_calendar %in% data_sets$crop_calendars)
  assert_that(is.character(irrigation), irrigation %in% c("irrigated", "rainfed", "combined"))
  assert_that(is.character(masked), masked %in% c("_masked", "masked", ""))
  assert_that(is.character(data_type), data_type %in% data_types$growing_season)
  assert_that(is.character(data_set), data_set %in% data_sets[[data_type]])
  assert_that(is.character(variable), variable %in% variables[[data_type]][[data_set]])
  assert_that(resolution %in% c("0.5", "2.5", "3.75", "1.5"))
  assert_that(is.number(lag))
  
  if (verbose) {
    print(sprintf("Detrending growing season data: %s", crop))
    print(sprintf("Crop calendar: %s, %s, %s", crop_calendar, irrigation, sub("_", "", masked)))
    print(sprintf("Data set: %s, %s", data_set, variable))
    print(sprintf("Resolution: %s", resolution))
  }
  
  # prepare data paths
  if (data_type == "yield") {
    input_dir = sprintf("%s/crop_yield_data/%s_crop_data_regridded",
                        data_path, data_set)
    output_dir = sprintf("%s/growing_season_data_detrended/%s_yield_data",
                         data_path, data_set)
    if (data_set == "deepak") {
      netcdf_file = sprintf("%s/%s_%s_yield_1961_2008_%sdeg.nc", 
                            input_dir, data_set, crop, resolution)
      output_file = sprintf("%s/%s_%s_yield_1961_2008_%sdeg.nc", 
                            output_dir, data_set, crop, resolution)
    } else if (data_set == "iizumi") {
      netcdf_file = sprintf("%s/%s_yield_estimate_50.0_%s_1982_2006_%sdeg.nc", 
                            input_dir, data_set, crop, resolution)
      output_file = sprintf("%s/%s_%s_yield_1982_2006_%sdeg.nc", 
                            output_dir, data_set, crop, resolution)
    }
    
  } else {
    if (masked == "masked") masked = "_masked" # replace, if there is a typo
    input_dir = sprintf("%s/growing_season_data_regridded/%s/%s/%s_%s_%s%s/%s",
                        data_path, data_type, data_set, crop_calendar, crop, irrigation, masked,
                        ifelse(lag == 0, "original", sprintf("lag_%s", lag)))
    
    netcdf_file = sprintf("%s/%s_%s_%s_%s_%s_%sdeg.nc", 
                          input_dir, crop, irrigation, data_set, variable, time_range, resolution)
    output_dir = sprintf("%s/growing_season_data_detrended/%s/%s/%s_%s_%s%s/%s",
                         data_path, data_type, data_set, crop_calendar, crop, irrigation, masked,
                         ifelse(lag == 0, "original", sprintf("lag_%s", lag)))
    
    output_file = sprintf("%s/%s_%s_%s_%s_%s_%sdeg.nc", 
                          output_dir, crop, irrigation, data_set, variable, time_range, resolution)
  }
  
  # some tests
  print(netcdf_file)
  assert_that(is.readable(netcdf_file))
  assert_that(length(output_file) == 1)
  
  if (skip_existing_files & file.exists(output_file)) {
    print(sprintf("Skipped, because file already exists: %s", output_file))
    return(output_file)
  }
  
  # prepare variable names
  if (variable == "yield") {
    var_names = "yield"
  } else {
    if (variable %in% c("tmp", "pre", "frs", "dtr", "DTR", "SMm", "TX90p", "TN10p")) {
      statistics = "mean"
    } else if (variable %in% c("TXx","Rx5day", "Rx1day", "tmx", "SMx")) {
      statistics = "max"
    } else if (variable %in% c("TNn", "tmn", "SMn")) {
      statistics = "min"
    } else if (variable %in% c("sc_PDSI_pm")) {
      statistics = c("mean", "max")
    } else {
      statistics = c("mean", "min", "max")
    }
    
    var_names = sprintf("%s_%s", variable, time_range)
    var_names = sprintf("%s_%s", var_names, statistics)
    var_names = sprintf("%s%s", var_names, ifelse(lag == 0, "", sprintf("_lag_%s", lag)))
  }
  
  if (data_type == "yield") {
    if (is.na(detrending_types))
      detrending_types = c("original", "lin", "quad", "ssa", "incr", "best_fit")
  } else {
    if (is.na(detrending_types))
      detrending_types = c("original", "lin", "quad", "ssa", "incr")
  }
  
  # do the calculations
  detrend_netcdf(netcdf_file = netcdf_file,
                 var_names = var_names,
                 detrending_types = detrending_types,
                 output_file = output_file,
                 verbose = verbose)
  
  if (verbose) print(sprintf("Detrended netcdf saved at: %s", output_file))
  return(output_file)
}





################################################################################
# Helper function: Fill in the missing values of the time series with a linear interpolation
seriesFill = function(y) {
  # Generate a generic x vector and a new y vector
  x = c(1:length(y))
  ynew = y
  # Find indices of missing values
  ii = which(is.na(y))
  # Approximate the missing values with a linear fit
  if (length(ii)>0) ynew[ii] = approx(x,y,xout=x[ii],rule=2)$y
  return(ynew)
} 

################################################################################
## Helper function: Compute a non-linear trend from a given time series (input) using SSA
## 	input:  1-dimensional vector
##	output: ssa-trend as a 1-dimensional vector

ssatrend = function(y, L = (length(y) - 1) %/% 2, fill = TRUE) {
  library(Rssa)
  # Compute the trend from SSA decomposition
  
  # If all (or only one) values in input are NA, then return NA
  if( sum(!is.na(y)) <= 1 ){
    trend = rep(NA, length(y))
    return(trend)
  }
  
  # If fill requested, make sure that NAs are filled in
  if (fill) y = seriesFill(y)
  
  # Get a new ssa object and
  # Reconstruct the first element
  ssa   = ssa(y, L = L)
  trend = reconstruct(ssa, groups = list(1))$F1
  
  # That's it! Return the trend...
  return(trend)
}


################################################################################
# SSA detrending - returning absolute residuals
detrending.ssa = function(x){
  x = as.ts(x)
  x.new = x - ssatrend(x,L=5)   
  return(as.numeric(x.new))
}

################################################################################
# SSA detrending - returning relative residuals
detrending.ssa.pct = function(x){    
  x = as.ts(x)
  model = ssatrend(x, L=5)
  x.new = (x - model) / model
  return(as.numeric(x.new))
}

################################################################################
# linear detrending - returning absolute residuals
detrending.linear = function(x){
  t = 1:length(x)
  lm = lm(x ~ t, na.action = na.exclude) #fit a linear model to the time series
  return(resid(lm)) #return a time series of residuals
}

################################################################################
# linear detrending - returning relative residuals
detrending.linear.pct = function(x){
  t = 1:length(x)
  lm = lm(x ~ t, na.action = na.exclude) #fit a linear model to the time series
  return(resid(lm)/predict(lm)*100) #return a time series of residuals relative to predicted value (in percent)
}

################################################################################
# quadratic detrending - returning absolute residuals
detrending.quadratic = function(x){
  t = 1:length(x)
  lm = lm(x ~ poly(t, 2), na.action = na.exclude) #fit a linear model to the time series
  return(resid(lm)) #return a time series of residuals
}

################################################################################
# quadratic detrending - returning relative residuals
detrending.quadratic.pct = function(x){
  t = 1:length(x)
  lm = lm(x ~ poly(t, 2), na.action = na.exclude) #fit a linear model to the time series
  return(resid(lm)/predict(lm)*100) #return a time series of residuals relative to predicted value (in percent)
}

################################################################################
# best fit detrending - returning absolute residuals
# best fit from 4 models: 1) de-meaning, 2) linear, 3) quadratic, 4) cubic
detrending.best.fit = function(x){
  t = 1:length(x)
  
  # fit four different models
  lm_fit = list()
  lm_fit[[1]] = lm(x ~ 1, na.action = na.exclude)
  lm_fit[[2]] = lm(x ~ t, na.action = na.exclude)
  lm_fit[[3]] = lm(x ~ poly(t, 2), na.action = na.exclude)
  lm_fit[[4]] = lm(x ~ poly(t, 3), na.action = na.exclude)
  
  # calculate AIC for each model
  aic = sapply(X = lm_fit, FUN = AIC)
  
  # select best model
  pos = which(aic == min(aic))
  if (length(pos) > 1) { # no model fits best
    return(rep(NA, length(x)))
  } else {
    lm_fit = lm_fit[[pos]]
    return(resid(lm_fit)) #return a time series of residuals
  }
}

################################################################################
# best fit detrending - returning relative residuals
# best fit from 4 models: 1) de-meaning, 2) linear, 3) quadratic, 4) cubic
detrending.best.fit.pct = function(x){
  t = 1:length(x)
  
  # fit four different models
  lm_fit = list()
  lm_fit[[1]] = lm(x ~ 1, na.action = na.exclude)
  lm_fit[[2]] = lm(x ~ t, na.action = na.exclude)
  lm_fit[[3]] = lm(x ~ poly(t, 2), na.action = na.exclude)
  lm_fit[[4]] = lm(x ~ poly(t, 3), na.action = na.exclude)
  
  # calculate AIC for each model
  aic = sapply(X = lm_fit, FUN = AIC)
  
  # select best model
  pos = which(aic == min(aic))
  if (length(pos) > 1) { # no model fits best
    return(rep(NA, length(x)))
  } else {
    lm_fit = lm_fit[[pos]]
    return(resid(lm_fit)/predict(lm_fit)*100) #return a time series of residuals relative to predicted value (in percent)
  }
}

################################################################################
# annual increments - returning absolute residuals
detrending.incr = function(x){
  y = NULL
  y[1] = NA
  y[2:length(x)] = x[2:length(x)] - x[1:(length(x)-1)]
  return(y)
}

################################################################################
# annual increments - returning relative residuals
detrending.incr.pct = function(x){
  y = NULL
  y[1] = NA
  y[2:length(x)] = (x[2:length(x)] - x[1:(length(x)-1)]) / x[1:(length(x)-1)] * 100
  return(y)
}

################################################################################
# wrapper function for all the different detrending types mentioned above
# options for detrending_type: lin, lin_pct, quad, quad_pct, best_fit, best_fit_pct, ssa, ssa_pct, incr, incr_pct

detrending = function(x, detrending_type, required_timesteps = 10, verbose = FALSE){
  
  if (verbose) print(sprintf("function: detrend"))
  
  # if time series has less than 'required_timesteps' values, return a vector of NAs
  if (sum(!is.na(x)) < required_timesteps) return(rep(NA, length(x)))
  
  if (detrending_type == "lin") {
    if (verbose) print(sprintf("linear detrending"))
    return(detrending.linear(x))
  } else if (detrending_type == "lin_pct") {
    if (verbose) print(sprintf("linear detrending, percentages"))
    return(detrending.linear.pct(x))
  } else if (detrending_type == "quad") {
    if (verbose) print(sprintf("quadratic detrending"))
    return(detrending.quadratic(x))
  } else if (detrending_type == "quad_pct") {
    if (verbose) print(sprintf("quadratic detrending, percentages"))
    return(detrending.quadratic.pct(x))
  } else if (detrending_type == "best_fit") {
    if (verbose) print(sprintf("best fit detrending"))
    return(detrending.best.fit(x))
  } else if (detrending_type == "best_fit_pct") {
    if (verbose) print(sprintf("best fit detrending, percentages"))
    return(detrending.best.fit.pct(x))
  } else if (detrending_type == "ssa") {
    if (verbose) print(sprintf("ssa detrending"))
    return(detrending.ssa(x))
  } else if (detrending_type == "ssa_pct") {
    if (verbose) print(sprintf("ssa detrending, percentages"))
    return(detrending.ssa.pct(x))
  } else if (detrending_type == "incr") {
    if (verbose) print(sprintf("annual increments"))
    return(detrending.incr(x))
  } else if (detrending_type == "incr_pct") {
    if (verbose) print(sprintf("annual increments, percentages"))
    return(detrending.incr.pct(x))
  }
}


# function: detrend_netcdf ----------------------------------------
# 
# This function detrends the data in a netcdf file
#
#   Reads in:
#     - netcdf_file: file name of the netcdf file
#     - var_names: variable name or names of the variable(s) to be detrended
#     - detrending_types: one or more detrending methods to be used for the detrending
#     - verbose: 
#   Returns:
#     - data_frame of new coordinates with names "lon" and "lat"

detrend_netcdf = function(netcdf_file,
                          var_names,
                          detrending_types = c("original", "lin", "lin_pct", "quad", "quad_pct",
                                               "ssa", "ssa_pct", "best_fit", "best_fit_pct",
                                               "incr", "incr_pct"),
                          output_file,
                          verbose = FALSE) {
  
  if (verbose) print("function: detrend_netcdf")
  
  # some tests
  assert_that(is.character(netcdf_file), is.readable(netcdf_file))
  assert_that(is.character(detrending_types))
  assert_that(all(detrending_types %in% 
                    c("original", "lin", "lin_pct", "quad", "quad_pct",
                      "ssa", "ssa_pct", "best_fit", "best_fit_pct",
                      "incr", "incr_pct")))
  assert_that(is.character(output_file))
  
  # open netcdf file
  data_nc = nc_open(netcdf_file)
  
  # retrieve meta data about file
  dims = infoNcdfDims(netcdf_file)
  # remove all dimensions with size = 1
  dims = dims %>% filter(length > 1)
  assert_that("time" %in% dims$name, "lat" %in% dims$name | "latitude" %in% dims$name,
              "lon" %in% dims$name | "longitude" %in% dims$name)
  lat_pos = which(dims$name == "lat" | dims$name == "latitude") 
  lon_pos = which(dims$name == "lon" | dims$name == "longitude") 
  time_pos = which(dims$name == "time") 
  time_size = dims[time_pos, 3]
  lat_size = dims[lat_pos, 3]
  lon_size = dims[lon_pos, 3]
  
  # coordinates = readNcdfCoordinates(netcdf_file)
  coordinates = list()
  dim_names = dims$name
  for (dim_name in dim_names) {
    coordinates[[dim_name]] = ncvar_get(data_nc, varid = dim_name)
  }
  lat_vals = coordinates[lat_pos][[1]]
  lon_vals = coordinates[lon_pos][[1]]
  time_vals = coordinates[time_pos][[1]]
  
  # loop through variables
  for (var_number in 1:length(var_names)) {
    var_name = var_names[var_number]
    if (verbose) print(sprintf("Variable: %s", var_name))
    var = ncvar_get(data_nc, var_name)
    # reorder data dimensions
    # 1: lon, 2: lat, 3: time
    var = aperm(var, c(lon_pos, lat_pos, time_pos))
    
    # create empty result array
    detrended_data = list()
    for (detrending_type in detrending_types) {
      detrended_data[[detrending_type]] = array(data = NA, dim = c(lon_size, lat_size, time_size))
    }
    
    # loop through coordinates and do the detrending for each grid cell
    for (lat_idx in 1:lat_size) {
      for (lon_idx in 1:lon_size) {
        
        # if (verbose) print(sprintf("lat: %.2f, lon: %.2f", lat_vals[lat_idx], lon_vals[lon_idx]))
        
        data_gridcell = var[lon_idx, lat_idx, ]
        
        # if all climate data are empty, go to next grid cell
        if (all(is.na(data_gridcell)))
          next
        
        for (detrending_type in detrending_types) {
          if (detrending_type == "original") {
            detrended_data[[detrending_type]][lon_idx, lat_idx, ] = data_gridcell
          } else {
            detrended_data[[detrending_type]][lon_idx, lat_idx, ] = 
              detrending(x = data_gridcell, detrending_type = detrending_type)
          }
        }
        
      } # lat_idx
    } # lon_idx
    
    # create netcdf for results
    if (var_number == 1) {
      output_lat = ncdim_def("lat", "latitude", lat_vals)
      output_lon = ncdim_def("lon", "longitude", lon_vals)
      output_time = data_nc$dim$time
      output_units = ncatt_get(data_nc, var_name, "units")$value
      if (output_units == 0) output_units = ""
    }
    
    detrended_var = list()
    for (detrending_type in detrending_types) {
      if (detrending_type == "original") {
        detrended_var[[detrending_type]] = ncvar_def(name = var_name,
                                                     units = output_units,
                                                     dim = list(output_lon, output_lat, output_time),
                                                     prec = "double")
      } else {
        detrended_var[[detrending_type]] = ncvar_def(name = sprintf("%s_detrended_%s", var_name, detrending_type),
                                                     units = output_units,
                                                     dim = list(output_lon, output_lat, output_time),
                                                     prec = "double")
      }
    }
    
    # create netcdf or add variables to existing netcdf
    if (var_number == 1) {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
      detrended_data_nc = nc_create(filename = output_file, vars = detrended_var)
    } else {
      for (detrending_type in detrending_types) {
        detrended_data_nc = ncvar_add(nc = detrended_data_nc, v = detrended_var[[detrending_type]])
      }
    }
    
    # put values into file
    for (detrending_type in detrending_types) {
      if (detrending_type == "original") {
        ncvar_put(detrended_data_nc, varid = var_name,
                  vals = detrended_data[[detrending_type]])
      } else {
        ncvar_put(detrended_data_nc, varid = sprintf("%s_detrended_%s", var_name, detrending_type),
                  vals = detrended_data[[detrending_type]])
      }
    }
    
  } # var_number
  nc_close(detrended_data_nc)
  nc_close(data_nc)
  
  return(output_file)
}


