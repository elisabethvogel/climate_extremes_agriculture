
# This script provides several functions for the statistical data analysis

# These functions require the package "netcdf.conversions".
# Installation:
# install.packages("devtools")
# devtools::install_github("elisabethvogel/netcdf.conversions", quiet = TRUE)

################################################################################

options(stringsasfactors = FALSE)

# load required libraries
# require(relaimpo)
require(raster) # load before plyr and dplyr
require(plyr)
require(dplyr)
require(reshape2)
require(assertthat)
require(lubridate)
# require(stringr)
require(netcdf.conversions)

# set project path
project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"
# read in required functions
source(sprintf("%s/code/other/file_information.R", project_path))

###############################################################################
# function: calculate_results_statistical_analysis
# 
#   This function is a wrapper function for the statistical analysis function: random_forest
#   
#   Reads in:
#   parameter_list (list), which has to include the following parameters:
#
#     - crop: "wheat", "maize", "rice" or "soybeans"
#     - yield_dataset: "deepak" or "iizumi"
#     - extremes_dataset: "hadex2" or "ghcndex"
#     - climate_dataset: "cru_ts_323"
#     - soil_moisture_dataset: "esa_soil_moisture"
#     - drought_indicator_dataset: "ncar_pdsi"
#     - crop_calendar: "agmip"
#     - irrigation: "combined", "irrigated", "rainfed"
#     - masked: "masked", ""
#     - detrending_yield: "lin", "quad", best_fit", "incr", "ssa"
#     - detrending_climate: "lin", "quad", "incr", "ssa"

#     - extreme_indicator_vars, climate_vars, soil_moisture_vars, drought_vars:
#       variable names to be used as predictors for the statistical analysis,
#       full name with growing season statistic, e.g. tmp_gs_mean, TXx_gs_max
#     - grouping_by: should the statistical analysis be carried out for each gridcell or for whole regions / globally
#       options are: coordinates (lat, lon), continent, un_region, country, global

#       or for each year (grouping_vars = "time") or for all years and grid cells together (grouping_vars = NA)
#     - include_time_as_predictor, include_location_as_predictor: 
#       logical, Should time or lat/lon be used as predictor variables? default = FALSE
#       Particularly needed for when the analysis is done with whole data set, not
#       for each grid cell individually.
#     - location_predictor: if location is included as predictor variable, which level of location
#       options are: coordinates (lat, lon), continent, un_region, country

#     - resolution: 0.5, 2.5, 3.75
#     - year_min, year_max: time range for the analysis, automatically it will be reduced to the years where all datasets have data
#     - stepwise_direction: "forward" or "backward", in case of stepwise regression
#     - need_all_vars: TRUE, FALSE: if TRUE, only gridcells (or whichever grouping variables are used), where all 
#       predictors have at least one datapoint are used
#       if FALSE, predictors with no data will be removed and the analysis is done with the rest of the predictors
#     - quadratic_terms: for which variables should quadratic terms be included in statistical analysis?
#     - comment: optional, will be used to create a subdirectory to distuinguish different results, default is the current date and time
#     - statistical_method: "random_forest"
#     - stepwise_direction: "forward" or "backward"
#     - save_files: should the output files (csv and netcdf) be saved?
#     - verbose
#   Returns:
#     - lots of stuff (--> in the [...]/results folder) plus a list with all filenames

calculate_results_statistical_analysis = function(parameter_list, data_all,
                                                  overwrite_data = FALSE) {
  
  # read in parameters ========================================
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  if (verbose) print("function: calculate_results_statistical_analysis")
  if (verbose) print(sprintf("Creating output directory: %s", output_dir))
  
  # create directories (if not done yet)
  sapply(c("data_used", "final", "csv", "netcdf", "model"), function(x) 
    dir.create(sprintf("%s/%s", output_dir, x), showWarnings = FALSE, recursive = TRUE))
  
  # Check if output file exists already (which means the whole analysis went through already)
  result_output_file = sprintf("%s/final/results_summary_dataframe.csv", output_dir)
  if (file.exists(result_output_file) && !overwrite_data) {
    if (verbose) print(sprintf("Results for this analysis have previously been computed. Return results from: %s", result_output_file))
    results = transpose_back_df(read.csv(result_output_file, row.names = NULL))
    return(results)
  }
  
  # Check if empty result files exists (meaning that NULL was returned previously, for example, because there is not enough data for this crop / region / irrigation combination)
  empty_results_file = "NULL_results.txt"
  if (file.exists(empty_results_file) && !overwrite_data) {
    if (verbose) print("Empty results previously computed. Returning NULL.")
    return(NULL)
  }
  
  # Start ========================================
  
  # No data available for this crop / region / irrigation combination
  if (is.null(data_all)) {
    if (save_files) {
      # create a text file that says this job is completed
      fn = "NULL_results.txt"
      write.csv(data.frame(empty = TRUE), fn, row.names = FALSE)
    }
    return(NULL)
  }
  
  # create 5 groups (each year is in one of those 5 groups) --> every cross_val_iteration in cross validation uses 80% of data
  if (!is.na(seed)) {
    set.seed(seed)
  }
  years_temp = unique(data_all$time)
  n_timesteps = length(years_temp)
  
  if (oos_calculation == TRUE) {
    years_groupings = data.frame(
      time = sample(x = years_temp, size = n_timesteps, replace = FALSE),
      # out-of-sample group:
      oos_group = ceiling(1:length(years_temp) / (n_timesteps / n_groups_out_of_bag_sampling))) 
    
    # test that there are n groups
    assert_that(length(unique(years_groupings$oos_group)) == n_groups_out_of_bag_sampling) 
    data_all = merge(data_all, years_groupings, by = "time")
    
  } else {
    data_all$oos_group = 0
  }
  
  
  # calculate model parameters and fitted values/residuals ========================================
  temp = calculate_model_parameters_fitted_values_residuals(data_all = data_all, parameter_list = parameter_list)
  model_parameters = temp$model_parameters
  fitted_values_residuals = temp$fitted_values_residuals
  
  # rename fitted_values_residuals
  fitted_values_residuals$train = fitted_values_residuals$train_using_all_data
  fitted_values_residuals$test = fitted_values_residuals$test_combined
  fitted_values_residuals$train_using_all_data = NULL
  fitted_values_residuals$test_combined = NULL
  
  # combine both datasets
  fitted_values_residuals = do.call(rbind, fitted_values_residuals)
  fitted_values_residuals$train = grepl(pattern = "train", x = row.names(fitted_values_residuals), fixed = TRUE)
  yield_data = data_all %>% select(lon, lat, time, starts_with("yield")) %>% unique()
  
  n_before = nrow(fitted_values_residuals)
  fitted_values_residuals = merge(fitted_values_residuals, yield_data, by = c("lon", "lat", "time"), all = FALSE)
  assert_that(n_before == nrow(fitted_values_residuals))
  
  # Calculate standardised statistical parameters from fitted values ========================================
  
  # Calculate standardised Rsq and RMSE (global, continent, un_region, country and grid cell)
  if (verbose) print("Calculate standardised R-squared and RMSE (for: global, continent, un_region, country, grid cell).")
  
  standardised_statistical_parameters = calculate_statistical_parameters_from_fitted_values(
    data_frame = fitted_values_residuals, 
    fitted_values_var = "fitted_values", 
    observed_values_var = "original",
    parameter_list = parameter_list,
    significant_only = FALSE)
  
  statistical_parameters_per_region = standardised_statistical_parameters$statistical_parameters_per_region
  output_global =  statistical_parameters_per_region[["global"]]
  output_continent = statistical_parameters_per_region[["continent"]] 
  
  output_continent = output_continent %>%
    select(train, continent, mean_r_squared, mean_rmse_rel_sd, mean_rmse_rel_mean, r_squared_all,
           rmse_rel_sd_all, rmse_rel_mean_all,
           r_squared_yield, r_squared_production) %>%
    recast(id.var = c("train", "continent"), measure.var = 
             c("mean_r_squared", "mean_rmse_rel_sd",
               "mean_rmse_rel_mean", "r_squared_all", 
               "rmse_rel_sd_all", "rmse_rel_mean_all",
               "r_squared_yield", "r_squared_production"),
           train ~ continent + variable)
  
  output_parameters_train = data.frame(output_global, output_continent) %>%
    filter(train == TRUE) %>% select(-train, -train.1)
  output_parameters_test = data.frame(output_global, output_continent) %>%
    filter(train == FALSE) %>% select(-train, -train.1)
  names(output_parameters_train) = sprintf("%s_train", names(output_parameters_train))
  names(output_parameters_test) = sprintf("%s_test", names(output_parameters_test))
  
  
  # prepare the results to return ========================================
  
  if (verbose) print("Prepare results to return")
  
  input_parameters = data.frame(statistical_method = statistical_method,
                                stepwise_direction = stepwise_direction,
                                crop = crop,
                                yield_dataset = yield_dataset,
                                extremes_dataset = extremes_dataset,
                                climate_dataset = climate_dataset,
                                soil_moisture_dataset = soil_moisture_dataset,
                                drought_indicator_dataset = drought_indicator_dataset,
                                crop_calendar = crop_calendar, 
                                irrigation = irrigation, 
                                masked = masked,
                                detrending_climate = detrending_climate, 
                                detrending_yield = detrending_yield,
                                extreme_indicator_vars = paste(extreme_indicator_vars, collapse = ","),
                                climate_vars = paste(climate_vars, collapse = ","),
                                soil_moisture_vars = paste(soil_moisture_vars, collapse = ", "),
                                grouping_by = grouping_by,
                                grouping_vars = paste(grouping_vars, collapse = ","),
                                use_time_as_predictor = as.character(use_time_as_predictor),
                                use_location_as_predictor = as.character(use_location_as_predictor),
                                location_predictor = location_predictor,
                                use_only_regions_with_frost_days = use_only_regions_with_frost_days,
                                standardise_variables = as.character(standardise_variables),
                                res = res,
                                year_min = min(data_all$time, na.rm = TRUE),
                                year_max = max(data_all$time, na.rm = TRUE),
                                time_range = time_range,
                                lag = lag,
                                shift_yield_time_by_x_years = shift_yield_time_by_x_years,
                                need_all_vars = as.character(need_all_vars),
                                quadratic_terms = paste(quadratic_terms, collapse = ","),
                                seed = seed,
                                oos_calculation = oos_calculation,
                                n_groups_out_of_bag_sampling = n_groups_out_of_bag_sampling,
                                crop_to_region = crop_to_region,
                                
                                source_area_harvested = source_area_harvested,
                                comment = comment,
                                number_of_groups = ifelse(exists("number_of_groups_tested"),
                                                          number_of_groups_tested, NA)
  )
  
  if (oos_calculation == TRUE) {
    results = data.frame(input_parameters, output_parameters_train, output_parameters_test)
  } else {
    results = data.frame(input_parameters, output_parameters_train)
  }
  
  
  if (save_files) {
    # write statistical parameters ========================================
    if (verbose) print("Save results.")
    write.csv(transpose_df(results), result_output_file, row.names = FALSE)
    
    # create a text file that says this job is completed
    fn = "job_completed.txt"
    write.csv(data.frame(completed = TRUE), fn, row.names = FALSE)
  }
  
  return(results)
}


###############################################################################
# helper function: checks for a grid cell / region, whether there are enough data points / observations
check_if_enough_variables_available = function(data, predictand, 
                                               predictors = 
                                                 names(data)[!names(data) %in% 
                                                               c("lon", "lat", "time", "crop",
                                                                 "crop_calendar", "irrigation",
                                                                 "masked", "resolution",
                                                                 predictand)],
                                               need_all_vars, 
                                               n_threshold = 20,
                                               verbose = FALSE) {
  
  
  # check that all variables are present in this time series ========================================
  if (all(is.na(data[, predictand]))) {
    return(-1000) # return -1000 when predictand is not available
  }
  
  # otherwise, check if any variable is missing
  col_is_NA = apply(data, 2, function(x) all(is.na(x)))
  
  if (any(col_is_NA) & need_all_vars == TRUE) {
    return(-999) # return -999 when some of the variables are not available but all variables are needed
  }
  
  # either need_all_vars = FALSE or all variables available
  
  # remove variables with no data
  data = data[, !col_is_NA, drop = FALSE]
  
  # look at remaining predictors
  predictors_new = intersect(predictors, names(data))
  if (length(predictors_new) == 0) {
    return(-998) # return -998 when only predictand, but no predictor variable is available
  }
  
  # check for variables that have too few data points
  col_is_NA = apply(data, 2, function(x) sum(!is.na(x)) < n_threshold)
  
  # remove variables with not enough data
  data = data[, !col_is_NA, drop = FALSE]
  
  # check if predictand still in data set
  if (!predictand %in% names(data)) {
    return(-1000) # return -1000, when not enough predictand data points
  }
  
  # look at remaining predictors
  predictors_new = intersect(predictors, names(data))
  if (length(predictors_new) == 0) {
    return(-998) # return -998 when no predictor variable is available
  }
  
  # check, if there are at least n_threshold years with overlapping observations for all variables
  if (nrow(na.exclude(data)) < n_threshold) {
    return(-997) # return -997, when overlapping time series are not long enough
  }
  
  # else, if everything went fine
  return(1)
  
}



###############################################################################
# function: random_forest ----------------------------------------
# 
#   This function returns the result of a stepwise forward regression for 
#   a given set of predictor variables and one predictor variable
#
#   Reads in:
#     - data: (data.frame) must contain the predictand and predictor variables
#     - predictand: (character) name of the predicted variable
#     - predictors: (character) name(s) of predictor variable(s)
#     - direction: (character) "forward" or "backward"
#     - need_all_vars: (logical) TRUE = if one variable is completely empty, -999 is returned
#       FALSE: if one variable is completely empty, it will be removed from the analysis
#     - return_data: "model_parameters" or "time_series"
#       model_parameters: returns a row-vector of parameters of the statistical model
#       (coefficient estimates, p-val, in_model, R-squared, adj-R-squared, RMSE)
#       time_series: returns fitted values and residuals
#       (not necessarily a time series, but usually it is dependent on time)
#     - time_var: (character) name of the time dimension, needed for return_data = "time_series"
#     - verbose
#   Returns:
#   output (list), including:
#       - model_parameters (list)
#         - importance_IncMSE (dataframe)
#         - importance_IncNodePurity (dataframe)
#         - r_squared_from_RF (numeric)
#         - oob_error_from_RF (numeric)
#         - n_total (numeric)
#         - n_train (numeric)
#         - n_test (numeric)
#       - time_series (list)
#         - train (dataframe)
#           - lon, lat, time, original, fitted_values, residuals
#         - test (dataframe)
#           - lon, lat, time, original, fitted_values, residuals

random_forest = function(data,
                         predictand,
                         predictors,
                         need_all_vars = TRUE,
                         time_var = "time", # when returning time series
                         train = NULL,
                         n_threshold = 20,
                         seed_random_forest = 123456,
                         verbose = FALSE,
                         manual_subsetting = FALSE,
                         sampling_fraction = NULL,
                         n_repeat = ifelse(manual_subsetting == TRUE, 10, 1),
                         extremes_dataset = NULL,
                         ntree = 500 / n_repeat,
                         res = NULL) {
  
  if (verbose) print("function: random_forest")
  
  library(randomForest)
  
  # some tests
  assert_that(is.data.frame(data))
  assert_that(predictand %in% names(data), all(predictors %in% names(data)))
  assert_that(is.logical(need_all_vars))
  assert_that(is.null(train) || (is.logical(train) && length(train) == nrow(data)))
  if (manual_subsetting == TRUE) {
    assert_that(!is.null(sampling_fraction))
    assert_that(!is.null(extremes_dataset), !is.null(res))
  }
  
  time = data[, time_var]
  data_original = data
  
  # create empty output dataframe for statistical parameters
  output = list()
  output$model_parameters = list()
  
  for (parameter in c("importance_IncMSE", "importance_IncNodePurity")) {
    output$model_parameters[[parameter]] = as.data.frame(t(rep(NA, length(predictors))))
    names(output$model_parameters[[parameter]]) = sprintf("%s_%s", predictors, parameter)
  }
  
  output$model_parameters$r_squared_from_RF = NA
  output$model_parameters$oob_error_from_RF = NA
  output$model_parameters$n = data.frame(n_total = NA, n_train = NA, n_test = NA)
  
  
  # create temporary dataframe for time series (will later be divided into training and test time series)
  time_series = data.frame(data_original[, intersect(names(data_original), c("time", "lon", "lat")), drop = FALSE], # include all dimensions that are available
                           original = data[, predictand])
  time_series$fitted_values = rep(NA, nrow(data))
  time_series$residuals = rep(NA, nrow(data))
  
  # check that all variables are present in this time series ========================================
  
  test = check_if_enough_variables_available(data = data, predictand = predictand,
                                             predictors = predictors,
                                             need_all_vars = need_all_vars,
                                             n_threshold = n_threshold,
                                             verbose = verbose)
  
  if (test < 0) { # just return empty output with NAs, not -999 etc anymore, because this can be found in "coverage" data
    
    # put time series into output
    output$time_series = list()
    if (is.null(train)) {
      output$time_series$train = time_series
      output$time_series$test = NULL
    } else {
      output$time_series$train = time_series[train == 1, ]
      output$time_series$test = time_series[train == 0, ]
    }
    
    output$model = NULL
    
    return(output)
  }
  
  # else, if test > 0, go ahead with the analysis ========================================
  
  # remove all variables that do not have enough values (check_if_enough_variables_available it was already tested
  # if this is consistent with the need_all_vars setting)
  remove_columns = apply(data, 2, function(x) {sum(!is.na(x)) < n_threshold})
  data = data[, !remove_columns, drop = FALSE]
  predictors_new = intersect(predictors, names(data))
  
  # subset the data into training and test data set 
  if (!is.null(train)) { #
    data_train = data[train, ]
    data_test = data[!train, ]
    if (nrow(data_train) == 0) data_train = NULL
    if (nrow(data_test) == 0) data_test = NULL
  } else { #  if all years are to be used 
    data_train = data
    data_test = NULL
  }
  
  data = na.exclude(data)
  data_train = na.exclude(data_train)
  data_test = na.exclude(data_test)
  n_total = nrow(data)
  n_train = nrow(data_train)
  n_test = nrow(data_test)
  assert_that(n_train >= n_train)
  
  if (is.null(n_train)) n_train = NA
  if (is.null(n_test)) n_test = NA
  
  # subset data if applicable (for random forest variable importance plots)
  if (manual_subsetting == TRUE && verbose)
    if (verbose) print("Repeating random forest training with subsamples of the training data that exclude grid cells that are too close together (to reduce spatial correlations due to coarser HadEX2 grid).")
  
  forest_list = list()
  
  for (repeat_idx in 1:n_repeat) {
    
    if (manual_subsetting == TRUE && verbose)
      if (verbose) print(sprintf("Training random forest: %s", repeat_idx))
    # subset 
    if (manual_subsetting == TRUE) {
      data_train_temp = create_manual_subset(data_frame = data_train, sampling_fraction = sampling_fraction, verbose = verbose, extremes_dataset = extremes_dataset, res = res)
    } else {
      data_train_temp = data_train
    }
    
    # do regression analysis using random forests using the train data set ========================================
      x = data_train_temp[, predictors_new, drop = FALSE]
      y = data_train_temp[, predictand]
      forest_temp = randomForest(x = x, y = y, 
                                 ntree = ntree, 
                                 importance = TRUE)
    forest_list[[repeat_idx]] = forest_temp
    
  } # repeat_idx
  
  # combine random forests
  if (verbose) print("Combining random forests...")
  forest = do.call(combine, forest_list)
  
  # calculate the mean variable importance
  
  if (manual_subsetting == TRUE) {
    output$model_parameters$oob_error_from_RF = NA # this does not exist anymore for combined forests
    output$model_parameters$r_squared_from_RF = NA # this does not exist anymore for combined forests
  } else {
    output$model_parameters$oob_error_from_RF = forest$mse[ntree]
    output$model_parameters$r_squared_from_RF = forest$rsq[ntree]
  }
  
  output$model_parameters$n$n_total = n_total
  output$model_parameters$n$n_train = n_train
  output$model_parameters$n$n_test = n_test
  
  # prepare the variable importance dataframe
  importance_temp = as.data.frame(forest$importance)
  names(importance_temp) = c("importance_IncMSE", "importance_IncNodePurity")
  importance_temp$predictor = row.names(importance_temp)
  importance_temp = merge(data.frame(predictor = predictors), importance_temp, all = TRUE)
  importance_temp = melt(importance_temp, id.vars = "predictor",
                         measure.vars = c("importance_IncMSE", "importance_IncNodePurity"),
                         na.rm = FALSE, value.name = "value")
  # convert into one line
  importance_temp = dcast(importance_temp, 0  ~ predictor + variable)[, -1]
  
  # add to output model parameters
  output$model_parameters$importance_IncMSE = select(importance_temp, ends_with("importance_IncMSE"))
  output$model_parameters$importance_IncNodePurity = select(importance_temp, ends_with("importance_IncNodePurity"))
  # if training data is not full dataset, then length of data is different from original length
  time_series$fitted_values = napredict(na.action(data), predict(object = forest, newdata = data[, predictors_new, drop = FALSE]))
  time_series$residuals = time_series$fitted_values - time_series$original
  
  # put time series into output
  output$time_series = list()
  if (is.null(train)) {
    output$time_series$train = time_series
    output$time_series$test = NULL
  } else {
    output$time_series$train = time_series[train == 1, ]
    output$time_series$test = time_series[train == 0, ]
  }
  
  # add model to output as well, in case it's needed
  output$model = forest
  
  return(output)
}


###############################################################################
# function: calculate_statistical_parameters_from_fitted_values ----------------------------------------
# 
#   This function calculates the R-squared value for any given set of fitted values (model values) and observations by creating a linear model between the two and returning the R-squared value of this linear regression.

calculate_statistical_parameters_from_fitted_values = function(
  data_frame = NULL,
  parameter_list,
  fitted_values_var = "fitted_values", 
  observed_values_var = "original",
  significant_only = FALSE, # should the statistical parameters  be averaged over significant areas only (analogue to what Ray did in his variability paper)
  significance_threshold = 0.1, # threshold for significance, only relevant if significant_only = TRUE
  return_data = TRUE,
  overwrite_data = FALSE) {
  
  # read in parameters
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  if (verbose) print("function: calculate_statistical_parameters_from_fitted_values")
  
  # some tests
  assert_that(is.logical(significant_only))
  if (significant_only == TRUE) {
    assert_that(is.numeric(significance_threshold), !is.null(significance_threshold), 
                significance_threshold >= 0, significance_threshold <= 1)
  }
  
  # check if standardised statistical parameters were calculated before ========================================
  
  files = c(sprintf("%s/csv/yields_and_production_observed_and_fitted_%s.csv", output_dir,
                    c("original", "global", "un_region", "country")),
            sprintf("%s/csv/statistical_parameters_per_region_%s.csv", output_dir,
                    c("global", "continent", "un_region", "country")),
            sprintf("%s/netcdf/parameters_calculated_from_observed_and_fitted_values_%s.nc", output_dir,
                    c("train", "test")))
  
  if (significant_only == TRUE) {
    files = gsub(x = files, ".csv", "_significant_only.csv")
    files = gsub(x = files, ".nc", "_significant_only.nc")
  }
  
  if (all(sapply(files, file.exists)) && !overwrite_data) {
    if (verbose) print("Standardised parameters were calculated before.")
    
    if (return_data == TRUE) {
      # data only returned for normal metrics, not the Ray metrics
      production_yield_timeseries_per_region = list()
      for (geographic_scale in c("original", "global", "continent", "un_region", "country")) {
        fn = sprintf("%s/csv/yields_and_production_observed_and_fitted_%s.csv", output_dir, geographic_scale)
        if (significant_only == TRUE) fn = gsub(x = fn, ".csv", "_significant_only.csv")
        production_yield_timeseries_per_region[[geographic_scale]] = read.csv(fn, row.names = NULL)
      }
      
      statistical_parameters_per_region = list()
      for (geographic_scale in c("global", "continent", "un_region", "country")) {
        fn = sprintf("%s/csv/statistical_parameters_per_region_%s.csv", output_dir, geographic_scale)
        if (significant_only == TRUE) fn = gsub(x = fn, ".csv", "_significant_only.csv")
        statistical_parameters_per_region[[geographic_scale]] = read.csv(fn, row.names = NULL)
      }
      
      fn1 = sprintf("%s/netcdf/parameters_calculated_from_observed_and_fitted_values_train.nc", output_dir)
      fn2 = sprintf("%s/netcdf/parameters_calculated_from_observed_and_fitted_values_test.nc", output_dir)
      if (significant_only == TRUE) {
        fn1 = gsub(x = fn1, ".nc", "_significant_only.nc")
        fn2 = gsub(x = fn2, ".nc", "_significant_only.nc")
      }
      parameters_gridcell = rbind(
        netcdf2dataframe(netcdf_file = fn1, verbose = verbose),
        netcdf2dataframe(netcdf_file = fn2, verbose = verbose))
      parameters_gridcell = add_location_information(parameters_gridcell, res = res, verbose = verbose)
      
      return(list(parameters_gridcell = parameters_gridcell,
                  statistical_parameters_per_region = statistical_parameters_per_region,
                  production_yield_timeseries_per_region = production_yield_timeseries_per_region))
      
    } else {
      return(NULL)
    }
    
  }
  
  # calculate standardised values ========================================
  
  # if no data is provided, prepare dataframe from information in parameter list
  if (is.null(data_frame)) {
    fitted_values_residuals_temp = list()
    fn = sprintf("%s/csv/fitted_values_residuals_0_of_5.csv", parameter_list$output_dir)
    fitted_values_residuals_temp$train = read.csv(fn, row.names = NULL)
    fn = sprintf("%s/csv/fitted_values_residuals_test_combined.csv", parameter_list$output_dir)
    fitted_values_residuals_temp$test = read.csv(fn, row.names = NULL)
    
    # combine both datasets
    fitted_values_residuals_temp = do.call(rbind, fitted_values_residuals_temp)
    fitted_values_residuals_temp$train = grepl(pattern = "train", x = row.names(fitted_values_residuals_temp), fixed = TRUE)
    
    fn = sprintf("%s/data_used/data_used.csv", parameter_list$output_dir)
    data_used = read.csv(fn, row.names = NULL)
    yield_data = data_used %>% select(lon, lat, time, starts_with("yield")) %>% unique()
    
    n_before = nrow(fitted_values_residuals_temp)
    fitted_values_residuals_temp = merge(fitted_values_residuals_temp, yield_data, by = c("lon", "lat", "time"), all = FALSE)
    assert_that(n_before == nrow(fitted_values_residuals_temp))
    
    data_frame = fitted_values_residuals_temp
  }
  
  # some tests
  assert_that(is.character(fitted_values_var), is.character(observed_values_var))
  assert_that(fitted_values_var %in% names(data_frame), observed_values_var %in% names(data_frame))
  assert_that(all(c("lon", "lat", "time", "train") %in% names(data_frame)))
  
  # define function for calculating stats per grid cell
  calculate_stats = function(observed_values, fitted_values, absolute_yield) {
    
    results = data.frame(r_squared = NA, rmse = NA, 
                         rmse_rel_sd = NA, rmse_rel_mean = NA,
                         p_val_model = NA, model_significant = NA)
    
    if (all(is.na(observed_values)) || all(is.na(fitted_values)))
      return(results)
    
    # bring both datasets to the correct length (needed for ANOVA later, otherwise row length is different)
    temp = data.frame(observed_values, fitted_values)
    observed_values = observed_values[complete.cases(temp)]
    fitted_values = fitted_values[complete.cases(temp)]
    assert_that(length(observed_values) == length(fitted_values))
    
    fit = lm(observed_values ~ fitted_values, data = data.frame(observed_values, fitted_values))
    if(!is.null(summary(fit)$r.squared)) {
      results$r_squared = summary(fit)$r.squared
    }
    rmse = sqrt(mean((fit$residuals)^2, na.rm = TRUE))
    if (!is.null(rmse / sd(observed_values, na.rm = TRUE))) {
      results$rmse_rel_sd = rmse / sd(observed_values, na.rm = TRUE)
    }
    if (!is.null(rmse / mean(absolute_yield, na.rm = TRUE))) {
      results$rmse_rel_mean = rmse / mean(absolute_yield, na.rm = TRUE)
    }
    # calculate p-value for whole model
    fit.0 = lm(observed_values ~ 1) # null model
    an = anova(fit.0, fit)
    p_val_model = an$`Pr(>F)`[2]
    if (!is.null(p_val_model)) {
      results$p_val_model = p_val_model
      results$model_significant[results$p_val_model < significance_threshold] = 1
      results$model_significant[results$p_val_model >= significance_threshold] = 0
    }
    
    return(results)
  }
  
  # retrieve location data for each grid cell
  if (any(! c("continent", "un_region", "country") %in% names(data_frame)))
    data_frame = add_location_information(data_frame, res = res, verbose = verbose)
  
  # calculate R-squared value and RMSE for each grid cell ========================================
  if (verbose) print("Calculate R2 and RMSE for each grid cell from fitted data.")
  parameters_gridcell = ddply(.data = data_frame, 
                              .variables = c("train", "lon", "lat"), 
                              .fun = function(x)
                                calculate_stats(
                                  observed_values = x[, observed_values_var], 
                                  fitted_values = x[, fitted_values_var],
                                  absolute_yield = x[, "yield"]))
  
  # add location data to parameters per gridcell
  parameters_gridcell = add_location_information(data_frame = parameters_gridcell,
                                                 res = res, verbose = verbose)
  
  if (significant_only == TRUE) {
    # reduced dataframe for only significant grid cells
    significant_gridcells = parameters_gridcell %>%
      filter(model_significant == 1) %>%
      select(lon, lat, train) %>%
      unique()
    
    # reduce grid cells and data frame to only those grid cells with significant values
    parameters_gridcell = parameters_gridcell %>% filter(model_significant == 1)
    data_frame = merge(significant_gridcells, data_frame, all = FALSE)
  }
  
  if (nrow(data_frame) == 0) {
    
    if (return_data == TRUE) {
      return(list(parameters_gridcell = NULL,
                  statistical_parameters_per_region = NULL,
                  production_yield_timeseries_per_region = NULL))
    } else {
      return(NULL)
    }
  }
  
  
  
  # calculate the mean R-squared value for different geographic aggregation levels
  if (verbose) print("Calculate regional mean R2 and RMSE.")
  mean_parameters_per_region = list()
  
  for (geographic_scale in c("global", "continent", "un_region", "country")) {
    
    if (geographic_scale == "global") {
      grouping = "train"
    } else {
      grouping = c("train", geographic_scale)
    }
    
    temp = ddply(
      .data = parameters_gridcell,
      .variables = grouping,
      .fun = function(x) {
        return(data.frame(mean_r_squared = mean(x$r_squared, na.rm = TRUE),
                          mean_rmse_rel_sd = mean(x$rmse_rel_sd, na.rm = TRUE),
                          mean_rmse_rel_mean = mean(x$rmse_rel_mean, na.rm = TRUE)))
      })
    
    mean_parameters_per_region[[geographic_scale]] = temp
  }
  
  # calculate the R-squared value and RMSE for each geographic aggregation level for all values
  if (verbose) print("Calculate R2 and RMSE for each region from all fitted data.")
  parameters_per_region_from_all_values = list()
  
  for (geographic_scale in c("global", "continent", "un_region", "country")) {
    
    if (geographic_scale == "global") {
      grouping = "train"
    } else {
      grouping = c("train", geographic_scale)
    }
    
    temp = ddply(.data = data_frame, 
                 .variables = grouping, 
                 .fun = function(x) {
                   temp = calculate_stats(
                     observed_values = x[, observed_values_var], 
                     fitted_values = x[, fitted_values_var],
                     absolute_yield = x[, "yield"])
                   return(data.frame(r_squared_all = temp$r_squared,
                                     rmse_rel_sd_all = temp$rmse_rel_sd,
                                     rmse_rel_mean_all = temp$rmse_rel_mean))
                 })
    
    # remove lines with empty region
    if (geographic_scale != "global") {
      temp = temp[!is.na(temp[, geographic_scale]), ]
    }
    
    parameters_per_region_from_all_values[[geographic_scale]] = temp
  }
  
  
  # Calculate yield and production from fitted values
  
  if (verbose) print("Calculate yield and production from fitted values")
  
  production_yield_timeseries = calculate_yields_production_from_fitted_values(
    data_frame = data_frame, parameter_list = parameter_list, 
    overwrite_data = overwrite_data, significant_only = significant_only)
  
  # calculate the R-squared for each region (for all years and test years)
  statistical_parameters_per_region = list()
  
  for (geographic_scale in c("global", "continent", "un_region", "country")) {
    
    if (geographic_scale == "global") {
      grouping = "train"
    } else {
      grouping = c("train", geographic_scale)
    }
    
    temp = ddply(
      .data = production_yield_timeseries[[geographic_scale]],
      .variables = grouping,
      .fun = function(x) {
        
        # calculate parameters from detrended yields
        temp = x %>%
          # select(yield_detrended_observed, yield_detrended_fitted) %>%
          select(yield_detrended_observed = yield_anomalies_observed, 
                 yield_detrended_fitted = yield_anomalies_fitted, # aggr. anomalies instead of detrending of mean yields
                 yield_observed) %>%
          na.omit()
        n = nrow(temp)
        if (n <= 3) {
          r_squared_yield = NA
          rmse_yield_rel_sd = NA
          rmse_yield_rel_mean = NA
        } else {
          fit = lm(yield_detrended_observed ~ yield_detrended_fitted + 0, data = temp)
          summary_fit = summary(fit)
          r_squared_yield = summary_fit$r.squared
          rmse_yield = sqrt(mean((summary_fit$residuals)^2, na.rm = T))
          rmse_yield_rel_sd = rmse_yield / sd(temp$yield_detrended_observed, na.rm = TRUE)
          rmse_yield_rel_mean = rmse_yield / mean(temp$yield_observed, na.rm = TRUE)
        }
        
        # calculate parameters from detrended production
        temp = x %>%
          # select(production_detrended_observed, production_detrended_fitted) %>%
          select(production_detrended_observed = production_anomalies_observed,
                 production_detrended_fitted = production_anomalies_fitted,
                 production_observed) %>%
          na.omit()
        n = nrow(temp)
        if (n <= 3) {
          r_squared_production = NA
          rmse_production_rel_sd = NA
          rmse_production_rel_mean = NA
        } else {
          fit = lm(production_detrended_observed ~ production_detrended_fitted + 0, data = temp)
          summary_fit = summary(fit)
          r_squared_production = summary_fit$r.squared
          rmse_production = sqrt(mean((summary_fit$residuals)^2, na.rm = T))
          rmse_production_rel_sd = rmse_production / sd(temp$production_detrended_observed, na.rm = TRUE)
          rmse_production_rel_mean = rmse_production / mean(temp$production_observed, na.rm = TRUE)
        }
        
        return(data.frame(r_squared_yield = r_squared_yield,
                          rmse_yield_rel_sd = rmse_yield_rel_sd,
                          rmse_yield_rel_mean = rmse_yield_rel_mean,
                          r_squared_production = r_squared_production,
                          rmse_production_rel_sd = rmse_production_rel_sd,
                          rmse_production_rel_mean = rmse_production_rel_mean))
      })
    
    temp = merge(
      parameters_per_region_from_all_values[[geographic_scale]], temp)
    
    temp = merge(mean_parameters_per_region[[geographic_scale]], temp)
    
    # remove lines with empty region
    if (geographic_scale != "global") {
      temp = temp[!is.na(temp[, geographic_scale]), ]
    }
    
    statistical_parameters_per_region[[geographic_scale]] = temp
  }
  
  
  # write output data ========================================
  
  if (save_files) {
    
    if (verbose) print("Saving results.")
    
    # statistical parameters per region
    for (geographic_scale in c("global", "continent", "un_region", "country")) {
      temp = statistical_parameters_per_region[[geographic_scale]]
      file_out = sprintf("%s/csv/statistical_parameters_per_region_%s.csv", 
                         output_dir, geographic_scale)
      if (significant_only) file_out = gsub(x = file_out, ".csv", "_significant_only.csv")
      write.csv(temp, file_out, row.names = FALSE)
    }
    
    # prepare lon, lats and fill the grid cell parameters (if significant only, because otherwise the dataframe2netcdf function
    # cannot guess the resolution properly
    
    if (significant_only == TRUE) {
      # fill lon / lat, because it's often not possible anymore to guess the correct resolution
      lon_res = seq(min(data_frame$lon, na.rm = TRUE), max(data_frame$lon, na.rm = TRUE), res)
      
      if (res == 3.75) {
        res_temp = 2.5
      } else {
        res_temp = res
      }
      lat_res = seq(min(data_frame$lat, na.rm = TRUE), max(data_frame$lat, na.rm = TRUE), res_temp)
      coordinates_temp = expand.grid(lon = lon_res, lat = lat_res) %>% arrange(lon, lat)
      assert_that(all(parameters_gridcell$lon %in% lon_res))
      assert_that(all(parameters_gridcell$lat %in% lat_res))
    } 
    
    # grid cell parameters
    file_out = sprintf("%s/netcdf/parameters_calculated_from_observed_and_fitted_values_train.nc", output_dir)
    if (significant_only == TRUE) file_out = gsub(x = file_out, ".nc", "_significant_only.nc")
    
    parameters_train = parameters_gridcell %>% 
      filter(is.na(train) | train == TRUE) %>% 
      select(-continent, -un_region, -country, -train)
    
    if (significant_only == TRUE) # merge both
      parameters_train = merge(coordinates_temp, parameters_train, all.x = TRUE, all.y = FALSE)
    
    dataframe2netcdf(data_frame = parameters_train,
                     netcdf_file = file_out, 
                     overwrite_existing = TRUE, verbose = verbose)
    
    if (oos_calculation == TRUE) {
      file_out = sprintf("%s/netcdf/parameters_calculated_from_observed_and_fitted_values_test.nc", output_dir)
      if (significant_only == TRUE) file_out = gsub(x = file_out, ".nc", "_significant_only.nc")
      parameters_test = parameters_gridcell %>%
        filter(is.na(train) | train == FALSE) %>% 
        select(-continent, -un_region, -country, -train)
      if (significant_only == TRUE) # merge both
        parameters_test = merge(coordinates_temp, parameters_test, all.x = TRUE, all.y = FALSE)
      dataframe2netcdf(data_frame = parameters_test,
                       netcdf_file = file_out, 
                       overwrite_existing = TRUE, verbose = verbose)
    }
  }
  
  #@ to do: why are yield and production time series for significant only the same as not only for significant time series?
  
  
  if (return_data == TRUE) {
    return(list(parameters_gridcell = parameters_gridcell,
                statistical_parameters_per_region = statistical_parameters_per_region,
                production_yield_timeseries_per_region = production_yield_timeseries))
  } else {
    return(NULL)
  }
  
}


################################################################################

# function: calculate_yields_production_from_fitted_values

calculate_yields_production_from_fitted_values = function(
  data_frame = NULL,
  file_path = NULL,
  parameter_list,
  return_data = TRUE,
  significant_only,
  overwrite_data = FALSE) {
  
  # read in parameters
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  # some tests
  assert_that(!is.null(significant_only), !is.na(significant_only), is.logical(significant_only))
  
  if (crop %in% c("winter_wheat", "spring_wheat")) {
    crop_standard = "wheat"
  } else {
    crop_standard = crop
  }
  
  if (verbose) print("function: calculate_yields_production_from_fitted_values")
  
  
  # check if yields and production have not been calculated yet ========================================
  
  files = sprintf("%s/csv/yields_and_production_observed_and_fitted_%s.csv", output_dir,
                  c("original", "global", "continent", "un_region", "country"))
  if (significant_only == TRUE) files = gsub(x = files, ".csv", "_significant_only.csv")
  
  if (all(sapply(files, file.exists)) && !overwrite_data) {
    if (verbose) print("Production and yield (from fitted values) were calculated before.")
    
    if (return_data == TRUE) {
      production_yield_timeseries = list()
      for (geographic_scale in c("original", "global", "continent", "un_region", "country")) {
        fn = sprintf("%s/csv/yields_and_production_observed_and_fitted_%s.csv", output_dir,
                     geographic_scale)
        if (significant_only == TRUE) fn = gsub(x = fn, ".csv", "_significant_only.csv")
        production_yield_timeseries[[geographic_scale]] = read.csv(fn, row.names = NULL)
      }
      return(production_yield_timeseries)
    } else {
      return(NULL)
    }
  }
  
  
  # if production and yield time series were not calculated before ========================================
  
  # some tests
  assert_that(!is.null(data_frame) || (!is.null(file_path) && is.readable(file_path)))
  
  if (is.null(data_frame)) {
    if (verbose) print("Reading in fitted and observed yield anomalies")
    data_frame = read.csv(file_path, row.names = NULL)
    data_frame = add_location_information(data_frame, res = res, verbose = verbose)
  }
  
  
  assert_that(all(c("lon", "lat", "time", "train",
                    "continent", "un_region", "country",
                    "yield", "yield_trend", "fitted_values") 
                  %in% names(data_frame)))
  assert_that(standardise_variables == FALSE || "yield_detrended_sd" %in% names(data_frame))
  
  
  # calculate yield from fitted values ========================================
  
  data_frame = data_frame %>%
    mutate(yield_observed_all = yield, # all original yield values
           yield_trend_all = yield_trend,
           yield_detrended_original_all = yield_detrended_original,
           yield_observed = ifelse(is.na(fitted_values), NA, yield_observed_all), # those yield values, where fitted values have values as well
           yield_trend = ifelse(is.na(fitted_values), NA, yield_trend),
           yield_detrended_original = ifelse(is.na(fitted_values), NA, yield_detrended_original_all)) %>% # remove yield values, where there are no fitted values, otherwise, the weighted means will be different
    select(-yield)
  
  if (standardise_variables == TRUE) {
    data_frame$yield_anomalies_fitted = data_frame$fitted_values * data_frame$yield_detrended_sd
  } else {
    data_frame$yield_anomalies_fitted = data_frame$fitted_values
  }
  data_frame$yield_fitted = data_frame$yield_trend + data_frame$yield_anomalies_fitted
  
  
  
  # calculate production for observed and fitted values ========================================
  
  if (source_area_harvested == "deepak") {
    
    if (verbose) print("Reading in area harvested data: Ray et al data")
    filename = sprintf("%s/crop_yield_data/deepak_crop_data_regridded/deepak_%s_area_harvested_1961_2008_%sdeg.nc", data_path, crop_standard, res)
    print(filename)
    area_harvested = netcdf2dataframe(filename, time_format = "%Y")
    area_harvested$level = NULL
    data_frame = merge(data_frame, area_harvested, by = c("lon", "lat", "time"), all.x = TRUE, all.y = FALSE)
  } else if (source_area_harvested == "mirca") {
    
    if (verbose) print("Reading in area harvested data: MIRCA2000 data")
    # read in MIRCA land use area
    filename = sprintf("%s/landuse_data/mirca2000_regridded/mirca_area_harvested_total_%s_ha_%sdeg.nc", data_path, crop_standard, res)
    area_harvested = netcdf2dataframe(filename)
    data_frame = merge(data_frame, area_harvested, by = c("lon", "lat"), all.x = TRUE, all.y = FALSE)
  }
  
  data_frame$production_observed = data_frame$yield_observed * data_frame$area_harvested
  data_frame$production_fitted = data_frame$yield_fitted * data_frame$area_harvested
  
  
  # calculate yield anomalies, yield and production time series for different geographic aggregation levels 
  
  production_yield_timeseries = list()
  
  if (verbose) print("Aggregate yield and production for regions and calculate regional R2 and RMSE.")
  
  for (geographic_scale in c("global", "continent", "un_region", "country")) {
    if (geographic_scale == "global") {
      grouping = c("train", "time")
    } else {
      grouping = c("train", geographic_scale, "time")
    }
    
    production_yield_timeseries[[geographic_scale]] = ddply(
      .data = data_frame,
      .variables = grouping,
      .fun = function(x) {
        yield_observed = weighted.mean(x$yield_observed, x$area_harvested, na.rm = TRUE)
        yield_fitted = weighted.mean(x$yield_fitted, x$area_harvested, na.rm = TRUE)
        production_observed = sum(x$production_observed, na.rm = TRUE)
        production_fitted = sum(x$production_fitted, na.rm = TRUE)
        
        yield_anomalies_observed = weighted.mean(x$yield_detrended_original,
                                                 x$area_harvested, na.rm = TRUE)
        yield_anomalies_fitted = weighted.mean(x$yield_anomalies_fitted, 
                                               x$area_harvested, na.rm = TRUE)
        
        production_anomalies_observed = sum(x$yield_detrended_original * x$area_harvested, 
                                            na.rm = TRUE)
        production_anomalies_fitted = sum(x$yield_anomalies_fitted * x$area_harvested, na.rm = TRUE)
        
        yield_trend = weighted.mean(x$yield_trend, x$area_harvested, na.rm = TRUE)
        production_trend = sum(x$yield_trend * x$area_harvested, na.rm = TRUE)
        
        yield_all_observed = weighted.mean(x$yield_observed_all, x$area_harvested, na.rm = TRUE)
        production_all_observed = sum(x$yield_observed_all * x$area_harvested, na.rm = TRUE)
        yield_trend_all_observed = weighted.mean(x$yield_trend_all, x$area_harvested, na.rm = TRUE)
        production_trend_all_observed = sum(x$yield_trend_all * x$area_harvested, na.rm = TRUE)
        
        return(data.frame(yield_observed = yield_observed,
                          yield_fitted = yield_fitted,
                          production_observed = production_observed,
                          production_fitted = production_fitted,
                          yield_anomalies_observed = yield_anomalies_observed,
                          yield_anomalies_fitted = yield_anomalies_fitted,
                          production_anomalies_observed = production_anomalies_observed,
                          production_anomalies_fitted = production_anomalies_fitted,
                          yield_trend = yield_trend,
                          production_trend = production_trend,
                          yield_all_observed = yield_all_observed,
                          production_all_observed = production_all_observed,
                          yield_trend_all_observed = yield_trend_all_observed,
                          production_trend_all_observed = production_trend_all_observed))
      })
    
    production_yield_timeseries[[geographic_scale]][, ".id"] = NULL
  }
  
  
  
  # I'm not detrending production and yield time series anymore, but rather use the aggregated (weighted mean) of yield and production anomalies per grid cell
  
  # detrend production and yield time series for different geographic aggregation levels ========================================
  
  production_yield_timeseries_detrended = list()
  
  for (geographic_scale in c("global", "continent", "un_region", "country")) {
    
    if (geographic_scale == "global") {
      grouping = c("train")
    } else {
      grouping = c("train", geographic_scale)
    }
    
    production_yield_timeseries_detrended[[geographic_scale]] = ddply(
      .data = production_yield_timeseries[[geographic_scale]],
      .variables = grouping,
      .fun = function(x) {
        yield_observed_detrended = detrending(x = x$yield_observed,
                                              detrending_type = detrending_yield)
        yield_fitted_detrended = detrending(x = x$yield_fitted,
                                            detrending_type = detrending_yield)
        production_observed_detrended = detrending(x = x$production_observed,
                                                   detrending_type = detrending_yield)
        production_fitted_detrended = detrending(x = x$production_fitted,
                                                 detrending_type = detrending_yield)
        
        return(data.frame(time = x$time,
                          yield_observed_detrended = yield_observed_detrended,
                          yield_fitted_detrended = yield_fitted_detrended,
                          production_observed_detrended = production_observed_detrended,
                          production_fitted_detrended = production_fitted_detrended
        ))
      })
    
    production_yield_timeseries_detrended[[geographic_scale]][, ".id"] = NULL
    
    production_yield_timeseries[[geographic_scale]] = merge(production_yield_timeseries[[geographic_scale]],
                                                            production_yield_timeseries_detrended[[geographic_scale]])
    
    # remove lines with empty region
    if (geographic_scale != "global") {
      production_yield_timeseries[[geographic_scale]] = production_yield_timeseries[[
        geographic_scale]][!is.na(production_yield_timeseries[[geographic_scale]][, geographic_scale]), ]
    }
  }
  
  
  #     # add train column back to data
  #     for (field in c("global", "continent", "un_region", "country")) {
  #       production_yield_timeseries[[field]] = merge(unique(select(data_frame, time, train)),
  #                                                    production_yield_timeseries[[field]],
  #                                                    by = "time", all.x = FALSE, all.y = TRUE)
  #     }
  
  # add original original grid-cell time series as well ========================================
  production_yield_timeseries[["original"]] = data_frame
  
  
  # save files ========================================
  
  if (save_files) {
    
    for (field in c("global", "continent", "un_region", "country", "original")) {
      
      # save as csv
      fn = sprintf("%s/csv/yields_and_production_observed_and_fitted_%s.csv", output_dir, field)
      if (significant_only == TRUE) fn = gsub(x = fn, ".csv", "_significant_only.csv")
      write.csv(production_yield_timeseries[[field]], fn, row.names = FALSE)
      
      # save the original data as netcdf
      if (field == "original") {
        fn = sprintf("%s/netcdf/yields_and_production_observed_and_fitted_%s.nc", output_dir, field)
        if (significant_only == TRUE) fn = gsub(x = fn, ".nc", "_significant_only.nc")
        temp = production_yield_timeseries[[field]]
        
        # select test data (original data is the same, but fitted data is different)
        if (oos_calculation == TRUE) {
          
          temp = temp[temp$train == 0, ]
          
          if (standardise_variables == TRUE) {
            dataframe2netcdf(select(temp, -continent, -un_region, -country, -yield_detrended_sd), 
                             netcdf_file = fn, overwrite_existing = TRUE)
          } else {
            dataframe2netcdf(select(temp, -continent, -un_region, -country), 
                             netcdf_file = fn, overwrite_existing = TRUE)
          }
        }
      }
    }
    
    if (return_data == TRUE) {
      return(production_yield_timeseries)
    } else {
      return(NULL)
    }
  }
  
}

################################################################################
# function: calculate_model_parameters_fitted_values_residuals
# This function calculates the model parameters and fitted values/residuals for both random forests and multiple regression,
# with or without grouping variables
# 
# Reads in:
# data_all: dataframe with all data, needs to include predictors and predictand
# file_path: or alternatively, the file path to the dataset to be used
# parameter_list: parameter list with all settings
# return_data: should results just be saved or also returned to the calling function?
# overwrite_data: if results were calculated before, should they be overwritten?


calculate_model_parameters_fitted_values_residuals = function(
  data_all = NULL,
  file_path = NULL,
  parameter_list,
  return_data = TRUE,
  overwrite_data = FALSE) {
  
  # read in parameters ========================================
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  if (verbose) print("function: calculate_model_parameters_fitted_values_residuals")
  
  # some tests
  assert_that(!is.null(data_all) || (!is.null(file_path) && is.readable(file_path)))
  assert_that(!is.null(data_all) || is.data.frame(data_all))
  
  ################################################################################
  
  # define all model parameter types to write for each cross_val_iteration as well as a mean over all cross_val_iterations
  all_types = c("importance_IncMSE", "importance_IncNodePurity",
                "n", "oob_error_from_RF", "r_squared_from_RF")
  
  ################################################################################
  
  # prepare result lists for each cross_val_iteration
  model_parameters = list()
  fitted_values_residuals = list()
  model = list()
  
  ################################################################################
  # go through cross_val_iterations
  
  for (cross_val_iteration in 0:n_groups_out_of_bag_sampling) {
    model_parameters_cross_val_iteration = list()
    fitted_values_residuals_cross_val_iteration = list()
    model_cross_val_iteration = list()
    
    if (verbose) {
      if (cross_val_iteration == 0)
        if (verbose) print(sprintf("cross_val_iteration 0/%s (all data)", n_groups_out_of_bag_sampling))
      else
        if (verbose) print(sprintf("cross_val_iteration %s/%s", cross_val_iteration, n_groups_out_of_bag_sampling))
    }
    
    # prepare the train vector
    if (cross_val_iteration > 0) {
      # those rows for which oos_group == cross_val_iteration, are in the test dataset
      data_all$train = data_all$oos_group != cross_val_iteration 
      cross_val_years = unique(data_all$time[data_all$train == FALSE])
    } else {
      data_all$train = NULL
      cross_val_years = "all years"
    }
    
    if (verbose) print(sprintf("years for cross validation (out-of-sample): %s", paste(cross_val_years, collapse = ", ")))
    
    ################################################################################
    
    # check if model parameters and fitted values/residuals already exist
    files = c(sprintf("%s/csv/results_%s_%s_of_%s.csv", output_dir, all_types, 
                      cross_val_iteration, n_groups_out_of_bag_sampling)) # add cross_val_iteration number to file name
    
    if (cross_val_iteration > 0) {
      files = c(files, 
                sprintf("%s/csv/fitted_values_residuals_train_%s_of_%s.csv",
                        output_dir, cross_val_iteration, n_groups_out_of_bag_sampling),
                sprintf("%s/csv/fitted_values_residuals_test_%s_of_%s.csv", 
                        output_dir, cross_val_iteration, n_groups_out_of_bag_sampling))
    } else {
      files = c(files, 
                sprintf("%s/csv/fitted_values_residuals_%s_of_%s.csv",
                        output_dir, cross_val_iteration, n_groups_out_of_bag_sampling))
    }
    
    if (all(sapply(files, file.exists)) && !overwrite_data) {
      
      if (verbose) print("Model parameters and fitted values/residuals were calculated before (skipping calculation).")
      
      # read in model parameters for this cross_val_iteration
      model_parameters_cross_val_iteration = list()
      for (type in all_types) { 
        fn = sprintf("%s/csv/results_%s_%s_of_%s.csv", output_dir, type, cross_val_iteration,
                     n_groups_out_of_bag_sampling)
        model_parameters_cross_val_iteration[[type]] = read.csv(fn, row.names = NULL)
      }
      
      # read in fitted values / residuals for this cross_val_iteration
      if (cross_val_iteration > 0) {
        fitted_values_residuals_cross_val_iteration = list()
        fn = sprintf("%s/csv/fitted_values_residuals_train_%s_of_%s.csv", output_dir, cross_val_iteration,
                     n_groups_out_of_bag_sampling)
        fitted_values_residuals_cross_val_iteration$train = read.csv(fn, row.names = NULL)
        fn = sprintf("%s/csv/fitted_values_residuals_test_%s_of_%s.csv", output_dir, cross_val_iteration,
                     n_groups_out_of_bag_sampling)
        fitted_values_residuals_cross_val_iteration$test = read.csv(fn, row.names = NULL)
      } else {
        fitted_values_residuals_cross_val_iteration = list()
        fn = sprintf("%s/csv/fitted_values_residuals_%s_of_%s.csv", output_dir, cross_val_iteration, 
                     n_groups_out_of_bag_sampling)
        fitted_values_residuals_cross_val_iteration$train = read.csv(fn, row.names = NULL)
        fitted_values_residuals_cross_val_iteration$test = NULL
      }
      
    } else { # if this cross_val_iteration has not been calculated before
      
      # if the model parameters and fitted values / residuals didn't exist before
      if (is.null(data_all)) {
        data_all = read.csv(file_path, row.names = NULL)
      }
      
      ################################################################################
      # Calculating model parameters and fitted values / residuals
      
      if (verbose) print("Calculating model parameters and fitted values/residuals")
      
      # if grouping_by == "global" ----------------------------------------
      if (is.null(grouping_vars)) {
        
        results_temp = random_forest(data_all,
                                     predictand = predictand,
                                     predictors = predictors,
                                     need_all_vars = need_all_vars, 
                                     train = data_all$train,
                                     verbose = ifelse(is.null(grouping_vars), TRUE, FALSE),
                                     manual_subsetting = manual_subsetting,
                                     sampling_fraction = sampling_fraction,
                                     extremes_dataset = extremes_dataset,
                                     res = res)
        
        # remove variable type from names
        for (type in c("importance_IncMSE", "importance_IncNodePurity")) {
          model_parameters_cross_val_iteration[[type]] = results_temp$model_parameters[[type]]
          names(model_parameters_cross_val_iteration[[type]]) =
            gsub(sprintf("_%s", type), "", names(model_parameters_cross_val_iteration[[type]]))
        }
        model_parameters_cross_val_iteration$oob_error_from_RF = data.frame(oob_error_from_RF = results_temp$model_parameters$oob_error_from_RF)
        model_parameters_cross_val_iteration$r_squared_from_RF = data.frame(r_squared_from_RF = results_temp$model_parameters$r_squared_from_RF)
        
        # grouping vars NULL & any type of analysis
        
        model_cross_val_iteration = results_temp$model
        model_parameters_cross_val_iteration$n = results_temp$model_parameters$n
        
        # fitted values and residuals
        fitted_values_residuals_cross_val_iteration = results_temp$time_series
        
      } else { # if grouping_vars non-empty ----------------------------------------
        
        # create strings from numbers, so that afterwards (after do.call(bind, ...)) they can be retrieved again from rownames
        if (grouping_by == "coordinates") {
          if (is.numeric(data_all$lat))
            data_all$lat = sprintf("%.2f", data_all$lat)
          if (is.numeric(data_all$lon))
            data_all$lon = sprintf("%.2f", data_all$lon)
        }
        
        lon_pos = which(grouping_vars == "lon")
        lat_pos = which(grouping_vars == "lat")
        
        results_temp = dlply(
          .data = data_all,
          .variables = grouping_vars,
          .progress = "text",
          .fun = function(x) {
            
            if (verbose) print(sprintf("- %s", unique(x[, grouping_vars])))
            
            temp = random_forest(
              x,
              predictand = predictand,
              predictors = predictors,
              need_all_vars = need_all_vars, 
              verbose = ifelse(is.null(grouping_vars), TRUE, FALSE),
              train = x$train,
              manual_subsetting = manual_subsetting,
              sampling_fraction = sampling_fraction,
              extremes_dataset = extremes_dataset,
              res = res)
            return(temp)
          })
        
        # combine all model parameters
        for (type in names(results_temp[[1]][["model_parameters"]])) {
          if (verbose) print(type)
          temp = do.call(rbind, lapply(results_temp, function(x) x$model_parameters[[type]]))
          assert_that(!is.null(temp))
          
          # for oob_error and r_squared --> create dataframe
          if (!is.data.frame(temp)) { 
            temp = as.data.frame(temp)
            assert_that(ncol(temp) == 1)
            names(temp) = type
          }
          
          # add names for grouping variables (e.g. continent names)
          if (length(grouping_vars) == 1) {
            temp[, grouping_vars] = row.names(temp)
          } else { # with lat/lon
            assert_that(all(grouping_vars %in% c("lon", "lat")))
            temp$lon = sapply(row.names(temp), function(x) strsplit(x, ".", fixed = TRUE)[[1]][which(grouping_vars == "lon")])
            temp$lat = sapply(row.names(temp), function(x) strsplit(x, ".", fixed = TRUE)[[1]][which(grouping_vars == "lat")])
          }
          
          row.names(temp) = NULL # has to be removed, so one can identify the cross_val_iterations number later (when merging all cross_val_iterations together)
          
          model_parameters_cross_val_iteration[[type]] = temp
        }
        
        for (type in c("importance_IncMSE", "importance_IncNodePurity")) {
          names(model_parameters_cross_val_iteration[[type]]) =
            gsub(sprintf("_%s", type), "", names(model_parameters_cross_val_iteration[[type]]))
        }
        
        # for which grid cells / groups do we have any results?
        model_parameters_available = model_parameters_cross_val_iteration$r_squared_from_RF$r_squared_from_RF
        model_parameters_available[model_parameters_available > - 900] = 1 # includes -1000, -999 etc. for no values and 1 for results available
        
        # how many groups were tested? (particularly relevant for grouping by grid cells)
        number_of_groups_tested = sum(model_parameters_available, na.rm = TRUE)
        
        
        # Combine fitted values / residuals from every cross_val_iteration / grouping var ----------------------------------------      
        if (verbose) print("Combining fitted values / residuals")
        
        # concatenate the training and out-of-sample test time series
        time_series_train = do.call(rbind, lapply(results_temp, function(x) x$time_series$train))
        time_series_test = do.call(rbind, lapply(results_temp, function(x) x$time_series$test))
        
        # add grouping_vars (only applicable if analysis is grouped by continents etc.)
        if (!is.null(time_series_train)) {
          if (length(grouping_vars) == 1) {
            time_series_train[, grouping_vars] = sapply(row.names(time_series_train), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
          } else {
            assert_that(all(grouping_vars %in% c("lon", "lat")))
            time_series_train$lon = as.numeric(sapply(row.names(time_series_train), function(x) paste(strsplit(x, ".", fixed = TRUE)[[1]][c((2*lon_pos)-1, 2*lon_pos)], collapse = ".")))
            time_series_train$lat = as.numeric(sapply(row.names(time_series_train), function(x) paste(strsplit(x, ".", fixed = TRUE)[[1]][c((2*lat_pos)-1, 2*lat_pos)], collapse = ".")))
          }
          row.names(time_series_train) = NULL # has to be removed, so one can identify the cross_val_iterations number later (when merging all cross_val_iterations together)
        }
        
        if (!is.null(time_series_test)) {
          if (length(grouping_vars) == 1) {
            time_series_test[, grouping_vars] = sapply(row.names(time_series_test), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
          } else {
            assert_that(all(grouping_vars %in% c("lon", "lat")))
            time_series_test$lon = as.numeric(sapply(row.names(time_series_test), function(x) paste(strsplit(x, ".", fixed = TRUE)[[1]][c((2*lon_pos)-1, 2*lon_pos)], collapse = ".")))
            time_series_test$lat = as.numeric(sapply(row.names(time_series_test), function(x) paste(strsplit(x, ".", fixed = TRUE)[[1]][c((2*lat_pos)-1, 2*lat_pos)], collapse = ".")))
          }
          row.names(time_series_test) = NULL
        }
        
        fitted_values_residuals_cross_val_iteration$train = time_series_train
        fitted_values_residuals_cross_val_iteration$test = time_series_test
        
        fitted_values_residuals_cross_val_iteration$test[, ".id"] = NULL
        fitted_values_residuals_cross_val_iteration$train[, ".id"] = NULL
        
        # Combine models together (not for coordinates) ----------------------------------------
        model_cross_val_iteration = lapply(results_temp, function(x) return(x$model))
        
      } # end of if grouping vars == "global"
      
      # put data into list
      cross_val_iteration_str = as.character(cross_val_iteration)
      model_parameters[[cross_val_iteration_str]] = model_parameters_cross_val_iteration
      fitted_values_residuals[[cross_val_iteration_str]] = fitted_values_residuals_cross_val_iteration
      # model[[cross_val_iteration_str]] = model_cross_val_iteration # not needed for now
      
      
      # write model parameters ========================================
      
      assert_that(all(all_types %in% names(model_parameters_cross_val_iteration)))
      
      if (save_files) {
        
        if (verbose) print("Write model model parameters, fitted values / residuals and model")
        
        for (type in all_types) {
          fn = sprintf("%s/csv/results_%s_%s_of_%s.csv", output_dir, 
                       type, cross_val_iteration, n_groups_out_of_bag_sampling)
          if (verbose) print(sprintf("Saving %s to: %s", type, fn))
          write.csv(model_parameters_cross_val_iteration[[type]], fn, row.names = FALSE)
        }
        
        if (!is.null(data_all$train)) {
          fn = sprintf("%s/csv/fitted_values_residuals_train_%s_of_%s.csv", output_dir, 
                       cross_val_iteration, n_groups_out_of_bag_sampling)
          if (verbose) print(sprintf("Saving fitted_values_residuals_train to: %s", fn))
          write.csv(fitted_values_residuals_cross_val_iteration$train, fn, row.names = FALSE)
          fn = sprintf("%s/ csv/fitted_values_residuals_test_%s_of_%s.csv", output_dir, 
                       cross_val_iteration, n_groups_out_of_bag_sampling)
          if (verbose) print(sprintf("Saving fitted_values_residuals_test to: %s", fn))
          write.csv(fitted_values_residuals_cross_val_iteration$test, fn, row.names = FALSE)
        } else {
          fn = sprintf("%s/csv/fitted_values_residuals_%s_of_%s.csv", output_dir, cross_val_iteration, n_groups_out_of_bag_sampling)
          if (verbose) print(sprintf("Saving fitted_values_residuals_train to: %s", fn))
          write.csv(fitted_values_residuals_cross_val_iteration$train, fn, row.names = FALSE)
        }
        
        # write model
        fn = sprintf("model/model_%s_of_%s.RData", cross_val_iteration, n_groups_out_of_bag_sampling)
        if (verbose) print(sprintf("Saving model to: %s", fn))
        assert_that(!is.null(model_cross_val_iteration)) # for some reason, the model was NULL sometimes
        saveRDS(object = model_cross_val_iteration, file = fn)
      }
      
    } # if all data existed
    
  } # cross_val_iterations of cross-validation
  
  ################################################################################
  # change lon, lat back to numeric values
  if (grouping_by == "coordinates") {
    if (is.character(data_all$lat))
      data_all$lat = as.numeric(data_all$lat)
    if (is.character(data_all$lon))
      data_all$lon = as.numeric(data_all$lon)
  }
  
  ################################################################################
  # combine all cross_val_iterations into one for both model parameters and fitted values/residuals
  
  if (verbose) print("Combine the model parameters from all cross_val_iterations.")
  
  # check if all cross_val_iterations are in the model_parameters / fitted values list (if some were skipped, those need to be read in)
  for (cross_val_iteration in 0:n_groups_out_of_bag_sampling) {
    cross_val_iteration_str = as.character(cross_val_iteration)
    
    # read in missing model parameters
    if (is.null(model_parameters[[cross_val_iteration_str]])) { # if this cross_val_iteration was not calculated here, but rather before, read in the data
      model_parameters[[cross_val_iteration_str]] = list()
      for (type in all_types) {
        fn = sprintf("%s/csv/results_%s_%s_of_%s.csv", output_dir, type, cross_val_iteration, n_groups_out_of_bag_sampling)
        model_parameters[[cross_val_iteration_str]][[type]] = read.csv(fn, row.names = NULL)
      }
    }
    
    # read in missing fitted values / residuals
    if (is.null(fitted_values_residuals[[cross_val_iteration_str]])) {
      if (cross_val_iteration > 0) { # only for test/train combination
        fn = sprintf("%s/csv/fitted_values_residuals_train_%s_of_%s.csv", output_dir, cross_val_iteration, n_groups_out_of_bag_sampling)
        fitted_values_residuals[[cross_val_iteration_str]]$train = read.csv(fn, row.names = NULL)
        fn = sprintf("%s/csv/fitted_values_residuals_test_%s_of_%s.csv", output_dir,
                     cross_val_iteration, n_groups_out_of_bag_sampling)
        fitted_values_residuals[[cross_val_iteration_str]]$test = read.csv(fn, row.names = NULL)
      } else {
        fn = sprintf("%s/csv/fitted_values_residuals_%s_of_%s.csv", output_dir, 
                     cross_val_iteration, n_groups_out_of_bag_sampling)
        fitted_values_residuals[[cross_val_iteration_str]]$train = read.csv(fn, row.names = NULL)
        fitted_values_residuals[[cross_val_iteration_str]]$test = NULL
      }
    }
    
    # model not needed for moment, so not read in again
  }
  
  
  ################################################################################
  # for model parameters: calculate mean of all model parameters
  
  # check if everything has been calculated before --> if yes, just read it in again
  files = sprintf("%s/csv/results_%s_all.csv", output_dir, all_types)
  files = c(files, sprintf("%s/csv/fitted_values_residuals_test_combined.csv", output_dir))
  
  if (all(sapply(files, file.exists))) {
    if (verbose) print("Combined model parameters and fitted values have been calculated before (reading them in).")
    
    for (type in all_types) {
      fn = sprintf("%s/csv/results_%s_all.csv", output_dir, type)
      model_parameters[["all"]][[type]] = read.csv(fn, row.names = NULL)
    }
    
    fn = sprintf("%s/csv/fitted_values_residuals_test_combined.csv", output_dir)
    fitted_values_residuals[["test_combined"]] = read.csv(fn, row.names = NULL)
    
    
  } else {
    
    # if not everything has been calculated before
    if (verbose) print("Calculating mean values of model parameters.")
    
    model_parameters[["all"]] = list()
    
    for (type in all_types) {
      
      if (verbose) print(sprintf("- %s", type))
      
      # combine all cross_val_iterations into one dataframe
      model_parameters_all = do.call(rbind, lapply(model_parameters, FUN = function(x) x[[type]]))
      
      
      # calculate the mean over each columns
      if (all(grouping_vars %in% names(model_parameters_all))) { # if each parameter was calculated for each group independently
        mean = ddply(.data = model_parameters_all, 
                     .variables = grouping_vars,
                     .fun = function(x) {
                       if (!is.null(grouping_vars))
                         x = x[, !names(x) %in% grouping_vars, drop = FALSE] # remove grouping vars (if there are any, from)
                       apply(x, 2, function(y) mean(y, na.rm = TRUE))
                     })
        mean[, ".id"] = NULL
        mean$cross_val_iteration = "mean"
        
        cv = ddply(.data = model_parameters_all, 
                   .variables = grouping_vars,
                   .fun = function(x) {
                     if (!is.null(grouping_vars))
                       x = x[, !names(x) %in% grouping_vars, drop = FALSE] # remove grouping vars (if there are any, from)
                     apply(x, 2, function(y) {sd(y) / mean(y, na.rm = TRUE) * 100})
                   })
        cv[, ".id"] = NULL
        cv$cross_val_iteration = "cv"
        
      } else { # for aggregated parameters, especially in_model_pct, in_model_aggregated_pct
        mean = ddply(.data = model_parameters_all, 
                     .variables = NULL,
                     .fun = function(x) {
                       if (!is.null(grouping_vars))
                         x = x[, !names(x) %in% grouping_vars, drop = FALSE] # remove grouping vars (if there are any, from)
                       apply(x, 2, function(y) mean(y, na.rm = TRUE))
                     })
        
        mean[, ".id"] = NULL
        mean$cross_val_iteration = "mean"
        
        cv = ddply(.data = model_parameters_all, 
                   .variables = NULL,
                   .fun = function(x) {
                     if (!is.null(grouping_vars))
                       x = x[, !names(x) %in% grouping_vars, drop = FALSE] # remove grouping vars (if there are any, from)
                     apply(x, 2, function(y) {sd(y) / mean(y, na.rm = TRUE) * 100})
                   })
        cv[, ".id"] = NULL
        cv$cross_val_iteration = "cv"
      }
      
      model_parameters_all$cross_val_iteration = sapply(row.names(model_parameters_all), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
      model_parameters_all = model_parameters_all[, names(mean)] # bring into same order
      model_parameters_all = rbind(model_parameters_all, mean = mean, cv = cv)
      model_parameters_all = select(model_parameters_all, cross_val_iteration, everything()) # put cross_val_iteration as 1st column
      
      model_parameters[["all"]][[type]] = model_parameters_all
    }
    
    # for fitted values / residuals: for test data, combine all time series
    temp = fitted_values_residuals
    temp[["0"]] = NULL
    
    if (oos_calculation == TRUE) {
      fitted_values_residuals_test_all = do.call(rbind, lapply(temp, FUN = function(x) x$test))
      # make comparable to data_all by using same time, lon, lat combinations
      fitted_values_residuals_test_all = merge(unique(data_all[, c("time", "lon", "lat")]), fitted_values_residuals_test_all, by = c("time", "lon", "lat"), all = TRUE)
      fitted_values_residuals_test_all = arrange(fitted_values_residuals_test_all, time, lon, lat)
      fitted_values_residuals[["test_combined"]] = fitted_values_residuals_test_all
      
      # compare original time series
      compareNA <- function(v1,v2) {
        same <- (v1 == v2) | (is.na(v1) & is.na(v2))
        same[is.na(same)] <- FALSE
        return(same)
      }
      data_all = arrange(data_all, time, lon, lat)
      assert_that(all(compareNA(data_all[, predictand], fitted_values_residuals_test_all$original)))
    }
    
    # saving combined model parameters and fitted values and residuals ========================================
    if (save_files) {
      
      if (verbose) print("Saving combined model parameters and fitted values/residuals.")
      
      for (type in all_types) {
        fn = sprintf("%/csv/results_%s_all.csv", output_dir, type)
        write.csv(model_parameters[["all"]][[type]], fn, row.names = FALSE)
      }
      
      if (oos_calculation == TRUE) {
        fn = sprintf("%s/csv/fitted_values_residuals_test_combined.csv", output_dir)
        write.csv(fitted_values_residuals[["test_combined"]], fn, row.names = FALSE)
        fn = sprintf("%s/netcdf/fitted_values_residuals_test_combined.nc", output_dir)
        dataframe2netcdf(data_frame = fitted_values_residuals[["test_combined"]], netcdf_file = fn, overwrite_existing = TRUE)
      }
      
    } # saving data
  } # has combined / mean model parameters and fitted values been calculated before?
  
  # return results
  fitted_values_residuals[["train_using_all_data"]] = add_location_information(fitted_values_residuals[["0"]]$train, res = res, verbose = verbose) # one dataframe using all data
  if (oos_calculation) {
    fitted_values_residuals[["test_combined"]] = add_location_information(fitted_values_residuals[["test_combined"]], res = res, verbose = verbose)
  }
  fitted_values_residuals[!names(fitted_values_residuals) %in% c("train_using_all_data", "test_combined")] = NULL
  
  # set back to old working directory
  if (return_data == TRUE) {
    return(list(model_parameters = model_parameters[["all"]], fitted_values_residuals = fitted_values_residuals))
  } else {
    return(NULL)
  }
}


###############################################################################
# Helper function: create_manual_subset
# This functions manually subsets the grid cells used for the statistical analysis in order
# to reduce spatial correlations due to the re-gridding of the coarse-resolution HadEX2 data.
create_manual_subset = function(data_frame, sampling_fraction, verbose, extremes_dataset, res) {
  
  if (verbose) print("function: create_manual_subset")
  
  grid_cells_old = select(data_frame, lon, lat) %>% unique()
  n_grid_cells_old = nrow(grid_cells_old)
  
  # prepare n_sample: how many grid cells should be in sub sample? max means, sample grid cells until no grid cells left
  if (sampling_fraction == "max") {
    n_sample = n_grid_cells_old
  } else {
    n_sample = round(n_grid_cells_old * sampling_fraction)
  }
  
  # define lat/lon range that is considered too close
  if (is.na(extremes_dataset)) {
    # just remove neighboring cells, so lon/lat range depending on resolution
    if (res == 3.75) {
      res_range = c(3.75, 2.5) # lon, lat
    } else if (res == 2.5) {
      res_range = c(2.5, 2.5)
    } else if (res == 1.5) {
      res_range = c(1.5, 1.5)
    } else {
      res_range = c(0.5, 0.5)
    }
  } else if (extremes_dataset == "hadex2") {
    res_range = c(3.75, 2.5)
  } else if (extremes_dataset == "ghcndex") {
    res_range = c(2.5, 2.5)
  } else if (extremes_dataset == "era_interim") {
    res_range = c(1.5, 1.5)
  }
  # sample grid cells and remove neighbouring grid cells
  grid_cells_subset = NULL
  for (i in 1:n_sample) {
    
    if (nrow(grid_cells_old) == 0 && sampling_fraction == "max") {
      next
    } else if (nrow(grid_cells_old) == 0 && sampling_fraction != "max") {
      stop() # shouldn't be the case but who knows
    }
    
    # select any grid cell
    grid_cell_temp = grid_cells_old[sample(1:nrow(grid_cells_old), 1), ]
    grid_cells_subset = rbind(grid_cells_subset, grid_cell_temp)
    
    # remove neighbouring gridcells
    grid_cells_old = grid_cells_old %>%
      # either lon or lat out of neighbourhood bounds
      filter(lon < grid_cell_temp$lon - res_range[1] | 
               lon > grid_cell_temp$lon + res_range[1] |
               lat < grid_cell_temp$lat - res_range[2] |
               lat > grid_cell_temp$lat + res_range[2])
    
    if (nrow(grid_cells_old) == 0 && i < n_sample && sampling_fraction != "max") {
      if (verbose) print(sprintf("Only %s (%.2f) instead of %s (%.2f) samples. Stopping data preparation. Consider smaller sampling_fraction.", i, i/n_grid_cells_old, n_sample, sampling_fraction))
      stop(sprintf("Only %s (%.2f) instead of %s (%.2f) samples. Stopping data preparation. Consider smaller sampling_fraction.", i, i/n_grid_cells_old, n_sample, sampling_fraction))
    }
  }
  
  n_grid_cells_new = nrow(grid_cells_subset)
  if (verbose) print(sprintf("Sampled %s (%.2f) out of %s grid cells.", 
                             n_grid_cells_new, n_grid_cells_new/n_grid_cells_old, n_grid_cells_old))
  
  # reduce to only subset grid cells
  data_frame_subset = merge(data_frame, grid_cells_subset, all = FALSE)
  
  return(data_frame_subset)
}


###############################################################################
# Helper functions: transpose and transpose back (needed to combine statistical outputs into dataframes)
transpose_df = function(data_frame) {
  transposed_df = as.data.frame(t(data_frame))
  transposed_df = data.frame(colnames(data_frame), transposed_df)
  colnames(transposed_df) = c("var", rownames(data_frame))
  rownames(transposed_df) = NULL
  return(transposed_df)
}

transpose_back_df = function(transposed_df) {
  data_frame = as.data.frame(t(transposed_df))
  colnames(data_frame) = transposed_df[, 1]
  data_frame = data_frame[-1, ]
  row.names(data_frame) = NULL
  return(data_frame)
}

