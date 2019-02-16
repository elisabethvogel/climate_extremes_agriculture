
project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"

# function: create_stasticial_analysis_parameter_list
create_stasticial_analysis_parameter_list = function(
  crop = "wheat",
  yield_dataset = "deepak", 
  extremes_dataset = NA,
  climate_dataset = NA,
  soil_moisture_dataset = NA,
  drought_indicator_dataset = NA,
  crop_calendar = "agmip",
  irrigation = "combined",
  mask_growing_season_areas = TRUE,
  detrending_climate = "ssa",
  detrending_yield = "best_fit",
  
  extreme_indicator_vars = c(
    "TXx_gs_max", "TX90p_gs_mean",
    "TNn_gs_min", "TN10p_gs_mean",
    "Rx5day_gs_max", "Rx1day_gs_max"),
  
  climate_vars = c(   
    "tmp_gs_mean", "pre_gs_mean", 
    "frs_gs_mean", "dtr_gs_mean"),
  
  soil_moisture_vars = c(
    "SMm_gs_mean", "SMn_gs_min", "SMx_gs_max"),
  
  drought_vars = c("spi_3_gs_mean", "spi_6_gs_mean"),
  
  grouping_by = "coordinates",
  use_time_as_predictor = FALSE,
  use_location_as_predictor = FALSE,
  location_predictor = NA,
  use_only_regions_with_frost_days = FALSE,
  
  res = 1.5,
  year_min = 1961, year_max = 2008,
  
  time_range = "gs",
  lag = 0,
  shift_yield_time_by_x_years = 0,
  crop_to_region = "global",
  need_all_vars = TRUE,
  quadratic_terms = NULL,
  comment = format(now(), "%Y-%m-%d-%H%M%S"),
  statistical_method,
  stepwise_direction = NA,
  oos_calculation = TRUE,
  seed = NA,
  n_groups_out_of_bag_sampling = 5,
  save_files = TRUE,
  source_area_harvested = "deepak", # deepak or mirca
  standardise_variables = FALSE,
  verbose = FALSE,
  type_of_analysis = NULL,
  manual_subsetting = FALSE,
  sampling_fraction = NULL,
  include_interaction_terms = FALSE) {
  
  if (oos_calculation == FALSE) n_groups_out_of_bag_sampling = 0
  
  
  # PREPARE PARAMETERS FOR DEFAULT ANALYSIS - FULL AND REDUCED STATISTICAL MODEL
  if (!is.null(type_of_analysis)) {
    
    quadratic = ifelse(statistical_method == "stepwise_regression", TRUE, FALSE)
    comment = sprintf("analysis_%s", type_of_analysis)
    
    
    # Full statistical model (previously "analysis 8")
    if (type_of_analysis == "full_model") {
      
      extremes_dataset = "hadex2"
      climate_dataset = "cru_ts_323"
      soil_moisture_dataset = NA
      drought_indicator_dataset = "cru_ts_323_spi"
      
      extreme_indicator_vars = c(
        "TXx_gs_max", "TX90p_gs_mean",
        "TNn_gs_min", "TN10p_gs_mean",
        "Rx5day_gs_max")
      climate_vars = c(   
        "pre_gs_mean",
        "tmp_gs_mean", 
        "dtr_gs_mean",
        "frs_gs_mean")
      soil_moisture_vars = NULL
      drought_vars = c("spi_6_gs_mean")
      if (quadratic == TRUE) {
        quadratic_terms = c("tmp_gs_mean", "pre_gs_mean")
      } else {
        quadratic_terms = NULL
      }
      
      # Reduced statistical model - using only tmp and precip (previously "analysis 9")
    } else if (type_of_analysis == "reduced_model") {
      
      extremes_dataset = NA
      climate_dataset = "cru_ts_323"
      soil_moisture_dataset = NA
      drought_indicator_dataset = NA
      
      extreme_indicator_vars = NULL
      climate_vars = c("pre_gs_mean", "tmp_gs_mean")
      soil_moisture_vars = NULL
      drought_vars = NULL
      if (quadratic == TRUE) {
        quadratic_terms = c("tmp_gs_mean", "pre_gs_mean")
      } else {
        quadratic_terms = NULL
      }   
      
      
    } # type of statistical analysis
    rm(quadratic)
  } 
  
  if (time_range == "growing_season")
    time_range = "gs"
  
  if (mask_growing_season_areas == TRUE) {
    masked = "masked"
  } else {
    masked = ""
  }
  
  if (is.na(extremes_dataset)) extreme_indicator_vars = NULL
  if (is.na(climate_dataset)) climate_vars = NULL
  if (is.na(soil_moisture_dataset)) soil_moisture_vars = NULL
  if (is.na(drought_indicator_dataset)) drought_vars = NULL
  
  # define the grouping_vars
  if (grouping_by == "coordinates") {
    grouping_vars = c("lon", "lat")
  } else if (grouping_by == "global") {
    grouping_vars = NULL
  } else {
    grouping_vars = grouping_by
  }
  
  # prepare predictand and predictor variables ========================================
  predictand = sprintf("yield_detrended_%s", detrending_yield)
  predictors = c(sprintf("%s_detrended_%s", 
                         c(extreme_indicator_vars, climate_vars, 
                           drought_vars, soil_moisture_vars),
                         detrending_climate), 
                 sprintf("%s_detrended_%s_2", quadratic_terms, detrending_climate))
  
  # add time and lat/lon to predictors, if needed
  if (use_time_as_predictor == TRUE)
    predictors = c("time", predictors)
  if (use_location_as_predictor == TRUE) {
    if (location_predictor == "coordinates") {
      predictors = c("lon", "lat", predictors)
    } else if (location_predictor == "continent") {
      predictors = c("continent", predictors)
    } else if (location_predictor == "un_region") {
      predictors = c("un_region", predictors)
    } else if (location_predictor == "country") {
      predictors = c("fao_country", predictors)
    }
  } else {
    location_predictor = NA
  }
  
  parameter_list = as.list(environment())
  parameter_list$output_dir = create_output_dir(parameter_list)
  parameter_list$output_dir_general = strsplit(parameter_list$output_dir, split = "deg")[[1]][1] # first part for general agricultural
  parameter_list$output_dir_general = sprintf("%sdeg", parameter_list$output_dir_general)
  
  # test the parameters
  test_parameters(parameter_list)
  
  # save parameter list
  fn = sprintf("%s/parameter_list/parameter_list.R", parameter_list$output_dir)
  dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
  save(parameter_list, file = fn)
  return(parameter_list)
}


###############################################################################
# helper function: test parameters

test_parameters = function(parameter_list) {
  
  source(sprintf("%s/code/other/file_information.R", project_path))
  
  # read in parameters
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  # some tests
  assert_that(length(parameter_list) == 49) # number of expected parameters
  assert_that(crop %in% c("wheat", "winter_wheat", "spring_wheat",
                          "maize", "rice", "soybeans"))
  if (crop %in% c("winter_wheat", "spring_wheat"))
    assert_that(use_only_regions_with_frost_days == FALSE) # frost regions is already some kind of distinction of wheat / spring wheat
  
  if (time_range == "growing_season") time_range = "gs"
  assert_that(time_range %in% c("gs", "annual"))
  
  assert_that(yield_dataset %in% data_sets$yield)
  assert_that(is.na(extremes_dataset) || extremes_dataset %in% data_sets$extreme_indicators)
  assert_that(is.na(climate_dataset) || climate_dataset %in% data_sets$climate_data)
  assert_that(is.na(soil_moisture_dataset) || soil_moisture_dataset %in% data_sets$soil_moisture)
  assert_that(is.na(drought_indicator_dataset) || drought_indicator_dataset %in% data_sets$drought_indicators)
  
  # create test variables (remove the time range and growing season statistic)
  test_extreme_indicator_vars = sapply(extreme_indicator_vars, function(x) 
    strsplit(x, sprintf("_%s", time_range))[[1]][1])
  test_climate_vars = sapply(climate_vars, function(x)
    strsplit(x, sprintf("_%s", time_range))[[1]][1])
  test_soil_moisture_vars = sapply(soil_moisture_vars, function(x)
    strsplit(x, sprintf("_%s", time_range))[[1]][1])
  test_drought_vars = sapply(drought_vars, function(x)
    strsplit(x, sprintf("_%s", time_range))[[1]][1])
  
  assert_that(all(test_extreme_indicator_vars %in% variables$extreme_indicators[[extremes_dataset]]))
  assert_that(all(test_climate_vars %in% variables$climate_data[[climate_dataset]]))
  assert_that(all(test_soil_moisture_vars %in% variables$soil_moisture[[soil_moisture_dataset]]))
  assert_that(all(test_drought_vars %in% variables$drought_indicators[[drought_indicator_dataset]]))
  
  assert_that(!is.na(crop_calendar %in% data_sets$crop_calendars))
  assert_that(irrigation %in% c("rainfed", "irrigated", "combined"))
  assert_that(is.logical(mask_growing_season_areas))
  assert_that(masked %in% c("masked", ""))
  
  assert_that(res %in% c("0.5", "2.5", "3.75", "1.5"))
  
  # interaction terms only when pre and tmp in predictors, and only for regression
  assert_that(include_interaction_terms == FALSE || 
                all(sprintf(c("pre_gs_mean_detrended_%s", "tmp_gs_mean_detrended_%s"), detrending_climate) %in% predictors))
  assert_that(include_interaction_terms == FALSE || statistical_method == "stepwise_regression") 
  
  if (is.na(crop_to_region)) 
    crop_to_region = "global"
  
  #   assert_that(crop_to_region %in% c("global", "Africa", "Asia", "Europe",
  #                                     "North America", "South America", "Australia"))
  assert_that(grouping_by %in% c("coordinates", "continent", "un_region", "country", "global"))
  assert_that(use_location_as_predictor == FALSE || !is.na(location_predictor))
  assert_that(is.na(location_predictor) || location_predictor %in%
                c("coordinates", "continent", "un_region", "country"))
  assert_that(is.na(location_predictor) || (location_predictor != grouping_by))
  
  
  assert_that(is.numeric(lag))
  assert_that(is.character(output_dir), is.writeable(output_dir))
  assert_that(is.character(output_dir_general), is.writeable(output_dir_general))
  
  assert_that(is.character(predictand), length(predictand) == 1)
  assert_that(is.character(predictors), length(predictors) >= 1)
  
  assert_that(is.logical(use_only_regions_with_frost_days))
  assert_that(is.logical(oos_calculation))
  if (oos_calculation == FALSE) assert_that(n_groups_out_of_bag_sampling == 0)
  assert_that(is.numeric(shift_yield_time_by_x_years), 
              floor(shift_yield_time_by_x_years) == shift_yield_time_by_x_years)
  assert_that(is.logical(manual_subsetting))
  
  
  if (manual_subsetting == TRUE) {
    assert_that(!is.null(sampling_fraction))
    assert_that((is.character(sampling_fraction) && sampling_fraction == "max") ||
                  (is.numeric(sampling_fraction) && (sampling_fraction > 0 || sampling_fraction < 1)))
    
    assert_that(statistical_method == "random_forest") # not needed for stepwise regression and thereforenot implemented
  }
  
  return(1)
}

###############################################################################
# helper function: create_output_dir
create_output_dir = function(parameter_list) {
  
  # read in parameters
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  # prepare predictor variables
  predictors = c(extreme_indicator_vars, climate_vars, drought_vars, soil_moisture_vars, 
                 sprintf("%s_2", quadratic_terms))
  predictors = predictors[!is.na(predictors)]
  
  # add interaction terms to predictors
  if (include_interaction_terms == TRUE) {
    # interaction term is written out as, for example, pre_gs_mean_detrended_lin:tmp_gs_mean_detrended_lin
    predictors = c(predictors, "pre-X-tmp_gs_mean")
  }
  predictors = sort(predictors)
  
  if (statistical_method == "random_forest") {
    method = "random_forest"
  } else if (statistical_method == "stepwise_regression") {
    method = sprintf("stepwise_%s_regression", stepwise_direction)
  }
  
  output_dir = sprintf(
    "%s/results/statistical_analysis/%s_%s_%s_%s_%s_%sdeg/%s_%s_%s_%s_%s/%s/%s_%s_%s_%s_lag_%s/%s/%s_%s/%s/%s",
    
    project_path,
    crop, crop_calendar, irrigation, masked, ifelse(use_only_regions_with_frost_days == TRUE, "frost_only", ""), res,
    ifelse(is.na(yield_dataset), "", yield_dataset), 
    ifelse(is.na(extremes_dataset), "", extremes_dataset), 
    ifelse(is.na(climate_dataset), "", climate_dataset), 
    ifelse(is.na(soil_moisture_dataset), "", soil_moisture_dataset), 
    ifelse(is.na(drought_indicator_dataset), "", drought_indicator_dataset),
    method,
    crop_to_region, year_min, year_max,
    time_range, lag,
    ifelse(shift_yield_time_by_x_years == 0, "", sprintf("shift_by_%s_yrs", shift_yield_time_by_x_years)),
    comment, ifelse(manual_subsetting == TRUE, sprintf("manual_subsample_%s", sampling_fraction), ""),
    
    
    paste(
      c("detrending_yield", detrending_yield, 
        "climate", detrending_climate,
        "std_vars", as.character(standardise_variables),
        "group_by", grouping_by,
        "time_pred", as.character(use_time_as_predictor), 
        "location_pred", as.character(use_location_as_predictor),
        ifelse(is.na(location_predictor), "", location_predictor),
        "need_all_vars", as.character(need_all_vars), 
        "seed", as.character(seed)), collapse = "_"),
    
    paste(predictors, collapse = "_"))
  
  # remove double underscores or underscores at the beginning/end of a folder name
  while (grepl("__", output_dir))
    output_dir = gsub("__", "_", output_dir)
  while (grepl("/_", output_dir))
    output_dir = sub("/_", "/", output_dir)
  while (grepl("_/", output_dir))
    output_dir = sub("_/", "/", output_dir)
  
  dir.create(sprintf("%s/csv", output_dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(sprintf("%s/netcdf", output_dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(sprintf("%s/model", output_dir), showWarnings = FALSE, recursive = TRUE)
  
  return(output_dir)
}
