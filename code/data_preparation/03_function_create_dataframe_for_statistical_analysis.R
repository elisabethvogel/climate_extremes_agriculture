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
data_path = sprintf("%s/data", project_path)

# function: create_dataframe_for_statistical_analysis ----------------------------------------
#
#   This function creates one large dataframe for the analysis
#
#   Reads in:
#     - parameter_list that contains information on
#       - crop: "wheat", "maize", "rice" or "soybeans"
#       - yield_dataset: "deepak" or "iizumi"
#       - extremes_dataset: "hadex2" or "ghcndex"
#       - climate_dataset: "cru_ts_323"
#       - soil_moisture_dataset: "esa_soil_moisture"
#       - drought_indicator_dataset: "ncar_pdsi"
#       - crop_calendar: "agmip"
#       - irrigation: "combined", "irrigated", "rainfed"
#       - masked: "masked", "_masked", ""
#       - detrending_yield: "lin", "quad", best_fit", "incr", "ssa"
#       - detrending_climate: "lin", "quad", "incr", "ssa"
#       - extreme_indicator_vars, climate_vars, soil_moisture_vars, drought_vars: variable
#         names to be used as predictors for the statistical analysis, full name
#         with growing season statistic, e.g. tmp_gs_mean, TXx_gs_max
#       - resolution: 0.5, 2.5, 3.75
#       - year_min, year_max: time range for the analysis, automatically it will be
#         reduced to the years where all datasets have data
#       - remove_empty_cells: logical, if TRUE, grid cells (i.e. lat/lon combinations)
#         with no data will be removed from the data frame
#
#     - remove_NA: logical, basically like na.omit
#     - return_format: "wide", "long"
#     - add_meta_data: TRUE, FALSE - if TRUE, adding information about crop, crop
#       calendar, resolution etc. to the data frame
#     - save_output: logical, if TRUE, dataframe will be saved to a csv file
#     - output_file: filename where to save output, optional, will be automatically
#       created otherwise
#     - verbose
#
#   Returns:
#     - data_frame with all data

create_dataframe_for_statistical_analysis =
  function(parameter_list, remove_empty_cells = TRUE,
           remove_NA = FALSE,  return_format = "wide",
           add_meta_data = TRUE, overwrite_existing = FALSE) {
    
    # read in parameters
    for (i in 1:length(parameter_list)) {
      assign(names(parameter_list[i]), parameter_list[[i]])
    }
    
    if (verbose)
      print("function: create_dataframe_for_statistical_analysis")
    
    resolution = as.character(res)
    
    # some tests
    assert_that(return_format %in% c("wide", "long"),
                is.logical(remove_empty_cells),
                is.logical(add_meta_data),
                is.logical(overwrite_existing))
    
    # check if file exists already ----------------------------------------
    fn = sprintf("%s/data_used/data_used.csv", output_dir)
    if (overwrite_existing == FALSE && file.exists(fn)) {
      if (verbose) print(sprintf("File already existed. Reading in: %s", fn))
      data = read.csv(fn, row.names = NULL)
      return(data)
    }
    
    # check if previously was confirmed that no data is available for this country/crop combination
    fn = sprintf("%s/data_used/data_used_NULL.csv", output_dir)
    if (overwrite_existing == FALSE && file.exists(fn)) {
      if (verbose) print("No data for this country available.")
      return(NULL)
    }
    
    # else start ----------------------------------------
    
    if (crop %in% c("winter_wheat", "spring_wheat")) {
      crop_standard = "wheat" 
    } else {
      crop_standard = crop
    }
    
    
    ###############################################################################
    # SUBSETTING
    
    # subset the data to a specific continent
    grid_cells = subset_to_region(crop_to_region, res)
    
    # read in landuse data and find those grid cells where irrigated/rainfed >= 80% of
    # area harvested (if applicable)
    grid_cells = subset_to_irrigation_pattern(irrigation, res, crop_standard, grid_cells = grid_cells)
    
    # reduce grid cells to those where winter / spring wheat is grown
    grid_cells = subset_to_spring_winter_wheat(crop, res, grid_cells = grid_cells)
    
    # if (use_only_regions_with_frost_days == TRUE) reduce grid cells to those where 
    # frost day frequency does not have a standard deviation of zero
    
    if (use_only_regions_with_frost_days) {
      # read in standard deviation data of CRU TS 323 FRS
      fn = sprintf("%s/growing_season_data_detrended/climate_data/cru_ts_323/%s_%s_%s_%s/%s/%s_%s_cru_ts_323_frs_%s_%sdeg_stdev.nc", 
                   data_path, crop_calendar, crop, irrigation, masked, 
                   ifelse(lag == 0, "original", sprintf("lag_%s", lag)),
                   crop, irrigation, time_range, res)
      var = sprintf("frs_%s_mean_detrended_%s%s", time_range, detrending_climate,
                    ifelse(lag == 0, "", sprintf("_lag_%s", lag)))
      frs_stdev = netcdf2dataframe(netcdf_file = fn, variables = var,
                                   grid_cells = grid_cells, remove_NA = remove_NA, verbose = verbose)
      
      # remove all grid cells where the standard deviation of FRS is zero, i.e. likely there are no frost days ever (technically this is not correct, but it's unlikely that FRS is non-zero and constant over all years)
      frs_stdev = frs_stdev[!is.na(frs_stdev[, var]) & frs_stdev[, var] != 0, ]
      grid_cells = select(frs_stdev, lon, lat)
      assert_that(nrow(grid_cells) > 0)
    }
    
    
    ###############################################################################
    # COMBINING EVERYTHING
    
    # put everything into a vector in order to better iterate through it
    data_set = list()
    data_set$yield = yield_dataset
    data_set$extreme_indicators = extremes_dataset
    data_set$climate_data = climate_dataset
    data_set$soil_moisture = soil_moisture_dataset
    data_set$drought_indicators = drought_indicator_dataset
    
    # create variable names
    variable_names = list()
    variable_names$yield = sprintf("yield_detrended_%s", detrending_yield)
    variable_names$extreme_indicators = sprintf("%s_detrended_%s", extreme_indicator_vars, detrending_climate)
    variable_names$climate_data = sprintf("%s_detrended_%s", climate_vars, detrending_climate)
    variable_names$soil_moisture = sprintf("%s_detrended_%s", soil_moisture_vars, detrending_climate)
    variable_names$drought_indicators = sprintf("%s_detrended_%s", drought_vars, detrending_climate)
    
    
    # to do: add correct variable names with gs and lag
    for (var_type in names(variable_names)) {
      variable_names[[var_type]] = sprintf("%s%s", variable_names[[var_type]],
                                           ifelse(lag == 0, "", sprintf("_lag_%s", lag)))
    }
    
    # prepare all file paths
    dir = prepare_climate_data_file_paths(data_set, yield_dataset, crop_calendar, crop_standard, irrigation, masked)
    
    ###############################################################################3
    # READING IN
    data_frame = data.frame()
    for (data_type in c("yield", "extreme_indicators", "climate_data",
                        "soil_moisture", "drought_indicators")) {
      
      if (is.na(data_set[[data_type]]))
        next
      
      for (variable in variable_names[[data_type]]) {
        
        # create filename
        if (data_type == "yield") {
          if (yield_dataset == "deepak") {
            filename = sprintf("%s/%s_%s_yield_1961_2008_%sdeg.nc",
                               dir[[data_type]], data_set[[data_type]], crop_standard, resolution)
          } else if (yield_dataset == "iizumi") {
            filename = sprintf("%s/%s_%s_yield_1982_2006_%sdeg.nc",
                               dir[[data_type]], data_set[[data_type]], crop_standard, resolution)
          }
        } else {
          temp_var = strsplit(variable, "_gs_") [[1]][1] # get rid of the "_gs_something" stuff, because that's not in the filename
          filename = sprintf(
            "%s/%s_%s_%s_%s_%s_%sdeg.nc",
            dir[[data_type]], crop_standard, irrigation, data_set[[data_type]], temp_var, time_range, resolution
          )
        }
        
        # if standardise data, adjust the filename
        if (standardise_variables == TRUE) {
          filename = gsub(x = filename, "deg.nc", "deg_standardised.nc")
        }
        
        if (verbose) print(sprintf("Reading in: %s", filename))
        
        data_temp = netcdf2dataframe(netcdf_file = filename, variables = variable,
                                     grid_cells = grid_cells, years = year_min:year_max, remove_NA = remove_NA, 
                                     time_format = "%Y", verbose = verbose)
        
        assert_that(nrow(data_temp) > 0)
        
        # add data to overall dataframe
        if (nrow(data_frame) == 0) {
          data_frame = data_temp
        } else {
          data_frame = merge_or_cbind(data_frame, data_temp, by = c("lon", "lat", "time"),
                                      verbose = verbose, all = FALSE)
        }
      }
    }
    
    # add:
    # - yield_detrended_sd: standard deviation of yield anomalies
    # - yield_detrended_original: absolute anomalies of yield (not standardised)
    # - yield_original: absolute yields (not detrended)
    # - yield_trend: trend of yield (absolute yields minus detrended yields)
    
    # yield_detrended_sd ----------------------------------------
    if (standardise_variables == TRUE) {
      
      if (yield_dataset == "deepak") {
        filename = sprintf("%s/%s_%s_yield_1961_2008_%sdeg_stdev.nc",
                           dir$yield, data_set$yield, crop_standard, resolution)
      } else if (yield_dataset == "iizumi") {
        filename = sprintf("%s/%s_%s_yield_1982_2006_%sdeg_stdev.nc",
                           dir$yield, data_set$yield, crop_standard, resolution)
      }
      
      if (verbose) print(sprintf("Reading in: %s", filename))
      
      data_temp = netcdf2dataframe(netcdf_file = filename, variables = variable_names$yield,
                                   grid_cells = grid_cells, remove_NA = remove_NA, 
                                   verbose = verbose)
      
      # rename the variable to "yield_detrended_sd"
      data_temp$yield_detrended_sd = data_temp[, variable_names$yield]
      data_temp[, variable_names$yield] = NULL
      
      assert_that(nrow(data_temp) > 0)
      data_frame = merge(data_frame, data_temp, all = FALSE)
    }
    
    # yield_detrended_original and yield_original ----------------------------------------
    if (yield_dataset == "deepak") {
      filename = sprintf("%s/%s_%s_yield_1961_2008_%sdeg.nc",
                         dir$yield, data_set$yield, crop_standard, resolution)
    } else if (yield_dataset == "iizumi") {
      filename = sprintf("%s/%s_%s_yield_1982_2006_%sdeg.nc",
                         dir$yield, data_set$yield, crop_standard, resolution)
    }
    
    if (verbose) print(sprintf("Reading in: %s", filename))
    
    data_temp = netcdf2dataframe(netcdf_file = filename, variables = c("yield", variable_names$yield),
                                 remove_NA = remove_NA, time_format = "%Y", years = year_min:year_max,
                                 grid_cells = grid_cells, verbose = verbose)
    
    # rename the detrended variable to "yield_detrended_original"
    data_temp$yield_detrended_original = data_temp[, variable_names$yield]
    data_temp[, variable_names$yield] = NULL
    
    assert_that(nrow(data_temp) > 0)
    data_frame = merge_or_cbind(data_frame, data_temp, all = FALSE)
    
    # yield_trend ----------------------------------------
    data_frame$yield_trend = data_frame$yield - data_frame$yield_detrended_original
    
    # add location information ========================================
    if (verbose) print("Adding location information.")
    data_frame = add_location_information(data_frame = data_frame, res = res, verbose = verbose)
    
    # create quadratic terms and add them to data set ========================================
    if (!is.null(quadratic_terms)) {
      for (var in sprintf("%s_detrended_%s", quadratic_terms, detrending_climate)) {
        assert_that(var %in% names(data_frame))
        data_frame[, sprintf("%s_2", var)] = (data_frame[, var]) * (data_frame[, var])
      }
    }
    
    # if the analysis is for a reduced statistical model, then subset the data to the 
    # same grid cells as for the equivalent complete statistical model
    data_frame = subset_gridcells_for_reduced_statistical_model(data_frame, parameter_list)
    
    ###############################################################################
    # PREPARE OUTPUT FORMAT (wide, remove empty grid cells)
    
    # melt dataframe
    data_frame_long = melt(data_frame, 
                           id.vars = c("lon", "lat", "time", "continent_code", "continent",
                                       "un_region_code", "un_region", "country_code", "country"),
                           variable.name = "variable")
    
    # remove gridcells that have no values at all
    test_na = function(x)
      all(is.na(x))
    
    data_frame_long = data_frame_long %>%
      group_by(lon, lat) %>%
      mutate(all_empty = test_na(value)) %>%
      ungroup() %>%
      filter(all_empty == FALSE) %>%
      select(-all_empty)
    
    if (nrow(data_frame_long) == 0) {
      print("No data for this country available.")
      fn = sprintf("%s/data_used/data_used_NULL.csv", output_dir)
      write.csv("", fn) # save it so that it isn't tried the next time again
      return(NULL)
    }
    
    # dcast into wide format again
    data_frame_wide = dcast(data_frame_long, lon + lat + time + continent_code +
                              continent + un_region_code + un_region + country_code + country ~ variable,
                            value.var = "value")
    return_data = data_frame_wide
    
    # order data by lat, lon and time
    return_data = return_data %>% arrange(lon, lat, time)
    
    # adding meta data to the data frame
    if (add_meta_data) { 
      return_data$crop = crop
      return_data$crop_calendar = crop_calendar
      return_data$irrigation = irrigation       
      return_data$masked = sub("_", "", masked)
      return_data$resolution = res
      return_data = select(return_data, lon, lat, time, crop, crop_calendar,
                           irrigation, masked, resolution, continent_code,
                           continent, un_region_code, un_region, country_code, country, everything())
      
    }
    
    # save data if needed
    
    if (save_files) {
      fn = sprintf("%s/data_used/data_used.csv", output_dir)
      dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
      write.csv(return_data, fn, row.names = FALSE)
      
      if (verbose) print(sprintf("Data frame saved as csv at: %s", fn))
      
      if (return_format == "wide") {
        fn = sprintf("%s/data_used/data_used.nc", output_dir)
        dataframe2netcdf(data = return_data, netcdf_file = fn, overwrite_existing = overwrite_existing)
      }
    }
    
    return(return_data)
  }


###############################################################################
# function: add_location_information
#
#   This function adds continent, un_region and country to dataframe that contains lat/lon coordinates.
#
#   Reads in:
#     - data_frame that has lat / lon coordinates
#     - res: resolution of the data
#   Returns:
#     - data_frame with lat, lon, continent, un_region, country information

add_location_information = function(data_frame, res, verbose = FALSE) {
  
  if (verbose) print("function: add_location_information")
  
  # some tests
  assert_that(is.data.frame(data_frame))
  assert_that(all(c("lat", "lon") %in% names(data_frame)))
  
  # add location information
  filename = sprintf("%s/location_metadata/location_metadata_%sdeg.csv", data_path, res)
  if (verbose) print(sprintf("Reading in: %s", filename))
  
  location = read.csv(filename, row.names = NULL)
  data_frame = data_frame[, !grepl("continent", names(data_frame)) &
                            !grepl("un_region", names(data_frame)) &
                            !grepl("country", names(data_frame)) ]
  data_frame = merge(data_frame, location, by = c("lon", "lat"), all.x = TRUE, all.y = FALSE)
  
  if ("time" %in% names(data_frame)) {
    data_frame = data_frame %>%
      arrange(lon, lat, time)
  } else {
    data_frame = data_frame %>%
      arrange(lon, lat)
  }
  
  
  return(data_frame)
  
}


###############################################################################
# Helper function - return grid cells for a given region
subset_to_region = function(crop_to_region, res) {
  
  if (crop_to_region != "global") {
    
    # check if it's a country or continent
    fn = sprintf("%s/data/fao_countries_nc/fao_country_codes_with_continents.csv", project_path)
    country_continents = read.csv(fn, row.names = NULL)
    
    # for a continent
    if (crop_to_region %in% country_continents$continent) {
      
      if (verbose) print(sprintf("Read in grid cells that belong to continent: %s", crop_to_region))
      fn = sprintf("%s/data/fao_un_region_nc/un_region_continent_%sdeg.nc", project_path, res)
      grid_cells = netcdf2dataframe(fn, variables = "continent_code")
      fn = sprintf("%s/data/fao_un_region_nc/fao_un_region_codes.csv", project_path)
      meta = read.csv(fn, row.names = NULL)
      meta = meta %>%
        select(continent_code = Continent_Code, continent = Continent) %>%
        unique()
      grid_cells = merge(grid_cells, meta)
      grid_cells = grid_cells %>%
        filter(continent == crop_to_region) %>%
        select(lon, lat)
      assert_that(nrow(grid_cells) > 0)
      
      # for a country
    } else if (crop_to_region %in% country_continents$country) {
      
      if (verbose) print(sprintf("Read in grid cells that belong to country: %s", crop_to_region))
      fn = sprintf("%s/data/fao_countries_nc/fao_countries_%sdeg.nc", project_path, res)
      grid_cells = netcdf2dataframe(fn, variables = "fao_country_code")
      fn = sprintf("%s/data/fao_countries_nc/fao_country_codes.csv", project_path)
      meta = read.csv(fn, row.names = NULL)
      meta = meta %>%
        select(fao_country_code = FAO_CODE, country = COUNTRY) %>%
        unique()
      grid_cells = merge(grid_cells, meta)
      grid_cells = grid_cells %>%
        filter(country == crop_to_region) %>%
        select(lon, lat)
      assert_that(nrow(grid_cells) > 0)
      
    } else {
      print("No country or continent found under this name. Returning NULL.")
      fn = sprintf("%s/data_used/data_used_NULL.csv", output_dir)
      write.csv("", fn)
      return(NULL)
    }
    
  } else {
    grid_cells = NULL # read in data for the whole world
  }
}


###############################################################################
# Helper function: subset to irrigated / rainfed grid cells
subset_to_irrigation_pattern = function(irrigation, res, crop_standard, threshold = 0.8,
                                        grid_cells = NULL) {
  
  if (irrigation %in% c("rainfed", "irrigated")) {
    # area fraction of this irrigation type
    fn = sprintf("%s/landuse_data/mirca2000_regridded/mirca_area_harvested_%s_%s_ha_%sdeg.nc", data_path, irrigation, crop_standard, res)
    area_irrigation = netcdf2dataframe(fn, grid_cells = grid_cells)
    
    # total area harvested
    fn = sprintf("%s/landuse_data/mirca2000_regridded/mirca_area_harvested_total_%s_ha_%sdeg.nc", data_path, crop_standard, res)
    area_total = netcdf2dataframe(fn, grid_cells = grid_cells)
    names(area_total)[names(area_total) == "area_harvested"] = "area_total"
    
    area_irrigation = merge_or_cbind(area_irrigation, area_total)
    area_irrigation = area_irrigation %>% 
      mutate(fraction = area_harvested / area_total) %>%
      filter(fraction >= threshold)
    
    grid_cells = select(area_irrigation, lon, lat)
    assert_that(nrow(grid_cells) > 0)
    assert_that(nrow(unique(grid_cells)) == nrow(grid_cells))
  }
  
  return(grid_cells)
}


###############################################################################
# Helper function: Subset to spring / winter wheat
subset_to_spring_winter_wheat = function(crop, res, grid_cells = NULL) {
  
  if (crop %in% c("winter_wheat", "spring_wheat")) {
    
    # classification by Sacks et al.
    fn = sprintf("%s/landuse_data/sacks_wheat_classification/sacks_wheat_classification_%sdeg.nc", data_path, res)
    grid_cells = netcdf2dataframe(fn, variables = crop, grid_cells = grid_cells)
    grid_cells$classification = grid_cells[, crop]
    grid_cells = grid_cells %>%
      filter(classification == 1) %>%
      select(lon, lat) %>%
      unique()
    
    if (nrow(grid_cells) == 0) {
      if (verbose) print("No grid cells for this crop type.")
      stop("No grid cells for this crop type.")
    }
    
    return(grid_cells)
  }
}


###############################################################################
# Helper function: Prepare all file paths to read in the climate data
# create empty object for file paths
prepare_climate_data_file_paths = function(data_set, yield_dataset,
                                           crop_calendar, crop_standard,
                                           irrigation, masked, lag=0) {
  
  # prepare empty list for file paths
  dir = list()
  
  if (masked == "masked") masked = "_masked"
  
  for (data_type in c(
    "yield", "extreme_indicators", "climate_data", "soil_moisture", "drought_indicators"
  )) {
    if (is.na(data_set[[data_type]])) {
      # if no dataset selected
      dir[[data_type]] = NA
    } else {
      if (data_type == "yield") {
        dir[[data_type]] = sprintf(
          "%s/growing_season_data_detrended/%s_yield_data",
          data_path, yield_dataset
        )
      } else {
        dir[[data_type]] = sprintf(
          "%s/growing_season_data_detrended/%s/%s/%s_%s_%s%s/%s",
          data_path, data_type, data_set[[data_type]], crop_calendar, crop_standard, irrigation, masked,
          ifelse(lag == 0, "original", sprintf("lag_%s", lag))
        )
      }
    }
  }
  
  return(dir)
}



###############################################################################
# Helper function: Subset the grid cells for the reduced statistical model

# Background: The reduced statistical model (using only mean temperature and precipitation)
# generally has more widespread coverage of gridcells globally, because the data is limited
# by fewer variables. The allow for a fair comparison between the reduced and full
# statistical model, this function can be used to subset the grid cells to the same ones
# that were used for a previously prepared data set for the complete statistical model.

subset_gridcells_for_reduced_statistical_model = function(data_frame, parameter_list) {
  
  # read in parameters
  for (i in 1:length(parameter_list)) {
    assign(names(parameter_list[i]), parameter_list[[i]])
  }
  
  # if the analysis is done for the reduced statistical model, ensure that the same grid cells are
  # used as for the full statistical model
  
  if (type_of_analysis == "reduced_model") {
    # create temporary parameter_list (because all the variables have to be consistent)
    plist_temp = create_stasticial_analysis_parameter_list(
      type_of_analysis = "full_model",
      statistical_method = statistical_method, 
      stepwise_direction = stepwise_direction,
      need_all_vars = need_all_vars,
      grouping_by = grouping_by,
      res = res, 
      standardise_variables = standardise_variables,
      seed = seed,
      oos_calculation = oos_calculation,
      n_groups_out_of_bag_sampling = n_groups_out_of_bag_sampling,
      save_files = save_files,
      crop = crop, 
      use_location_as_predictor = use_location_as_predictor,
      location_predictor = location_predictor,
      crop_calendar = crop_calendar,
      irrigation = irrigation, 
      mask_growing_season_areas = mask_growing_season_areas, 
      detrending_climate = detrending_climate,
      detrending_yield = detrending_yield,
      year_min = year_min, year_max = year_max, 
      verbose = verbose,
      yield_dataset = yield_dataset, 
      source_area_harvested = source_area_harvested,
      crop_to_region = crop_to_region,
      use_only_regions_with_frost_days = use_only_regions_with_frost_days,
      manual_subsetting = manual_subsetting,
      sampling_fraction = sampling_fraction,
      include_interaction_terms = include_interaction_terms
    )
    
    # test if results have been calculated before
    fn = sprintf("%s/job_completed.txt", plist_temp$output_dir)
    assert_that(file.exists(fn))
    
    # read in data_used
    fn = sprintf("%s/data_used/data_used.nc", plist_temp$output_dir)
    data_used_temp = netcdf2dataframe(fn, variables = c(predictand, plist_temp$predictors),
                                      time_format = "%Y")
    
    # reduce to only those years/grid cells were all values are present
    data_used_temp = na.omit(data_used_temp)
    data_used_temp = select(data_used_temp, lon, lat, time)
    
    # merge with data use to reduce the grid cells and time
    data_frame = merge(data_used_temp, data_frame, by = c("lon", "lat", "time"), all = FALSE)
  }
  
  return(data_frame)
}


###############################################################################
# Helper function: test if dimensions in a data frame are equal
# function: dimensions_are_equal ----------------------------------------
# 
# This function tests if two data frames have same dimensions (lon, lat, time have the same values in the same order)
#
#   Reads in:
#     - data_frame_1, data_frame_2: two data_frames with same dimension names
#     - dimension names, default: lon, lat, time
#     - verbose
#   Returns:
#     - logical: TRUE if dimensions are the same, FALSE if dimensions are different

dimensions_are_equal = function(data_frame_1, data_frame_2, 
                                dimensions = c("lon", "lat", "time"),
                                verbose = FALSE) {
  
  if (verbose) print("function: dimensions_are_equal")
  
  # some tests
  assert_that(is.data.frame(data_frame_1), is.data.frame(data_frame_2))
  assert_that(all(dimensions %in% names(data_frame_1)))
  assert_that(all(dimensions %in% names(data_frame_2)))
  
  for (dimension in dimensions) {
    if (nrow(data_frame_1[dimension]) != nrow(data_frame_2[dimension])) {
      if (verbose) print("Dimensions do not have same length.")
      return(FALSE)
    }
    
    if (any(data_frame_1[dimension] != data_frame_2[dimension])) {
      if (verbose) print("Dimensions do not have same values.")
      return(FALSE)
    }
  }
  
  if (verbose) print("Dimensions are equal.")
  return(TRUE)
}



merge_or_cbind = function(data_frame1, data_frame2,
                          by = intersect(names(data_frame1), names(data_frame2)),
                          all = FALSE,
                          all.x = all, all.y = all,
                          verbose = FALSE) {
  
  assert_that(is.data.frame(data_frame1), is.data.frame(data_frame2))
  assert_that(all(by %in% intersect(names(data_frame1), names(data_frame2))))
  
  # test if time, lon, lat are the same as in data_frame
  if (dimensions_are_equal(data_frame1, data_frame2, dimensions = by, verbose = verbose)) {
    data_frame = cbind(data_frame1, data_frame2[, !names(data_frame2) %in% by, drop = FALSE])
  } else {
    # merge instead of cbind
    data_frame = merge(data_frame1, data_frame2, by = by, all.x = all.x, all.y = all.y)
  }
}
