# This script calculates standardised time series of all climate variables (i.e. division by standard deviation).

project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"
source(sprintf("%s/code/data_preparation/02_function_calculate_standardised_timeseries.R", project_path), local = TRUE)

source(sprintf("%s/code/other/file_information.R", project_path))

parallel_calculation = FALSE
redirecting_output_to_a_file = FALSE # for parallel
n_cluster = 6
skip_existing_files = TRUE
verbose = TRUE

################################################################################
# set parameters
################################################################################

res = 1.5 # detrending at resolution of 1.5
crop_calendar = "agmip"
irrigation = "combined" # both irrigated and rainfed
masked = "masked" # only use grid cells where Sacks et al. and SAGE calendar have growing season dates

# crops = c("wheat", "maize", "rice", "soybeans")
crops = c("maize")
data_types_temp = c("yield", "extreme_indicators", "climate_data", "drought_indicators")
detrending_types = c("original", "ssa")

if (parallel_calculation) {
  require(foreach)
  require(doParallel)
  # register cores for computation
  cl = makeCluster(n_cluster)
  registerDoParallel(cl)
}

################################################################################
if (parallel_calculation == TRUE) {
  foreach (data_type_temp = data_types_temp) %:%
    foreach (data_set = data_sets[[data_type_temp]]) %:%
    foreach (variable = variables[[data_type_temp]][[data_set]]) %:%
    foreach (crop = crops) %:%
    
    source(sprintf("%s/code/data_preparation/02_function_calculate_standardised_timeseries.R", project_path), local = TRUE)
  
  if (data_type == "yield" && (irrigation != "combined" || crop_calendar != "agmip" || masked != "masked"))
    # only calculate yield data once
    return(NULL)
  
  # don't standardise HADEX2 in GHCNDEX resolution and vice versa
  if ((data_set == "hadex2" && res == 2.5) || (data_set == "ghcndex" && res == 3.75) ||
      (data_set == "era_interim" && res != 1.5))
    return(NULL)
  
  ################################################################################
  # redirect output to a file instead of to display so one can still read what's happening (and for easier debugging)
  if (redirecting_output_to_a_file) {
    dir.create(sprintf("%s/code/cluster_jobs", project_path), showWarnings = FALSE, recursive = TRUE)
    fn_out = sprintf("%s/code/cluster_jobs/calculating_standardised_ts_%s_%s_%s_%s_%s_%s_%s_%s_%s_output.txt", 
                     project_path, data_type_temp, data_set, variable, lag, crop, crop_calendar, irrigation, masked, res) 
    file_out = file(fn_out, open = "wt")
    sink(fn_out, type = "output") # only captures print messages and other output, error messages seem to be "caught" first by foreach function --> couldn't find a way to catch them with sink
  }
  
  ################################################################################
  calculate_standardised_timeseries(crop, crop_calendar, irrigation, masked, 
                                    data_type_temp, data_set, variable, res, 
                                    skip_existing_files, verbose)
  
  ################################################################################
  # stop the sink connection
  if (redirecting_output_to_a_file)
    sink(type = "output")
  
} else {
  for(data_type_temp in data_types_temp) {
    for(data_set in data_sets[[data_type_temp]]) {
      for(variable in variables[[data_type_temp]][[data_set]]) {
        for(crop in crops) {
          
          print(sprintf("%s, %s, %s, %s, %s, %s, %s, %s", 
                        crop, crop_calendar, irrigation,
                        masked, data_type_temp, data_set, variable, res))
          
          if (data_type == "yield" && (irrigation != "combined" || crop_calendar != "agmip" || masked != "masked") )
            # only calculate yield data once
            return(NULL)
          
          # don't standardise HADEX2 in GHCNDEX resolution and vice versa
          if ((data_set == "hadex2" && res == 2.5) || (data_set == "ghcndex" && res == 3.75) ||
              (data_set == "era_interim" && res != 1.5))
            return(NULL)
          
          ################################################################################
          calculate_standardised_timeseries(crop, crop_calendar, irrigation, masked, 
                                            data_type_temp, data_set, variable, res, 
                                            skip_existing_files, verbose)
          
        }
      }
    }
  }
}

print("Completed!")
