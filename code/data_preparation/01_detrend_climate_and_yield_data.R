###### Detrending of all growing season data

project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"

source(file.path(project_path, "code/data_preparation/01_functions_detrending.R"), local = TRUE)
source(file.path(project_path, "code/other/file_information.R"))

################################################################################
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

if (parallel_calculation == TRUE) {
  
  foreach (crop = crops) %:%
    foreach (data_type = data_types_temp) %:% # without yield
    foreach (data_set = data_sets[[data_type]]) %:%
    foreach (variable = variables[[data_type]][[data_set]]) %dopar% {
      
      source(sprintf("%s/code/01_functions_detrending.R", project_path), local = TRUE) # includes detrending function
      
      if (data_type == "yield" && (irrigation != "combined" || crop_calendar != "agmip" || masked != "masked") )
        # only calculate yield data once
        return(NULL)
      
      # don't detrend HADEX2 in GHCNDEX resolution and vice versa
      if ((data_set == "hadex2" && res == 2.5) || (data_set == "ghcndex" && res == 3.75) ||
          (data_set == "era_interim" && res != 1.5))
        return(NULL)
      
      ################################################################################
      # redirect output to a file instead of to display so one can still read what's happening (and for easier debugging)
      if (redirecting_output_to_a_file) {
        dir.create(sprintf("%s/code/cluster_jobs", project_path), showWarnings = FALSE, recursive = TRUE)
        
        fn_out = sprintf("%s/code/cluster_jobs/detrending_gs_data_%s_%s_%s_%s_%s_%s_%s_%s_output.txt", 
                         project_path, lag, crop, crop_calendar, irrigation, masked, data_type, data_set, variable) 
        file_out = file(fn_out, open = "wt")
        sink(fn_out, type = "output") # only captures print messages and other output, error messages seem to be "caught" first by foreach function --> couldn't find a way to catch them with sink
      }
      
      detrend_growing_season_data(crop, crop_calendar,
                                  irrigation, masked, 
                                  data_type, data_set, variable, 
                                  res, 
                                  detrending_types = detrending_types,
                                  skip_existing_files = skip_existing_files,
                                  verbose = verbose)
      
      # stop the sink connection
      if (redirecting_output_to_a_file)
        sink(type = "output")
    }
  
} else {
  
  for (crop in crops) {
    for (data_type in data_types_temp) { # without yield
      for (data_set in data_sets[[data_type]]) {
        for (variable in variables[[data_type]][[data_set]]) {
          
          if (data_type == "yield" && (irrigation != "combined" || crop_calendar != "agmip" || masked != "masked") )
            # only calculate yield data once
            next
          
          # don't detrend HADEX2 in GHCNDEX resolution and vice versa
          if ((data_set == "hadex2" && res == 2.5) || (data_set == "ghcndex" && res == 3.75) ||
              (data_set == "era_interim" && res != 1.5))
            next
          
          detrend_growing_season_data(crop, crop_calendar,
                                      irrigation, masked, 
                                      data_type, data_set, variable, 
                                      res, 
                                      detrending_types = detrending_types,
                                      skip_existing_files = skip_existing_files,
                                      verbose = verbose)
          
        }
      }
    }
  }
}

print("Completed!")
