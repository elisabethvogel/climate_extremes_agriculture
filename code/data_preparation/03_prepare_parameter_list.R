# This script calculates standardised time series of all climate variables (i.e. division by standard deviation).

project_path = "/Volumes/Seagate/PhD/Programming_Analysis/03_extremes_and_crop_yields/github"

parameter_list = create_stasticial_analysis_parameter_list(
  
  crop = crop,
  statistical_method = statistical_method,
  type_of_analysis = type_of_analysis,
  
  
  yield_dataset = "deepak", 
  crop_calendar = "agmip",
  irrigation = "combined",
  mask_growing_season_areas = TRUE,
  detrending_climate = "ssa",
  detrending_yield = "ssa",
  grouping_by = "global",
  
  use_time_as_predictor = FALSE,
  use_location_as_predictor = FALSE,
  use_only_regions_with_frost_days = FALSE,
  res = 1.5,
  year_min = 1961,
  year_max = 2008,
  oos_calculation = TRUE,
  seed = 100,
  n_groups_out_of_bag_sampling = 5,
  save_files = TRUE,
  source_area_harvested = "deepak",
  standardise_variables = TRUE,
  verbose = TRUE,
  manual_subsetting = TRUE,
  sampling_fraction = "max",
  
  include_interaction_terms = FALSE
  
)
