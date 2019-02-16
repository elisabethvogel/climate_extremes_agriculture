# This script contains some meta-data of all data, so that it is easier to reference to data in a dynamic way.

# This script includes:
#
# Data set names: data_sets[[data_type]]
# * data_type = yield / extreme_indicators / climate_data / soil_moisture / drought_indicators
#
# Variable information - how are variables named for each data_type: variables[[data_type]][[data_source]]
# * data_type = yield / extreme_indicators / climate_data / soil_moisture / drought_indicators
# * data_source = depends on data_type, but for example "deepak", "cru_ts_323", etc...
#
# Data path information: data_paths[[time_frame]][[data_type]][[data_source]]
# * time_frame = growing_season / complete / static
# * data_type = yield / extreme_indicators / climate_data / soil_moisture / drought_indicators
# * data_source = depends on data_type, but for example "deepak", "cru_ts_323", etc...

################################################################################
# crops
crops = c("wheat", "winter_wheat", "spring_wheat", "rice", "maize", "soybeans")

################################################################################
# data types

data_types = list()

data_types$growing_season = c("yield", "extreme_indicators", 
                              "climate_data", "drought_indicators")

data_types$complete = c("extreme_indicators", "climate_data", "drought_indicators")
data_types$static = c("crop_calendars")

################################################################################
# data sets

data_sets = list()
data_sets$yield = c("deepak")
data_sets$extreme_indicators = c("hadex2")
data_sets$climate_data = c("cru_ts_323")
data_sets$drought_indicators = c("cru_ts_323_spi")
data_sets$crop_calendars = c("agmip")

################################################################################
# variables

variables = list()
variables$yield$deepak = c("yield")
variables$extreme_indicators$hadex2 = 
  c("TXx", "TX90p", "TNn", "TN10p", "Rx5day", "Rx1day", "DTR")
variables$climate_data$cru_ts_323 = c("tmp", "pre", "frs", "dtr")
variables$drought_indicators$cru_ts_323_spi = c("spi_3", "spi_6")
variables$crop_calendars$agmip = c("plant", "harvest")

################################################################################
# resolution

resolution = list()

resolution$yield$deepak = c(0.5)
# resolution$extreme_indicators$hadex2 = c("3.75")
resolution$extreme_indicators$hadex2 = c(0.5) # 1.5?
# resolution$extreme_indicators$ghcndex = c("2.5")
resolution$climate_data$cru_ts_323 = c(0.5)
resolution$drought_indicators$cru_ts_323_spi = c(0.5)

################################################################################
# data paths 
base_data_path = file.path(project_path, "data")

data_paths = list()
# data paths - complete climate data
data_paths$complete$extreme_indicators$hadex2 =
  sprintf("%s/climate_data/extreme_indicators/hadex2_gridded_data", base_data_path)
data_paths$complete$climate_data$cru_ts_323 = 
  sprintf("%s/climate_data/cru_ts_323", base_data_path)
data_paths$complete$drought_indicators$cru_ts_323_spi = 
  sprintf("%s/hydrological_data/drought_indicators/cru_ts_323_spi",
          base_data_path)

# data paths - growing season
data_paths$growing_season$yield$deepak = 
  sprintf("%s/yield_data_detrended/deepak_crop_data", base_data_path)
data_paths$growing_season$extreme_indicators$hadex2 =
  sprintf("%s/growing_season_data_detrended/hadex2", base_data_path)
data_paths$growing_season$climate_data$cru_ts_323 = 
  sprintf("%s/growing_season_data_detrended/climate_data/cru_ts_323", base_data_path)
data_paths$growing_season$drought_indicators$cru_ts_323_spi = 
  sprintf("%s/growing_season_data_detrended/drought_indicators/cru_ts_323_spi", base_data_path)

# data paths - crop calendars
data_paths$static$crop_calendars$agmip = sprintf("%s/crop_calendar/agmip_crop_calendar_month_regridded", base_data_path)
