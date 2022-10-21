library(tidyverse)
library(lubridate)
library(furrr)
library(agclimtools)
library(agclimdata)
library(tictoc)

location_index = as.numeric(commandArgs(trailingOnly=TRUE))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

source("sunburn/scripts/growth_model.R")

locations <- cdl_grid |> 
  filter(cdl_code == 68) |> 
  left_join(grid_stats) |>   
  left_join(elev_grid) |> 
  select(location, n, lat, lon, tz_long, elevation = min_elev) |> 
  filter(!is.na(elevation), !is.na(lat)) |> 
  filter(!location %in% read_lines("sunburn/output/completed_locations.txt"))

sel_location <- locations$location[location_index]

fst_single_location <- function(file_path, lat, solar_fraction, elevation){
  tic("process file")
  print(file_path)
  df <- read_binary(file_path, num_vars = 8) |> 
    # filter only simulated years
#    mutate(year = year(date)) |> 
#    filter(year %in% c(1975:2005, 
#                      seq(2040 - 15, 2040 + 15),
#                      seq(2060 - 15, 2060 + 15),
#                       seq(2080 - 15, 2080 + 15))) |> 
    mutate(doy = yday(date)) |> 
    add_hourly_temps(lat = lat, hour = 0:23) |> 
    add_growth_model() |> 
    # remove all temperatures less than 0 as boundary 
    # conductance can't be calculated below 0
    filter(temp > 0, diameter > 0) |>
    filter(between(month(date), 3, 9),
           between(hour, 4, 22)) |> 
    #  join with fractional solar data
    left_join(solar_fraction, by = c("doy", "hour"))|> 
    # calculate daily dew point temp
    mutate(tdew = dew_point(Rmax, tmin), 
           fruit_windspeed = canopy_wind_speed(windspeed),
           # calculate hourly solar radiation 
           # factor of 24 daily to hourly
           solar = solar_frac*SRAD*24,
           year = year(date)) |> 
    expand_grid(solar_red = c(1, 0.8, 0.7)) |>
    mutate(windspeed_red = if_else(solar_red == 1, 1, 0.6)) |> 
    # calculate fruit surface temp
    mutate(
      fst =
        fruit_surface_temp(temp,
                           windspeed*windspeed_red,
                           solar*solar_red,
                           tdew,
                           fruit_diameter = diameter,
                           fruit_reflectance = 0.325, 
                           fruit_sunlit_prop = 0.5)) 
    df_out <- df |>
      mutate(fst = if_else(cultivar == "Honeycrisp" & month(date) >= 9, 0, fst),
             sunburn = if_else(cultivar == "Honeycrisp",
                               fst > 43,
                               fst > 46)) |>
      group_by(year, cultivar, solar_red, date) |> 
      summarise(sunburn_day = max(sunburn),
                sunburn_hours = sum(sunburn)) |> 
      group_by(year, cultivar, solar_red) |>
      summarise(sunburn_hours = sum(sunburn_hours, na.rm = TRUE),
                sunburn_days = sum(sunburn_day))
    toc()
    
    return(df_out)
}

calc_sunburn <- function(location, lat, lon, tz_long, elevation){
  
  # Calculate fractional solar radiation
  solar_fraction <- expand_grid(hour = 0:23, doy = 1:366) |> 
    mutate(solar = clear_sky(hour, doy, lat, lon, time_zone_long = tz_long, elevation)) |> 
    group_by(doy) |> 
    mutate(solar_frac = solar/sum(solar)) |> 
    select(doy, hour, solar_frac)
  
  expand_grid(model = maca_models,
              climate_proj = climate_projections,
              location = location) |> 
    mutate(data = map(file.path(kamiak_maca_path(location),
                                model, climate_proj, location),
                      fst_single_location,
                      lat = lat,
                      solar_fraction = solar_fraction,
                      elevation = elevation)) |>
    unnest(data) 
    # mutate(period = case_when(between(year, 1975, 2005) ~ "Historical",
    #                           between(year, 2040 - 15, 2040 + 15) ~ "2040s",
    #                           between(year, 2060 - 15, 2060 + 15) ~ "2060s",
    #                           between(year, 2080 - 15, 2080 + 15) ~ "2080s")) |>
    # group_by(climate_proj, period, solar_red, cultivar) |> 
    # summarise(sunburn_hours = median(sunburn_hours),
    #           sunburn_days = median(sunburn_days))
}

sunburn_hours <- locations |> 
  filter(location == sel_location) |> 
  mutate(data = pmap(list(location, lat, lon, tz_long, elevation), calc_sunburn)) |> 
  select(-location, -tz_long, -elevation, -n) |> 
  unnest(data)

write_csv(sunburn_hours, paste0("sunburn/data/sunburn_", sel_location, ".csv"))
