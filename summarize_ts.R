library(tidyverse)
library(furrr)
library(vroom)
library(agclimdata)
library(agclimtools)

index <- as.numeric(commandArgs(trailingOnly=TRUE))
sel_model <- maca_models[index]

options(dplyr.summarise.inform = FALSE)

summarize_map <- function(file_path){
  
  print(file_path)
  
  vroom(file_path, progress = FALSE, col_types = cols()) |> 
    filter(model == sel_model)
}

plan(multicore)

files <- list.files("sunburn/data", pattern = "sunburn_data*") |> 
  str_remove_all("sunburn_|.csv")


df_loc <- grid_stats |> 
  mutate(region = case_when(state == "WA" ~ "Washington",
                            state %in% c("ME", "NH", "VT",
                                         "NY", "OH", "PA",
                                         "MA", "CT", "RI", 
                                         "NJ", "DE", "MD", 
                                         "MI", "WV") ~ "North East")) |> 
  filter(!is.na(region)) |> 
  inner_join(cdl_grid |> filter(cdl_code == 68)) |> 
  filter(location %in% files)

df <- future_map_dfr(paste0("sunburn/data/sunburn_", df_loc$location, ".csv"), summarize_map) |> 
  left_join(df_loc) |> 
  group_by(climate_proj, model, solar_red, cultivar, region, year) |> 
  summarize(mean_sunburn_hours = weighted.mean(sunburn_hours, n),
            mean_sunburn_days = weighted.mean(sunburn_days, n))

write_csv(df, paste0("sunburn/output/sunburn_ts", sel_model, ".csv"))


