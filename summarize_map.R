library(tidyverse)
library(furrr)
library(vroom)

options(dplyr.summarise.inform = FALSE)

summarize_map <- function(file_path){
  
  print(file_path)
  
  vroom(file_path, progress = FALSE, col_types = cols()) |> 
    mutate(period = case_when(between(year, 1975, 2005) ~ "Historical",
                              between(year, 2040 - 15, 2040 + 15) ~ "2040s",
                              between(year, 2060 - 15, 2060 + 15) ~ "2060s",
                              between(year, 2080 - 15, 2080 + 15) ~ "2080s")) |>
    group_by(location, climate_proj, period, solar_red, cultivar) |>
    summarise(sunburn_hours = median(sunburn_hours),
              sunburn_days = median(sunburn_days))
}

files <- list.files("sunburn/data", pattern = "sunburn_data_*", recursive = TRUE, full.names = TRUE)


plan(multicore)
df <- future_map_dfr(files, summarize_map)

write_csv(df, "sunburn/output/sunburn_map.csv")


