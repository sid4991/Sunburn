add_growth_model <- function(df, sel_cultivars = c("Cripps Pink", "Honeycrisp")) {
  # tictoc::tic("growth model")
  df_out <- df |> 
    # calculate degree hours
    mutate(temp_f = (temp*9/5) + 32,
           DH = if_else(temp_f >= 42 & temp_f <= 77.64, temp_f - 42, 0)) |> 
    # sum degree hours per day
    group_by(date) |> 
    summarise(DD = sum(DH)/24) |> 
    mutate(year = year(date)) |>
    group_by(year) |>
    mutate(DD = cumsum(DD),
           honeycrisp_prop = -0.5854134 + 0.001121*DD + -1.92e-07*DD^2,
           cripps_pink_prop = -0.257021 + 0.0007148*DD + -1.01e-07*DD^2,
           honeycrisp_prop = pmin(pmax(honeycrisp_prop, 0), 1),
           cripps_pink_prop = pmin(pmax(cripps_pink_prop, 0), 1),
           `Honeycrisp` = honeycrisp_prop*.085,
           `Cripps Pink` = cripps_pink_prop*.075) |>
    ungroup() |> 
    select(date, `Honeycrisp`, `Cripps Pink`) |> 
    pivot_longer(c(`Honeycrisp`, `Cripps Pink`), names_to = "cultivar", 
                 values_to = "diameter")|> 
    filter(cultivar %in% sel_cultivars) |> 
    #  join back to original data frame
    full_join(df, by = "date")
  # tictoc::toc()
  return(df_out)
}

  