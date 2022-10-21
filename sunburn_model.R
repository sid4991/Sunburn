calc_windspeed_fruit <- function(u, canopy_height=3, fruit_height=2, a = 0.4){
  u_star <- (0.4*u)/log((10 - .65*canopy_height)/(.1*canopy_height))
  windspeed_canopy <- (u_star/0.4)*log(3.5)
  windspeed_fruit <- windspeed_canopy*exp(a*((fruit_height/canopy_height)-1))
  return(windspeed_fruit)
}

calc_atmospheric_emissivity <- function(e_a, T_a){
  # Air temp cannot be negative
  epsilon_a <- 1.24*(((e_a*10)/(T_a+273))^(1/7))
  return(epsilon_a)
}

calc_vapour_pressure <- function(temp){
  e_a = 0.6108*exp((17.27*temp)/(temp +237.3))
  return(e_a)
}

calc_air_density_over_pressure <- function(P, T_a){
  rho_air <- 3.486*(1/(1.01*(T_a + 273)))
  return(rho_air)
}

calc_delta <- function(T_a){
  delta <- (2503*exp((17.27*T_a)/(T_a + 237.3)))/((T_a + 237.3)^2)
  return(delta)
}

calc_boundary_conductance <- function(u, fruit_diameter){
  g_a <- 1.4*.135*sqrt(u/(0.84*fruit_diameter))
  return(g_a)
}

calc_dew_point <- function(relative_humidity, temp){
  # https://en.wikipedia.org/wiki/Dew_point
  a = 6.1121 #mbar
  b = 18.678
  c = 257.14 #°C
  d = 234.5 #°C.
  
  gamma <- log(relative_humidity/100) + 
    (b*temp)/(c+temp)
  
  T_dew <- (c*gamma)/(b-gamma)
  return(T_dew)
}

calc_longwave <- function(T_a, T_ground, atmospheric_emissivity, 
                          ground_emissivity, fruit_groundlit_prop){
  stefan_boltzmann <- 5.67E-8 #Wm^-2T^-4

  longwave <- atmospheric_emissivity*stefan_boltzmann*(T_a + 273)^4 +
    fruit_groundlit_prop*ground_emissivity*stefan_boltzmann*(T_ground + 273)^4
  return(longwave)
}

calc_shortwave <- function(s_rad, fruit_sunlit_prop, fruit_reflectance){
  shortwave <- fruit_sunlit_prop*(1-fruit_reflectance)*s_rad
  return(shortwave)
}

calc_fruit_surface_temp <- 
  function(T_a, 
           u, 
           s_rad,
           T_dew,
           T_ground = T_a,
           ground_emissivity = 0.97,
           fruit_diameter = 0.08, 
           fruit_emissivity = 0.95,
           fruit_reflectance = 0.6,
           fruit_sunlit_prop = 1,
           fruit_groundlit_prop = 0,
           fruit_surface_conductance = 5E-5,
           print_all_vars = FALSE, 
           possibly = TRUE){
    
  # Constants #####
  stefan_boltzmann <- 5.67E-8 #Wm^-2T^-4
  specific_heat_air <- 29.3 #Jmol^-1C^-1
  latent_heat_of_vaporisation <- 2.429E6 #J/Kg
  
  g_a <- calc_boundary_conductance(u, fruit_diameter)
  delta <- calc_delta(T_a)
  air_density_over_pressure <- calc_air_density_over_pressure(pressure, T_a)
  vapour_pressure_air <- calc_vapour_pressure(T_a)
  vapour_pressure_dew <- calc_vapour_pressure(T_dew)
  atmospheric_emissivity <- 
    calc_atmospheric_emissivity(vapour_pressure_dew, T_a)
  
  # Radiation 
  shortwave <- calc_shortwave(s_rad, fruit_sunlit_prop, fruit_reflectance)
  longwave <- calc_longwave(T_a, T_ground, atmospheric_emissivity,
                            fruit_groundlit_prop, ground_emissivity)
  
  incoming_radiation <- longwave + shortwave
  
  denominator <- (specific_heat_air*g_a) +
    (latent_heat_of_vaporisation*
        fruit_surface_conductance*
        air_density_over_pressure*delta)
  
  # a_0
  a_0_numerator <- incoming_radiation +
    specific_heat_air*g_a*T_a +
    (latent_heat_of_vaporisation*fruit_surface_conductance*
        air_density_over_pressure*delta*T_a) -
    (latent_heat_of_vaporisation*fruit_surface_conductance*
        air_density_over_pressure*(vapour_pressure_air - vapour_pressure_dew))
  
  a_0 <- -((a_0_numerator/denominator)+273)
  
  # a_4
  a_4 <- (fruit_emissivity*stefan_boltzmann)/denominator
  
  if (print_all_vars) {
    cat(paste(
      "\n---------------------------------",
      "\nAir Temp:", round(T_a),
      "\nWind Speed:", round(u, 2), 
      "\nSolar:", round(s_rad),
      "\nDew Point Temp:", T_dew,
      "\nBoundary Conductance:", round(g_a, 2),
      "\nDelta:", round(delta, 2),
      "\nAir Density Over Pressure:", round(air_density_over_pressure, 2),
      "\nVapour Pressure (Air):", round(vapour_pressure_air, 2),
      "\nVapour Pressure (Dew):", round(vapour_pressure_dew, 2),
      "\nAtmospheric Emissivity:", round(atmospheric_emissivity, 2),
      "\nShortwave:", round(shortwave),
      "\nLongwave:", round(longwave),
      "\na_0:", a_0,
      "\na_4:", a_4))
  }
  
  solve <- function(a_0, a_4){
    
    sol <- polyroot(c(a_0,1,0,0,a_4))

    # filter solution set to only positive real solution
    sol_real <- Re(sol[abs(Im(sol)) < 1E-5])
    sol_final <- sol_real[sol_real > 0] - 273
    return(sol_final)
  }
  
  # solve quartic
  
  if (possibly) {
    sol <- purrr::map2_dbl(a_0, a_4, possibly(solve, otherwise = NA))
  } else {
    sol <- purrr::map2_dbl(a_0, a_4, solve)
  }
  
  return(sol)

}

