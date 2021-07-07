#### Calculate derived variables ####

# calculate cstar
b_cstar <- function(tco2, phosphate, talk){

  cstar = tco2  - (params_local$rCP * phosphate)  - 0.5 * (talk - (params_local$rNP * phosphate))
  return(cstar)

  }

# calculate phosphate star
b_phosphate_star <- function(phosphate, oxygen){

  phosphate_star = phosphate + (oxygen / params_local$rPO) - params_local$rPO_offset
  return(phosphate_star)

  }


# calculate silicate star
b_silicate_star <- function(silicate, nitrate){
  
  silicate_star = silicate - nitrate
  return(silicate_star)
  
}


# calculate apparent oxygen utilization
# supplied oxygen data must be in mol kg-1

b_aou <- function(sal, tem, depth, oxygen) {

  oxygen_sat_m3 <- gas_satconc(
    S = sal,
    t = tem,
    P = 1.013253,
    species = "O2"
  )

  rho <-
    gsw_pot_rho_t_exact(SA = sal,
                        t = tem,
                        p = depth,
                        p_ref = 10.1325)

  oxygen_sat_kg = oxygen_sat_m3 * (1000 / rho)

  aou = oxygen_sat_kg - oxygen

  return(aou)

}


#### Map variables from MLR coefficients ####

# map cant from MLR coefficients and predictor variables

b_cant <- function(df) {
  
  df <- df %>%
    mutate(cant_intercept = `delta_coeff_(Intercept)`)
  
  vars = params_local$MLR_predictors
  
  for (i_var in vars) {
    df <- df %>%
      mutate(!!sym(paste("cant_", i_var, sep = "")) :=
               !!sym(i_var) *
               !!sym(paste("delta_coeff_", i_var, sep = "")))
  }
  
  df <- df %>%
    select(-contains("delta_coeff_"))
  
  df <- df %>%
    mutate(cant = reduce(select(., starts_with("cant_")), `+`))
  
  df <- df %>%
    mutate(cant_pos = if_else(cant < 0, 0, cant))

  return(df)

}


# map target variable from MLR coefficients and predictor variables

b_target_model <- function(df) {
  
  df <- df %>%
    mutate(!!sym(paste(params_local$MLR_target, "intercept", sep = "_")) :=
                   `coeff_(Intercept)`)
  
  vars = params_local$MLR_predictors
  
  for (i_var in vars) {
    df <- df %>%
      mutate(!!sym(paste(params_local$MLR_target, i_var, sep = "_")) :=
               !!sym(i_var) *
               !!sym(paste("coeff_", i_var, sep = "")))
  }
  
  df <- df %>%
    select(-contains("coeff_"))
  
  df <- df %>%
    mutate(!!sym(params_local$MLR_target) :=
             reduce(select(., starts_with(paste(params_local$MLR_target, "_", sep = ""))), `+`))

  return(df)

}
