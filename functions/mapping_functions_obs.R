#### Averaging of mapped fields ####

m_cant_model_average <- function(df) {

  df <- df %>%
    select(lon, lat, depth, eras, basin, basin_AIP,
           starts_with("cant"),
           gamma)
  
  df_average <- df %>%
    fgroup_by(lon, lat, depth, eras, basin, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df_average)

}

m_cant_model_average_data_source <- function(df) {

  df <- df %>%
    select(lon, lat, depth, eras, basin, basin_AIP, data_source,
           starts_with("cant"),
           gamma)
  
  df_average <- df %>%
    fgroup_by(lon, lat, depth, eras, basin, basin_AIP, data_source) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df_average)

}

m_target_model_average <- function(df) {

  df <- df %>%
    select(lon, lat, depth, era, eras, basin, basin_AIP, gamma, 
           params_local$MLR_target)

  df <- df %>%
    fgroup_by(lon, lat, depth, era, eras, basin, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df)

}

m_target_model_average_data_source <- function(df) {

  df <- df %>%
    select(lon, lat, depth, era, eras, basin, basin_AIP, gamma, data_source,
           params_local$MLR_target)

  df <- df %>%
    fgroup_by(lon, lat, depth, era, eras, basin, basin_AIP, data_source) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df)

}

m_cant_predictor_zonal_mean <- function(df) {

  df <- df %>%
    fselect(lat, depth, eras, basin, basin_AIP,
            cant_intercept:gamma) %>%
    fgroup_by(lat, depth, eras, basin, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE))
    }

  return(df)

}

m_cant_predictor_zonal_mean_data_source <- function(df) {

  df <- df %>%
    fselect(lat, depth, eras, basin, basin_AIP, data_source,
            cant_intercept:gamma) %>%
    fgroup_by(lat, depth, eras, basin, basin_AIP, data_source) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE))
    }

  return(df)

}

m_cant_zonal_mean <- function(df) {

  df <- df %>%
    fselect(lat, depth, eras, basin, basin_AIP,
            cant, cant_pos, gamma, cant_sd, cant_pos_sd, gamma_sd) %>%
    fgroup_by(lat, depth, eras, basin, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_mean"),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df)

}

m_cant_zonal_mean_data_source <- function(df) {

  df <- df %>%
    fselect(lat, depth, eras, basin, basin_AIP, data_source,
            cant, cant_pos, gamma, cant_sd, cant_pos_sd, gamma_sd) %>%
    fgroup_by(lat, depth, eras, basin, basin_AIP, data_source) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_mean"),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df)

}


m_target_zonal_mean <- function(df) {

  df <- df %>%
    select(lat, depth, era, eras, basin, basin_AIP,
          gamma, gamma_sd,
          params_local$MLR_target, paste(params_local$MLR_target, "sd", sep = "_"))

  df <- df %>%
    fgroup_by(lat, depth, era, eras, basin, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_mean"),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df)

}


# calculate dcant column inventory [mol m-2] from dcant concentration [umol kg-1]
# inventories are calculated for a range of predefined inventory depth
m_dcant_inv <- function(df) {

  for (i_inventory_depth in params_global$inventory_depths) {

  # filter integration depth
  df_sub <- df %>%
    filter(depth <= i_inventory_depth)

  depth_level_volume <- tibble(
    depth = unique(df_sub$depth)) %>%
    arrange(depth)

  # determine depth level volume of each depth layer
  depth_level_volume <- depth_level_volume %>%
    mutate(layer_thickness_above = replace_na((depth - lag(depth)) / 2, 0),
           layer_thickness_below = replace_na((lead(depth) - depth) / 2, 0),
           layer_thickness = layer_thickness_above + layer_thickness_below) %>%
    select(-c(layer_thickness_above,
              layer_thickness_below))

  df_sub <- full_join(df_sub, depth_level_volume)

  # calculate cant layer inventory
  df_sub <- df_sub %>%
    mutate(dcant_layer_inv = dcant * layer_thickness * 1.03,
           dcant_pos_layer_inv = dcant_pos * layer_thickness * 1.03) %>%
    select(-layer_thickness)

  # sum up layer inventories to column inventories
  df_sub_inv <- df_sub %>%
    group_by(lon, lat, basin_AIP) %>%
    summarise(
      dcant     = sum(dcant_layer_inv, na.rm = TRUE) / 1000,
      dcant_pos = sum(dcant_pos_layer_inv, na.rm = TRUE) / 1000
    ) %>%
    ungroup()

  df_sub_inv <- df_sub_inv %>%
    mutate(inv_depth = i_inventory_depth)

  if (exists("df_inv")) {
    df_inv <- bind_rows(df_inv, df_sub_inv)
  }

  if (!exists("df_inv")) {
    df_inv <- df_sub_inv
  }


  }

  return(df_inv)

}


# calculate tcant column inventory [mol m-2] from tcant concentration [umol kg-1]
# inventories are calculated for a range of predefined inventory depth
m_tcant_inv <- function(df) {

  for (i_inventory_depth in params_global$inventory_depths) {

  # filter integration depth
  df_sub <- df %>%
    filter(depth <= i_inventory_depth)

  depth_level_volume <- tibble(
    depth = unique(df_sub$depth)) %>%
    arrange(depth)

  # determine depth level volume of each depth layer
  depth_level_volume <- depth_level_volume %>%
    mutate(layer_thickness_above = replace_na((depth - lag(depth)) / 2, 0),
           layer_thickness_below = replace_na((lead(depth) - depth) / 2, 0),
           layer_thickness = layer_thickness_above + layer_thickness_below) %>%
    select(-c(layer_thickness_above,
              layer_thickness_below))

  df_sub <- full_join(df_sub, depth_level_volume)

  # calculate cant layer inventory
  df_sub <- df_sub %>%
    mutate(tcant_layer_inv = tcant * layer_thickness * 1.03,
           tcant_pos_layer_inv = tcant_pos * layer_thickness * 1.03) %>%
    select(-layer_thickness)

  # sum up layer inventories to column inventories
  df_sub_inv <- df_sub %>%
    group_by(lon, lat, basin_AIP) %>%
    summarise(
      tcant     = sum(tcant_layer_inv, na.rm = TRUE) / 1000,
      tcant_pos = sum(tcant_pos_layer_inv, na.rm = TRUE) / 1000
    ) %>%
    ungroup()

  df_sub_inv <- df_sub_inv %>%
    mutate(inv_depth = i_inventory_depth)

  if (exists("df_inv")) {
    df_inv <- bind_rows(df_inv, df_sub_inv)
  }

  if (!exists("df_inv")) {
    df_inv <- df_sub_inv
  }


  }

  return(df_inv)

}

# calculate cant inventory within each density slab
m_cant_slab_inv_data_source <- function(df) {
  
    depth_level_volume <- df %>% 
      distinct(data_source, depth) %>% 
      arrange(depth)
    
    # determine depth level volume of each depth layer
    depth_level_volume <- depth_level_volume %>%
      group_by(data_source) %>% 
      mutate(layer_thickness_above = replace_na((depth - lag(depth)) / 2, 0),
             layer_thickness_below = replace_na((lead(depth) - depth) / 2, 0),
             layer_thickness = layer_thickness_above + layer_thickness_below) %>%
      ungroup() %>% 
      select(-c(layer_thickness_above,
                layer_thickness_below))
    
    df <- full_join(df, depth_level_volume)
    
    # calculate cant grid cell inventory
    df <- df %>%
      mutate(surface_area = earth_surf(lat, lon),
             cant_grid_inv = cant * layer_thickness * 1.03 * surface_area,
             cant_pos_grid_inv = cant_pos * layer_thickness * 1.03 * surface_area) %>%
      select(-c(layer_thickness, surface_area))
    

    df_slab_inv <- df %>%
      group_by(basin_AIP, data_source, gamma_slab) %>% 
      summarise(cant_total = sum(cant_grid_inv)*12*1e-18,
                cant_pos_total = sum(cant_pos_grid_inv)*12*1e-18) %>% 
      ungroup()

  return(df_slab_inv)
  
}

# calculate mean cant concentration within each grid cell of density slab
m_cant_slab <- function(df) {

  df_group <- df %>%
    group_by(lat, lon, gamma_slab, eras) %>%
    summarise(cant_pos = mean(cant_pos, na.rm = TRUE),
              cant = mean(cant, na.rm = TRUE),
              depth_max = max(depth, na.rm = TRUE)) %>%
    ungroup()

  return(df_group)

}

m_cant_slab_data_source <- function(df) {

  df_group <- df %>%
    group_by(lat, lon, gamma_slab, data_source) %>%
    summarise(cant_pos = mean(cant_pos, na.rm = TRUE),
              cant = mean(cant, na.rm = TRUE),
              depth_max = max(depth, na.rm = TRUE)) %>%
    ungroup()

  return(df_group)

}


# calculate zonal mean section
# old function that requires eras
m_zonal_mean_section <- function(df) {

  zonal_mean_section <- df %>%
    select(-lon) %>%
    fgroup_by(lat, depth, eras, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_mean"),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(zonal_mean_section)

}
# new function that uses all numeric variables and only required grouping variables
m_zonal_mean_sd <- function(df) {

  zonal_mean_section <- df %>%
    select(c(lat, depth, basin_AIP), where(is.numeric), -lon) %>%
    fgroup_by(lat, depth, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_mean"),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(zonal_mean_section)

}


#### Horizontal gridding ####

# cut lat and lon to a 1 x 1 deg horizontal grid
m_grid_horizontal <- function(df) {

  df <- df %>%
    mutate(
      lat = cut(lat, seq(-90, 90, 1), seq(-89.5, 89.5, 1)),
      lat = as.numeric(as.character(lat)),
      lon = cut(lon, seq(20, 380, 1), seq(20.5, 379.5, 1)),
      lon = as.numeric(as.character(lon))
    )

  return(df)

}

# cut lat and lon to a 5 x 5 deg horizontal grid
m_grid_horizontal_coarse <- function(df) {

  df <- df %>%
    mutate(
      lat_grid = cut(lat, seq(-90, 90, 5), seq(-87.5, 87.5, 5)),
      lat_grid = as.numeric(as.character(lat_grid)),
      lon_grid = cut(lon, seq(20, 380, 5), seq(22.5, 377.5, 5)),
      lon_grid = as.numeric(as.character(lon_grid))
    )

  return(df)

}



#### Neutral density slab assignment ####
# cut neutral density gamma into specific slabs for basins
m_cut_gamma <- function(df, var) {

  var <- sym(var)

  df_Atl <- df %>%
    filter(basin_AIP == "Atlantic") %>%
    mutate(gamma_slab = cut(!!var, params_local$slabs_Atl))

  df_Ind_Pac <- df %>%
    filter(basin_AIP %in% c("Indian", "Pacific")) %>%
    mutate(gamma_slab = cut(!!var, params_local$slabs_Ind_Pac))

  df <- bind_rows(df_Atl, df_Ind_Pac)

  return(df)

}
