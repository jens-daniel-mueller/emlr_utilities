#### Averaging of mapped dcant fields ####

m_dcant_3d_average <- function(df) {

  df <- df %>%
    select(data_source, lon, lat, depth, basin_AIP,
           starts_with("dcant"),
           gamma)
  
  df_average <- df %>%
    fgroup_by(data_source, lon, lat, depth, basin_AIP) %>% {
      add_vars(fgroup_vars(.,"unique"),
               fmean(., keep.group_vars = FALSE),
               fsd(., keep.group_vars = FALSE) %>% add_stub(pre = FALSE, "_sd"))
    }

  return(df_average)

}



# calculate layer thickness for each grid cell
# to be used in inventory and slab inventory calculations
m_layer_thickness <- function(df) {
  
  depth_level_volume <- df %>%
    distinct(depth) %>%
    arrange(depth) %>% 
    drop_na()
  
  # determine depth level volume of each depth layer
  depth_level_volume <- depth_level_volume %>%
    mutate(
      layer_thickness_above = (depth - lag(depth, default = 0)) / 2,
      layer_thickness_below = (lead(depth, default = max(depth)) - depth) / 2,
      layer_thickness = layer_thickness_above + layer_thickness_below
    ) %>%
    ungroup() %>%
    select(-c(layer_thickness_above,
              layer_thickness_below))
  
  df <- full_join(df, depth_level_volume)
  
  return(df)
  
}


# calculate dcant column inventory [mol m-2] from dcant concentration [umol kg-1]
# inventories are calculated for a range of predefined inventory depth
m_dcant_inv <- function(df) {
  
  for (i_inventory_depth in params_global$inventory_depths) {

  # filter integration depth
  df_sub <- df %>%
    filter(depth <= i_inventory_depth)

  df_sub <- m_layer_thickness(df_sub)

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


# calculate dcant column inventory [mol m-2] from dcant concentration [umol kg-1]
# inventories are calculated within predefined depth layers
m_dcant_inv_layer <- function(df) {
  
  layer <- params_global$inventory_depth_layer
  df <- m_layer_thickness(df)
  
  for (i in 1:(length(layer) - 1)) {
    # i <- 19
    # filter integration depth
    df_sub <- df %>%
      filter(depth >= layer[i],
             depth < layer[i + 1])
    
    if (nrow(df_sub) > 0) {
      # calculate cant layer inventory
      df_sub <- df_sub %>%
        mutate(
          dcant_layer_inv = dcant * layer_thickness * 1.03,
          dcant_pos_layer_inv = dcant_pos * layer_thickness * 1.03
        ) %>%
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
        mutate(inv_depth = (layer[i] + layer[i + 1]) / 2)
        # mutate(inv_depth = paste(layer[i], "-", layer[i + 1]))
      
      if (exists("df_inv")) {
        df_inv <- bind_rows(df_inv, df_sub_inv)
      }
      
      if (!exists("df_inv")) {
        df_inv <- df_sub_inv
      }
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

  df_sub <- m_layer_thickness(df_sub)

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

# calculate dcant budgets
m_dcant_budget <- function(df) {
  
  molC_to_PgC <- 12*1e-15
  
  df <- df %>% 
    mutate(surface_area = earth_surf(lat, lon),
           dcant_grid = dcant*surface_area*molC_to_PgC,
           dcant_pos_grid = dcant_pos*surface_area*molC_to_PgC)
  
    # calculate cant grid cell inventory
  df_budget <- df %>%
    group_by(data_source, inv_depth, method) %>% 
    summarise(dcant = sum(dcant_grid, na.rm = TRUE),
              dcant = round(dcant,3),
              dcant_pos = sum(dcant_pos_grid, na.rm = TRUE),
              dcant_pos = round(dcant_pos,3)) %>% 
    ungroup() %>% 
    pivot_longer(cols = dcant:dcant_pos,
                 names_to = "estimate",
                 values_to = "value")
  
  return(df_budget)
  
}

# calculate cant inventory within each density slab
# requires to run m_layer_thickness()
m_dcant_slab_budget <- function(df) {
  
  df <- m_layer_thickness(df)
  
  # calculate cant grid cell inventory
  df <- df %>%
    mutate(
      surface_area = earth_surf(lat, lon),
      volume_grid = layer_thickness * surface_area,
      dcant_grid = dcant * 1.03 * volume_grid,
      dcant_pos_grid = dcant_pos * 1.03 * volume_grid
    ) %>%
    select(-c(surface_area))
  
  
  df_slab_inv <- df %>%
    group_by(basin_AIP, gamma_slab) %>%
    summarise(
      dcant = sum(dcant_grid) * 12 * 1e-18,
      dcant_pos = sum(dcant_pos_grid) * 12 * 1e-18,
      slab_volume = sum(volume_grid)
    ) %>%
    ungroup()
  
  return(df_slab_inv)
  
}

# calculate mean dcant concentration within each grid cell of density slab
# requires to run m_layer_thickness()
m_dcant_slab_concentration <- function(df) {

  df <- m_layer_thickness(df)
  
  # calculate cant grid cell inventory
  df <- df %>%
    mutate(surface_area = earth_surf(lat, lon),
           layer_volume = layer_thickness * surface_area) %>%
    select(-c(surface_area))
  
  
  df_group <- df %>%
    group_by(lat, lon, gamma_slab) %>%
    summarise(dcant_pos = mean(dcant_pos, na.rm = TRUE),
              dcant = mean(dcant, na.rm = TRUE),
              dcant_pos_sd = mean(dcant_pos_sd, na.rm = TRUE),
              dcant_sd = mean(dcant_sd, na.rm = TRUE),
              layer_thickness = sum(layer_thickness),
              layer_volume = sum(layer_volume),
              n_layer = n()) %>%
    ungroup()

  return(df_group)

}



# zonal mean section of gridded data
# calculated for all numeric variables
# uses only mandatory gridding variables
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


# zonal mean section of bottle data with preceding lat/depth gridding of observations
# calculated for all numeric variables
# uses only mandatory gridding variables
m_zonal_mean_sd_bottle <- function(df) {
  
  df <- df %>%
    mutate(
      depth = cut(depth,
                       seq(0, 1e4, 200),
                       seq(100, 1e4, 200)),
      depth = as.numeric(as.character(depth)),
      lat = cut(lat,
                     seq(-100, 100, 2),
                     seq(-99, 100, 2)),
      lat = as.numeric(as.character(lat))
    )

  zonal_mean_section <- m_zonal_mean_sd(df)

  return(zonal_mean_section)

}



#### Temporal gridding ####

# cut years into decades

m_grid_decade <- function(year) {

  decade <- cut(
      year,
      c(1989, 2000, 2010, 2021),
      right = FALSE,
      labels = c("1989-1999", "2000-2009", "2010-2020")
    )
  
  return(decade)

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
  
  if (params_local$MLR_basins != "1") {
    df_Atl <- df %>%
      filter(basin_AIP == "Atlantic") %>%
      mutate(gamma_slab = cut(!!var, params_local$slabs_Atl))
    
    df_Ind_Pac <- df %>%
      filter(basin_AIP %in% c("Indian", "Pacific")) %>%
      mutate(gamma_slab = cut(!!var, params_local$slabs_Ind_Pac))
    
    df <- bind_rows(df_Atl, df_Ind_Pac)
  } else {
    df <- df %>%
      mutate(gamma_slab = cut(!!var, params_local$slabs_Ind_Pac))
  }
  
  return(df)
  
}
