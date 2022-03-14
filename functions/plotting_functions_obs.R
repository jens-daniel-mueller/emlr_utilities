#### Section plots  ####

# Global section along transect
# Color scale continuous (default) or divergent
p_section_global <-
  function(df,
           var,
           var_name = var,
           title_text = "Global section",
           subtitle_text = "N-Atl -> SO -> N-Pac",
           col = "continuous") {

    var <- sym(var)

    # subset data along section
    df_sec <- left_join(section_global_coordinates, df)

    # prepare base section plot
    section_base <- df_sec %>%
      ggplot(aes(dist, depth, z = !!var)) +
      scale_y_reverse() +
      labs(y = "Depth (m)") +
      guides(fill = guide_colorsteps(barheight = unit(8, "cm")))

    # add chose color scale (default continuous)
    if (col == "continuous") {

      section <- section_base +
        geom_contour_filled() +
        geom_vline(data = section_global_coordinates %>% filter(lat == 0.5),
                   aes(xintercept = dist),
                   col = "white") +
        scale_fill_viridis_d(name = var_name)

    } else {

      title_text <- "Global section | Color range 99th percentile"

      max <- df_sec %>%
        select(!!var) %>%
        pull %>%
        abs() %>%
        quantile(0.99, na.rm = TRUE)

      breaks = seq(-max, max, length.out = 20)

      section <- section_base +
        geom_contour_filled(breaks = breaks) +
        geom_vline(data = section_global_coordinates %>% filter(lat == 0.5),
                   aes(xintercept = dist),
                   col = "white") +
        scale_fill_scico_d(palette = "vik", drop = FALSE,
                           name = var_name)

    }

    # cut out surface water section
    surface <-
      section +
      coord_cartesian(expand = 0,
                      ylim = c(500, 0)) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      labs(y = "Depth (m)")

    # cut out deep water section
    deep <-
      section +
      coord_cartesian(expand = 0,
                      ylim = c(params_global$plotting_depth, 500)) +
      labs(x = "Distance (Mm)", y = "Depth (m)")

    # combine surface and deep water section
    surface / deep +
      plot_layout(guides = "collect") +
      plot_annotation(title = title_text,
                      subtitle = subtitle_text)

  }



# plot sections at regular lon intervals
p_section_climatology_regular <-
  function(df,
           var,
           var_name = var,
           surface = "n",
           col = "continuous",
           title_text = "Latitudinal sections",
           subtitle_text = "at predefined longitude intervals") {

  var <- sym(var)

  # plot base section
  section <- df %>%
    filter(lon %in% params_global$longitude_sections_regular) %>%
    ggplot(aes(lat, depth, z = !!var)) +
    guides(fill = guide_colorsteps(barheight = unit(7, "cm"))) +
    scale_x_continuous(breaks = seq(-100, 100, 40),
                       limits = c(-85,85)) +
    facet_wrap( ~ lon, ncol = 3, labeller = label_both) +
    theme(axis.title.x = element_blank()) +
    labs(subtitle = subtitle_text)

  # plot layer for chose color scale (default continuous)
  if (col == "continuous") {

    section <- section +
      geom_contour_filled() +
      scale_fill_viridis_d(name = var_name) +
      labs(title = title_text)

  } else {

    title_text <- "Latitudinal sections | Color range 99th percentile"

    max <- df %>%
      select(!!var) %>%
      pull %>%
      abs() %>%
      quantile(0.99)

    breaks = seq(-max, max, length.out = 20)

    section <- section +
      geom_contour_filled(breaks = breaks) +
      scale_fill_scico_d(palette = "vik",
                         drop = FALSE,
                         name = var_name) +
      labs(title = title_text)

  }
  
  
  # select surface or deep water y range
  if (surface == "n") {

    section <- section +
      scale_y_continuous(trans = trans_reverser("sqrt"),
                         breaks = c(100,500,seq(1000,5000,1000))) +
      coord_cartesian(expand = 0,
                      ylim = c(params_global$plotting_depth, 0))

  } else {
    
    section <- section +
      coord_cartesian(expand = 0,
                      ylim = c(params_local$depth_min, 0))

  }

  section

}


# Zonal mean section of dcant estimates
p_section_zonal <-
  function(df,
           var = "dcant",
           var_name = var,
           col = "continuous",
           gamma = "gamma_mean",
           plot_slabs = "y",
           drop_slabs = 1,
           breaks = c(-Inf, seq(0,16,2), Inf),
           legend_title = expression(atop(Delta * C[ant],
                                          (mu * mol ~ kg ^ {-1}))),
           title_text = "Zonal mean section",
           subtitle_text = "") {

    var <- sym(var)
    gamma <- sym(gamma)
    
    if (var == "dcant_mean"){
      legend_title <- expression(atop(Delta * C["ant"],
                                      (mol ~ m ^ {-2})))
    }
    
    if (var == "dcant_pos_mean"){
      legend_title <- expression(atop(Delta * C["ant,pos"],
                                      (mol ~ m ^ {-2})))
    }
    
    # plot base section
    section <- df %>%
      ggplot() +
      guides(fill = guide_colorsteps(barheight = unit(8, "cm"),
                                     show.limits = TRUE)) +
      scale_y_reverse() +
      scale_x_continuous(breaks = seq(-100, 100, 20),
                         limits = c(-85,85))

    # plot layer for chose color scale (default continuous)
    if (col == "continuous") {

      breaks_n <- length(breaks) - 1

      section <- section +
        geom_contour_filled(aes(lat, depth, z = !!var),
                            breaks = breaks) +
        scale_fill_viridis_d(drop = FALSE,
                             name = legend_title)
    } else {
      
      breaks <- c(-Inf, seq(-6,6,1), Inf)
      
      if (var == "dcant_bias"){
        legend_title <- expression(atop(Delta * Delta * C["ant"],
                                        (mol ~ m ^ {-2})))
      }
      
      if (var == "dcant_pos_bias"){
        legend_title <- expression(atop(Delta * Delta * C["ant,pos"],
                                        (mol ~ m ^ {-2})))
      }

      section <- section +
        geom_contour_filled(aes(lat, depth, z = !!var),
                            breaks = breaks) +
        scale_fill_scico_d(palette = "vik",
                           drop = FALSE,
                           name = legend_title)

    }


    # plot isoneutral density lines if chosen (default yes)
    if (plot_slabs == "y") {

      # select slab breaks for plotted basin
      if (unique(df$basin_AIP) == "Atlantic") {
        slab_breaks <- params_local$slabs_Atl
      } else {
        slab_breaks <- params_local$slabs_Ind_Pac
      }


      section <- section  +
        geom_hline(yintercept = params_local$depth_min,
                   col = "white",
                   linetype = 2) +
        geom_contour(aes(lat, depth, z = !!gamma),
                     breaks = slab_breaks,
                     col = "black") +
        geom_text_contour(
          aes(lat, depth, z = !!gamma),
          breaks = slab_breaks,
          col = "black",
          skip = drop_slabs
        )

    }

    # cut surface water section
    surface <-
      section +
      coord_cartesian(
        expand = 0,
        ylim = c(500, 0)
      ) +
      labs(y = "Depth (m)",
           title = title_text,
           subtitle = subtitle_text) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

    # cut deep water section
    deep <-
      section +
      coord_cartesian(
        expand = 0,
        ylim = c(params_global$plotting_depth, 500)
      ) +
      labs(x = expression(latitude~(degree*N)), y = "Depth (m)")


    # combine surface and deep water section
    surface / deep +
      plot_layout(guides = "collect")

  }


# Zonal mean section of dcant estimates
p_section_zonal_continous_depth <-
  function(df,
           var = "dcant",
           var_name = var,
           col = "continuous",
           gamma = "gamma_mean",
           plot_slabs = "y",
           drop_slabs = 1,
           breaks = c(-Inf, seq(0,16,2), Inf),
           legend_title = expression(atop(Delta * C[ant],
                                          (mu * mol ~ kg ^ {-1}))),
           title_text = "Zonal mean section",
           subtitle_text = "") {

    var <- sym(var)
    gamma <- sym(gamma)
    
    if (var == "dcant_mean"){
      legend_title <- expression(atop(Delta * C["ant"],
                                      (mol ~ m ^ {-2})))
    }
    
    if (var == "dcant_pos_mean"){
      legend_title <- expression(atop(Delta * C["ant,pos"],
                                      (mol ~ m ^ {-2})))
    }
    
    # plot base section
    section <- df %>%
      ggplot() +
      guides(fill = guide_colorsteps(barheight = unit(8, "cm"),
                                     show.limits = TRUE)) +
      scale_y_continuous(
        trans = trans_reverser("sqrt"),
        breaks = c(0,100, 500, seq(1500, 5000, 1000)),
        limits = c(params_global$inventory_depth_standard, 0),
        # breaks = seq(0,4900, 1000),
        name = "Depth (m)"
      )+
      scale_x_continuous(breaks = seq(-100, 100, 20),
                         limits = c(-80,65),
                         name = "Latitude (°N)") +
      coord_cartesian(expand = 0) +
      labs(title = title_text,
           subtitle = subtitle_text)

    # plot layer for chose color scale (default continuous)
    if (col == "continuous") {

      breaks_n <- length(breaks) - 1

      section <- section +
        geom_contour_filled(aes(lat, depth, z = !!var),
                            breaks = breaks) +
        colorspace::scale_fill_discrete_sequential(palette = "Rocket",
                                                   drop = FALSE,
                                                   name = legend_title)
        # scale_fill_viridis_d(drop = FALSE,
        #                      name = legend_title)
    } else {
      
      breaks <- c(-Inf, seq(-6,6,1), Inf)
      
      if (var == "dcant_bias"){
        legend_title <- expression(atop(Delta * Delta * C["ant"],
                                        (mol ~ m ^ {-2})))
      }
      
      if (var == "dcant_pos_bias"){
        legend_title <- expression(atop(Delta * Delta * C["ant,pos"],
                                        (mol ~ m ^ {-2})))
      }

      section <- section +
        geom_contour_filled(aes(lat, depth, z = !!var),
                            breaks = breaks) +
        colorspace::scale_fill_discrete_divergingx(palette = "RdBu",
                                                   drop = FALSE,
                                                   name = legend_title)
        # scale_fill_scico_d(palette = "cork",
        #                    drop = FALSE,
        #                    name = legend_title)

    }


    # plot isoneutral density lines if chosen (default yes)
    if (plot_slabs == "y") {

      # select slab breaks for plotted basin
      if (unique(df$basin_AIP) == "Atlantic") {
        slab_breaks <- params_local$slabs_Atl
      } else {
        slab_breaks <- params_local$slabs_Ind_Pac
      }


      section <- section  +
        geom_hline(yintercept = params_local$depth_min,
                   col = "white",
                   linetype = 2) +
        geom_contour(aes(lat, depth, z = !!gamma),
                     breaks = slab_breaks,
                     col = "black") +
        geom_text_contour(
          aes(lat, depth, z = !!gamma),
          breaks = slab_breaks,
          col = "black",
          skip = drop_slabs
        )

    }
    
    section
}


#### Map plots  ####


# plot column inventory map of cant estimate
# Color scale continuous (default) for pos cant or divergent for all cant
p_map_cant_inv <-
  function(df,
           var = "dcant",
           col = "continuous",
           breaks = c(-Inf, seq(0,16,2), Inf),
           title_text = "Column inventory map",
           subtitle_text = NULL) {
    
    var <- sym(var)
    legend_title <- var
    
    if (var == "dcant"){
      legend_title <- expression(atop(Delta * C["ant"],
                                      (mol ~ m ^ {-2})))
    }
    
    if (var == "dcant_scaled"){
      legend_title <- expression(atop(beta,
                                      (mol ~ m ^ {-2} ~ µatm ^ {-1})))
    }
    
    if (var == "dcant_pos"){
      legend_title <- expression(atop(Delta * C["ant,pos"],
                                      (mol ~ m ^ {-2})))
    }

    if (var == "dcant_mean"){
      legend_title <- expression(atop(Delta * C["ant,mean"],
                                      (mol ~ m ^ {-2})))
    }

    if (var == "dcant_mean_bias"){
      legend_title <- expression(atop(Delta * C["ant,mean"] ~ bias,
                                      (mol ~ m ^ {-2})))
    }

    if (var == "dcant_sd"){
      legend_title <- expression(atop(Delta * C["ant,sd"],
                                      (mol ~ m ^ {-2})))
    }
    
    if (var == "tcant"){
      legend_title <- expression(atop(C["ant"],
                                      (mol ~ m ^ {-2})))
      breaks = c(-Inf, seq(0,60,10), Inf)
    }
    
    if (var == "tcant_pos"){
      legend_title <- expression(atop(C["ant,pos"],
                                      (mol ~ m ^ {-2})))
      breaks = c(-Inf, seq(0,60,10), Inf)
    }
    
    if (col == "continuous") {
      
      breaks_n <- length(breaks) - 1
      
      df <- df %>%
        mutate(var_int = cut(!!var,
                             breaks,
                             right = FALSE))
      map +
        geom_tile(data = df,
                    aes(lon, lat, fill = var_int)) +
        colorspace::scale_fill_discrete_sequential(palette = "Rocket",
                                       drop = FALSE,
                                       name = legend_title)+
        # scale_fill_viridis_d(drop = FALSE,
        #                      name = legend_title) +
        guides(fill = guide_colorsteps(barheight = unit(6, "cm"))) +
        labs(title = title_text,
             subtitle = subtitle_text)
      
    } else if (col == "divergent") {
      breaks = params_global$breaks_cant_inv
      
      map +
        geom_tile(data = df,
                    aes(lon, lat, fill = cut(!!var, breaks))) +
        colorspace::scale_fill_discrete_divergingx(palette = "RdBu",
                           drop = FALSE,
                           name = expression(atop(Delta * C[ant],
                                                  (mu * mol ~ kg ^ {-1})))) +
        # scale_fill_scico_d(palette = "cork",
        #                    drop = FALSE,
        #                    name = expression(atop(Delta * C[ant],
        #                                           (mu * mol ~ kg ^ {-1})))) +
        guides(fill = guide_colorsteps(barheight = unit(6, "cm")))  +
        labs(title = title_text,
             subtitle = subtitle_text)
      
    } else if (col == "bias") {
      breaks = c(-Inf, seq(-6,6,1),Inf)
      
      if (var == "dcant_bias"){
        legend_title <- expression(atop(Delta * C["ant"] ~ bias,
                                        (mol ~ m ^ {-2})))
      }
      
      if (var == "dcant_pos_bias"){
        legend_title <- expression(atop(Delta * C["ant,pos"] ~ bias,
                                        (mol ~ m ^ {-2})))
      }
      
      map +
        geom_tile(data = df,
                    aes(lon, lat, fill = cut(!!var, breaks))) +
        scale_fill_scico_d(palette = "vik",
                           drop = FALSE,
                           name = legend_title) +
        guides(fill = guide_colorsteps(barheight = unit(6, "cm")))  +
        labs(title = title_text,
             subtitle = subtitle_text)
    }
  }




# plot map of mean cant within density slab
# Color scale continuous (default) for pos cant or divergent for all cant
p_map_dcant_slab <-
  function(df,
           var = "dcant_pos",
           col = "continuous",
           breaks = c(-Inf, seq(0,14,2),Inf),
           legend_title = NULL,
           title_text = "Density slab average",
           subtitle_text = NULL) {
    var <- sym(var)
    
    slab_map <-
      map +
      geom_tile(data = df,
                  aes(lon, lat, fill = cut(!!var,
                                           breaks,
                                           right = FALSE))) +
      guides(fill = guide_colorsteps(barheight = unit(6, "cm"))) +
      labs(title = title_text,
           subtitle = subtitle_text)
    
    
    # plot map for chose color scale (default continuous)
    if (col == "continuous") {
      
      if (is.null(legend_title)) {
        legend_title = expression(atop(Delta * C[ant],
                                       (mu * mol ~ kg ^ {-1})))
      }

      if (var == sym("dcant_sd")) {
        legend_title = expression(atop(Delta * C[ant] ~ SD,
                                       (mu * mol ~ kg ^ {-1})))
      }
      
      slab_map +
        scale_fill_viridis_d(name = legend_title,
                             drop = FALSE)
      
    } else {
      breaks <- c(-Inf, seq(-8,8,2), Inf)
      
      if (is.null(legend_title)){
      legend_title <- expression(atop(Delta * C[ant],
                                      (mu * mol ~ kg ^ {-1})))
      }
      
      slab_map +
        scale_fill_scico_d(palette = "vik",
                           drop = FALSE,
                           name = legend_title)
      
    }
  }




# Maps at predefined depth layers
# Color scale continuous (default) or divergent
p_map_climatology <-
  function(df,
           var,
           title_text = "Distribution maps",
           subtitle_text = "at predefined depth levels",
           col = "continuous") {

    var <- sym(var)

    # filter depth levels
    df <- df %>%
      filter(depth %in% params_global$depth_levels)

    # prepare map
    map_base <-
      map +
      geom_tile(data = df,
                  aes(lon, lat, fill = !!var)) +
      geom_tile(data = section_global_coordinates,
                  aes(lon, lat), fill = "white") +
      facet_wrap( ~ depth, labeller = label_both) +
      labs(title = title_text,
           subtitle = subtitle_text)

    # add chose color scale (default continuous)
    if (col == "continuous") {
    map_base +
      scale_fill_viridis_c()

    } else {

      max <- df %>%
        select(!!var) %>%
        pull %>%
        abs() %>%
        max()

      limits <- c(-1, 1) * max

      map_base +
        scale_fill_scico(
          palette = "vik",
          limit = limits
        )

    }


  }


#### Other plots  ####


# plot column inventory map of cant estimate
# Color scale continuous (default) for pos cant or divergent for all cant
p_prop_prop <-
  function(df,
           var1,
           var2,
           limit_percentile = 1,
           bins = 50) {
    var1 <- sym(var1)
    var2 <- sym(var2)
    
    
    df <- df %>% 
      filter(!!var1 <= quantile(!!var1, limit_percentile, na.rm = TRUE),
             !!var1 >= quantile(!!var1, 1-limit_percentile, na.rm = TRUE),
             !!var2 <= quantile(!!var2, limit_percentile, na.rm = TRUE),
             !!var2 >= quantile(!!var2, 1-limit_percentile, na.rm = TRUE))
    
    
    # calculate equal axis limits and binwidth
    
    axis_lims <- df %>%
      drop_na() %>%
      summarise(max_value = max(c(max(!!var1),
                                  max(!!var2))),
                min_value = min(c(min(!!var1),
                                  min(!!var2))))
    
    binwidth_value <-
      (axis_lims$max_value - axis_lims$min_value) / bins
    axis_lims <- c(axis_lims$min_value, axis_lims$max_value)
    
    ggplot(df, aes(x = !!var1,
                   y = !!var2)) +
      geom_bin2d(binwidth = binwidth_value) +
      scale_fill_viridis_c(trans = "log10") +
      geom_abline(slope = 1, col = 'red') +
      coord_equal(xlim = axis_lims,
                  ylim = axis_lims) +
      facet_wrap(~ basin_AIP)
    
  }

