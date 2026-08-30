plot_grouped_bar <- function(table_list, value_col,
                                group_names = NULL,
                                color_vec = NULL,
                                p_col = NULL,
                                p_thresh = 0.05,
                                var_names = NULL,
                                name_mapping = NULL,
                                use_alpha = TRUE){
  
  df <- dplyr::bind_rows(table_list, .id = "group")
  
  # optional variable filtering
  if(!is.null(var_names)){
    df <- df[df$phenotype %in% var_names, ]
  }
  
  # optional name mapping
  if(!is.null(name_mapping)){
    df <- merge(df, name_mapping, by.x = "phenotype", by.y = "key", all.x = TRUE)
    df$plot_label <- ifelse(is.na(df$plot_names), df$phenotype, df$plot_names)
  } else {
    df$plot_label <- df$phenotype
  }
  
  if(!is.null(group_names)){
    df$group <- factor(df$group, labels = group_names)
  }
  
  df$signif <- TRUE
  if(!is.null(p_col)){
    df$signif <- df[[p_col]] < p_thresh
  }
  
  ggplot(df, aes(x = reorder(plot_label, .data[[value_col]]),
                 y = .data[[value_col]])) +
    
    geom_col(
      aes(
        fill = group,
        alpha = signif
      ),
      position = "dodge"
    ) +
    
    coord_flip() +
    
    scale_fill_manual(values = color_vec) +
    
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5)) +
    
    theme_minimal() +
    labs(x = NULL, y = value_col)
}

plot_lollipop <- function(table_list,
                          value_col,
                          label_col = "phenotype",
                          group_names = NULL,
                          color_vec = NULL,
                          color_col = NULL,
                          color = "steelblue",
                          p_col = NULL,
                          p_thresh = 0.05,
                          var_names = NULL,
                          name_mapping = NULL,
                          group_col = NULL,
                          position_dodge_width = 0.6,
                          line_width = 1,
                          point_size = 5) {

  # Handle single dataframe -- group_col is optional here; omitting it means
  # "one implicit group" (no dodge), covering what used to be a separate,
  # ungrouped plot_lollipop() with a single fixed color.
  if (is.data.frame(table_list)) {

    df <- table_list
    df$group <- if (!is.null(group_col)) df[[group_col]] else "1"

  } else {

    # Handle list of tables
    df <- dplyr::bind_rows(table_list, .id = "group")
  }

  # Optional filtering
  if (!is.null(var_names)) {
    df <- df[df[[label_col]] %in% var_names, ]
  }

  # Optional name mapping
  if (!is.null(name_mapping)) {

    df <- merge(
      df,
      name_mapping,
      by.x = label_col,
      by.y = "key",
      all.x = TRUE
    )

    df$plot_label <- ifelse(
      is.na(df$plot_names),
      df[[label_col]],
      df$plot_names
    )

  } else {

    df$plot_label <- df[[label_col]]
  }

  # Optional group names
  if (!is.null(group_names)) {
    df$group <- factor(
      df$group,
      levels = unique(df$group),
      labels = group_names
    )
  }

  # Significance
  df$signif <- TRUE

  if (!is.null(p_col)) {
    df$signif <- df[[p_col]] < p_thresh
  }

  # Color can come from a different column than the dodge/legend group --
  # defaults to `group` (color = which group a point belongs to). With
  # neither color_col supplied nor more than one group, there's nothing
  # meaningful to color-map, so fall back to a single fixed color and skip
  # the legend entirely -- this is the old single-series plot_lollipop() case.
  df$color_val <- if (!is.null(color_col)) df[[color_col]] else df$group
  mapped_color <- !is.null(color_col) || length(unique(df$group)) > 1

  p <- ggplot(
    df,
    aes(
      y = reorder(plot_label, .data[[value_col]]),
      x = .data[[value_col]],
      group = group
    )
  )

  if (mapped_color) {

    p <- p +
      geom_linerange(
        aes(
          xmin = 0,
          xmax = .data[[value_col]],
          color = color_val,
          alpha = signif
        ),
        position = position_dodge(width = position_dodge_width),
        linewidth = line_width
      ) +
      geom_point(
        aes(
          color = color_val,
          alpha = signif
        ),
        position = position_dodge(width = position_dodge_width),
        size = point_size
      ) +
      scale_color_manual(values = color_vec)

  } else {

    p <- p +
      geom_linerange(
        aes(
          xmin = 0,
          xmax = .data[[value_col]],
          alpha = signif
        ),
        linewidth = line_width,
        color = color
      ) +
      geom_point(
        aes(alpha = signif),
        size = point_size,
        color = color
      )
  }

  p +
    scale_alpha_manual(
      values = c(`TRUE` = 1, `FALSE` = 0.3)
    ) +
    theme_minimal() +
    labs(
      x = value_col,
      y = NULL
    )
}

plot_dumbbell <- function(table_list, value_col,
                             group_names = NULL,
                             color_vec = NULL,
                             p_col = NULL,
                             p_filter = FALSE,
                             p_thresh = 0.05,
                             var_names = NULL,
                             name_mapping = NULL,
                             connect = TRUE){
  
  df <- dplyr::bind_rows(table_list, .id = "group")
  
  if(!is.null(var_names)){
    df <- df[df$phenotype %in% var_names, ]
  }
  
  if(!is.null(name_mapping)){
    df <- merge(df, name_mapping, by.x = "phenotype", by.y = "key", all.x = TRUE)
    df$plot_label <- ifelse(is.na(df$plot_names), df$phenotype, df$plot_names)
  } else {
    df$plot_label <- df$phenotype
  }
  
  if(!is.null(group_names)){
    df$group <- factor(df$group, labels = group_names)
  }
  
  if(p_filter && !is.null(p_col)){
    df <- df[df[[p_col]] < p_thresh, ]
  }
  
  ggplot(df, aes(x = .data[[value_col]],
                 y = reorder(plot_label, .data[[value_col]]))) +
    geom_point(aes(color = group), size = 2) +
    {if(connect) geom_line(aes(group = phenotype), color = "grey80")} +
    scale_color_manual(values = color_vec) +
    theme_minimal() +
    labs(x = value_col, y = NULL)
}

plot_heatmap <- function(table_list, value_col,
                            group_names = NULL,
                            low_col = "blue",
                            mid_col = "white",
                            high_col = "red",
                            p_col = NULL,
                            p_filter = FALSE,
                            p_thresh = 0.05,
                            var_names = NULL,
                            name_mapping = NULL){
  
  df <- dplyr::bind_rows(table_list, .id = "group")
  
  if(!is.null(var_names)){
    df <- df[df$phenotype %in% var_names, ]
  }
  
  if(!is.null(name_mapping)){
    df <- merge(df, name_mapping, by.x = "phenotype", by.y = "key", all.x = TRUE)
    df$plot_label <- ifelse(is.na(df$plot_names), df$phenotype, df$plot_names)
  } else {
    df$plot_label <- df$phenotype
  }
  
  if(!is.null(group_names)){
    df$group <- factor(df$group, labels = group_names)
  }
  
  if(p_filter && !is.null(p_col)){
    df <- df[df[[p_col]] < p_thresh, ]
  }
  
  ggplot(df, aes(x = group,
                 y = plot_label,
                 fill = .data[[value_col]])) +
    geom_tile() +
    scale_fill_gradient2(low = low_col,
                         mid = mid_col,
                         high = high_col,
                         midpoint = median(df[[value_col]], na.rm = TRUE)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = value_col)
}

plot_brain_map <- function(table,
                           value_col,
                           atlas = c("dk", "aseg"),
                           p_col = NULL,
                           p_filter = FALSE,
                           p_thresh = 0.05,
                           fill_name = "value",
                           title = "",
                           limits = NULL,
                           low_color  = "#F2F0EC",
                           mid_color  = "#D99A9A",
                           high_color = "#D95C72",
                           na_color = "grey90",
                           legend_width_cm = 6){
  
  atlas <- match.arg(atlas)
  
  df <- table
  
  if(p_filter && !is.null(p_col)){
    df <- df[df[[p_col]] < p_thresh, ]
  }
  
  fill_scale <- if(is.null(limits)){
    
    scale_fill_gradientn(
      colours = c(low_color, mid_color, high_color),
      na.value = na_color,
      guide = guide_colorbar(barwidth = unit(legend_width_cm, "cm"))
    )
    
  } else {
    
    scale_fill_gradientn(
      colours = c(low_color, mid_color, high_color),
      values = scales::rescale(c(limits[1],
                                 mean(limits),
                                 limits[2])),
      limits = limits,
      oob = scales::squish,
      na.value = na_color,
      guide = guide_colorbar(barwidth = unit(legend_width_cm, "cm"))
    )
  }
  
  if(atlas == "dk"){
    
    plot_data <- merge(df, ggseg::dk, by = "label")
    
    ggplot(plot_data) +
      ggseg::geom_brain(
        atlas = ggseg::dk,
        position = ggseg::position_brain(hemi ~ side),
        aes(fill = .data[[value_col]])
      ) +
      fill_scale +
      theme_void() +
      labs(title = title, fill = fill_name)
    
  } else {
    
    plot_data <- merge(df, ggseg::aseg, by = "label")
    
    ggplot(plot_data) +
      ggseg::geom_brain(
        atlas = ggseg::aseg,
        aes(fill = .data[[value_col]])
      ) +
      fill_scale +
      theme_void() +
      labs(title = title, fill = fill_name)
  }
}

# 
# plot_brain_map <- function(table,
#                            value_col,
#                            atlas = c("dk", "aseg"),
#                            p_col = NULL,
#                            p_filter = FALSE,
#                            p_thresh = 0.05,
#                            fill_name = "value",
#                            title = "",
#                            limits = NULL){
#   
#   atlas <- match.arg(atlas)
#   
#   df <- table
#   
#   if(p_filter && !is.null(p_col)){
#     df <- df[df[[p_col]] < p_thresh, ]
#   }
#   
#   fill_scale <- if(is.null(limits)){
#     scale_fill_gradient(
#       low = "grey90",
#       high = "red",
#       na.value = "grey90"
#     )
#   } else {
#     scale_fill_gradient(
#       low = "grey90",
#       high = "red",
#       na.value = "grey90",
#       limits = limits
#     )
#   }
#   
#   if(atlas == "dk"){
#     
#     plot_data <- merge(df, ggseg::dk, by = "label")
#     
#     ggplot(plot_data) +
#       ggseg::geom_brain(
#         atlas = ggseg::dk,
#         position = ggseg::position_brain(hemi ~ side),
#         aes(fill = .data[[value_col]])
#       ) +
#       fill_scale +
#       theme_void() +
#       labs(title = title, fill = fill_name)
#     
#   } else {
#     
#     plot_data <- merge(df, ggseg::aseg, by = "label")
#     
#     ggplot(plot_data) +
#       ggseg::geom_brain(
#         atlas = ggseg::aseg,
#         aes(fill = .data[[value_col]])
#       ) +
#       fill_scale +
#       theme_void() +
#       labs(title = title, fill = fill_name)
#   }
# }


plot_region_of_change <- function(table,
                                  var_names = NULL,
                                  name_mapping = NULL,
                                  alpha_range = c(0.2, 1),
                                  positive_color = "#C95A5A",
                                  negative_color = "#5A7BC9",
                                  segment_width = 6,
                                  r2_limits = NULL,
                                  alpha_order = FALSE){
  
  df <- table
  
  df$phenotype <- rownames(df)
  
  # optional variable filtering
  if(!is.null(var_names)){
    df <- df[df$phenotype %in% var_names, ]
  }
  
  # optional phenotype labels
  if(!is.null(name_mapping)){
    
    df <- merge(
      df,
      name_mapping[, c("key", "plot_names")],
      by.x = "phenotype",
      by.y = "key",
      all.x = TRUE
    )
    
    df$plot_label <- ifelse(
      is.na(df$plot_names),
      df$phenotype,
      df$plot_names
    )
    
  } else {
    
    df$plot_label <- df$phenotype
  }
  
  # remove rows with no detectable region
  if(sum(is.na(df$region_start) & is.na(df$region_end)) > 0){
    message(glue("These phenotypes have no region start/region end: {paste(df[(is.na(df$region_start) & is.na(df$region_end)),]$plot_label, collapse = \",\")}"))
  }
  df <- df[!is.na(df$region_start) &
             !is.na(df$region_end), ]
  
  # derivative direction
  df$direction <- ifelse(
    df$avg_first_deriv >= 0,
    "Positive",
    "Negative"
  )
  
  # alphabetical ordering, keeping L/R pairs together
  if (alpha_order == TRUE) {
    base_label <- gsub("^[LR][-_ ]*", "", df$plot_label)
    
    ord <- order(
      base_label,
      ifelse(grepl("^L[-_ ]*", df$plot_label), 0,
             ifelse(grepl("^R[-_ ]*", df$plot_label), 1, 0))
    )
    
    df$plot_label <- factor(
      df$plot_label,
      levels = rev(unique(df$plot_label[ord]))
    )
  }
  
  if(is.null(r2_limits)){
    
    alpha_scale <- scale_alpha_continuous(
      range = alpha_range,
      name = expression("Partial " * R^2)
    )
    
  } else {
    
    alpha_scale <- scale_alpha_continuous(
      range = alpha_range,
      limits = r2_limits,
      oob = scales::squish,
      name = expression("Partial " * R^2)
    )
    
  }
  
  ggplot(
    df,
    aes(
      y = plot_label,
      x = region_start,
      xend = region_end,
      yend = plot_label,
      color = direction,
      alpha = r2_full_vs_none
    )
  ) +
    geom_segment(
      linewidth = segment_width,
      lineend = "butt"
    ) +
    scale_color_manual(
      values = c(
        Positive = positive_color,
        Negative = negative_color
      )
    ) +
    alpha_scale +
    labs(
      x = "Region of Change",
      y = NULL,
      color = "Direction"
    ) +
    labs(
      x = "Region of Change",
      y = NULL,
      color = "Direction"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}


