visualize_spline_model <- function(modobj, 
                                   smooth_var,
                                   plot_og_scale = FALSE,
                                   mean_lookup = NULL,
                                   sd_lookup = NULL,
                                   pnorm_smooth_var = FALSE,
                                   plabels = NULL, 
                                   point_size = 2,
                                   line_size = 2,
                                   font_size = 16, 
                                   color_line = "black",
                                   color_point = "black") {
  
  if (any(class(modobj) == "gam")) {
    model <- modobj
  } else if (class(modobj$gam) == "gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  
  np <- 10000
  df <- model$model
  
  theseVars <- attr(model$terms, "term.labels")
  varClasses <- attr(model$terms, "dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  thisPred <- data.frame(init = rep(0, np))
  
  for (v in seq_along(theseVars)) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    
    if (thisVar == smooth_var) {
      thisPred[, smooth_var] <- seq(
        min(df[, smooth_var], na.rm = TRUE),
        max(df[, smooth_var], na.rm = TRUE),
        length.out = np
      )
    } else {
      switch(
        thisClass,
        "numeric" = {
          thisPred[, thisVar] <- median(df[, thisVar], na.rm = TRUE)
        },
        "factor" = {
          thisPred[, thisVar] <- levels(df[, thisVar])[[1]]
        },
        "ordered" = {
          thisPred[, thisVar] <- levels(df[, thisVar])[[1]]
        }
      )
    }
  }
  
  pred <- thisPred %>% dplyr::select(-init)
  p <- data.frame(predict(model, pred, se.fit = TRUE))
  pred <- cbind(pred, p)
  pred$selo <- pred$fit - 2 * pred$se.fit
  pred$sehi <- pred$fit + 2 * pred$se.fit
  if(plot_og_scale){
    if(pnorm_smooth_var){
      pred[[paste0(smooth_var, "_plot")]] <- pnorm(pred[[smooth_var]])
      df[[paste0(smooth_var, "_plot")]] <- pnorm(df[[smooth_var]])
    }
    else {
      pred[[paste0(smooth_var, "_plot")]] <- map2_dbl(pred[[smooth_var]], smooth_var, ~ (.x*sd_lookup[[.y]]) + mean_lookup[[.y]] )
      df[[paste0(smooth_var, "_plot")]] <- map2_dbl(df[[smooth_var]], smooth_var, ~ (.x*sd_lookup[[.y]]) + mean_lookup[[.y]] )
    }
    pred[["fit_plot"]] <- pnorm(pred$fit)
    pred[["se.fit_plot"]] <- pnorm(pred$se.fit)
    pred[["selo_plot"]] <- pnorm(pred$selo)
    pred[["sehi_plot"]] <- pnorm(pred$sehi)
    
    df[[paste0(thisResp, "_plot")]] <- pnorm(df[[thisResp]])
  } else{
    pred[[paste0(smooth_var, "_plot")]] <- pred[[smooth_var]]
    pred[["fit_plot"]] <- pred$fit
    pred[["se.fit_plot"]] <- pred$se.fit
    pred[["selo_plot"]] <- pred$selo
    pred[["sehi_plot"]] <- pred$sehi
    
    df[[paste0(smooth_var, "_plot")]] <- df[[smooth_var]]
    df[[paste0(thisResp, "_plot")]] <- df[[thisResp]]
  }
  
  ggplot(df, aes(x = .data[[paste0(smooth_var, "_plot")]], y = .data[[paste0(thisResp, "_plot")]])) +
    geom_point(alpha = 0.3, stroke = 0, size = point_size, color = color_point) +
    geom_ribbon(
      data = pred,
      aes(
        x = .data[[paste0(smooth_var,"_plot")]],
        ymin = selo_plot,
        ymax = sehi_plot
      ),
      inherit.aes = FALSE,
      alpha = 0.35,
      linetype = 0,
      color = color_line
    ) +
    geom_line(
      data = pred,
      aes(
        x = .data[[paste0(smooth_var,"_plot")]],
        y = fit_plot
      ),
      inherit.aes = FALSE,
      linewidth = line_size,
      color = color_line
    ) +
    labs(title = plabels,
         ylab = thisResp,
         xlab = smooth_var) +
    theme(
      text = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
}


# visualize_spline_model <- function(modobj, 
#                                    smooth_var, 
#                                    plabels = NULL, 
#                                    point_size = 2,
#                                    line_size = 2,
#                                    font_size = 16, 
#                                    color_line = "black",
#                                    color_point = "black") {
#   
#   if (any(class(modobj) == "gam")) {
#     model <- modobj
#   } else if (class(modobj$gam) == "gam") {
#     model <- modobj$gam
#   } else {
#     stop("Can't find a gam object to plot")
#   }
#   
#   np <- 10000
#   df <- model$model
#   
#   theseVars <- attr(model$terms, "term.labels")
#   varClasses <- attr(model$terms, "dataClasses")
#   thisResp <- as.character(model$terms[[2]])
#   
#   thisPred <- data.frame(init = rep(0, np))
#   
#   for (v in seq_along(theseVars)) {
#     thisVar <- theseVars[[v]]
#     thisClass <- varClasses[thisVar]
#     
#     if (thisVar == smooth_var) {
#       thisPred[, smooth_var] <- seq(
#         min(df[, smooth_var], na.rm = TRUE),
#         max(df[, smooth_var], na.rm = TRUE),
#         length.out = np
#       )
#     } else {
#       switch(
#         thisClass,
#         "numeric" = {
#           thisPred[, thisVar] <- median(df[, thisVar], na.rm = TRUE)
#         },
#         "factor" = {
#           thisPred[, thisVar] <- levels(df[, thisVar])[[1]]
#         },
#         "ordered" = {
#           thisPred[, thisVar] <- levels(df[, thisVar])[[1]]
#         }
#       )
#     }
#   }
#   
#   pred <- thisPred %>% dplyr::select(-init)
#   
#   p <- data.frame(predict(model, pred, se.fit = TRUE))
#   pred <- cbind(pred, p)
#   pred$selo <- pred$fit - 2 * pred$se.fit
#   pred$sehi <- pred$fit + 2 * pred$se.fit
#   
#   ggplot(df, aes(x = .data[[smooth_var]], y = .data[[thisResp]])) +
#     geom_point(alpha = 0.3, stroke = 0, size = point_size, color = color_point) +
#     geom_ribbon(
#       data = pred,
#       aes(
#         x = .data[[smooth_var]],
#         ymin = selo,
#         ymax = sehi
#       ),
#       inherit.aes = FALSE,
#       alpha = 0.35,
#       linetype = 0,
#       color = color_line
#     ) +
#     geom_line(
#       data = pred,
#       aes(
#         x = .data[[smooth_var]],
#         y = fit
#       ),
#       inherit.aes = FALSE,
#       linewidth = line_size,
#       color = color_line
#     ) +
#     labs(title = plabels) +
#     theme(
#       text = element_text(size = font_size),
#       axis.text = element_text(size = font_size),
#       legend.text = element_text(size = font_size),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.background = element_blank(),
#       plot.background = element_blank()
#     )
# }