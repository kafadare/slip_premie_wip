## The modifications in this wp.taki script has the option to return the % of points outside of the 95CI, instead of a plot.
## specified via the plot = TRUE/FALSE flag.
wp.taki_EK <- function (object = NULL, xvar = NULL, resid = NULL, n.inter = 4, df_tmp = NULL,
                        xlim.worm = 3.5, ylim.worm = 12 * sqrt(n.inter/length(resid)), plot = TRUE,
                        ...) 
{
  get.wp.x <- function(y.in) {
    qq <- as.data.frame(qqnorm(y.in, plot = FALSE))
    return(qq$x)
  }
  get.wp.y <- function(y.in) {
    qq <- as.data.frame(qqnorm(y.in, plot = FALSE))
    return(qq$y - qq$x)
  }
  get.lims <- function(resid) {
    lims.df <- data.frame(zval = seq(-xlim.worm, xlim.worm, 
                                     0.1))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * 
                                                    (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se
    lims.df$high <- -lims.df$low
    return(lims.df)
  }
  deparen <- function(expr) {
    while (is.language(expr) && !is.name(expr) && deparse(expr[[1L]])[1L] == 
           "(") expr <- expr[[2L]]
    expr
  }
  if (is.null(object) && is.null(resid)) 
    stop(paste("A fitted object with resid() method or the argument resid should be used ", 
               "\n", ""))
  resid <- if (is.null(object)) {
    resid
  }
  else {
    resid(object)
  }
  DataExist <- FALSE
  if (!is.null(object) && any(grepl("data", names(object$call))) && 
      !(object$call["data"] == "sys.parent()()")) {
    DaTa <- eval(object$call[["data"]])
    DataExist <- TRUE
  }
  if (is.null(xvar)) {
    lims.df <- data.frame(zval = seq(-xlim.worm, xlim.worm, 
                                     0.1))
    lims.df$p <- pnorm(lims.df$zval)
    lims.df$se <- (1/dnorm(lims.df$zval)) * (sqrt(lims.df$p * 
                                                    (1 - lims.df$p)/length(resid)))
    lims.df$low <- qnorm((1 - 0.95)/2) * lims.df$se
    lims.df$high <- -lims.df$low
    wp.df <- data.frame(y = resid %>% get.wp.y, x = resid %>% 
                          get.wp.x)
    wp.dt <- wp.df %>% arrange(x) %>% as.data.table()
    lims.dt <- as.data.table(lims.df)
    combo.dt <- lims.dt[wp.dt, on = .(zval == x), roll = TRUE]
    if (sum(is.na(combo.dt)) > 0) {
      warning("missing some CI values, try increasing xlim.worm")
    }
    n_outer <- combo.dt %>% mutate(outer = ifelse((y < low | 
                                                     y > high), 1, 0)) %>% summarise(n = n(), n_outer = sum(outer)) %>% 
      mutate(pcnt = n_outer/n, x = xlim.worm * 0.95, y = ylim.worm * 
               0.95)
    p <- ggplot(wp.df, aes(x = x, y = y)) + geom_smooth(method = lm, 
                                                        formula = y ~ poly(x, 3)) + geom_point() + theme_classic() + 
      xlab("Unit Normal Quantile") + ylab("Deviation") + 
      {
        if (is.finite(xlim.worm)) 
          xlim(c(-xlim.worm, xlim.worm))
      } + {
        if (is.finite(ylim.worm)) {
          coord_cartesian(ylim = c(-ylim.worm, ylim.worm))
        }
        else {
          coord_cartesian(ylim = c(-1, 1))
        }
      } + geom_line(data = lims.df, aes(x = zval, y = low), 
                    linetype = "dashed") + geom_line(data = lims.df, 
                                                     aes(x = zval, y = high), linetype = "dashed") + geom_hline(yintercept = 0, 
                                                                                                                linetype = "dashed") + geom_text(data = n_outer, 
                                                                                                                                                 mapping = aes(x = x, y = y, label = scales::percent(pcnt)), 
                                                                                                                                                 inherit.aes = FALSE, color = "blue")
    if (plot == TRUE) {
      return(p)
    }  else{
      return(n_outer)
    }    
    
  }
  else {
    if (!is(xvar, "formula")) {
      if (length(resid) != length(xvar)) 
        stop("Error - incorrect length of predictor vector...")
      if (is.factor(xvar)) {
        wp.df <- data.frame(z = xvar)
      }
      else {
        wp.df <- data.frame(z = cut(xvar, breaks = quantile(xvar, 
                                                            probs = seq(0, 1, length.out = n.inter + 1)), 
                                    include.lowest = TRUE))
      }
    }
    if (is(xvar, "formula")) {
      if (DataExist) {
        xvar.vec <- eval(deparen(deparen(xvar)[[2L]]), 
                         envir = as.environment(DaTa))
        if (is.factor(xvar.vec)) {
          wp.df <- data.frame(z = xvar.vec)
        }
        else {
          wp.df <- data.frame(z = cut(xvar.vec, breaks = quantile(xvar.vec, 
                                                                  probs = seq(0, 1, length.out = n.inter + 
                                                                                1)), include.lowest = TRUE))
        }
      }
      else {
        stop("Dataframe missing, exiting...")
      }
    }
    wp.df$resid <- resid
    wp.df <- wp.df %>% group_by(z) %>% mutate(y = resid %>% 
                                                get.wp.y, x = resid %>% get.wp.x)
    lims.df <- data.frame(z = wp.df$z, resid = wp.df$resid)
    lims.df <- lims.df %>% group_by(z) %>% reframe(lims.df = get.lims(resid)) %>% 
      tidyr::unnest(cols = c(lims.df))
    wp.dt <- wp.df %>% group_by(z) %>% arrange(x, .by_group = TRUE) %>% 
      as.data.table()
    lims.dt <- as.data.table(lims.df)
    combo.dt <- lims.dt[wp.dt, on = .(z == z, zval == x), 
                        roll = TRUE]
    if (sum(is.na(combo.dt)) > 0) {
      warning("missing some CI values, try increasing xlim.worm")
    }
    n_outer <- combo.dt %>% mutate(outer = ifelse((y < low | 
                                                     y > high), 1, 0)) %>% group_by(z) %>% summarise(n = n(), 
                                                                                                     n_outer = sum(outer)) %>% mutate(pcnt = n_outer/n)
    p <- ggplot(wp.df, aes(x = x, y = y)) + geom_smooth(method = lm, 
                                                        formula = y ~ poly(x, 3)) + geom_point() + facet_wrap(~z) + 
      theme_classic() + xlab("Unit Normal Quantile") + 
      ylab("Deviation") + {
        if (is.finite(xlim.worm)) 
          xlim(c(-xlim.worm, xlim.worm))
      } + {
        if (is.finite(ylim.worm)) {
          ylim(c(-ylim.worm, ylim.worm))
        }
        else {
          ylim(c(-1, 1))
        }
      } + geom_line(data = lims.df, aes(x = zval, y = low), 
                    linetype = "dashed") + geom_line(data = lims.df, 
                                                     aes(x = zval, y = high), linetype = "dashed") + geom_hline(yintercept = 0, 
                                                                                                                linetype = "dashed") + geom_text(data = n_outer, 
                                                                                                                                                 mapping = aes(x = Inf, y = Inf, label = scales::percent(pcnt)), 
                                                                                                                                                 hjust = 1.5, vjust = 1.5, label.size = 0.15, color = "blue")
    
    
    if (plot == TRUE) {
      return(p)
    }  else{
      return(n_outer)
    }    
    
    }
}