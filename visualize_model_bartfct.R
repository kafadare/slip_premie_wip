visualize_model <- function(modobj,smooth_var, int_var ,
                            #group_var, 
                            plabels = NULL,check_diagnostics = F,derivative_plot = F){
  
  this_font_size = font_size*1.25
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  s<-summary(model)
  
  ## Generate custom line plot
  np <- 10000 #number of predicted values
  df = model$model
  
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  if (!is.null(int_var)) {
    # We will produce and interaction plot
    if (!any(grepl(x=as.character(model$formula),pattern = int_var))) {
      warning("int_var not recognized in model formula!")
      return()
    }
    switch (varClasses[int_var],
            "numeric" = {
              q <- quantile(df[,int_var],probs = c(.10,0.50,.90)) #pick 10%, 50%, and 90% to plot
              bigq <- q[[3]]
              midq <- q[[2]]
              smallq <- q[[1]]
              values <- c(bigq,midq, smallq)
              labs <- c(sprintf("high (%1.2f)",bigq),sprintf("mid (%1.2f)",midq), sprintf("low (%1.2f)",smallq))
              
              q <-quantile(rescale(df[,int_var],c(0,1)),probs = c(0,.5,1))
              limit_values <- c(q[[1]],q[[length(q)]])
              midpoint_val <- unname(q[[2]])
              cbar_vals <- unname(q)
              
              theseLabs = rep(values,each = np)
              grad_fill = T
            },
            "factor" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = rep(values,each = np)
              grad_fill = F
            },
            "ordered" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = ordered(rep(values,each = np),levels = values)
              grad_fill = F
            }
    )
    
    labPred <- data.frame(init = rep(0,np*length(labs)))
    labPred[,int_var] = theseLabs
    labPred$lab = rep(labs,each = np)
    labPred <- labPred[,names(labPred) !="init"]
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else if (v == int_var) {
        next
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    
    thisPred <- thisPred %>% select(-init)
    thisPred <- do.call("rbind", replicate(length(labs), thisPred, simplify = FALSE))
    
    pred <- cbind(labPred,thisPred)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    #pred[,group_var] = NA #these columns have to exist in the dataframe for plotting
    pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
    
    low_color  = "#91bfdb"   # muted blue (keep or slightly deepen)
    mid_color  = "#b39ddb"   # muted purple (balanced midpoint)
    high_color = "#fc8d59"   # muted orange-red (keep)
    
    low_line   = "#4575b4"   # darker blue
    mid_line   = "#7b6ba8"   # darker muted purple
    high_line  = "#f46d43"   # darker orange-red
    
    if (grad_fill == T) {
      p1 <- ggplot(data = df, aes(x = get(smooth_var),y = get(thisResp), color = get(int_var))) +
        geom_point(alpha = 0.55, stroke = 0, size = point_size) + 
        #geom_line(aes_string(group = group_var),alpha = .5) +
        scale_color_gradientn(colors = c(low_color, mid_color, high_color), values = cbar_vals,name = "") +
        geom_ribbon(data = pred,aes(x = get(smooth_var), ymin = get("selo"),ymax = get("sehi"), fill = get("lab")),
                    alpha = .18, linetype = 0) +
        scale_fill_manual(values = c(high_color,mid_color, low_color)) +
        geom_line(data = pred,aes(x = get(smooth_var), y = get("fit"),group = get("lab")),linewidth = line_size) +
        labs(title = plabels)
    } else {
      
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) +
        geom_point(alpha = .35,stroke = 0, size = point_size) + 
        #geom_line(aes_string(group = group_var),alpha = .3) +
        scale_color_brewer(type = "qual",palette = "Set1",direction = 1) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = int_var),alpha = .5, linetype = 0) +
        scale_fill_brewer(type = "qual",palette = "Set1",direction = 1) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",color = int_var),size = line_size) +
        labs(title = plabels)
    }
  } else {
    
    # No interaction variable, just produce a single line plot
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    pred <- thisPred %>% select(-init)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    #pred[,group_var] = NA
    pred[,thisResp] = 1
    
    p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp)) +
      geom_point(alpha = .3,stroke = 0, size = point_size) +
      #geom_line(aes_string(group = group_var),alpha = .3) +
      geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi"),alpha = .5, linetype = 0) +
      geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),size = line_size) +
      labs(title = plabels)
  }
  
  if (derivative_plot == T) {
    # We will add a bar that shows where the derivative is significant.
    # First make some adjustments to the line plot.
    p1<- p1+theme(text = element_text(size=this_font_size),
                  axis.text = element_text(size = this_font_size),
                  axis.title.y = element_text(size = this_font_size),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.text = element_text(size = this_font_size),
                  legend.title = element_text(size = this_font_size),
                  axis.title = element_text(size = this_font_size),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(.2, .2, 0, .2), "cm")) #Top, left,Bottom, right
    scatter <- list(p1)
    
    # Now add the plots using the derivative plotting function
    if (any(grepl(x = row.names(s$s.table),pattern =  ":") & grepl(x=row.names(s$s.table),pattern = int_var))) {
      # Factor levels separately if there is an interaction in the model.
      f<-formula(model) # current formula
      fterms <- terms(f)
      fac <- attr(fterms, "factors")
      idx <- which(as.logical(colSums(fac[grep(x=row.names(fac),pattern = int_var),])))
      new_terms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
      new_formula <- formula(new_terms) # Formula without any interaction terms in the model.
      
      #add derivative gradients for each level of the factor
      num_levels <- length(levels(df[,int_var]))
      level_colors <- suppressWarnings(RColorBrewer::brewer.pal(num_levels,"Set1")) #use the same palette as the line plot
      plotlist = vector(mode = "list",length = num_levels+1) # we will be building a list of plots
      plotlist[1] = scatter # first the scatter plot
      
      for (fcount in 1:num_levels) {
        this_level <- levels(df[,int_var])[fcount]
        df$subset <- df[,int_var] == this_level
        #df$group_var <- df[,group_var]
        #this_mod <- gamm(formula = new_formula,data = df,subset = subset,random=list(group_var=~1))
        this_mod <- gamm(formula = new_formula,data = df,subset = subset)
        # this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        
        if (fcount != num_levels & fcount != 1){
          # get rid of redundant junk
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          this_d$theme$legend.text = element_blank()
        }
        if (fcount == 1) {
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.text = element_blank()
        }
        if (fcount == num_levels) {
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
        }
        this_d$labels$fill=NULL
        plotlist[fcount+1] = list(this_d)
      }
      pg<-plot_grid(rel_heights = c(16*num_levels,rep(num_levels,num_levels-1),3*num_levels),plotlist = plotlist,align = "v",axis = "lr",ncol = 1)
      final_plot <- pg
      #print(final_plot)
    } else {
      # No need to split
      d1 <- get_derivs_and_plot(modobj = modobj,smooth_var = smooth_var)
      scatter <- list(p1)
      bar <- list(d1)
      allplots <- c(scatter,bar)
      pg<-plot_grid(rel_heights = c(16,3),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
      final_plot <- pg
      #print(final_plot)
    }
    
  }    else {
    # No derivative plot
    p1<- p1+theme(text = element_text(size=font_size),
                  axis.text = element_text(size = font_size),
                  legend.text = element_text(size = font_size),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank())
    final_plot <- p1
    #print(final_plot)
  }
  
  if (check_diagnostics == T) {
    cp <- check(b,
                a.qq = list(method = "tnorm",
                            a.cipoly = list(fill = "light blue")),
                a.respoi = list(size = 0.5),
                a.hist = list(bins = 10))
    print(cp)
  }
  return(final_plot)
}