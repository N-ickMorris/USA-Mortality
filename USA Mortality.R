# -----------------------------------------------------------------------------------
# ---- Input Information ------------------------------------------------------------
# -----------------------------------------------------------------------------------

# set the path of where the input files are
mywd = "C:/ ... /USA Mortality Reports"

# -----------------------------------------------------------------------------------
# ---- Packages ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# data handling
require(data.table)
require(tm)
require(gtools)
require(stringdist)
require(zoo)
require(missRanger)
require(TTR)

# plotting
require(ggplot2)
require(gridExtra)
require(GGally)
require(scales)
require(scatterplot3d)
require(gridGraphics)
require(corrplot)
require(VIM)

# modeling
require(fpc)
require(caret)
require(ranger)
require(cluster)
require(car)
require(nortest)
require(neuralnet)
require(h2o)

# parallel computing
require(foreach)
require(parallel)
require(doSNOW)

}

# -----------------------------------------------------------------------------------
# ---- Functions --------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# these are functions i like to use

# ---- prints the data types of each column in a data frame -------------------------

types = function(dat)
{
  # make dat into a data.frame
  dat = data.frame(dat)
  
  # get the column names
  column = colnames(dat)
  
  # get the class of the columns
  data.type = sapply(1:ncol(dat), function(i) class(dat[,i]))
  
  # compute the number of levels for each column
  levels = sapply(1:ncol(dat), function(i) ifelse(data.type[i] == "factor", length(levels(droplevels(dat[,i]))), 0))
  
  return(data.frame(column, data.type, levels))
}

# ---- converts all columns to a character data type --------------------------------

tochar = function(dat)
{
  # make dat into a data.frame
  dat = data.frame(dat)
  
  # get the column names
  column = colnames(dat)
  
  # get the values in the columns and convert them to character data types
  values = lapply(1:ncol(dat), function(i) as.character(dat[,i]))
  
  # combine the values back into a data.frame
  dat = data.frame(do.call("cbind", values), stringsAsFactors = FALSE)
  
  # give dat its column names
  colnames(dat) = column
  
  return(dat)
}

# ---- a qualitative color scheme ---------------------------------------------------

mycolors = function(n)
{
  require(grDevices)
  return(colorRampPalette(c("#e41a1c", "#0099ff", "#4daf4a", "#984ea3", "#ff7f00", "#ff96ca", "#a65628"))(n))
}

# ---- emulates the default ggplot2 color scheme ------------------------------------

ggcolor = function(n, alpha = 1)
{
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# ---- plots various cluster plots --------------------------------------------------

plot.clusters = function(dat, cluster.column.name = NULL, distance.matrix = NULL, DC.title = "Discriminant Coordinate Cluster Plot", pairs.title = "Cluster Pairs Plot", silhouette.title = "Silhouette Width", font.size = 20, pairs.plot.font.size = 12, rotate = 0, cor.size = 3)
{
  # load packages we need
  require(data.table)
  require(ggplot2)
  require(GGally)
  require(scatterplot3d)
  require(gridGraphics)
  require(grid)
  require(scales)
  require(cluster)
  require(fpc)
  
  # this function emulates the default ggplot2 color scheme
  ggcolor = function(n)
  {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # error check
  if(is.null(cluster.column.name))
  {
    print("you must specify a value for the parameter: cluster.column.name")
    
  } else
  {
    # ---- computing discriminant coordiantes -------------------------------
    
    # make dat into a data table
    dat = data.table(dat)
    
    # extract the cluster column from dat
    clusters = as.numeric(unname(unlist(dat[, cluster.column.name, with = FALSE])))
    
    # remove the cluster column from dat
    dat = dat[, !cluster.column.name, with = FALSE]
    
    # compute the discriminant coordinates, and extract the first 3
    dat.dc = data.table(discrproj(x = dat, 
                                  clvecd = clusters,
                                  method = "dc")$proj[,1:3])
    
    # rename the columns appropriately
    setnames(dat.dc, paste0("DC", 1:3))
    
    # give dat.dc the cluster column, and make it into a factor for plotting purposes
    dat.dc[, Cluster := factor(clusters, levels = sort(unique(clusters)))]
    
    # ---- plotting 2D discriminant coordiantes -----------------------------
    
    # create a 2D cluster plot across the first 2 discriminant coordinates
    plot.2D = ggplot(dat.dc, aes(x = DC1, y = DC2, fill = Cluster)) +
      stat_density_2d(geom = "polygon", color = NA, alpha = 1/3) +
      theme_bw(font.size) +
      ggtitle(DC.title) +
      theme(legend.position = "top", 
            legend.key.size = unit(.25, "in"), 
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(fill = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1))
    
    # store this plot in a list
    output = list("plot.2D" = plot.2D)
    
    # ---- plotting 3D discriminant coordiantes -----------------------------
    
    # create a table indicating which cluster should get which ggplot color
    color.dat.dc = data.table(Cluster = levels(dat.dc$Cluster),
                              Color = ggcolor(length(levels(dat.dc$Cluster))))
    
    # set Cluster as the key column in color.dat.dc and dat.dc
    setkey(dat.dc, Cluster)
    setkey(color.dat.dc, Cluster)
    
    # join color.dat.dc onto dat.dc
    dat.dc = color.dat.dc[dat.dc]
    
    # convert Cluster back into a factor data type
    dat.dc[, Cluster := factor(Cluster, levels = sort(unique(Cluster)))]
    
    # here is my default font size for base R plotting
    font.default = 20
    
    # compute the desired adjustment according to the specified value of font.size
    font.adjust = font.size / font.default
    
    # adjust the font of the title, axis titles, axis labels, and legend
    font.title = 2 * font.adjust
    font.axis.title = 1.125 * font.adjust
    font.axis.label = 1.125 * font.adjust
    font.legend = 1.75 * font.adjust
    
    # here are my 4 default angles for viewing a 3D scatterplot
    angles = c(45, 135, 225, 315)
    
    # apply the specified rotation
    angles = angles + rotate
    
    # set up 4 plot ID numbers so each plot angle has a position in the plot window
    plot.id = 2:5
    
    # set up a legend ID number so the legend has a postion across the top of the plot window
    legend.id = c(1, 1)
    
    # set up a matrix that defines the plot layout
    plot.layout = matrix(c(legend.id, plot.id), nrow = 3, ncol = 2, byrow = TRUE)
    
    # create a new plot window
    windows()
    plot.new()
    
    # define plot margins
    par(mar = c(0, 0, 3.5, 0))
    
    # apply the layout to the plot window
    layout(mat = plot.layout, heights = c(1, 1.5, 1.5))
    
    # produce a dummy plot as a place holder for the legend and title
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    
    # produce the title
    title(main = DC.title, cex.main = font.title)
    
    # produce the legend
    legend("top", inset = 0, bty = "n", cex = font.legend, horiz = TRUE,
           title = "Clusters", legend = levels(dat.dc$Cluster), fill = ggcolor(length(levels(dat.dc$Cluster))))
    
    # create 3D cluster plots across the first 3 discriminant coordinates
    for(angle in angles)
    {
      # build the scatterplot
      scatterplot3d(x = dat.dc$DC1, y = dat.dc$DC2, z = dat.dc$DC3, color = alpha(dat.dc$Color, 1/3), 
                    xlab = "DC1", ylab = "DC2", zlab = "DC3", mar = c(3, 3, 0, 3), cex = 1.5,
                    pch = 16, cex.lab = font.axis.title, cex.axis = font.axis.label, angle = angle)
    }
    
    # save this plot as a grid object
    grid.echo()
    plot.3D = grid.grab()
    
    # add this plot to our output list
    output$plot.3D = plot.3D
    
    # close the graphics window
    graphics.off()
    
    # ---- plotting clusters across variable pairs ---------------------------
    
    # give dat the cluster column, and make it into a factor for plotting purposes
    dat[, Cluster := factor(clusters, levels = sort(unique(clusters)))]
    
    # plot the clusters across all variable pairs
    plot.pairs = ggpairs(dat,
                         mapping = aes(color = Cluster, fill = Cluster),
                         columns = which(names(dat) != "Cluster"),
                         lower = list(continuous = wrap(ggally_points, size = 1.5, alpha = 1/3)), 
                         upper = list(continuous = wrap(ggally_cor, size = cor.size)),
                         diag = list(continuous = wrap(ggally_densityDiag, alpha = 1/3)),
                         title = pairs.title,
                         legend = grab_legend(plot.2D + theme_classic(base_size = pairs.plot.font.size) + theme(legend.position = "top"))) + 
      theme_classic(base_size = pairs.plot.font.size) +
      theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
    
    # remove the cluster column from dat
    dat[, Cluster := NULL]
    
    # add this plot to our output list
    output$plot.pairs = plot.pairs
    
    # ---- plotting silhouette widths ----------------------------------------
    
    # if the user gave a distance matrix then lets compute silhouette widths
    if(!is.null(distance.matrix))
    {
      # compute the silhouette widths
      dat.sil = silhouette(x = clusters,
                           dist = distance.matrix)
      
      # compute the summary of dat.sil
      dat.sil.sum = summary(dat.sil)
      
      # extract the avg widths from dat.sil.sum
      dat.sil.avg = data.table(Cluster = as.numeric(names(dat.sil.sum$clus.avg.widths)),
                               Average_Width = round(unname(dat.sil.sum$clus.avg.widths), 2))
      
      # order dat.sil.avg by Cluster
      dat.sil.avg = dat.sil.avg[order(Cluster)]
      
      # extract the cluster sizes from dat.sil.sum
      dat.sil.size = data.table(Cluster = as.numeric(names(dat.sil.sum$clus.sizes)),
                                Size = as.numeric(unname(dat.sil.sum$clus.sizes)))
      
      # order dat.sil.size by Cluster
      dat.sil.size = dat.sil.size[order(Cluster)]
      
      # combine dat.sil.avg and dat.sil.size into a table that will go on our plot
      dat.sil.tab = cbind(dat.sil.avg, dat.sil.size[,!"Cluster"])
      
      # convert dat.sil into a data table
      dat.sil = data.table(dat.sil[1:nrow(dat.sil), 1:ncol(dat.sil)])
      
      # sort dat.sil by cluster and sil_width
      dat.sil = dat.sil[order(cluster, -sil_width)]
      
      # give dat.sil an ID column for plotting purposes
      dat.sil[, ID := 1:nrow(dat.sil)]
      
      # convert cluster to a factor for plotting purposes
      dat.sil[, cluster := factor(cluster, levels = sort(unique(cluster)))]
      
      # aggregate sil_width by cluster in dat.sil to determine where to place dat.sil.tab in the plot 
      dat.agg = dat.sil[, .(sil_width.min = min(sil_width),
                            sil_width.max = max(sil_width)),
                        by = cluster]
      
      # build the four corners of the dat.sil.tab to place it in the plot
      # find the cluster with the smallest peak and set the peak's sil_width as the ymin
      ymin = as.numeric(min(dat.agg$sil_width.max))
      
      # find the sil_width of the max peak and set it as the ymax
      ymax = as.numeric(max(dat.sil$sil_width))
      
      # extract the cluster with the smallest peak from dat.agg
      small.peak = dat.agg[which.min(sil_width.max), cluster]
      
      # find the first ID number in dat.sil for the cluster with the smallest peak, and set that as xmin
      xmin = min(dat.sil[cluster == small.peak, ID])
      
      # find the last ID number in dat.sil for the cluster with the smallest peak, and set that as xmax
      xmax = max(dat.sil[cluster == small.peak, ID])
      
      # plot the silhouette width and add the dat.sil.tab to it
      plot.sil.width = ggplot(dat.sil, aes(x = ID, y = sil_width, fill = cluster, color = cluster)) + 
        geom_bar(stat = "identity", position = "dodge") +
        # annotation_custom(tableGrob(as.matrix(dat.sil.tab), rows = NULL, 
        #                             theme = ttheme_default(base_size = font.size,
        #                                                    colhead = list(fg_params = list(col = "black"), bg_params = list(fill = "lightgray", col = "black")),
        #                                                    core = list(fg_params = list(hjust = 0.5), bg_params = list(fill = c("white"), col = "black")))),
        #                   xmin = xmin, 
        #                   xmax = xmax, 
        #                   ymin = ymin, 
        #                   ymax = ymax) +
        ggtitle(silhouette.title) +
        labs(x = "Observation", y = "Silhouette Width", fill = "Cluster", color = "Cluster") +
        theme_bw(base_size = font.size) +
        theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
        guides(fill = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1),
               color = guide_legend(nrow = 1))
      
      # add this plot to our output list
      output$plot.sil.width = plot.sil.width
    }
    
    # add the DC data to the output list
    output$dat.dc = dat.dc
    
    # helpful message
    print("use grid.draw() to see the 3D cluster plot, for example: grid.draw(my.plot.clusters$plot.3D)")
    
    return(output)
  }
}

# ---- prints out a dat file object in ampl syntax ----------------------------------

ampl = function(dat, object = "param", name = "c")
{
  # converts all columns to a character data type 
  tochar = function(dat)
  {
    # make dat into a data.frame
    dat = data.frame(dat)
    
    # get the column names
    column = colnames(dat)
    
    # get the values in the columns and convert them to character data types
    values = lapply(1:ncol(dat), function(i) as.character(dat[,i]))
    
    # combine the values back into a data.frame
    dat = data.frame(do.call("cbind", values), stringsAsFactors = FALSE)
    
    # give dat its column names
    colnames(dat) = column
    
    return(dat)
  }
  
  # make sure the data is a data frame object
  dat = tochar(dat)
  
  # every parameter/set object in an ampl dat file must end with a semicolon
  # so set up 1 semicolon to give to dat
  semicolon = c(";", rep(" ", ncol(dat) - 1))
  
  # add this semicolon as the last row of the data frame
  result = data.frame(rbind(dat, semicolon))
  
  # every parameter/set object in an ample dat file must begin with the name of the object and what it equals
  # for example: param c := 
  # so set up a header to give to dat
  header = c(paste(object, name, ":="), rep(" ", ncol(dat) - 1))
  
  # update the column names of dat to be the header we created
  colnames(result) = header
  
  # print out the result without any row names
  # print out the result left adjusted
  # print(result, right = FALSE, row.names = FALSE)
  
  return(result)	
}

# ---- compares the quantiles of emprical data against the quantiles of any statistical distribution 

ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), alpha = 0.33, basefont = 20, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
{
  require(ggplot2)
  
  # compute the sample quantiles and theoretical quantiles
  q.function = eval(parse(text = paste0("q", distribution)))
  d.function = eval(parse(text = paste0("d", distribution)))
  x = na.omit(x)
  ord = order(x)
  n = length(x)
  P = ppoints(length(x))
  df = data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  # compute the quantile line
  Q.x = quantile(df$ord.x, c(probs[1], probs[2]))
  Q.z = q.function(c(probs[1], probs[2]), ...)
  b = diff(Q.x) / diff(Q.z)
  coef = c(Q.x[1] - (b * Q.z[1]), b)
  
  # compute the confidence interval band
  zz = qnorm(1 - (1 - conf) / 2)
  SE = (coef[2] / d.function(df$z, ...)) * sqrt(P * (1 - P) / n)
  fit.value = coef[1] + (coef[2] * df$z)
  df$upper = fit.value + (zz * SE)
  df$lower = fit.value - (zz * SE)
  
  # plot the qqplot
  p = ggplot(df, aes(x = z, y = ord.x)) + 
    geom_point(color = "blue", alpha = alpha) +
    geom_abline(intercept = coef[1], slope = coef[2], size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    coord_cartesian(ylim = c(min(df$ord.x), max(df$ord.x))) + 
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # conditional additions
  if(main != "")(p = p + ggtitle(main))
  
  return(p)
}

# ---- plots 6 residual plots -------------------------------------------------------

residplots = function(actual, fit, binwidth = NULL, from = NULL, to = NULL, by = NULL, histlabel.y = -10, n = NULL, basefont = 20)
{
  require(ggplot2)
  
  residual = actual - fit 
  DF = data.frame("actual" = actual, "fit" = fit, "residual" = residual)
  
  rvfPlot = ggplot(DF, aes(x = fit, y = residual)) + 
    geom_point(na.rm = TRUE) +
    stat_smooth(method = "loess", se = FALSE, na.rm = TRUE, color = "blue") +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Fitted values") +
    ylab("Residuals") +
    ggtitle("Residual vs Fitted Plot") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), note = TRUE, alpha = 0.33, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
  {
    # compute the sample quantiles and theoretical quantiles
    q.function = eval(parse(text = paste0("q", distribution)))
    d.function = eval(parse(text = paste0("d", distribution)))
    x = na.omit(x)
    ord = order(x)
    n = length(x)
    P = ppoints(length(x))
    DF = data.frame(ord.x = x[ord], z = q.function(P, ...))
    
    # compute the quantile line
    Q.x = quantile(DF$ord.x, c(probs[1], probs[2]))
    Q.z = q.function(c(probs[1], probs[2]), ...)
    b = diff(Q.x) / diff(Q.z)
    coef = c(Q.x[1] - (b * Q.z[1]), b)
    
    # compute the confidence interval band
    zz = qnorm(1 - (1 - conf) / 2)
    SE = (coef[2] / d.function(DF$z, ...)) * sqrt(P * (1 - P) / n)
    fit.value = coef[1] + (coef[2] * DF$z)
    DF$upper = fit.value + (zz * SE)
    DF$lower = fit.value - (zz * SE)
    
    # plot the qqplot
    p = ggplot(DF, aes(x = z, y = ord.x)) + 
      geom_point(color = "black", alpha = alpha) +
      geom_abline(intercept = coef[1], slope = coef[2], size = 1, color = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
      coord_cartesian(ylim = c(min(DF$ord.x), max(DF$ord.x))) + 
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # conditional additions
    if(main != "")(p = p + ggtitle(main))
    
    return(p)
  }
  
  qqPlot = ggqq(residual, 
                alpha = 1,				  
                main = "Normal Q-Q Plot", 
                xlab = "Theoretical Quantiles", 
                ylab = "Residuals")
  
  rvtPlot = ggplot(data.frame("x" = 1:length(DF$residual), "y" = DF$residual), aes(x = x, y = y)) + 
    geom_line(na.rm = TRUE) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Obs. Number") +
    ylab("Residuals") +
    ggtitle("Residual Time Series") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  variogramDF = function(x)
  {
    n = length(x) - 2
    
    num = sapply(1:n, function(k)
      sapply(1:(length(x) - k), function(i)
        x[i + k] - x[i]))
    
    num = sapply(1:length(num), function(j)
      var(num[[j]]))
    
    den = var(sapply(1:(length(x) - 1), function(i)
      x[i + 1] - x[i]))
    
    val = num / den
    
    DF = data.frame("Lag" = 1:n, "Variogram" = val)
    
    return(DF)
  }
  
  DFv = variogramDF(x = DF$residual)
  
  varioPlot = ggplot(DFv, aes(x = Lag, y = Variogram)) + 
    geom_point() +
    geom_line(color = "blue") +
    xlab("Lag") +
    ylab("Variogram") +
    ggtitle("Variogram of Residuals") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  test = t.test(DF$residual)
  
  CI = data.frame("x" = test$estimate, 
                  "LCB" = test$conf.int[1], 
                  "UCB" = test$conf.int[2], 
                  row.names = 1)
  
  histPlot = ggplot(DF, aes(x = residual)) +
    geom_histogram(color = "white", fill = "black", binwidth = binwidth) +
    geom_segment(data = CI, aes(x = LCB, xend = LCB, y = 0, yend = Inf), color = "blue") +
    geom_segment(data = CI, aes(x = UCB, xend = UCB, y = 0, yend = Inf), color = "blue") +
    annotate("text", x = CI$x, y = histlabel.y, 
             label = "T-Test C.I.", size = 5, 
             color = "blue", fontface = 2) + 
    ggtitle("Residual Histogram") +
    labs(x = "Residuals", y = "Frequency") +
    theme_bw(base_size = basefont) +
    theme(legend.key.size = unit(.25, "in"),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(class(from) != "NULL" & class(to) != "NULL" & class(by) != "NULL") (histPlot = histPlot + scale_x_continuous(breaks = seq(from = from, to = to, by = by)))
  
  ggacf = function(x, n = NULL, conf.level = 0.95, main = "ACF Plot", xlab = "Lag", ylab = "Autocorrelation", basefont = 20) 
  {
    if(class(n) == "NULL")
    {
      n = length(x) - 2
    }
    
    ciline = qnorm((1 - conf.level) / 2) / sqrt(length(x))
    bacf = acf(x, lag.max = n, plot = FALSE)
    bacfdf = with(bacf, data.frame(lag, acf))
    bacfdf = bacfdf[-1,]
    
    p = ggplot(bacfdf, aes(x = lag, y = acf)) + 
      geom_bar(stat = "identity", position = "dodge", fill = "black") +
      geom_hline(yintercept = -ciline, color = "blue", size = 1) +
      geom_hline(yintercept = ciline, color = "blue", size = 1) +
      geom_hline(yintercept = 0, color = "red", size = 1) +
      labs(x = xlab, y = ylab) +
      ggtitle(main) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    return(p)
  }
  
  acfPlot = ggacf(x = DF$residual, main = "ACF Plot of Residuals", basefont = basefont, n = n)
  
  return(list("rvfPlot" = rvfPlot, 
              "qqPlot" = qqPlot, 
              "rvtPlot" = rvtPlot, 
              "varioPlot" = varioPlot, 
              "histPlot" = histPlot, 
              "acfPlot" = acfPlot))
}

}

# -----------------------------------------------------------------------------------
# ---- Prepare the Data -------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# set the work directory
setwd(mywd)

# import all of the data
dat = data.table(read.csv("USA Mortality.csv", stringsAsFactors = FALSE))
va.pop =  data.table(read.csv("VA-Population.csv", stringsAsFactors = FALSE))
va.dem =  data.table(read.csv("VA-Dementia.csv", stringsAsFactors = FALSE))

# update population to be numeric
dat[, Population := as.numeric(Population)]

# remove missing values
dat = na.omit(dat)

# update indicators that should be factors
dat[, Cause := as.factor(Cause)]
dat[, Age := factor(Age, levels = unique(Age))]
dat[, Gender := as.factor(Gender)]
dat[, Race := as.factor(Race)]
dat[, Hispanic := as.factor(Hispanic)]

# combine VA data
# get the VA population numbers for 65+
va.pop = melt(data = va.pop[,.(Year, Age_65.69, Age_70.74, Age_75.79, Age_80.84, Age_85.)], 
              id.vars = "Year", 
              variable.name = "Age", 
              value.name = "Population")

# add up population numbers across all age groups by year
va.pop = va.pop[,.(Age = "65+",
                   Population = sum(Population)),
                by = .(Year)]

# add up dementia prevelence across all age groups by year
va.dem = va.dem[,.(Age = "65+",
                   Veterans_LowEstimate = sum(Veterans_LowEstimate),
                   Veterans_MiddleEstimate = sum(Veterans_MiddleEstimate),
                   Veterans_HighEstimate = sum(Veterans_HighEstimate)), 
                by = .(Year)]

# inner join va.pop onto va.dem
setkey(va.pop, Year, Age)
setkey(va.dem, Year, Age)
va = data.table(va.pop[va.dem, nomatch = 0])

# lets combine USA Alzhiemers Mortality data with the VA data
# this allows us to contrast the prevelence of dementia in US citizens vs. US military veterans
# note: Alzhiemers is the most prevelent of all dementias in the USA
alz = data.table(dat[Cause %in% c("All causes of death", "Alzheimer's disease") &
                     Age %in% c("65-74 years", "75-84 years", "85+ years"),
                     .(Age = "65+",
                       Deaths = sum(Deaths),
                       Population = sum(Population)),
                     by = .(Cause, Year)])

# update alz to compare Alzheimer's deaths to All causes of death
alz = cbind(alz[Cause %in% "Alzheimer's disease", !"Population"], 
            alz[Cause %in% "All causes of death", .(allDeaths = Deaths)])

# update alz such that it can be combined with va
alz = alz[,.(Group = "Alzheimers Prevelence in Deceased USA Citizens",
             Year = Year, 
             Age = Age,
             Population = allDeaths,
             LowEstimate = Deaths,
             MiddleEstimate = Deaths,
             HighEstimate = Deaths)]

# update va such that it can be combined with alz
va = va[,.(Group = "Expected Dementia Prevelence in Living USA Military Veterans",
           Year = Year, 
           Age = Age,
           Population = Population,
           LowEstimate = Veterans_LowEstimate,
           MiddleEstimate = Veterans_MiddleEstimate,
           HighEstimate = Veterans_HighEstimate)]

# combine alz and va
va = rbind(va, alz)

# lets update Population in dat to take on the value of "All causes of death"
# get a subset of "All causes of death"
all.deaths = data.table(dat[Cause == "All causes of death", .(Age, Gender, Hispanic, Race, Year, allDeaths = Deaths)])

# remove "All causes of death" from dat
dat = dat[Cause != "All causes of death"]

# give dat an id column preserve current row order
dat[, id := 1:nrow(dat)]

# join allDeaths onto dat
setkey(all.deaths, Age, Gender, Hispanic, Race, Year)
setkey(dat, Age, Gender, Hispanic, Race, Year)
dat = all.deaths[dat]

# order dat by id then remove id
dat = dat[order(id)]
dat[, id := NULL]

}

# -----------------------------------------------------------------------------------
# ---- Study USA Mortality ----------------------------------------------------------
# -----------------------------------------------------------------------------------

{
  
# lets compare all possible population pairs across Gender, Race, and Hispanic for older people (age: 65+)
# lets ensure that when we compare two populations they share the same Gender or Ethnicity (Race & Hispanic)
# set up all possible population subsets
DT = data.table(expand.grid(Gender = c("Male", "Female", "All"),
                            Ethnicity = c("American Indian or Alaska Native", "Asian or Pacific Islander",       
                                          "Black or African American", "White", "Hispanic or Latino", 
                                          "Not Hispanic or Latino", "All")))

# give DT an id number for each population subset
DT[, id := 1:nrow(DT)]

# set up all possible population pairs based on id from DT
comb = expand.grid(pop1 = DT$id, pop2 = DT$id)

# row sort comb
comb = setnames(data.table(t(apply(comb, 1, sort))), c("pop1", "pop2"))

# remove all duplicates
comb = comb[!duplicated(comb)]

# remove instances where a population is compared to itself
comb = comb[!(pop1 == pop2)]

# join DT onto comb by id for pop1
setnames(comb, c("id", "pop2"))
setkey(comb, id)
setkey(DT, id)
comb = DT[comb]

# join DT onto comb by id for pop2
setnames(comb, c("Gender1", "Ethnicity1", "pop1", "id"))
setkey(comb, id)
setkey(DT, id)
comb = DT[comb]

# update column names in comb and remove id and pop1
comb[, id := NULL]
comb[, pop1 := NULL]
setnames(comb, c("Gender2", "Ethnicity2", "Gender1", "Ethnicity1"))

# remove instances where All is compared to a subset of Gender
comb = comb[!(Gender1 == "All" & Gender2 != "All")]
comb = comb[!(Gender2 == "All" & Gender1 != "All")]

# remove instances where All is compared to a subset of Ethnicity
comb = comb[!(Ethnicity1 == "All" & Ethnicity2 != "All")]
comb = comb[!(Ethnicity2 == "All" & Ethnicity1 != "All")]

# remove instances where Not Hispanic or Latino is not compared to Hispanic or Latino
comb = comb[!(Ethnicity1 == "Not Hispanic or Latino" & Ethnicity2 != "Hispanic or Latino")]
comb = comb[!(Ethnicity2 == "Not Hispanic or Latino" & Ethnicity1 != "Hispanic or Latino")]

# only keep instances where the pair of populations share one demographic indicator
comb = comb[Gender1 == Gender2 | Ethnicity1 == Ethnicity2]

# convert all columns into character data types
comb = data.table(tochar(comb))

# choose the number of workers and tasks for parallel processing
workers = max(1, floor((2/3) * detectCores()))
# workers = 1
tasks = nrow(comb)

# create a log file
myfile = "log.txt"
file.create(myfile)

# set up a cluster if workers > 1, otherwise don't set up a cluster
if(workers > 1)
{
  # setup parallel processing
  cl = makeCluster(workers, type = "SOCK", outfile = "")
  registerDoSNOW(cl)
  
  # define %dopar%
  `%fun%` = `%dopar%`
  
  # write out start time to log file
  sink(myfile, append = TRUE)
  cat("\n------------------------------------------------\n")
  cat("fisher's test\n")
  cat(paste(workers, "workers started at", Sys.time(), "\n"))
  sink()
  
} else
{
  # define %do%
  `%fun%` = `%do%`
  
  # write out start time to log file
  sink(myfile, append = TRUE)
  cat("\n------------------------------------------------\n")
  cat("fisher's test\n")
  cat(paste("task 1 started at", Sys.time(), "\n"))
  sink()
}

# run a fishers test on each pair of populations in comb
fishers = foreach(i = 1:tasks) %fun%
{
  # load packages we need for our tasks
  require(data.table)
  require(ggplot2)
  require(stats)
  
  # lets extract older Hispanics and older whites who have died
  test = data.table(dat[Age %in% c("65-74 years", "75-84 years", "85+ years"),
                        .(Deaths = sum(Deaths),
                          allDeaths = sum(allDeaths)),
                        by = .(Cause, Race, Hispanic, Gender, Year)])
  
  # aggregate test by ignoring Gender if Gender1 equals All
  if(comb$Gender1[i] == "All")
  {
    test = test[,.(Deaths = sum(Deaths),
                  allDeaths = sum(allDeaths)),
                by = .(Cause, Race, Hispanic, Year)]
    
    # aggregate test by ignoring Race if Ethnicity2 equals Not Hispanic or Latino
    if(comb$Ethnicity2[i] == "Not Hispanic or Latino")
    {
      test = test[,.(Deaths = sum(Deaths),
                     allDeaths = sum(allDeaths)),
                  by = .(Cause, Hispanic, Year)]
      
      # relabel Hispanic as Ethnicity
      test[, Ethnicity := Hispanic]
      test[, Hispanic := NULL]
      
      # aggregate test by ignoring Hispanic if Ethnicity2 doesn't equal Hispanic or Latino
    } else if(comb$Ethnicity2[i] != "Hispanic or Latino")
    {
      test = test[,.(Deaths = sum(Deaths),
                     allDeaths = sum(allDeaths)),
                  by = .(Cause, Race, Year)]
      
      # relabel Race as Ethnicity
      test[, Ethnicity := Race]
      test[, Race := NULL]
      
    } else
    {
      # remove instances that aren't Ethnicity1 and aren't Hispanic for Race and Hispanic respectively
      test = test[!(Race != comb$Ethnicity1[i] & Hispanic == "Not Hispanic or Latino")]
      
      # remove instances that are both Ethnicity1 and Hispanic for Race and Hispanic respectively
      test = test[!(Race == comb$Ethnicity1[i] & Hispanic == comb$Ethnicity2[i])]
      
      # create an ethnicity column for Ethnicity1 and Hispanic
      test[, Ethnicity := ifelse(Race == comb$Ethnicity1[i], comb$Ethnicity1[i], "Hispanic or Latino")]
      test[, Ethnicity := as.factor(Ethnicity)]
      
      # lets aggregate on Ethnicity
      test = data.table(test[, .(Deaths = sum(Deaths), allDeaths = sum(allDeaths)),
                             by = .(Cause, Ethnicity, Year)])
      
    }
    
    # get the timeline for all causes of death
    timeline = expand.grid(Cause = unique(test$Cause), Year = unique(test$Year))
    
    # build a matrix showing the death numbers and living numbers
    fisher.matrix = lapply(1:nrow(timeline), function(j)
      matrix(c(unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity1[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity1[i], .(allDeaths)]) - unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity1[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity2[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity2[i], .(allDeaths)]) - unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity2[i], .(Deaths)])), 
             ncol = 2))
    
    # lets run a fisher test to compare population mortality rates
    fisher.results = lapply(1:nrow(timeline), function(j) fisher.test(fisher.matrix[[j]]))
    names(fisher.results) = paste(timeline$Cause, timeline$Year)
    
    # create a table of annual odds ratios for each fishers test
    DT = data.table(oddsRatio = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$estimate)),
                    CI95LB = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$conf.int)[1]),
                    CI95UB = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$conf.int)[2]))
    
    # add more problem data to DT
    DT[, Ratio := paste0("[", comb$Ethnicity1[i], " Mortality] / [", comb$Ethnicity2[i], " Mortality]")]
    DT[, Cause := timeline$Cause]
    DT[, Year := timeline$Year]
    
    # open up a graphics window
    # windows()
    
    # lets plot the odds ratio
    DT.plot = ggplot(DT, aes(x = Year, y = oddsRatio, color = Cause, fill = Cause, group = Cause)) + 
      geom_point(size = 3, na.rm = TRUE) +
      geom_line(size = 1, na.rm = TRUE) +
      geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "black") + 
      geom_ribbon(aes(ymin = CI95LB, ymax = CI95UB), alpha = 1/3, color = NA) + 
      # scale_color_manual(values = color.set) +
      # scale_y_continuous(labels = dollar) + 
      ggtitle("Multiple Causes of Death in the USA") + 
      labs(x = "Year", y = DT$Ratio[1], color = "", fill = "") + 
      theme_bw(base_size = 30) +
      theme(legend.position = "top", 
            legend.key.size = unit(.25, "in"),
            plot.title = element_text(hjust = 0.5),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
    
    # aggregate test by ignoring Race and Hispanic if Ethnicity1 equals All
  } else if(comb$Ethnicity1[i] == "All")
  {
    test = test[,.(Deaths = sum(Deaths),
                   allDeaths = sum(allDeaths)),
                by = .(Cause, Gender, Year)]
    
    # get the timeline for all causes of death
    timeline = expand.grid(Cause = unique(test$Cause), Year = unique(test$Year))
    
    # build a matrix showing the death numbers and living numbers
    fisher.matrix = lapply(1:nrow(timeline), function(j)
      matrix(c(unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Gender == comb$Gender1[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Gender == comb$Gender1[i], .(allDeaths)]) - unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Gender == comb$Gender1[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Gender == comb$Gender2[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Gender == comb$Gender2[i], .(allDeaths)]) - unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Gender == comb$Gender2[i], .(Deaths)])), 
             ncol = 2))
    
    # lets run a fisher test to compare population mortality rates
    fisher.results = lapply(1:nrow(timeline), function(j) fisher.test(fisher.matrix[[j]]))
    names(fisher.results) = paste(timeline$Cause, timeline$Year)
    
    # create a table of annual odds ratios for each fishers test
    DT = data.table(oddsRatio = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$estimate)),
                    CI95LB = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$conf.int)[1]),
                    CI95UB = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$conf.int)[2]))
    
    # add more problem data to DT
    DT[, Ratio := paste0("[", comb$Gender1[i], " Mortality] / [", comb$Gender2[i], " Mortality]")]
    DT[, Cause := timeline$Cause]
    DT[, Year := timeline$Year]
    
    # open up a graphics window
    # windows()
    
    # lets plot the odds ratio
    DT.plot = ggplot(DT, aes(x = Year, y = oddsRatio, color = Cause, fill = Cause, group = Cause)) + 
      geom_point(size = 3, na.rm = TRUE) +
      geom_line(size = 1, na.rm = TRUE) +
      geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "black") + 
      geom_ribbon(aes(ymin = CI95LB, ymax = CI95UB), alpha = 1/3, color = NA) + 
      # scale_color_manual(values = color.set) +
      # scale_y_continuous(labels = dollar) + 
      ggtitle("Multiple Causes of Death in the USA") + 
      labs(x = "Year", y = DT$Ratio[1], color = "", fill = "") + 
      theme_bw(base_size = 30) +
      theme(legend.position = "top", 
            legend.key.size = unit(.25, "in"),
            plot.title = element_text(hjust = 0.5),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
    
  } else
  {
    # aggregate test by ignoring Race if Ethnicity2 equals Not Hispanic or Latino
    if(comb$Ethnicity2[i] == "Not Hispanic or Latino")
    {
      test = test[,.(Deaths = sum(Deaths),
                     allDeaths = sum(allDeaths)),
                  by = .(Cause, Gender, Hispanic, Year)]
      
      # relabel Hispanic as Ethnicity
      test[, Ethnicity := Hispanic]
      test[, Hispanic := NULL]
      
      # aggregate test by ignoring Hispanic if Ethnicity2 doesn't equal Hispanic or Latino
    } else if(comb$Ethnicity2[i] != "Hispanic or Latino")
    {
      test = test[,.(Deaths = sum(Deaths),
                     allDeaths = sum(allDeaths)),
                  by = .(Cause, Gender, Race, Year)]
      
      # relabel Race as Ethnicity
      test[, Ethnicity := Race]
      test[, Race := NULL]
      
    } else
    {
      # remove instances that aren't Ethnicity1 and aren't Hispanic for Race and Hispanic respectively
      test = test[!(Race != comb$Ethnicity1[i] & Hispanic == "Not Hispanic or Latino")]
      
      # remove instances that are both Ethnicity1 and Hispanic for Race and Hispanic respectively
      test = test[!(Race == comb$Ethnicity1[i] & Hispanic == comb$Ethnicity2[i])]
      
      # create an ethnicity column for Ethnicity1 and Hispanic
      test[, Ethnicity := ifelse(Race == comb$Ethnicity1[i], comb$Ethnicity1[i], "Hispanic or Latino")]
      test[, Ethnicity := as.factor(Ethnicity)]
      
      # lets aggregate on Ethnicity
      test = data.table(test[, .(Deaths = sum(Deaths), allDeaths = sum(allDeaths)),
                             by = .(Cause, Gender, Ethnicity, Year)])
      
    }
    
    # get the timeline for all causes of death
    timeline = expand.grid(Cause = unique(test$Cause), Year = unique(test$Year))
    
    # build a matrix showing the death numbers and living numbers
    fisher.matrix = lapply(1:nrow(timeline), function(j)
      matrix(c(unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity1[i] & Gender == comb$Gender1[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity1[i] & Gender == comb$Gender1[i], .(allDeaths)]) - unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity1[i] & Gender == comb$Gender1[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity2[i] & Gender == comb$Gender2[i], .(Deaths)]),
               unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity2[i] & Gender == comb$Gender2[i], .(allDeaths)]) - unlist(test[Cause == timeline$Cause[j] & Year == timeline$Year[j] & Ethnicity == comb$Ethnicity2[i] & Gender == comb$Gender2[i], .(Deaths)])), 
             ncol = 2))
    
    # lets run a fisher test to compare population mortality rates
    fisher.results = lapply(1:nrow(timeline), function(j) fisher.test(fisher.matrix[[j]]))
    names(fisher.results) = paste(timeline$Cause, timeline$Year)
    
    # create a table of annual odds ratios for each fishers test
    DT = data.table(oddsRatio = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$estimate)),
                    CI95LB = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$conf.int)[1]),
                    CI95UB = sapply(1:nrow(timeline), function(j) as.numeric(fisher.results[[j]]$conf.int)[2]))
    
    # add more problem data to DT
    DT[, Ratio := paste0("[", comb$Gender1[i], " ", comb$Ethnicity1[i], " Mortality] / [", comb$Gender2[i], " ", comb$Ethnicity2[i], " Mortality]")]
    DT[, Cause := timeline$Cause]
    DT[, Year := timeline$Year]
    
    # open up a graphics window
    # windows()
    
    # lets plot the odds ratio
    DT.plot = ggplot(DT, aes(x = Year, y = oddsRatio, color = Cause, fill = Cause, group = Cause)) + 
      geom_point(size = 3, na.rm = TRUE) +
      geom_line(size = 1, na.rm = TRUE) +
      geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "black") + 
      geom_ribbon(aes(ymin = CI95LB, ymax = CI95UB), alpha = 1/3, color = NA) + 
      # scale_color_manual(values = color.set) +
      # scale_y_continuous(labels = dollar) + 
      ggtitle("Multiple Causes of Death in the USA") + 
      labs(x = "Year", y = DT$Ratio[1], color = "", fill = "") + 
      theme_bw(base_size = 30) +
      theme(legend.position = "top", 
            legend.key.size = unit(.25, "in"),
            plot.title = element_text(hjust = 0.5),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  }
  
  # combine the test, data, and plot into a list as an output
  output = list("tests" = fisher.results,
                "ratio.data" = DT,
                "ratio.plot" = DT.plot)
  
  # clean up RAM
  gc()
  
  # export progress information
  sink(myfile, append = TRUE)
  cat(paste("task", i, "of", tasks, "finished at", Sys.time(), "\n"))
  sink()
  
  return(output)
}

# write out end time to log file
sink(myfile, append = TRUE)
cat(paste(tasks, "tasks finished at", Sys.time(), "\n"))
sink()

# end the cluster if it was set up
if(workers > 1)
{
  stopCluster(cl)
}

windows()
fishers[[1]]$ratio.plot

# lets test the following claims made by the alzheimer's association:
# (1) older African-Americans are about two times more likely than older whites to have Alzheimer's disease and other dementias.
fishers[[22]]$ratio.plot

# (2) older Hispanics are about one and one-half times more likely than older whites to have Alzheimer's disease and other dementias.
fishers[[35]]$ratio.plot

# (3) conditions such as high blood pressure and diabetes, both of which are risk factors for Alzheimer's and other dementias, are more common in African-Americans and Hispanics than in whites
fishers[[22]]$ratio.plot
fishers[[35]]$ratio.plot

# extract all of the alzheimers odds ratios
alz.odds = rbindlist(lapply(1:length(fishers), function(i) 
{
  # extract the ratio data for comparison i
  output = data.table(fishers[[i]]$ratio.data)
  
  # extract the alzheimers subset of output
  output = output[Cause == "Alzheimer's disease"]
  
  return(output)
}))

# compute the distance from 1 for each ratio
alz.odds[, distance := oddsRatio - 1]

# lets pick a span (in terms of years) to exponentially average distance (this gives more recent years more weight)
exp.span = 4

# exponentially average the distance together by Ratio
alz.odds.exp = data.table(alz.odds[, .(distance = tail(EMA(distance, n = 1,  ratio = 2 / (exp.span + 1)), 1)), 
                                   by = .(Ratio)])

# give alz.odds.exp an id number
alz.odds.exp[, id := 1:nrow(alz.odds.exp)]

# order alz.odds.exp by distance
alz.odds.exp = alz.odds.exp[order(abs(distance), decreasing = TRUE)]

# assign a ranking to each ratio
alz.odds.exp[, rank := 1:nrow(alz.odds.exp)]

# identify the target odds ratios that show which demographic is most likely to develop alziheimers
alz.odds.target = data.table(alz.odds.exp[rank %in% c(4, 8, 11, 15, 18, 20, 24, 25, 29, 32)])

}

# -----------------------------------------------------------------------------------
# ---- Study VA Dementia ------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# open up a graphics window
windows()

# lets compare the projected VA dementia prevelence with the known citizen dementia prevelence
VA.plot = ggplot(va, aes(x = Year, y = MiddleEstimate / Population, color = Group, fill = Group, group = Group)) + 
  geom_point(size = 3, na.rm = TRUE) +
  geom_line(size = 1, na.rm = TRUE) +
  geom_ribbon(aes(ymin = LowEstimate / Population, 
                  ymax = HighEstimate / Population), alpha = 1/3, color = NA) + 
  scale_y_continuous(labels = percent) + 
  ggtitle("Veteran Affairs Dementia Prevelence") + 
  labs(x = "Year", y = "Population with Dementia", color = "", fill = "") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "top", 
        legend.key.size = unit(.25, "in"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))

VA.plot

}

#####################################################################################
# the next section is useful for seeing the difference between two population tests #
# but it doesn't need to be run becuase the results above are the final results #####
#####################################################################################

# -----------------------------------------------------------------------------------
# ---- Proportion Test vs. Fisher Test ----------------------------------------------
# -----------------------------------------------------------------------------------

{

# lets extract older African-Americans and older whites who have died
dat.test = data.table(dat[Race %in% c("Black or African American", "White") & 
                          Age %in% c("65-74 years", "75-84 years", "85+ years"),
                          .(Deaths = sum(Deaths),
                            Population = sum(Population)),
                          by = .(Cause, Race, Year)])

# build a matrix showing the death numbers and population numbers
comb = expand.grid(Cause = unique(dat.test$Cause), Year = unique(dat.test$Year))
prop.matrix = lapply(1:nrow(comb), function(i)
  matrix(c(unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "Black or African American", .(Deaths)]), 
           unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "White", .(Deaths)]),
           unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "Black or African American", .(Population)]), 
           unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "White", .(Population)])), 
         ncol = 2))

# lets run a proportion test to test (1)
prop.results = lapply(1:nrow(comb), function(i) prop.test(prop.matrix[[i]], correct = FALSE))
names(prop.results) = paste(comb$Cause, comb$Year)

# we can estimate the proportion values in prop.results for any test k
k = 1
prop.stats = c(prop.matrix[[k]][1,1] / prop.matrix[[k]][1,2], prop.matrix[[k]][2,1] / prop.matrix[[k]][2,2])
prop.stats

# build a matrix showing the death numbers and living numbers
fisher.matrix = lapply(1:nrow(comb), function(i)
  matrix(c(unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "Black or African American", .(Deaths)]),
           unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "Black or African American", .(Population)]) - unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "Black or African American", .(Deaths)]),
           unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "White", .(Deaths)]),
           unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "White", .(Population)]) - unlist(dat.test[Cause == comb$Cause[i] & Year == comb$Year[i] & Race == "White", .(Deaths)])), 
         ncol = 2))

# lets run a fisher test to test (1)
fisher.results = lapply(1:nrow(comb), function(i) fisher.test(fisher.matrix[[i]]))
names(fisher.results) = paste(comb$Cause, comb$Year)

# we can estimate the odds ratios in fisher.results for any test k
k = 1
odds.ratio = (fisher.matrix[[k]][1,1] / fisher.matrix[[k]][1,2]) / (fisher.matrix[[k]][2,1] / fisher.matrix[[k]][2,2])
odds.ratio

# if we subtract the odds ratio by 1 then we get the percent difference between these proportions
k = 1
as.numeric(fisher.results[[k]]$estimate) - 1

# create a table of annual odds ratios for each fishers test
DT1 = data.table(oddsRatio = sapply(1:nrow(comb), function(i) as.numeric(fisher.results[[i]]$estimate)),
                 CI95LB = sapply(1:nrow(comb), function(i) as.numeric(fisher.results[[i]]$conf.int)[1]),
                 CI95UB = sapply(1:nrow(comb), function(i) as.numeric(fisher.results[[i]]$conf.int)[2]))

# add more problem data to DT1
DT1[, Ratio := "[Black or African American Mortality] / [White Mortality]"]
DT1[, Cause := comb$Cause]
DT1[, Year := comb$Year]

# open up a graphics window
windows()

# lets plot the odds ratio
DT1.plot = ggplot(DT1, aes(x = Year, y = oddsRatio, color = Cause, fill = Cause, group = Cause)) + 
  geom_point(size = 3, na.rm = TRUE) +
  geom_line(size = 1, na.rm = TRUE) +
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "black") + 
  geom_ribbon(aes(ymin = CI95LB, ymax = CI95UB), alpha = 1/3, color = NA) + 
  # scale_color_manual(values = color.set) +
  # scale_y_continuous(labels = dollar) + 
  ggtitle("Multiple Causes of Death in the USA") + 
  labs(x = "Year", y = "[Black or African American Mortality] / [White Mortality]", color = "", fill = "") + 
  theme_bw(base_size = 20) +
  theme(legend.position = "top", 
        legend.key.size = unit(.25, "in"),
        plot.title = element_text(hjust = 0.5),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))

DT1.plot

}











