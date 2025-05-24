# Get axis attributes
GetAxisAttributes <- function(x,
                              y,
                              bin.size) {
  # Bin x and y axis values
  x.bin <- discretise(x = x,
                      by = bin.size)
  y.bin <- discretise(x = y,
                      by = bin.size)
  
  # Get levels of x and y bins
  x.bin.levels <- levels(x = x.bin)
  y.bin.levels <- levels(x = y.bin)
  
  # Get x and y bin values
  x.bin.values <- as.numeric(x = unlist(x = stringr::str_extract_all(string = x.bin.levels,
                                                                     pattern = "-?[0-9.]+",
                                                                     simplify = FALSE)))
  y.bin.values <- as.numeric(x = unlist(x = stringr::str_extract_all(string = y.bin.levels,
                                                                     pattern = "-?[0-9.]+",
                                                                     simplify = FALSE)))
  
  # Get min and max of x and y axis values
  min.x <- min(x.bin.values)
  min.y <- min(y.bin.values)
  max.x <- max(x.bin.values)
  max.y <- max(y.bin.values)
  
  # Get range of x and y axis values
  range.x <- abs(x = max.x - min.x) # get range of x axis
  range.y <- abs(x = max.y - min.y) # get range of y axis
  
  # Get ratio of x and y ranges
  range.ratio <- range.x / range.y
  
  # Get maximum range value from both axes
  range.max <- max(range.x, range.y)
  
  # Return data
  return(list(min.x = min.x,
              min.y = min.y,
              max.x = max.x,
              max.y = max.y,
              range.x = range.x,
              range.y = range.y,
              range.ratio = range.ratio,
              range.max = range.max,
              x.bin.levels = x.bin.levels,
              y.bin.levels = y.bin.levels))
}

# Get bin attributes
GetBinAttributes <- function(row.data,
                             clusters.of.interest,
                             cluster.labels,
                             cluster.colors,
                             n.pad) {
  row.data <- unlist(x = row.data)
  if (any(row.data > 0)) {
    # Order clusters by their percentage of cells in each bin
    clusters.order <- order(row.data, decreasing = TRUE)
    row.data <- row.data[clusters.order]
    clusters.of.interest <- clusters.of.interest[clusters.order]
    
    # Keep only clusters represented by at least 1 cell in the corresponding bin
    clusters.to.keep <- row.data > 0
    clusters.order <- clusters.order[clusters.to.keep]
    row.data <- row.data[clusters.to.keep]
    clusters.of.interest <- clusters.of.interest[clusters.to.keep]
    cluster.labels <- cluster.labels[clusters.of.interest]
    cluster.colors <- cluster.colors[clusters.of.interest]
    
    # Get ID and label attributes for corresponding bin
    bin.attributes <- c(clusters.of.interest[1],
                        cluster.labels[1])
    if (length(x = row.data) > 1) {
      for (i in 2:(length(x = row.data))) {
        bin.attributes <- paste(bin.attributes,
                                c(clusters.of.interest[i],
                                  cluster.labels[i]),
                                sep = {
                                  if (row.data[i-1] > row.data[i]) {
                                    ">"
                                  } else {
                                    "="
                                  }
                                })
      }
    }
    
    # Get color attribute for corresponding bin
    if (length(x = row.data) > 1) { # decreases time for computing as this operation is costly
      # Create sRGB colors
      cluster.colors.rgb <- colorspace::sRGB(t(x = col2rgb(col = cluster.colors)))
      # Define iterative color name
      color1 <- cluster.colors.rgb[1,]
      # Iteratively mix colors
      for (i in 2:length(x = row.data)) {
        color2 <- cluster.colors.rgb[i,]
        alpha <- row.data[i] / sum(row.data[1:i])
        color1 <- colorspace::mixcolor(color1 = color1,
                                       color2 = color2,
                                       alpha = alpha)
      }
      # Get hex color code
      bin.color <- rgb(red = color1@coords[1,"R"],
                       green = color1@coords[1,"G"],
                       blue = color1@coords[1,"B"],
                       maxColorValue = 255)
      # Add color to bin attributes
      bin.attributes <- c(bin.attributes,
                          bin.color)
    } else {
      # Add color to bin attributes
      bin.attributes <- c(bin.attributes,
                          cluster.colors)
    }
    
    # Right pad bin order for cluster breaks ordering
    bin.order <- stringr::str_pad(string = paste0(clusters.order,
                                                  collapse = ""),
                                  width = n.pad,
                                  side = "right",
                                  pad = 0)
    
    # Add decimal value to order bins with one major cluster first
    bin.order <- paste(bin.order,
                       stringr::str_count(string = bin.attributes[1],
                                          pattern = "="),
                       sep = ".")
    
    # Add bin order to attributes
    bin.attributes <- c(bin.attributes,
                        bin.order)
    
    # Return bin attributes
    return(as.list(x = bin.attributes))
  } else {
    # Return NA
    return(as.list(rep(x = NA, times = 4)))
  }
}

# Mix bin cluster colors with white if cluster is imputed (decrease alpha values with increasing imputation iteration)
MixBinColors <- function(row.data) {
  # Get alpha value
  alpha <- unlist(row.data[,"imputed_bin_cluster_color_alpha"])
  
  # Get bin cluster color
  bin.color <- unlist(row.data[,"imputed_bin_cluster_color"])
  
  if (alpha < 1) {
    # Create sRGB color for bin cluster
    bin.color.rgb <- colorspace::sRGB(t(x = col2rgb(col = bin.color)))
    
    # Create sRGB color for white
    white.rgb <- colorspace::sRGB(R = 255, G = 255, B = 255)
    
    # Mix color with white at 1 - alpha color contribution
    bin.color.rgb <- colorspace::mixcolor(color1 = bin.color.rgb,
                                          color2 = white.rgb,
                                          alpha = 1 - alpha)
    # Get hex color code
    bin.color <- rgb(red = bin.color.rgb@coords[1,"R"],
                     green = bin.color.rgb@coords[1,"G"],
                     blue = bin.color.rgb@coords[1,"B"],
                     maxColorValue = 255)
  }
  
  # Return bin color
  return(bin.color)
}
