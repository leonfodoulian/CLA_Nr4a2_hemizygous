# Find statistically significant differences in cell type proportions and plot data
CellTypeProportions <- function(
    x = NULL,
    clusters = NULL,
    sample = NULL,
    group = NULL,
    trend = FALSE,
    robust = TRUE,
    transform = "logit",
    return.plot = TRUE,
    markdown = FALSE,
    reference.group,
    pseudocount = 0,
    cluster.colors,
    bar.border.stroke = 1,
    label.size = 15,
    label.position = c("top", "bottom"),
    point.size = 2,
    sd.line.size = 1,
    y.axis.title = "log<sub>2</sub> ratio",
    font.family = "Arial",
    font.face = "plain",
    axis.text.size = 12,
    axis.title.size = 12,
    axis.ticks.length = 3,
    axis.line.size = rel(x = 1)
){
  # Find statistically significant differences in cell type proportions
  cell.type.prop <- speckle::propeller(
    x = x,
    clusters = clusters,
    sample = sample,
    group = group,
    trend = trend,
    robust = robust,
    transform = transform
  )
  
  # log2 transform ratio
  cell.type.prop$log2PropRatio <- log2(x = cell.type.prop$PropRatio)
  
  # Add * label for significant proportion differences
  cell.type.prop$label <- ifelse(test = cell.type.prop$FDR < 0.05,
                                 yes = "*",
                                 no = "")
  
  if (isTRUE(x = return.plot)) {
    # Verify pseudocount value
    if (!is.numeric(x = pseudocount) || pseudocount < 0) {
      pseudocount <- 0
    }
    
    # Verify label position argument
    label.position <- match.arg(arg = label.position,
                                choices = c("top", "bottom"))
    
    # Create colors if cluster.colors is missing
    if (missing(x = cluster.colors)) {
      cluster.colors <- viridis::viridis(n = length(x = unique(x = clusters)))
    }
    
    if (!missing(x = reference.group)) {
      # Get sample x group pairs
      group <- setNames(object = as.character(x = group),
                        nm = as.character(x = sample))
      group <- group[!duplicated(x = names(x = group))]
      
      # Prepare data for points plotting
      sample.prop <- table(sample = sample,
                           clusters = clusters) # compute ncells per cluster and sample
      sample.prop <- sample.prop + pseudocount # add pseudocount (if > 0, samples with 0 cells will be plotted)
      sample.prop <- sample.prop / rowSums(x = sample.prop) # normalize by sample
      sample.prop <- as.data.frame(x = sample.prop,
                                   responseName = "proportion",
                                   stringsAsFactors = FALSE) # transform to data frame
      if (sum(sample.prop$proportion == 0) > 0) {
        message("Removing ", sum(sample.prop$proportion == 0), " data point(s): to plot the data, set pseudocount > 0")
      }
      sample.prop <- sample.prop[sample.prop$proportion != 0,] # remove pairs with 0 cells
      sample.prop$group <- group[sample.prop$sample] # add group
      
      # Compute log2 ratio of proportions relative to reference group
      sample.prop$reference_prop_mean <- cell.type.prop[sample.prop$clusters, paste0("PropMean.", reference.group)]
      sample.prop$log2_prop_ratio <- log2(x = sample.prop$proportion / sample.prop$reference_prop_mean)
      
      # Compute mean +/- sd of log2 proportion ratios for each cluster
      mean.sd <- unname(
        obj = unlist(
          x = lapply(
            X = split(x = sample.prop[sample.prop$group != reference.group,]$log2_prop_ratio,
                      f = sample.prop[sample.prop$group != reference.group,]$clusters),
            FUN = function(x) {
              return(c(mean(x = x) - sd(x = x),
                       mean(x = x) + sd(x = x)))
            }
          )
        )
      )
      
      # Get log2 proportion ratios for test group
      log2.prop.ratio <- sample.prop[sample.prop$group != reference.group,]$log2_prop_ratio
      
      # Compute min and max y
      min.y <- 0.5 * floor(x = min(log2.prop.ratio, mean.sd) / 0.5)
      max.y <- 0.5 * ceiling(x = max(log2.prop.ratio, mean.sd) / 0.5)
      
      # Get label y value
      # label.y <- mean(x = c(min.y, min(sample.prop$log2_prop_ratio)))
    } else {
      # Set sample.prop to NULL
      sample.prop <- NULL
      
      # Compute min and max y
      min.y <-  0.5 * floor(x = min(cell.type.prop$log2PropRatio) / 0.5)
      max.y <-  0.5 * ceiling(x = max(cell.type.prop$log2PropRatio) / 0.5)
      
      # Get label y value
      # label.y <- mean(x = c(min.y, min(cell.type.prop$log2PropRatio)))
    }
    
    # Get label y value
    if (label.position == "top") {
      label.y <- max.y - 0.05
    } else {
      label.y <- min.y + 0.05
    }
    
    # Set axis text theme element
    if (isTRUE(x = markdown)) {
      theme_element_text <- function(...) { ggtext::element_markdown(...) }
    } else {
      theme_element_text <- function(...) { ggplot2::element_text(...) }
    }
    
    # Plot data
    prop.plot <- ggplot2::ggplot(data = cell.type.prop,
                                 mapping = ggplot2::aes(x = BaselineProp.clusters,
                                                        y = log2PropRatio,
                                                        color = BaselineProp.clusters)) +
      ggplot2::geom_col(size = bar.border.stroke,
                        fill = NA,
                        show.legend = FALSE) +
      ggplot2::geom_text(mapping = ggplot2::aes(label = label),
                         y = label.y,
                         size = label.size * (1/72 * 25.4),
                         colour = "black",
                         family = font.family,
                         fontface = "bold") +
      ggplot2::scale_y_continuous(limits = c(min.y, max.y),
                                  breaks = seq(from = min.y,
                                               to = max.y,
                                               by = 0.5),
                                  expand = c(0,0)) +
      ggplot2::scale_color_manual(values = cluster.colors) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = "cell types",
                    y = y.axis.title,
                    colour = "cell types") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = theme_element_text(size = axis.text.size,
                                                      angle = 45,
                                                      hjust = 1,
                                                      vjust = 1,
                                                      family = font.family,
                                                      face = font.face,
                                                      colour = cluster.colors),
                     axis.text.y = ggplot2::element_text(size = axis.text.size,
                                                         family = font.family,
                                                         face = font.face,
                                                         colour = "black"),
                     axis.title = ggtext::element_markdown(size = axis.title.size,
                                                           family = font.family,
                                                           face = font.face,
                                                           colour = "black"),
                     legend.text = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(size = axis.line.size,
                                                       colour = "black",
                                                       lineend = "square"),
                     axis.ticks = ggplot2::element_line(size = axis.line.size,
                                                        colour = "black",
                                                        lineend = "square"),
                     axis.ticks.length = ggplot2::unit(x = axis.ticks.length,
                                                       units = "mm"),
                     plot.title = ggplot2::element_blank(),
                     legend.position = "none",
                     strip.text = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     legend.background = ggplot2::element_blank())
    
    # Add points and sd range to plot
    if (!missing(x = reference.group)) {
      prop.plot <- prop.plot +
        ggplot2::stat_summary(data = sample.prop[sample.prop$group != reference.group,],
                              mapping = ggplot2::aes(x = clusters,
                                                     y = log2_prop_ratio,
                                                     colour = clusters),
                              geom = "segment",
                              size = sd.line.size,
                              fun.data = function(x) {
                                return(c(y = mean(x = x) - sd(x = x),
                                         yend = mean(x = x) + sd(x = x)))
                              }) +
        ggplot2::geom_point(data = sample.prop[sample.prop$group != reference.group,],
                            mapping = ggplot2::aes(x = clusters,
                                                   y = log2_prop_ratio,
                                                   colour = clusters),
                            size = point.size,
                            shape = 19,
                            stroke = 0)
    }
  } else {
    # Return NULL plot
    prop.plot <- NULL
  }
  
  # Return list of data and plot
  return(list(cell.type.prop = cell.type.prop,
              sample.prop = sample.prop,
              prop.plot = prop.plot))
}
