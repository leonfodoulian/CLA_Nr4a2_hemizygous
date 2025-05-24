# Expand cluster attributes with labels and alpha for each label
ExpandClusterAttributes <- function(
    cluster.attributes,
    labels,
    alpha,
    prefix = TRUE
){
  # Check if 'labels' and 'alpha' have the same length
  if (length(x = labels) != length(x = alpha)) {
    stop(paste("'labels' and 'alpha' should have the same length", sep = ""))
  }
  
  # Get number of replications
  nrep <- length(x = labels)
  
  # Add labels to cluster breaks
  if (prefix) {
    cluster.breaks <- paste(rep(x = cluster.attributes$cluster.breaks,
                                each = nrep),
                            labels)
  } else {
    cluster.breaks <- paste(labels,
                            rep(x = cluster.attributes$cluster.breaks,
                                each = nrep)) 
  }
  
  # Expand cluster colors
  cluster.colors <- setNames(object = rep(x = cluster.attributes$cluster.colors,
                                          each = nrep),
                             nm = cluster.breaks)
  
  # Create alpha for each label x cluster level
  cluster.alpha <- setNames(object = rep(x = alpha,
                                         length.out = length(x = cluster.breaks)),
                            nm = cluster.breaks)
  
  # Blend cluster colors with alpha for each label x cluster level
  cluster.blended.colors <- setNames(object = BlendColorAlpha(color = cluster.colors,
                                                              alpha = cluster.alpha),
                                     nm = cluster.breaks)
  
  # Return list of attributes
  return(list(cluster.breaks = cluster.breaks,
              cluster.colors = cluster.colors,
              cluster.blended.colors = cluster.blended.colors,
              cluster.alpha = cluster.alpha))
}
