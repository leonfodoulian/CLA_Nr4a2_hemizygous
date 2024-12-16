# Cluster RNAscope data
ClusterRNAscopeData <- function(
    rnascope.data,
    channels,
    n.clusters,
    log.transform = TRUE,
    pseudocount = 1,
    use.fastcluster = TRUE
){
  # Prepare data for clustering
  rownames(x = rnascope.data) <- rnascope.data$cell
  counts.data <- rnascope.data[,channels] # keep counts columns
  cells.to.keep <- rowSums(x = counts.data) > 0 # remove cells with zero counts
  counts.data <- counts.data[cells.to.keep,] # remove cells with zero counts
  
  if (isTRUE(x = log.transform)) {
    # Log transform the data
    counts.data <- log(x = counts.data + pseudocount)
  }
  
  # Scale data
  counts.data.scaled <- scale(x = counts.data,
                              center = TRUE,
                              scale = TRUE)
  colnames(x = counts.data.scaled) <- paste(colnames(x = counts.data.scaled),
                                            "_scaled",
                                            sep = "")
  
  # Compute euclidean distances between pairs of cells
  dist.data <- dist(x = counts.data.scaled,
                    method = "euclidean",
                    diag = FALSE,
                    upper = FALSE)
  
  # Compute ward.D2 hierarchical clustering
  if (isTRUE(x = use.fastcluster)) {
    hc.ward.D2 <- fastcluster::hclust(d = dist.data,
                                      method = "ward.D2",
                                      members = NULL)
  } else {
    hc.ward.D2 <- stats::hclust(d = dist.data,
                                method = "ward.D2")
  }
  
  # Cut tree at defined numbers of clusters (n.clusters)
  names(x = n.clusters) <- n.clusters
  clusters.ward.D2.l <- lapply(X = n.clusters,
                               FUN = function(k,
                                              hc.ward.D2) {
                                 # Cut tree at k clusters
                                 clusters.ward.D2 <- as.data.frame(x = cutree(tree = hc.ward.D2,
                                                                              k = k))
                                 # Rename column with number of clusters
                                 colnames(x = clusters.ward.D2) <- paste("k", k, "_clusters",
                                                                         sep = "")
                                 # Add cell names to data
                                 clusters.ward.D2$cell <- rownames(x = clusters.ward.D2)
                                 # Return data
                                 return(clusters.ward.D2)
                               },
                               hc.ward.D2 = hc.ward.D2)
  
  # Merge all clustering results
  clusters.ward.D2 <- Reduce(f = function(x, y) { merge(x = x,
                                                        y = y,
                                                        by = "cell",
                                                        sort = FALSE)},
                             x = clusters.ward.D2.l)
  
  # Add cell names to scaled counts data
  counts.data.scaled <- as.data.frame(x = counts.data.scaled)
  counts.data.scaled$cell <- rownames(x = counts.data.scaled)
  
  # Merge RNAscope data with scaled counts data and clustering data
  full.data <- Reduce(f = function(x, y) { merge(x = x,
                                                 y = y,
                                                 by = "cell",
                                                 all = TRUE,
                                                 sort = FALSE)},
                      x = list(rnascope.data,
                               counts.data.scaled,
                               clusters.ward.D2))
  
  # Re-arrange data as in rnascope.data
  full.data <- full.data[match(x = rnascope.data$cell, table = full.data$cell),]
  rownames(x = full.data) <- full.data$cell
  
  # Add tag to recognize clustered cells
  full.data$is_clustered <- ifelse(test = full.data$cell %in% clusters.ward.D2$cell,
                                   yes = TRUE,
                                   no = FALSE)
  
  # Return list of data
  return(list(full.data = full.data,
              # dist.data = dist.data,
              hc.ward.D2 = hc.ward.D2,
              clusters.ward.D2 = clusters.ward.D2))
}
