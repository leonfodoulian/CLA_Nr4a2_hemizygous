# Compute VI between a pair of clustering results
ClusteringPairVI <- function(
    meta.data,
    clustering.pair
){
  # Calculate VI between clustering pairs
  vi <- igraph::compare(
    comm1 = meta.data[[clustering.pair$comm1]],
    comm2 = meta.data[[clustering.pair$comm2]],
    method = "vi"
  )
  
  # Get clustering pairs
  pair1 <- clustering.pair$comm1
  pair2 <- clustering.pair$comm2
  pair <- paste(pair1,
                pair2,
                sep = " - ")
  
  # Get number of clusters in each clustering result
  n.clust.pair1 <- length(x = unique(x = meta.data[[clustering.pair$comm1]]))
  n.clust.pair2 <- length(x = unique(x = meta.data[[clustering.pair$comm2]]))
  
  # Get change in number of clusters between clustering pairs
  n.clust.diff <- n.clust.pair2 - n.clust.pair1
  
  # Calculate normalized VI
  normvi <- vi / (abs(x = n.clust.diff) + 1)
  
  # Return data
  return(data.frame(vi = vi,
                    normvi = normvi,
                    pair = pair,
                    pair1 = pair1,
                    pair2 = pair2,
                    n.clust.pair1 = n.clust.pair1,
                    n.clust.pair2 = n.clust.pair2,
                    n.clust.diff = n.clust.diff))
}
