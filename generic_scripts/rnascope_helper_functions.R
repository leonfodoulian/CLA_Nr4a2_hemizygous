# Get attributes from file name
GetFileAttributes <- function(file,
                              attributes,
                              channel.attributes = 5:7, 
                              split.pattern = "_",
                              prefix = "",
                              suffix = ".xlsx",
                              path = ".",
                              subdirectory = "") {
  # Remove prefix and suffix from file name
  file.name <- gsub(pattern = paste(paste0("^", prefix, "/"),
                                    paste0(suffix, "$"),
                                    sep = "|"),
                    replacement = "",
                    x = file)
  
  # Get meta data from file name
  meta.data <- unlist(x = stringr::str_split(string = file.name,
                                             pattern = split.pattern))
  names(x = meta.data) <- attributes
  
  # Separate channel data from meta data
  channels <- meta.data[channel.attributes]
  meta.data <- meta.data[-channel.attributes]
  
  # Define output directory path for results
  path <- file.path(path,
                    paste(channels,
                          collapse = "_"),
                    subdirectory,
                    paste(meta.data,
                          collapse = "_"))
  
  # Return list of data
  return(list(file = file,
              meta.data = meta.data,
              channels = channels,
              path = path))
}

# Load and prepare RNAscope data
LoadRNAscopeData <- function(file,
                             objects.to.keep = "PathCellObject",
                             meta.data,
                             channels) {
  # Rename meta.data if not named
  if (is.null(x = names(x = meta.data))) {
    names(x = meta.data) <- paste("meta_data",
                                  1:length(x = meta.data),
                                  sep = "_")
  }
  # Import data
  rnascope.data <- readxl::read_xlsx(path = file,
                                     na = "NaN")
  # Convert data to data frame
  rnascope.data <- as.data.frame(x = rnascope.data)
  # Keep only rows corresponding to a given cell
  rnascope.data <- rnascope.data[rnascope.data$Name %in% objects.to.keep,]
  # Add image information to data
  rnascope.data[names(x = meta.data)] <- as.list(x = meta.data)
  # Create cell column
  rnascope.data$cell <- paste(paste(meta.data,
                                    collapse = "_"),
                              "cell",
                              1:nrow(x = rnascope.data),
                              sep = "_")
  # Get names of counts columns
  counts.columns <- grep(pattern = "\\:\\sNum spots estimated$",
                         x = colnames(x = rnascope.data),
                         value = TRUE)
  # Subset dataset to keep only relevant columns for downstream analysis
  rnascope.data <- rnascope.data[,c(names(x = meta.data), "cell", counts.columns, "Centroid X µm", "Centroid Y µm", "Nucleus: Area")]
  # Rename counts columns with gene name only and remove some characters from column names
  colnames(x = rnascope.data) <- gsub(pattern = "(Subcellular\\:\\s)|(\\:\\sNum\\sspots\\sestimated)",
                                      replacement = "",
                                      x = colnames(x = rnascope.data))
  channel.columns <- grep(pattern = "Channel",
                          x = colnames(x = rnascope.data),
                          value = TRUE)
  colnames(x = rnascope.data) <- c(names(x = meta.data), "cell", unname(obj = channels[channel.columns]), "centroid_x", "centroid_y", "nucleus_area")
  # Return data
  return(as.data.frame(x = rnascope.data))
}

# Cluster RNAscope data
ClusterRNAscopeData <- function(rnascope.data,
                                channels,
                                n.clusters,
                                log.transform = TRUE,
                                pseudocount = 1,
                                use.fastcluster = TRUE) {
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

# Merge clusters using NMI score computed on complete vs individual clustering of cells from related clusters
MergeClustersNMI <- function(rnascope.data,
                             hc,
                             hcdata,
                             k,
                             channels,
                             log.transform = TRUE,
                             pseudocount = 1,
                             use.fastcluster = TRUE,
                             nmi.thresh) {
  # Get dendrogram height for chosen k
  hpk <- sort(x = hc$height,
              decreasing = TRUE)[k]
  
  # Create segments and labels data from clustering result (used to identify related clusters)
  hc.data <- dendro_data_k(hc = hc,
                           hcdata = hcdata,
                           k = k)
  
  # Depth of merging process
  merging.depth <- 1
  
  # Prepare an empty list to append data to return
  out.data <- list()
  
  repeat {
    # Indicate progress of merging process
    cat.verbose(x = paste("Merging depth", merging.depth, "started."),
                verbose = TRUE)
    
    # Get number of clusters from hc.data
    clusters <- sort(x = unique(x = hc.data$labels$clust))
    
    # Check if number of clusters is > 1
    if (length(x = clusters) == 1) {
      cat.verbose(x = "Only 1 cluster left: cannot proceed.",
                  verbose = TRUE)
      break()
    }
    
    # Subset hc.data$segments at height for chosen k
    hc.data.hpk <- hc.data$segments[hc.data$segments$y > hpk,]
    hc.data.hpk$order <- 1:nrow(x = hc.data.hpk)
    
    # Get names of related clusters and prepare hc.data.hpk for plotting
    related.clusters.l <- unname(
      obj = lapply(
        X = split(x = hc.data.hpk,
                  f = hc.data.hpk$y),
        FUN = function(hc.data.hpk.y) {
          # Get unique clusters
          hpk.clust <- sort(x = unique(x = hc.data.hpk.y[hc.data.hpk.y$line == 0,]$clust))
          
          # Get names of related clusters
          if (length(x = hpk.clust) == 0) {
            related.clusters <- NULL
          } else {
            related.clusters <- paste(hpk.clust,
                                      collapse = "-") 
          }
          
          # Prepare hc.data.hpk for plotting
          if (length(x = hpk.clust) > 1) {
            hc.data.hpk.y$child_clust <- ifelse(test = hc.data.hpk.y$clust %in% hpk.clust[1],
                                                yes = "child cluster 1",
                                                no = "child cluster 2")
          } else {
            hc.data.hpk.y$child_clust <- NA
          }
          
          return(list(hc.data.hpk = hc.data.hpk.y,
                      related.clusters = related.clusters))
        }))
    
    # Get names of related clusters
    related.clusters <- unique(x = unlist(x = lapply(X = related.clusters.l,
                                                     FUN = '[[',
                                                     "related.clusters")))
    
    # Bind hc.data.hpk by row
    hc.data.hpk <- dplyr::bind_rows(lapply(X = related.clusters.l,
                                           FUN = '[[',
                                           "hc.data.hpk"))
    
    # Replace NA with child cluster attribute for each cluster if one exists and reorder data (for plotting)
    hc.data.hpk <- dplyr::bind_rows(lapply(
      X = split(x = hc.data.hpk,
                f = hc.data.hpk$clust),
      FUN = function(hc.data.hpk.clust) {
        if (!all(is.na(x = hc.data.hpk.clust$child_clust))) {
          hc.data.hpk.clust$child_clust <- as.character(x = na.omit(object = unique(x = hc.data.hpk.clust$child_clust)))
        }
        return(hc.data.hpk.clust)
      }))
    hc.data.hpk <- hc.data.hpk[order(hc.data.hpk$order),]
    
    # Get names of related clusters to test for robustness of clustering
    clusters.to.test <- grep(pattern = "-",
                             x = related.clusters,
                             value = TRUE)
    names(x = clusters.to.test) <- clusters.to.test
    
    # Get NMI values between complete vs individual clustering of cells from related clusters
    nmi.data.l <- lapply(
      X = strsplit(x = clusters.to.test,
                   split = "-"),
      FUN = function(clust.comb,
                     rnascope.data,
                     channels,
                     log.transform,
                     pseudocount,
                     use.fastcluster) {
        # Subset RNAscope data for cluster combination
        rnascope.data.sub <- rnascope.data[rnascope.data$merged_clusters %in% unlist(x = clust.comb),]
        
        # Individually cluster cells from related clusters
        clustered.data.sub <- ClusterRNAscopeData(rnascope.data = rnascope.data.sub,
                                                  channels = channels,
                                                  n.clusters = 2,
                                                  log.transform = log.transform,
                                                  pseudocount = pseudocount,
                                                  use.fastcluster = use.fastcluster)
        
        # Get NMI between the two clusters
        clust.comb.nmi <- igraph::compare(comm1 = as.numeric(x = as.factor(x = clustered.data.sub$full.data$merged_clusters)),
                                          comm2 = clustered.data.sub$full.data$k2_clusters,
                                          method = "nmi")
        
        # Return data
        return(clust.comb.nmi = data.frame(nmi = clust.comb.nmi,
                                           cluster_pair = c("cluster_pair_1", "cluster_pair_2"),
                                           old_clusters = as.character(x = clust.comb),
                                           stringsAsFactors = FALSE))
      },
      rnascope.data = rnascope.data,
      channels = channels,
      log.transform = log.transform,
      pseudocount = pseudocount,
      use.fastcluster = use.fastcluster)
    
    # Bind data by rows
    nmi.data <- dplyr::bind_rows(nmi.data.l,
                                 .id = "related_clusters")
    
    # Replave "-" with "" (characters are problematic for igraph::compare())
    nmi.data$related_clusters <- gsub(pattern = "-",
                                      replacement = ".",
                                      x = nmi.data$related_clusters)
    
    # Determine if related clusters should be merged
    nmi.data$cluster_to_merge <- ifelse(test = nmi.data$nmi < nmi.thresh,
                                        yes = TRUE,
                                        no = FALSE)
    
    # Get names of clusters not tested
    other.clusters <- setdiff(x = related.clusters,
                              y = c(clusters.to.test, nmi.data$old_clusters))
    
    # Add data related to clusters not tested to nmi data
    if (length(x = other.clusters) > 0) {
      nmi.data <- rbind(nmi.data,
                        data.frame(related_clusters = other.clusters,
                                   nmi = NA,
                                   cluster_pair = "no_pair",
                                   old_clusters = other.clusters,
                                   cluster_to_merge = FALSE))
    }
    
    # Get new (merged) cluster identities
    nmi.data$merged_clusters <- ifelse(test = nmi.data$cluster_to_merge,
                                       yes = nmi.data$related_clusters,
                                       no = nmi.data$old_clusters)
    
    # Add nmi.data and hc.data.hpk to the list to return
    out.data[[paste("merging_depth", merging.depth, sep = "_")]] <- list(nmi.data = nmi.data,
                                                                         hc.data.hpk = hc.data.hpk)
    
    # Check if there are clusters to merge
    if (!any(nmi.data$cluster_to_merge)) {
      cat.verbose(x = "No clusters to merge: halting merging process.",
                  verbose = TRUE)
      break()
    } else {
      silent <- lapply(
        X = split(x = nmi.data,
                  f = nmi.data$related_clusters),
        FUN = function(rc.nmi.data) {
          if (all(rc.nmi.data$cluster_to_merge)) {
            cat.verbose(x = paste("   Merging clusters ",
                                  paste(rc.nmi.data$old_clusters,
                                        collapse = " and "),
                                  " in cluster ",
                                  unique(x = rc.nmi.data$related_clusters),
                                  " (nmi=",
                                  signif(x = unique(x = rc.nmi.data$nmi),
                                         digits = 2),
                                  "; <",
                                  nmi.thresh,
                                  ")",
                                  sep = ""))
          }
        })
    }
    
    # Prepare new (merged) cluster identities to append to data
    merged_clusters <- nmi.data$merged_clusters
    names(x = merged_clusters) <- nmi.data$old_clusters
    merged_clusters <- merged_clusters[order(names(x = merged_clusters))] # sort for easier manipulation
    
    # Replace merged cluster identities in RNAscope data
    rnascope.data$merged_clusters <- merged_clusters[as.character(x = rnascope.data$merged_clusters)]
    
    # Replace cluster identities in hc data
    hc.data$labels$clust <- merged_clusters[as.character(x = hc.data$labels$clust)] # labels data
    hc.data$segments$clust <- c("0" = "0", merged_clusters)[as.character(x = hc.data$segments$clust)] # segments data
    
    # Change identity of parent branch of a pair of clusters to merge
    clusters.to.merge <- unique(x = nmi.data$merged_clusters[nmi.data$cluster_to_merge])
    for (i in clusters.to.merge) {
      # Get y of a pair of clusters to merge
      y.ctm <- max(hc.data$segments$y[hc.data$segments$clust == i])
      # Get y of parent branch of pair of clusters to merge
      y.pb <- hc.data$segments$y[hc.data$segments$line == 1 & hc.data$segments$yend == y.ctm]
      # Get xend of parent branch of pair of clusters to merge
      xend.pb <- hc.data$segments$xend[hc.data$segments$line == 1 & hc.data$segments$yend == y.ctm]
      # Change cluster identity of parent branch
      hc.data$segments$clust[hc.data$segments$line == 1 & hc.data$segments$y == y.pb & hc.data$segments$xend == xend.pb] <- i
    }
    hc.data$segments$line[hc.data$segments$clust != 0] <- 0
    
    # Indicate progress of merging process
    cat.verbose(x = "Merging completed.",
                verbose = TRUE)
    
    # Add rnascope.data and hc.data to the list to return
    out.data[[paste("merging_depth", merging.depth, sep = "_")]] <- append(
      x = out.data[[paste("merging_depth", merging.depth, sep = "_")]],
      values = list(rnascope.data = rnascope.data,
                    hc.data = hc.data))
    
    # Increase depth of merging process
    merging.depth <- merging.depth + 1
  }
  
  # Add parameter indicating that clusters got merged
  out.data$clusters_merged <- ifelse(test = merging.depth > 1,
                                     yes = TRUE,
                                     no = FALSE)
  
  # Return list of datas
  return(out.data)
}

# Plot and save violin plots
PlotStackedViolinPlots <- function(rnascope.data,
                                   clust.res.col,
                                   channels,
                                   cluster.colors,
                                   cluster.labels,
                                   plot.title,
                                   y.scale.trans = "log10",
                                   path) {
  # Set clusters as character
  rnascope.data[[clust.res.col]] <- as.character(x = rnascope.data[[clust.res.col]])
  
  # Get names of clusters
  clusters <- na.omit(object = stringr::str_sort(x = unique(x = rnascope.data[[clust.res.col]]),
                                                 numeric = TRUE))
  
  # Re-define cluster colors and labels (in case is needed)
  cluster.colors <- cluster.colors[as.character(x = clusters)]
  cluster.labels <- cluster.labels[as.character(x = clusters)]
  
  # Plot violin plots
  vln.plots <- StackViolinPlots(
    data.use = t(x = rnascope.data[,channels]),
    genes.use = channels,
    cluster.data = rnascope.data[,clust.res.col, drop = FALSE],
    cellnames.as.rownames.in.cluster.data = TRUE,
    cluster.ident.name = clust.res.col,
    cluster.colors = cluster.colors,
    cluster.breaks = cluster.labels,
    cluster.labels = cluster.labels,
    is.log.transformed = FALSE,
    log.scale = "log",
    pseudocount.use = 1,
    y.scale.trans = y.scale.trans,
    point.size = 2,
    alpha.use = 1,
    vln.border.colour = "black",
    vln.border.stroke = 0,
    plot.title.size = 6,
    plot.title.angle = 45,
    plot.title.face = "italic",
    font.family = "Arial",
    import.font = FALSE,
    hjust.use = 0,
    vjust.use = 0,
    x.axis.title = "max count ",
    axis.text.x.size = 5,
    axis.text.y.size = 5,
    axis.text.face = "plain",
    axis.ticks.length = 0.5,
    axis.line.size = rel(0.5),
    round.to.ceiling = TRUE,
    verbose = TRUE)
  
  # Add mouse and image info to plot title
  vln.plots <- vln.plots +
    plot_annotation(title = plot.title,
                    theme = theme(plot.title = element_text(family = "Arial",
                                                            face = "plain",
                                                            colour = "black",
                                                            size = 6,
                                                            hjust = 0.5,
                                                            vjust = 0.5)))
  
  # Save plot in pdf format
  FixSizeAndSave(plot = vln.plots,
                 filename = file.path(path,
                                      paste(gsub(pattern = "\\s[\n]*",
                                                 replacement = "_",
                                                 x = plot.title),
                                            "hclust_wardD2_violin_plots",
                                            y.scale.trans,
                                            "scale.pdf",
                                            sep = "_")),
                 is.ggassemble = TRUE,
                 panel.width = 0.5,
                 panel.height = 0.325 * length(x = clusters),
                 margin = 5,
                 unit.use = "cm",
                 use.ggsave = TRUE,
                 useDingbats = FALSE)
  
  # Return plot
  return(vln.plots)
}

# Plot and save dendrograms
PlotDendrograms <- function(hc.data,
                            cluster.colors,
                            cluster.labels,
                            plot.title,
                            path) {
  # Render cluster labels in hc.data to factor
  hc.data$segments$cluster.label <- ifelse(test = hc.data$segments$clust == 0,
                                           yes = NA,
                                           no = hc.data$segments$clust)
  hc.data$segments$cluster.label <- factor(x = hc.data$segments$cluster.label)
  
  # Re-define cluster colors and labels (in case is needed)
  cluster.colors <- cluster.colors[as.character(x = sort(x = unique(x = hc.data$labels$clust)))]
  cluster.labels <- cluster.labels[as.character(x = sort(x = unique(x = hc.data$labels$clust)))]
  
  # Plot dendrogram
  dendrogram.plot <- ggplot(data = hc.data$segments) + 
    geom_segment(mapping = aes(x = x,
                               y = y,
                               xend = xend,
                               yend = yend,
                               colour = cluster.label),
                 lineend = "square",
                 size = 0.1) +
    scale_x_continuous(expand = c(0,0)) +
    # scale_y_reverse(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_colour_manual(values = cluster.colors,
                        labels = cluster.labels,
                        na.value = "black",
                        name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2,
                                                     shape = 15))) +
    labs(title = paste(plot.title,
                       " (n=", nrow(x = hc.data$labels), ")",
                       sep = "")) +
    # coord_flip(clip = "off") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.title = element_text(size = 6,
                                    family = "Arial",
                                    colour = "black",
                                    hjust = 0.5,
                                    vjust = 0,
                                    margin = margin(t = 0,
                                                    r = 0,
                                                    b = 1,
                                                    l = 0,
                                                    unit = "mm")),
          legend.text = element_text(size = 5, 
                                     colour = "black",
                                     family = "Arial"),
          legend.margin = margin(t = 0,
                                 r = 0,
                                 b = 0,
                                 l = 0,
                                 unit = "mm"),
          legend.key.size = unit(x = 0,
                                 units = "mm"),
          plot.margin = margin(t = 0,
                               r = 1,
                               b = 0,
                               l = 1,
                               unit = "mm"))
  
  # Save plot in pdf format
  FixSizeAndSave(plot = dendrogram.plot,
                 filename = file.path(path,
                                      paste(gsub(pattern = "\\s[\n]*",
                                                 replacement = "_",
                                                 x = plot.title),
                                            "hclust_wardD2_clustering_dendrogram.pdf",
                                            sep = "_")),
                 is.ggassemble = FALSE,
                 panel.width = 7,
                 panel.height = 3,
                 unit.use = "cm",
                 margin = 0.25,
                 use.ggsave = TRUE,
                 useDingbats = FALSE)
  
  # Return plot
  return(dendrogram.plot)
}

# Plot and save spatial maps
PlotSpatialMaps <- function(rnascope.data,
                            clust.res.col,
                            composite = TRUE,
                            map.coordinates,
                            cluster.colors,
                            cluster.labels,
                            plot.title,
                            path) {
  # Set clusters as character
  rnascope.data[[clust.res.col]] <- as.character(x = rnascope.data[[clust.res.col]])
  
  if (any(is.na(x = rnascope.data[[clust.res.col]]))) {
    # Attribute "no expression" as cluster identity for non expressing cells
    rnascope.data[[clust.res.col]] <- ifelse(test = is.na(x = rnascope.data[[clust.res.col]]),
                                             yes = "no expression",
                                             no = rnascope.data[[clust.res.col]])
    
    # Reorder RNAscope data and place non expressing cells above
    rnascope.data$order <- ifelse(test = rnascope.data[[clust.res.col]] == "no expression",
                                  yes = 1,
                                  no = 2)
    rnascope.data <- rnascope.data[order(rnascope.data$order),]
    
    # Define color and label to non expressing cells
    cluster.colors["no expression"] <- "lightgrey"
    cluster.labels["no expression"] <- "no expression"
  }
  
  # Create a list of RNAscope datas for spatial maps generation
  rnascope.data.l <- list(complete.data = rnascope.data)
  if (isTRUE(x = composite)) {
    rnascope.data.l <- append(x = rnascope.data.l,
                              values = split(x = rnascope.data,
                                             f = rnascope.data[[clust.res.col]]))
    # Get names of clusters
    clusters <- na.omit(object = stringr::str_sort(x = unique(x = rnascope.data[[clust.res.col]]),
                                                   numeric = TRUE))
    # Reorder slots of rnascope.data.l
    rnascope.data.l <- rnascope.data.l[c("complete.data", clusters)]
  }
  
  # Get min and max of x and y centroid values 
  min.x <- min(rnascope.data$centroid_x)
  min.y <- min(rnascope.data$centroid_y)
  max.x <- max(rnascope.data$centroid_x)
  max.y <- max(rnascope.data$centroid_y)
  
  # Get range of x and y axis values
  range.x <- abs(x = max.x - min.x) # get range of x axis
  range.y <- abs(x = max.y - min.y) # get range of y axis
  
  # Get ratio of x and y ranges
  range.ratio <- range.x / range.y
  
  # Create a list with x and y centroid values
  axis.info <- list(min.x = min.x,
                    min.y = min.y,
                    max.x = max.x,
                    max.y = max.y,
                    range.x = range.x,
                    range.y = range.y,
                    range.ratio = range.ratio)
  
  # Generate spatial maps
  spatial.map.plot.l <- mapply(
    rnascope.data = rnascope.data.l,
    dataset.name = names(x = rnascope.data.l),
    FUN = function(rnascope.data,
                   dataset.name,
                   clust.res.col,
                   map.coordinates,
                   axis.info,
                   cluster.colors,
                   cluster.labels) {
      # Get names of clusters
      clusters <- na.omit(object = stringr::str_sort(x = unique(x = rnascope.data[[clust.res.col]]),
                                                     numeric = TRUE))
      
      # Re-define cluster colors and labels (in case is needed)
      cluster.colors <- cluster.colors[as.character(x = clusters)]
      cluster.labels <- cluster.labels[as.character(x = clusters)]
      
      # Get maximum range value from both axes
      range.max <- max(axis.info$range.x, axis.info$range.y)
      
      # Plot spatial map
      spatial.map.plot <- Plot2DEmbedding(
        embedding.data = rnascope.data[,c("centroid_x","centroid_y")],
        cluster.data = rnascope.data[,clust.res.col, drop = FALSE],
        cellnames.as.rownames.in.embedding.data = TRUE,
        cellnames.as.rownames.in.cluster.data = TRUE,
        embedding.x.name = "centroid_x",
        embedding.y.name = "centroid_y",
        embedding.x.label = map.coordinates$x.direction,
        embedding.y.label = map.coordinates$y.direction,
        cluster.ident.name = clust.res.col,
        cluster.colors = cluster.colors,
        cluster.breaks = cluster.labels,
        cluster.labels = cluster.labels,
        add.border.to.points = FALSE,
        point.size = 0.5,
        border.stroke = NULL,
        border.colour = NULL,
        alpha.use = 1,
        legend.title.face = "plain",
        font.family = "Arial",
        import.font = FALSE,
        fix.aspect.ratio = FALSE,
        legend.point.size = 2,
        axis.title.size = 5,
        legend.title.size = 6,
        legend.text.size = 5,
        legend.text.space = unit(x = 3, units = "mm"),
        legend.ncol = 1,
        axis.line.size = rel(0.5),
        arrow.size = rel(2),
        range.scale = 0.05
      )
      
      # Adjust axes to standardize across all plots
      spatial.map.plot <- spatial.map.plot +
        scale_x_continuous(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                           breaks = c(axis.info$min.x,
                                      axis.info$min.x + 0.05 * range.max)) + # equal arrow size for x and y axes after plot resizing
        scale_y_reverse(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                        breaks = c(axis.info$max.y - 0.05 * range.max, # equal arrow size for x and y axes after plot resizing
                                   axis.info$max.y)) +
        expand_limits(x = c(axis.info$min.x - 0.025 * range.max,
                            axis.info$max.x),
                      y = c(axis.info$min.y,
                            axis.info$max.y + 0.025 * range.max)) +
        theme(aspect.ratio = 1/axis.info$range.ratio)
      
      # Add plot title if complete data is plotted, else remove legends
      if (dataset.name == "complete.data") {
        spatial.map.plot <- spatial.map.plot +
          labs(title = plot.title) +
          theme(plot.title = element_text(size = 6,
                                          family = "Arial",
                                          colour = "black",
                                          hjust = 0.5,
                                          vjust = 0,
                                          margin = margin(t = 0,
                                                          r = 0,
                                                          b = 1,
                                                          l = 0,
                                                          unit = "mm")))
      } else {
        spatial.map.plot <- spatial.map.plot +
          theme(legend.position = "none")
      }
      
      # Reorder data in plot and match the order in RNAscope data
      spatial.map.plot$data <- spatial.map.plot$data[match(x = rnascope.data$cell, table = spatial.map.plot$data$Cell),]
      
      # Return plot
      return(spatial.map.plot)
    },
    MoreArgs = list(clust.res.col = clust.res.col,
                    map.coordinates = map.coordinates,
                    axis.info = axis.info,
                    cluster.colors = cluster.colors,
                    cluster.labels = cluster.labels),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE)
  
  # Create a patchwork plot with complete or composite spatial maps 
  spatial.map.plot <- patchwork::wrap_plots(spatial.map.plot.l,
                                            nrow = 1,
                                            guides = "collect")
  
  # Save plot in pdf format
  FixSizeAndSave(plot = spatial.map.plot,
                 filename = file.path(path,
                                      paste(gsub(pattern = "\\s[\n]*",
                                                 replacement = "_",
                                                 x = plot.title),
                                            "hclust_wardD2_spatial_maps.pdf",
                                            sep = "_")),
                 is.ggassemble = TRUE,
                 panel.width = min(7, 7 * range.ratio),
                 panel.height = min(7, 7 / range.ratio),
                 unit.use = "cm",
                 margin = 5,
                 use.ggsave = TRUE,
                 useDingbats = FALSE)
  
  # Return plot
  return(spatial.map.plot)
}

# Add coloured leafs to dendrogram
ColoredLeafs <- function(n,
                         clustering.results,
                         clustering.k,
                         cluster.colors) {
  if (is.leaf(object = n)) {
    # Get attributes of leaf
    leaf.attribute <- attributes(x = n)
    # Get label from attributes
    leaf.label <- leaf.attribute$label
    # Get cluster of label
    leaf.cluster <- clustering.results[clustering.results$cell == leaf.label, clustering.k]
    # Get color of cluster
    leaf.color <- cluster.colors[leaf.cluster]
    # Add data to attributes
    attr(x = n, which = "nodePar") <- c(leaf.attribute$nodePar,
                                        list(cex = 0.5,
                                             pch = 15,
                                             col = unname(obj = leaf.color)))
  }
  # Return data
  return(n)
}
