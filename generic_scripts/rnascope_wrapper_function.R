# RNAscope wrapper function
rnascope_wrapper_function <- function(files,
                                      rnascope.data,
                                      objects.to.keep = "PathCellObject",
                                      meta.data.l,
                                      channels,
                                      n.clusters = 2:12,
                                      log.transform = TRUE,
                                      pseudocount = 1,
                                      use.fastcluster = TRUE,
                                      compute.optimal.k = FALSE,
                                      nmi.thresh = 0.55,
                                      plot.violin.plots = TRUE,
                                      plot.dendrograms = TRUE,
                                      plot.spatial.maps = TRUE,
                                      split.spatial.maps.by,
                                      cluster.labels,
                                      cluster.colors,
                                      map.coordinates.l,
                                      file.prefix,
                                      main.path,
                                      spatial.maps.path.l) {
  # Mention in output files/file names if data was log transformed for clustering
  if (isTRUE(x = log.transform)) {
    file.prefix <- paste(file.prefix,
                         "log_transformed",
                         sep = "_")
  }
  
  if (!missing(x = files)) {
    # Load and prepare RNAscope data
    rnascope.data.l <- mapply(FUN = LoadRNAscopeData,
                              file = files,
                              meta.data = meta.data.l,
                              MoreArgs = list(objects.to.keep = objects.to.keep,
                                              channels = channels),
                              SIMPLIFY = FALSE,
                              USE.NAMES = TRUE)
    
    # Bind RNAscope data by row
    rnascope.data <- dplyr::bind_rows(rnascope.data.l,
                                      .id = NULL)
  } else if (missing(x = rnascope.data) && missing(x = files)) {
    stop("Either a list of files or an RNAscope dataset should be provided to the function")
  }
  
  # Cluster RNAscope data
  clustered.data <- ClusterRNAscopeData(rnascope.data = rnascope.data,
                                        channels = channels,
                                        n.clusters = n.clusters,
                                        log.transform = log.transform,
                                        pseudocount = pseudocount,
                                        use.fastcluster = use.fastcluster)
  
  # Save dendrogram plot
  pdf(file = file.path(main.path,
                       paste(file.prefix,
                             "hclust_wardD2_clustering_dendrogram.pdf",
                             sep = "_")),
      width = 15 / 2.54,
      height = 10 / 2.54)
  plot(x = clustered.data$hc.ward.D2,
       labels = FALSE,
       frame.plot = FALSE,
       hang = -1)
  dev.off()
  
  if (isTRUE(x = plot.dendrograms)) {
    # Create segments and labels data for dendrogram plot from clustering result
    hc.data <- ggdendro::dendro_data(model = clustered.data$hc.ward.D2,
                                     type = "rectangle")
    
    # Add cluster attributes to segments and labels data for various numbers of clusters
    hc.data.l <- lapply(X = n.clusters,
                        FUN = dendro_data_k,
                        hc = clustered.data$hc.ward.D2,
                        hcdata = hc.data)
  }
  
  if (isTRUE(x = plot.dendrograms) && isTRUE(x = compute.optimal.k)) {
    # Detect clustering k above which cluster merging using NMI score occurs
    k.optimal.s <- sapply(X = n.clusters,
                          FUN = function(k,
                                         clustered.data,
                                         hc.data,
                                         channels,
                                         log.transform,
                                         pseudocount,
                                         use.fastcluster,
                                         nmi.thresh) {
                            # Prepare data for cluster merging based on NMI score
                            full.data.nmi <- clustered.data$full.data # get full RNAscope data
                            full.data.nmi$merged_clusters <- full.data.nmi[[paste("k", k, "_clusters", sep = "")]] # select k for cluster merging
                            cols.to.remove <- grepl(pattern = "^k|_scaled|is_clustered",
                                                    x = colnames(x = full.data.nmi)) # select columns that might interfere with MergeClustersNMI function
                            full.data.nmi[,cols.to.remove] <- NULL # remove columns that might interfere with MergeClustersNMI function
                            
                            # Merge clusters using NMI score
                            clustered.data.merged <- MergeClustersNMI(rnascope.data = full.data.nmi,
                                                                      hc = clustered.data$hc.ward.D2,
                                                                      hcdata = hc.data,
                                                                      k = k,
                                                                      channels = channels,
                                                                      log.transform = log.transform,
                                                                      pseudocount = pseudocount,
                                                                      use.fastcluster = use.fastcluster,
                                                                      nmi.thresh = nmi.thresh)
                            
                            # Check if cluster merging occured
                            if (!clustered.data.merged$clusters_merged) {
                              return(k)
                            } else {
                              return(NA)
                            }
                          },
                          clustered.data = clustered.data,
                          hc.data = hc.data,
                          channels = channels,
                          log.transform = log.transform,
                          pseudocount = pseudocount,
                          use.fastcluster = use.fastcluster,
                          nmi.thresh = nmi.thresh,
                          simplify = TRUE,
                          USE.NAMES = TRUE)
    k.optimal <- min(as.numeric(x = names(x = k.optimal.s[is.na(x = k.optimal.s)]))) - 1
    
    # Add "optimal" to column name of k.optimal clustering result
    colnames(x = clustered.data$full.data) <- gsub(pattern = paste("k", k.optimal, "_clusters",
                                                                   sep = ""),
                                                   replacement = paste("k", k.optimal, "_clusters_optimal",
                                                                       sep = ""),
                                                   x = colnames(x = clustered.data$full.data))
  }
  
  # Save clustered data as RDS file
  saveRDS(object = clustered.data,
          file = file.path(main.path,
                           paste(file.prefix,
                                 "hclust_wardD2_clustered_data.rds",
                                 sep = "_")))
  
  # Get column name of clustering results
  clust.res.cols <- grep(pattern = "k[0-9]+_clusters",
                         x = colnames(x = clustered.data$full.data),
                         value = TRUE)
  clust.res.cols <- stringr::str_sort(x = clust.res.cols,
                                      numeric = TRUE) # sort in case is needed
  names(x = clust.res.cols) <- n.clusters
  
  # Get various plot titles for each k for violin and dendrogram plots and spatial maps
  plot.titles <- paste(file.prefix,
                       clust.res.cols, # various clusters as defined above
                       sep = "\n")
  names(x = plot.titles) <- n.clusters
  
  if (isTRUE(x = plot.violin.plots)) {
    # Plot and save violin plots for each clustering result (explorative)
    vln.plots.l <- lapply(X = list("identity",
                                   "log10"),
                          FUN = function(y.scale.trans,
                                         clust.res.cols,
                                         plot.titles,
                                         rnascope.data,
                                         channels,
                                         cluster.colors,
                                         cluster.labels,
                                         main.path) {
                            # Define path where to save violin plots
                            violin.plots.path <- file.path(main.path,
                                                           paste("violin_plots",
                                                                 y.scale.trans,
                                                                 "scale",
                                                                 sep = "_"))
                            
                            # Create directory where to save violin plots
                            dir.create(path = violin.plots.path)
                            
                            # Plot and save violin plots for each clustering result
                            mapply(clust.res.col = clust.res.cols,
                                   plot.title = plot.titles,
                                   FUN = PlotStackedViolinPlots,
                                   MoreArgs = list(rnascope.data = rnascope.data,
                                                   channels = channels,
                                                   cluster.colors = cluster.colors,
                                                   cluster.labels = cluster.labels,
                                                   y.scale.trans = y.scale.trans,
                                                   path = violin.plots.path),
                                   SIMPLIFY = FALSE,
                                   USE.NAMES = TRUE)
                          },
                          clust.res.cols = clust.res.cols,
                          plot.titles = plot.titles,
                          rnascope.data = clustered.data$full.data[clustered.data$full.data$is_clustered,],
                          channels = channels,
                          cluster.colors = cluster.colors,
                          cluster.labels = cluster.labels,
                          main.path = main.path)
  }
  
  if (isTRUE(x = plot.dendrograms)) {
    # Define path where to save dendrograms
    clustering.dendrograms.path <- file.path(main.path,
                                             "clustering_dendrograms")
    
    # Create directory where to save dendrograms
    dir.create(path = clustering.dendrograms.path)
    
    # Plot and save dendrograms for each k (explorative)
    dendrogram.plots.l <- mapply(hc.data = hc.data.l,
                                 plot.title = plot.titles,
                                 FUN = PlotDendrograms,
                                 MoreArgs = list(cluster.colors = cluster.colors,
                                                 cluster.labels = cluster.labels,
                                                 path = clustering.dendrograms.path),
                                 SIMPLIFY = FALSE,
                                 USE.NAMES = TRUE)
  }
  
  # # Get dendrogram data from clustering results
  # dendrogram.data <- as.dendrogram(object = clustered.data$hc.ward.D2)
  # 
  # # Plot and save dendrograms for each k (explorative)
  # dendrogram.plots.base.l <- lapply(
  #   X = gsub(pattern = "_optimal",
  #            replacement = "",
  #            x = clust.res.cols),
  #   FUN = function(clust.res.col,
  #                  dendrogram.data,
  #                  clustered.data,
  #                  cluster.colors,
  #                  main.path,
  #                  file.prefix) {
  #     # Add coloured leafs to dendrogram
  #     dendrogram.data <- dendrapply(X = dendrogram.data,
  #                                   FUN = ColoredLeafs,
  #                                   clustering.results = clustered.data$clusters.ward.D2,
  #                                   clustering.k = clust.res.col,
  #                                   cluster.colors = cluster.colors)
  # 
  #     # Save dendrogram plot
  #     pdf(file = file.path(main.path,
  #                          paste(file.prefix,
  #                                clust.res.col,
  #                                "hclust_wardD2_clustering_dendrogram_baseR.pdf",
  #                                sep = "_")),
  #         width = 15 / 2.54,
  #         height = 10 / 2.54)
  #     plot(x = dendrogram.data,
  #          leaflab = "none",
  #          main = clust.res.col,
  #          xlab = clustered.data$hc.ward.D2$call,
  #          ylab = "Height")
  #     dev.off()
  # 
  #     # Return "DONE"
  #     return("DONE")
  #   },
  #   dendrogram.data = dendrogram.data,
  #   clustered.data = clustered.data,
  #   cluster.colors = cluster.colors,
  #   main.path = main.path,
  #   file.prefix = file.prefix)
  
  if (isTRUE(x = plot.spatial.maps)) {
    # Plot and save spatial maps
    spatial.maps.l <- mapply(
      clust.res.col = clust.res.cols,
      plot.title = plot.titles,
      FUN = function(rnascope.data,
                     split.spatial.maps.by,
                     clust.res.col,
                     map.coordinates.l,
                     cluster.colors,
                     cluster.labels,
                     plot.title,
                     spatial.maps.path.l) {
        # Split RNAscope data by slides
        rnascope.data.l <- split(x = rnascope.data,
                                 f = rnascope.data[,split.spatial.maps.by],
                                 sep = "_")
        
        # Remove empty datasets
        empty.data <- lapply(X = rnascope.data.l, FUN = nrow) == 0
        rnascope.data.l <- rnascope.data.l[!empty.data]
        
        # Sort datasets by names
        rnascope.data.l <- rnascope.data.l[sort(x = names(x = rnascope.data.l))]
        
        # Create as many plot titles as there are slides
        plot.titles.l <- list()
        for (n in names(x = rnascope.data.l)) {
          plot.titles.l[[n]] <- gsub(pattern = "\\\n",
                                     replacement = paste0("\n",n,"\n"),
                                     x = plot.title)
        }
        
        # Plot and save spatial maps 
        mapply(rnascope.data = rnascope.data.l,
               map.coordinates = map.coordinates.l,
               plot.title = plot.titles.l,
               path = spatial.maps.path.l,
               FUN = PlotSpatialMaps,
               MoreArgs = list(clust.res.col = clust.res.col,
                               composite = TRUE,
                               cluster.colors = cluster.colors,
                               cluster.labels = cluster.labels),
               SIMPLIFY = FALSE,
               USE.NAMES = TRUE)
      },
      MoreArgs = list(rnascope.data = clustered.data$full.data,
                      split.spatial.maps.by = split.spatial.maps.by,
                      map.coordinates.l = map.coordinates.l,
                      cluster.colors = cluster.colors,
                      cluster.labels = cluster.labels,
                      spatial.maps.path.l = spatial.maps.path.l),
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE)
  }
  
  # Return "DONE"
  return("DONE")
}
