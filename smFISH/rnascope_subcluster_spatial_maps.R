# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/")

# Set stringsAsFactors to FALSE
options(stringsAsFactors = FALSE)

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/spatial_registration_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/rnascope_claustro_insular_maps_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/discretise.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(patchwork)
require(ggplot2)
require(extrafont)
require(ggtext)
require(ggnewscale)
require(ggrastr)
require(uniqtag)

# Create directory where to save files
dir <- "./figures/subcluster_spatial_maps"
dir.create(path = dir,
           recursive = TRUE)

# Define specific parameters for each file
file.attributes.l <- list(
  "C2.Nr4a2_C3.Rprm_C4.Col6a1" = list(rds.file = "./results/e20220913_QuPath/C2.Nr4a2_C3.Rprm_C4.Col6a1/all_cells/e20220913_QuPath_C2.Nr4a2_C3.Rprm_C4.Col6a1_log_transformed_hclust_wardD2_clustered_data.rds",
                                      probe.set = "C2.Nr4a2_C3.Rprm_C4.Col6a1",
                                      clusters.to.subset = c(1,3,5),
                                      slide.to.subset = "Slide4",
                                      zone.to.subset = "zone3",
                                      clust.res.col = "k5_clusters"),
  "C2.Nr4a2_C3.Rprm_C4.Col6a1_2" = list(rds.file = "./results/e20220913_QuPath/C2.Nr4a2_C3.Rprm_C4.Col6a1/all_cells/e20220913_QuPath_C2.Nr4a2_C3.Rprm_C4.Col6a1_log_transformed_hclust_wardD2_clustered_data.rds",
                                        probe.set = "C2.Nr4a2_C3.Rprm_C4.Col6a1",
                                        clusters.to.subset = c(1,3,4,6),
                                        slide.to.subset = "Slide4",
                                        zone.to.subset = "zone1",
                                        clust.res.col = "k6_clusters"),
  "C2.Egr2_C3.Ndufs4_C4.Nr4a2" = list(rds.file = "./results/e20220929_QuPath/C2.Egr2_C3.Ndufs4_C4.Nr4a2/k3_cluster3_cells/e20220929_QuPath_C2.Egr2_C3.Ndufs4_C4.Nr4a2_log_transformed_hclust_wardD2_clustered_data.rds",
                                      probe.set = "C2.Egr2_C3.Ndufs4_C4.Nr4a2",
                                      clusters.to.subset = c(4,8,10),
                                      slide.to.subset = "Slide6",
                                      zone.to.subset = "zone2",
                                      clust.res.col = "k11_clusters"),
  "C2.Lxn_C3.Ndufs4_C4.Nr4a2" = list(rds.file = "./results/e20220929_QuPath/C2.Lxn_C3.Ndufs4_C4.Nr4a2/k2_cluster2_cells/e20220929_QuPath_C2.Lxn_C3.Ndufs4_C4.Nr4a2_log_transformed_hclust_wardD2_clustered_data.rds",
                                     probe.set = "C2.Lxn_C3.Ndufs4_C4.Nr4a2",
                                     clusters.to.subset = c(4,6),
                                     slide.to.subset = "Slide5",
                                     zone.to.subset = "zone3",
                                     clust.res.col = "k7_clusters"),
  "C2.Egr2_C3.Slc17a6_C4.Nnat" = list(rds.file = "./results/e20220929_QuPath/C2.Egr2_C3.Slc17a6_C4.Nnat/k2_cluster2_cells/e20220929_QuPath_C2.Egr2_C3.Slc17a6_C4.Nnat_log_transformed_hclust_wardD2_clustered_data.rds",
                                      probe.set = "C2.Egr2_C3.Slc17a6_C4.Nnat",
                                      clusters.to.subset = c(5,7),
                                      slide.to.subset = "Slide1",
                                      zone.to.subset = "zone3",
                                      clust.res.col = "k11_clusters"),
  "C2.Npsr1_C3.Rorb_C4.Nr4a2" = list(rds.file = "./results/e20220929_QuPath/C2.Npsr1_C3.Rorb_C4.Nr4a2/all_cells/e20220929_QuPath_C2.Npsr1_C3.Rorb_C4.Nr4a2_log_transformed_hclust_wardD2_clustered_data.rds",
                                     probe.set = "C2.Npsr1_C3.Rorb_C4.Nr4a2",
                                     clusters.to.subset = c(4,8),
                                     slide.to.subset = "Slide3",
                                     zone.to.subset = "zone1",
                                     clust.res.col = "k8_clusters"),
  "C2.Npy_C3.Nr4a2_C4.Gad1" = list(rds.file = "./results/e20220929_QuPath/C2.Npy_C3.Nr4a2_C4.Gad1/all_cells/e20220929_QuPath_C2.Npy_C3.Nr4a2_C4.Gad1_log_transformed_hclust_wardD2_clustered_data.rds",
                                   probe.set = "C2.Npy_C3.Nr4a2_C4.Gad1",
                                   clusters.to.subset = c(1,12),
                                   slide.to.subset = "Slide2",
                                   zone.to.subset = "zone1",
                                   clust.res.col = "k12_clusters"),
  "C2.Lxn_C3.Oprk1_C4.Fosl2" = list(rds.file = "./results/e20230412_QuPath/C2.Lxn_C3.Oprk1_C4.Fosl2/k2_cluster2_cells/e20230412_QuPath_C2.Lxn_C3.Oprk1_C4.Fosl2_log_transformed_hclust_wardD2_clustered_data.rds",
                                    probe.set = "C2.Lxn_C3.Oprk1_C4.Fosl2",
                                    clusters.to.subset = c(4,7),
                                    slide.to.subset = "Slide2",
                                    zone.to.subset = "zone3",
                                    clust.res.col = "k7_clusters"),
  "C2.Lxn_C3.Oprk1_C4.Fosl2_2" = list(rds.file = "./results/e20230412_QuPath/C2.Lxn_C3.Oprk1_C4.Fosl2/k2_cluster2_cells/e20230412_QuPath_C2.Lxn_C3.Oprk1_C4.Fosl2_log_transformed_hclust_wardD2_clustered_data.rds",
                                      probe.set = "C2.Lxn_C3.Oprk1_C4.Fosl2",
                                      clusters.to.subset = c(2,4),
                                      slide.to.subset = "Slide2",
                                      zone.to.subset = "zone4",
                                      clust.res.col = "k4_clusters"),
  "C2.Nr4a2_C3.Nnat_C4.Rprm" = list(rds.file = "./results/e20241217_QuPath/C2.Nr4a2_C3.Nnat_C4.Rprm/all_cells/e20241217_QuPath_C2.Nr4a2_C3.Nnat_C4.Rprm_log_transformed_hclust_wardD2_clustered_data.rds",
                                    probe.set = "C2.Nr4a2_C3.Nnat_C4.Rprm",
                                    clusters.to.subset = c(2,4,5,6,7),
                                    slide.to.subset = "Slide1",
                                    zone.to.subset = "zone3",
                                    clust.res.col = "k7_clusters"))

# Load and preprocess RNAscope data
rnascope.data.l <- unlist(x = lapply(X = unname(obj = file.attributes.l),
                                     FUN = function(file.attributes) {
                                       # Load RNAscope data
                                       rnascope.data <- readRDS(file = file.attributes$rds.file)$full.data
                                       
                                       # Add clustering data to a generic column
                                       rnascope.data$clustering_column <- rnascope.data[[file.attributes$clust.res.col]]
                                       
                                       # Subset data for cluster of interest
                                       cells.to.keep <- rnascope.data[[file.attributes$clust.res.col]] %in% file.attributes$clusters.to.subset &
                                         !is.na(x = rnascope.data[[file.attributes$clust.res.col]]) &
                                         rnascope.data$slide == file.attributes$slide.to.subset &
                                         rnascope.data$zone == file.attributes$zone.to.subset
                                       rnascope.data <- rnascope.data[cells.to.keep,]
                                       
                                       # Get data attributes (used for plot title and name later)
                                       data.attributes <- paste(gsub(pattern = ".*\\/|_hclust.*",
                                                                     replacement = "",
                                                                     x = file.attributes$rds.file),
                                                                paste(unique(x = rnascope.data$slide),
                                                                      unique(x = rnascope.data$zone),
                                                                      unique(x = rnascope.data$mouse),
                                                                      unique(x = rnascope.data$genotype),
                                                                      unique(x = rnascope.data$level),
                                                                      sep = "_"),
                                                                file.attributes$clust.res.col,
                                                                sep = "\n")
                                       
                                       # Split data by cluster
                                       rnascope.data.l <- split(x = rnascope.data,
                                                                f = rnascope.data[[file.attributes$clust.res.col]])
                                       names(x = rnascope.data.l) <- paste0(data.attributes,
                                                                            "_hclust_wardD2_cluster",
                                                                            names(x = rnascope.data.l),
                                                                            "_spatial_maps_uniform_scale.pdf")
                                       
                                       # Return list of data
                                       return(rnascope.data.l)
                                     }),
                          recursive = FALSE,
                          use.names = TRUE)

# Get axis attributes
axis.attributes.l <- lapply(X = rnascope.data.l,
                            FUN = function(rnascope.data) {
                              GetAxisAttributes(x = rnascope.data$centroid_x,
                                                y = rnascope.data$centroid_y,
                                                bin.size = 10)
                            })

# Get maximum ranges of x and y coordinates
max.range.x <- max(unlist(x = lapply(X = axis.attributes.l,
                                     FUN = "[[",
                                     "range.x")))
max.range.y <- max(unlist(x = lapply(X = axis.attributes.l,
                                     FUN = "[[",
                                     "range.y")))

# Expand x and y coordinates for each data
axis.attributes.l <- lapply(X = axis.attributes.l,
                            FUN = function(axis.attributes,
                                           max.range.x,
                                           max.range.y) {
                              # Get difference in ranges of x and y coordinates
                              diff.range.x <- max.range.x - axis.attributes$range.x 
                              diff.range.y <- max.range.y - axis.attributes$range.y
                              
                              # Adjust minimum and maximum x and y coordinate values
                              axis.attributes$min.x <- axis.attributes$min.x - diff.range.x/2
                              axis.attributes$max.x <- axis.attributes$max.x + diff.range.x/2
                              axis.attributes$min.y <- axis.attributes$min.y - diff.range.y/2
                              axis.attributes$max.y <- axis.attributes$max.y + diff.range.y/2
                              
                              # Adjust range of x and y coordinate values
                              axis.attributes$range.x <- max.range.x
                              axis.attributes$range.y <- max.range.y
                              
                              # Adjust ratio of x and y ranges
                              axis.attributes$range.ratio <- axis.attributes$range.x / axis.attributes$range.y
                              
                              # Adjust maximum range value from both axes
                              axis.attributes$range.max <- max(axis.attributes$range.x, axis.attributes$range.y)
                              
                              # Return data
                              return(axis.attributes)
                            },
                            max.range.x = max.range.x,
                            max.range.y = max.range.y)

# Plot spatial maps
spatial.map.plots <- mapply(
  plot.data = rnascope.data.l,
  plot.title = names(x = rnascope.data.l),
  axis.attributes = axis.attributes.l,
  FUN = function(plot.data,
                 plot.title,
                 axis.attributes,
                 path) {
    # Set clusters as character
    plot.data[["clustering_column"]] <- as.character(x = plot.data[["clustering_column"]])
    
    # Define cluster colors, breaks and labels
    cluster.breaks <- unique(x =  plot.data[["clustering_column"]])
    cluster.colors <- setNames(object = "#A9A9A9",
                               nm = cluster.breaks)
    cluster.labels <- setNames(object = cluster.breaks,
                               nm = cluster.breaks)
    
    # Plot spatial map
    spatial.map.plot <- Plot2DEmbedding(
      embedding.data = plot.data[,c("centroid_x", "centroid_y")],
      cluster.data = plot.data[,"clustering_column", drop = FALSE],
      cellnames.as.rownames.in.embedding.data = TRUE,
      cellnames.as.rownames.in.cluster.data = TRUE,
      embedding.x.name = "centroid_x",
      embedding.y.name = "centroid_y",
      embedding.x.label = "lateral",
      embedding.y.label = "dorsal",
      cluster.ident.name = "clustering_column",
      cluster.colors = cluster.colors,
      cluster.breaks = cluster.breaks,
      cluster.labels = cluster.labels,
      # add.border.to.points = FALSE,
      add.border.to.points = TRUE,
      point.size = 1.5,
      # border.stroke = NULL,
      border.stroke = NA,
      # border.colour = NULL,
      border.colour = "transparent",
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
      arrow.type = "closed",
      axis.lineend = "butt",
      axis.gap = 0,
      range.scale = 0.05
    )
    
    # Adjust axes to standardize across all plots
    spatial.map.plot <- spatial.map.plot +
      {
        if (unique(x = plot.data$zone) %in% c("zone1", "zone3")) {
          scale_x_continuous(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                             breaks = c(axis.attributes$min.x,
                                        axis.attributes$min.x + 0.1 * axis.attributes$range.max)) # equal arrow size for x and y axes after plot resizing
        } else if (unique(x = plot.data$zone) %in% c("zone2", "zone4")) {
          scale_x_reverse(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                          breaks = c(axis.attributes$max.x - 0.1 * axis.attributes$range.max, # equal arrow size for x and y axes after plot resizing
                                     axis.attributes$max.x))
        }
      } +
      scale_y_reverse(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                      breaks = c(axis.attributes$max.y - 0.1 * axis.attributes$range.max, # equal arrow size for x and y axes after plot resizing
                                 axis.attributes$max.y)) +
      {
        if (unique(x = plot.data$zone) %in% c("zone1", "zone3")) {
          expand_limits(x = c(axis.attributes$min.x - 0.025 * axis.attributes$range.max,
                              axis.attributes$max.x),
                        y = c(axis.attributes$min.y,
                              axis.attributes$max.y + 0.025 * axis.attributes$range.max))
        } else if (unique(x = plot.data$zone) %in% c("zone2", "zone4")) {
          expand_limits(x = c(axis.attributes$min.x,
                              axis.attributes$max.x + 0.025 * axis.attributes$range.max),
                        y = c(axis.attributes$min.y,
                              axis.attributes$max.y + 0.025 * axis.attributes$range.max))
        }
      } +
      labs(title = gsub(pattern = "_hclust.*",
                        replacement = "",
                        x = plot.title)) +
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
            legend.text = element_markdown(size = 6,
                                           family = "Arial",
                                           colour = "black"),
            aspect.ratio = 1/axis.attributes$range.ratio)
    
    # Remove clipping of data
    spatial.map.plot$coordinates$clip <- "off"
    
    # Save plot in pdf format
    FixSizeAndSave(plot = spatial.map.plot,
                   filename = file.path(path,
                                        gsub(pattern = "\n",
                                             replacement = "_",
                                             x = plot.title)),
                   is.ggassemble = FALSE,
                   # panel.width = min(7, 7 * range.ratio),
                   # panel.height = min(7, 7 / range.ratio),
                   panel.width = 4.2, # figure limit for width
                   panel.height = 4.2 / axis.attributes$range.ratio,
                   unit.use = "cm",
                   margin = 2,
                   use.ggsave = TRUE,
                   useDingbats = FALSE)
    
    # Return plot
    return(spatial.map.plot)
  },
  MoreArgs = list(path = dir),
  USE.NAMES = TRUE,
  SIMPLIFY = FALSE)
