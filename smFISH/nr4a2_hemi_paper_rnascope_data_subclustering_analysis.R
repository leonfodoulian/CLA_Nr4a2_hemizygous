# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/results/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/rnascope_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/rnascope_wrapper_function.R")
source(file = "/Users/leonfodoulian/scData/dendro_data_k.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/element_textbox_highlight_5.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")
source(file = "/Users/leonfodoulian/scData/cat.verbose.R")

# Required packages
require(readxl)
require(ggplot2)
require(reshape2)
require(ggdendro)
require(extrafont)
require(paletteer)
require(spsUtil)

# Define data attributes for subclustering
subclust.attributes <- list("C2.Egr2_C3.Ndufs4_C4.Nr4a2" = list(dir = "e20220929_QuPath/C2.Egr2_C3.Ndufs4_C4.Nr4a2/all_cells",
                                                                probe.set = "C2.Egr2_C3.Ndufs4_C4.Nr4a2",
                                                                cluster.to.subset = 3,
                                                                clust.res.col = "k3_clusters"),
                            "C2.Egr2_C3.Slc17a6_C4.Nnat" = list(dir = "e20220929_QuPath/C2.Egr2_C3.Slc17a6_C4.Nnat/all_cells",
                                                                probe.set = "C2.Egr2_C3.Slc17a6_C4.Nnat",
                                                                cluster.to.subset = 2,
                                                                clust.res.col = "k2_clusters"),
                            "C2.Lxn_C3.Ndufs4_C4.Nr4a2" = list(dir = "e20220929_QuPath/C2.Lxn_C3.Ndufs4_C4.Nr4a2/all_cells",
                                                               probe.set = "C2.Lxn_C3.Ndufs4_C4.Nr4a2",
                                                               cluster.to.subset = 2,
                                                               clust.res.col = "k2_clusters"),
                            "C2.Lxn_C3.Oprk1_C4.Fosl2" = list(dir = "e20230412_QuPath/C2.Lxn_C3.Oprk1_C4.Fosl2/all_cells",
                                                              probe.set = "C2.Lxn_C3.Oprk1_C4.Fosl2",
                                                              cluster.to.subset = 2,
                                                              clust.res.col = "k2_clusters"))

# Define zone to map coordinate correspondence
zone.map.coordinate.correspondence <- c("zone1" = "lateral",
                                        "zone2" = "medial",
                                        "zone3" = "lateral",
                                        "zone4" = "medial")

# Define attributes for spatial map generation
split.spatial.maps.by <- c("slide",
                           "zone",
                           "mouse",
                           "genotype",
                           "level")

# Define various numbers of clusters to compute
n.clusters <- 2:12
names(x = n.clusters) <- n.clusters

# Define cluster labels
cluster.labels <- as.character(x = 1:12)
names(x = cluster.labels) <- 1:12

# Define cluster colors
cluster.colors <- as.character(x = paletteer::paletteer_d(palette = "LaCroixColoR::paired",
                                                          n = 12,
                                                          type = "discrete"))
names(x = cluster.colors) <- names(x = cluster.labels)

# Run wrapper function for each experiment
lapply(X = subclust.attributes,
       FUN = function(subclust.attribute,
                      zone.map.coordinate.correspondence,
                      n.clusters,
                      split.spatial.maps.by,
                      cluster.labels,
                      cluster.colors) {
         # Get names of channels from probe set
         channels <- gsub(pattern = ".*\\.",
                          replacement = "",
                          x = unlist(x = stringr::str_split(string = subclust.attribute$probe.set,
                                                            pattern = "_")))
         
         # List directories from parent folder
         parent.dirs <- list.dirs(path = subclust.attribute$dir,
                                  full.names = TRUE,
                                  recursive = TRUE)
         
         # Define new subdirectory where to save files
         subdirectory <- gsub(pattern = "s$",
                              replacement = paste0(subclust.attribute$cluster.to.subset,
                                                   "_cells"),
                              x = subclust.attribute$clust.res.col)
         
         # Get list of new directories by replacing "all_cells" with new subdirectory 
         child.dirs <- gsub(pattern = "all_cells",
                            replacement = subdirectory,
                            x = parent.dirs)
         
         # Create directories where to save files
         lapply(X = child.dirs,
                FUN = function(path) {
                  dir.create(path = path,
                             recursive = TRUE)
                })
         
         # Get list of directories where spatial maps are saved
         spatial.maps.path.l <- child.dirs[grepl(pattern = "Slide", x = child.dirs)]
         names(x = spatial.maps.path.l) <- spatial.maps.path.l
         
         # Create data with appropriate spatial map coordinates for each image
         map.coordinates.l <- lapply(X = strsplit(x = spatial.maps.path.l,
                                                  split = "_"),
                                     FUN = function(spatial.maps.path,
                                                    zone.map.coordinate.correspondence) {
                                       zones <- names(x = zone.map.coordinate.correspondence)
                                       zone <- zones[zones %in% spatial.maps.path]
                                       data.frame(x.direction = unname(obj = zone.map.coordinate.correspondence[zone]),
                                                  y.direction = "dorsal")
                                     },
                                     zone.map.coordinate.correspondence = zone.map.coordinate.correspondence)
         
         # Adjust spatial map coordinate of one image with mismatch between zone and map coordinate correspondence
         if ("e20230412_QuPath/C2.Lxn_C3.Oprk1_C4.Fosl2/k2_cluster2_cells/Slide1_zone3_B6_wt.wt_level44" %in% names(x = map.coordinates.l)) {
           map.coordinates.l[["e20230412_QuPath/C2.Lxn_C3.Oprk1_C4.Fosl2/k2_cluster2_cells/Slide1_zone3_B6_wt.wt_level44"]]$x.direction <- "medial"
         }
         
         # Define RNAscope rds file
         rds.file <- list.files(path = subclust.attribute$dir,
                                pattern = ".rds",
                                full.names = TRUE)
         
         # Load RNAscope data
         rnascope.data <- readRDS(file = rds.file)$full.data
         
         # Get exact clustering result column name
         clust.res.col <- grep(pattern = subclust.attribute$clust.res.col,
                               x = colnames(x = rnascope.data),
                               value = TRUE)
         
         # Subset data for cluster of interest
         cells.to.keep <- rnascope.data[[clust.res.col]] == subclust.attribute$cluster.to.subset & !is.na(x = rnascope.data[[clust.res.col]])
         rnascope.data <- rnascope.data[cells.to.keep,]
         
         # Subset data for columns of interest
         columns.to.keep <- c("cell", split.spatial.maps.by, channels, "centroid_x", "centroid_y", "nucleus_area")
         rnascope.data <- rnascope.data[,columns.to.keep]
         
         # Run wrapper function on RNAscope data
         quiet(x = rnascope_wrapper_function(rnascope.data = rnascope.data,
                                             channels = channels,
                                             n.clusters = n.clusters,
                                             log.transform = TRUE,
                                             pseudocount = 1,
                                             use.fastcluster = TRUE,
                                             compute.optimal.k = TRUE,
                                             nmi.thresh = 0.55,
                                             plot.violin.plots = TRUE,
                                             plot.dendrograms = TRUE,
                                             plot.spatial.maps = TRUE,
                                             split.spatial.maps.by = split.spatial.maps.by,
                                             cluster.labels = cluster.labels,
                                             cluster.colors = cluster.colors,
                                             map.coordinates.l = map.coordinates.l,
                                             file.prefix = paste(gsub(pattern = "\\/.*",
                                                                      replacement = "",
                                                                      x =  subclust.attribute$dir),
                                                                 subclust.attribute$probe.set,
                                                                 sep = "_"),
                                             main.path = unique(x = gsub(pattern = "/Slide.*",
                                                                         replacement = "",
                                                                         x = spatial.maps.path.l)),
                                             spatial.maps.path.l = spatial.maps.path.l),
               print_cat = TRUE,
               message = TRUE,
               warning = TRUE)
       },
       zone.map.coordinate.correspondence = zone.map.coordinate.correspondence,
       n.clusters = n.clusters,
       split.spatial.maps.by = split.spatial.maps.by,
       cluster.labels = cluster.labels,
       cluster.colors = cluster.colors)
