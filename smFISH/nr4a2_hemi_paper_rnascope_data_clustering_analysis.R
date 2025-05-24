# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/data/")

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

# List experiment directories on which function will be executed
dirs <- list.dirs(path = ".",
                  full.names = FALSE,
                  recursive = FALSE)
dirs <- grep(pattern = "_QuPath",
             x = dirs,
             value = TRUE)
names(x = dirs) <- dirs

# Define attributes to extract from file names
attributes <- c("slide",
                "zone",
                "mouse",
                "genotype",
                "level",
                "Channel 2",
                "Channel 3",
                "Channel 4")

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

# Define whether optimal k should be computed
compute.optimal.k.l <- list(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE)
names(x = compute.optimal.k.l) <- dirs

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
mapply(dir = dirs,
       compute.optimal.k = compute.optimal.k.l,
       FUN = function(dir,
                      compute.optimal.k,
                      attributes,
                      zone.map.coordinate.correspondence,
                      n.clusters,
                      split.spatial.maps.by,
                      cluster.labels,
                      cluster.colors) {
         # List files to load
         files <- list.files(path = dir,
                             pattern = ".xlsx",
                             full.names = TRUE)
         names(x = files) <- files
         
         # Keep only RNAscope data
         files <- files[grepl(pattern = "Slide", x = files)]
         
         # Extract file attributes from file names
         file.attributes.l <- lapply(X = files,
                                     FUN = GetFileAttributes,
                                     attributes = attributes,
                                     channel.attributes = 6:8,
                                     split.pattern = "_",
                                     prefix = dir,
                                     suffix = ".xlsx",
                                     path = file.path("../results", dir),
                                     subdirectory = "all_cells")
         
         # Get list of directories where to save files
         path.l <- lapply(X = file.attributes.l,
                          FUN = "[[",
                          "path")
         
         # Create directories where to save files
         lapply(X = path.l,
                FUN = function(path) {
                  dir.create(path = path,
                             recursive = TRUE)
                })
         
         # Split file attributes by probe sets
         file.attributes.l <- split(x = file.attributes.l,
                                    f = gsub(pattern = ".*QuPath/|/Slide.*|/all_cells.*",
                                             replacement = "",
                                             x = path.l))
         
         # Run wrapper function for each probe set
         mapply(file.attributes = file.attributes.l,
                probe.set = names(x = file.attributes.l),
                FUN = function(file.attributes,
                               probe.set,
                               dir,
                               zone.map.coordinate.correspondence,
                               n.clusters,
                               compute.optimal.k,
                               split.spatial.maps.by,
                               cluster.labels,
                               cluster.colors) {
                  # Get files to load
                  files <- lapply(X = file.attributes,
                                  FUN = "[[",
                                  "file")
                  
                  # Get meta data
                  meta.data.l <- lapply(X = file.attributes,
                                        FUN = "[[",
                                        "meta.data")
                  
                  # Get list of channels (probe set)
                  channels <- gsub(pattern = ".*\\.",
                                   replacement = "",
                                   x = unlist(x = unique(x = lapply(X = file.attributes,
                                                                    FUN = "[[",
                                                                    "channels"))))
                  
                  # Get list of directories where spatial maps are saved
                  spatial.maps.path.l <- lapply(X = file.attributes,
                                                FUN = "[[",
                                                "path")
                  
                  # Create data with appropriate spatial map coordinates for each image
                  map.coordinates.l <- lapply(X = file.attributes,
                                              FUN = function(file.attributes,
                                                             zone.map.coordinate.correspondence) {
                                                data.frame(x.direction = unname(obj = zone.map.coordinate.correspondence[file.attributes$meta.data["zone"]]),
                                                           y.direction = "dorsal")
                                              },
                                              zone.map.coordinate.correspondence = zone.map.coordinate.correspondence)
                  
                  # Adjust spatial map coordinate of one image with mismatch between zone and map coordinate correspondence
                  if ("e20230412_QuPath/Slide1_zone3_B6_wt.wt_level44_C2.Lxn_C3.Oprk1_C4.Fosl2.xlsx" %in% names(x = map.coordinates.l)) {
                    map.coordinates.l[["e20230412_QuPath/Slide1_zone3_B6_wt.wt_level44_C2.Lxn_C3.Oprk1_C4.Fosl2.xlsx"]]$x.direction <- "medial"
                  }
                  
                  # Run wrapper function on RNAscope data
                  quiet(x = rnascope_wrapper_function(files = files,
                                                      objects.to.keep = "PathCellObject",
                                                      meta.data.l = meta.data.l,
                                                      channels = channels,
                                                      n.clusters = n.clusters,
                                                      log.transform = TRUE,
                                                      pseudocount = 1,
                                                      use.fastcluster = TRUE,
                                                      compute.optimal.k = compute.optimal.k,
                                                      nmi.thresh = 0.55,
                                                      plot.violin.plots = TRUE,
                                                      plot.dendrograms = TRUE,
                                                      plot.spatial.maps = TRUE,
                                                      split.spatial.maps.by = split.spatial.maps.by,
                                                      cluster.labels = cluster.labels,
                                                      cluster.colors = cluster.colors,
                                                      map.coordinates.l = map.coordinates.l,
                                                      file.prefix = paste(dir,
                                                                          probe.set,
                                                                          sep = "_"),
                                                      main.path = unique(x = gsub(pattern = "/Slide.*",
                                                                                  replacement = "",
                                                                                  x = spatial.maps.path.l)),
                                                      spatial.maps.path.l = spatial.maps.path.l),
                        print_cat = TRUE,
                        message = TRUE,
                        warning = TRUE)
                  
                },
                MoreArgs = list(dir = dir,
                                zone.map.coordinate.correspondence = zone.map.coordinate.correspondence,
                                n.clusters = n.clusters,
                                compute.optimal.k = compute.optimal.k,
                                split.spatial.maps.by = split.spatial.maps.by,
                                cluster.labels = cluster.labels,
                                cluster.colors = cluster.colors),
                SIMPLIFY = FALSE,
                USE.NAMES = TRUE)
       },
       MoreArgs = list(attributes = attributes,
                       zone.map.coordinate.correspondence = zone.map.coordinate.correspondence,
                       n.clusters = n.clusters,
                       split.spatial.maps.by = split.spatial.maps.by,
                       cluster.labels = cluster.labels,
                       cluster.colors = cluster.colors),
       SIMPLIFY = FALSE,
       USE.NAMES = TRUE)
