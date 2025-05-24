# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/")

# Set stringsAsFactors to FALSE
options(stringsAsFactors = FALSE)

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/rcpp_ImputeNearestNeighbors.R")
source(file = "/Users/leonfodoulian/scData/spatial_registration_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/rnascope_claustro_insular_maps_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/discretise.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Compile Rcpp function
Rcpp::sourceCpp(code = rcpp_ImputeNearestNeighbors)

# Required packages
require(patchwork)
require(ggplot2)
require(extrafont)
require(ggtext)
require(ggnewscale)
require(ggrastr)
require(uniqtag)

# Create directory where to save files
dirs <- c("./figures/CLA_map",
          "./figures/CLA_map/violin_plots",
          "./figures/CLA_map/raw_maps",
          "./figures/CLA_map/rotated_maps",
          "./figures/CLA_map/clustered_maps",
          "./figures/CLA_map/2d_density_maps",
          "./figures/CLA_map/hdr_maps",
          "./figures/CLA_map/binned_maps")
lapply(X = dirs,
       FUN = dir.create,
       recursive = TRUE)

# List RDS files
rds.files <- c("results/e20230417_QuPath/C2.Nr4a2_C3.Syt17_C4.Vglut2/all_cells/e20230417_QuPath_C2.Nr4a2_C3.Syt17_C4.Vglut2_log_transformed_hclust_wardD2_clustered_data.rds",
               "results/e20230425_QuPath/C2.Nr4a2_C3.Rprm_C4.Ctgf/all_cells/e20230425_QuPath_C2.Nr4a2_C3.Rprm_C4.Ctgf_log_transformed_hclust_wardD2_clustered_data.rds")
names(x = rds.files) <- gsub(pattern = ".*QuPath_|_log_transformed.*",
                             replacement = "",
                             x = rds.files)

# Define reference points file
# rp.file <- "./e20230417_CLA_map/reference_points/e20230425_Reference_points_20230831_16h43.xlsx"
# rp.file <- "./e20230417_CLA_map/reference_points/e20230425_Reference_points_20230908_18h04.xlsx"
rp.file <- "reference_points/e20230417_e20230425_reference_points_20240222_08h59.xlsx"

# Define cluster names
cluster.names.l <- list("C2.Nr4a2_C3.Syt17_C4.Vglut2" = c("1" = "Vglut2+",
                                                          "2" = "noise",
                                                          "3" = "Nr4a2+",
                                                          "4" = "noise",
                                                          "5" = "Syt17+",
                                                          "6" = "Syt17+;Vglut2+",
                                                          "7" = "Shell",
                                                          "8" = "CLA",
                                                          "not clustered" = "not clustered"),
                        "C2.Nr4a2_C3.Rprm_C4.Ctgf" = c("1" = "Rprm+",
                                                       "2" = "noise",
                                                       "3" = "Ctgf+",
                                                       "4" = "CLA",
                                                       "not clustered" = "not clustered"))

# Define reference images
reference.images <- c("Slide1_zone2_B6_wt.wt_level40_C2.Nr4a2_C3.Rprm_C4.Ctgf.tif",
                      "Slide3_zone1_B6_wt.wt_level44_C2.Nr4a2_C3.Rprm_C4.Ctgf.tif",
                      "Slide5_zone2_B6_wt.wt_level53_C2.Nr4a2_C3.Rprm_C4.Ctgf.tif",
                      "Slide8_zone2_B6_wt.wt_level68_C2.Nr4a2_C3.Rprm_C4.Ctgf.tif")

# Images to exclude
images.to.exclude <- c("Slide8_zone1_B6_wt.wt_level68_C2.Nr4a2_C3.Rprm_C4.Ctgf.tif", # Ctgf+ cells in Rprm+ cells in the medial bottom of the slide (not the case for other slides)
                       "Slide8_zone1_B6_wt.wt_level68_C2.Nr4a2_C3.Syt17_C4.Vglut2.tif") # slide slightly anterior with curved piriform cortex distorting the other slides at level 68

# Define clusters of interest
clusters.of.interest <- c("CLA", "Shell", "Ctgf+", "Rprm+", "Syt17+", "Vglut2+", "Syt17+;Vglut2+")

# Load RNAscope data
rnascope.data.l <- lapply(X = rds.files,
                          readRDS)

# Preprocess RNAscope data
rnascope.data.l <- mapply(
  rnascope.data = rnascope.data.l,
  cluster.names = cluster.names.l,
  cluster.col = list("C2.Nr4a2_C3.Syt17_C4.Vglut2" = "k8_clusters",
                     "C2.Nr4a2_C3.Rprm_C4.Ctgf" = "k4_clusters"),
  cell.suffix = list("C2.Nr4a2_C3.Syt17_C4.Vglut2" = "_1",
                     "C2.Nr4a2_C3.Rprm_C4.Ctgf" = "_2"),
  cluster.to.remove = list("C2.Nr4a2_C3.Syt17_C4.Vglut2" = "Nr4a2+", # remove Nr4a2+ cells from the C2.Nr4a2_C3.Syt17_C4.Vglut2 probe set as they correspond to a mixed population
                           "C2.Nr4a2_C3.Rprm_C4.Ctgf" = NULL),
  FUN = function(rnascope.data,
                 cluster.names,
                 cluster.col,
                 cell.suffix,
                 cluster.to.remove) {
    # Prepare data
    rnascope.data <- rnascope.data$full.data
    # Give unique name to cells
    rnascope.data$cell <- paste0(rnascope.data$cell, cell.suffix)
    # Rename row names of data
    rownames(x = rnascope.data) <- rnascope.data$cell
    # Create a new column corresponding to clustering results
    rnascope.data$clustering_results <- ifelse(test = is.na(x = rnascope.data[[cluster.col]]),
                                               yes = "not clustered",
                                               no = rnascope.data[[cluster.col]])
    # Rename clusters
    rnascope.data$clustering_results <- cluster.names[rnascope.data$clustering_results]
    # Remove clusters from data if necessary
    if (!is.null(x = cluster.to.remove)) {
      rnascope.data <- rnascope.data[!(rnascope.data$clustering_results %in% cluster.to.remove),]
    }
    # Set x coordinates to the right side
    rnascope.data$centroid_x <- ifelse(test = rnascope.data$zone == "zone1",
                                       yes = rnascope.data$centroid_x,
                                       no = -rnascope.data$centroid_x)
    # Return data
    return(list(all.data = rnascope.data,
                spatial.data = rnascope.data[c("cell", "slide", "zone", "mouse", "genotype", "level", "centroid_x", "centroid_y", "is_clustered", "clustering_results")]))
  },
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE)

# Store all data in another object for gene expression figure generation
all.data <- lapply(X = rnascope.data.l,
                   FUN = "[[",
                   "all.data")

# Bind spatial data by rows
rnascope.data <- dplyr::bind_rows(lapply(X = rnascope.data.l,
                                         FUN = "[[",
                                         "spatial.data"),
                                  .id = "probe_set")

# Define image name for each cell
rnascope.data$image <- paste(rnascope.data$slide,
                             rnascope.data$zone,
                             rnascope.data$mouse,
                             rnascope.data$genotype,
                             rnascope.data$level,
                             paste0(rnascope.data$probe_set,
                                    ".tif"),
                             sep = "_")

# Load reference points data
reference.points <- readxl::read_xlsx(path = rp.file,
                                      na = "NaN")

# Get meta data for each reference point
rp.meta.data <- dplyr::bind_rows(lapply(
  X = strsplit(x = reference.points$Image,
               split = "_"),
  FUN = function(meta.data) {
    # Group probe set info
    meta.data <- c(meta.data[1:5],
                   paste(meta.data[6:8],
                         collapse = "_"))
    # Remove ".tif" from probe set info
    meta.data <- gsub(pattern = "\\.tif",
                      replacement = "",
                      x = meta.data)
    # Name meta data entries
    names(x = meta.data) <- c("slide",
                              "zone",
                              "mouse",
                              "genotype",
                              "level",
                              "probe_set")
    # Return meta data
    return(meta.data)
  }),
  .id = NULL)

# Bind reference points data to meta data
reference.points <- dplyr::bind_cols(reference.points,
                                     rp.meta.data)

# Subset reference points data to columns of interest
reference.points <- as.data.frame(x = reference.points[c("Image", "slide", "zone", "mouse", "genotype", "level", "probe_set", "Name", "Centroid X µm", "Centroid Y µm")])
colnames(x = reference.points) <- c("image", "slide", "zone", "mouse", "genotype", "level", "probe_set", "reference_points", "centroid_x", "centroid_y")

# Re-structure reference points data in wide format
reference.points <- data.table::dcast(data = data.table::setDT(reference.points),
                                      formula = image + slide + zone + mouse + genotype + level + probe_set ~ reference_points,
                                      value.var = c("centroid_x", "centroid_y"))

# Define whether image corresponds to reference image
reference.points$is_reference_image <- reference.points$image %in% reference.images

# Define dorsal and ventral coordinates for image rotation
reference.points$centroid_x_dorsal <- ifelse(test = reference.points$level != "level68",
                                             yes = reference.points$centroid_x_EPd,
                                             no = reference.points$centroid_x_CLAd)
reference.points$centroid_y_dorsal <- ifelse(test = reference.points$level != "level68",
                                             yes = reference.points$centroid_y_EPd,
                                             no = reference.points$centroid_y_CLAd)
reference.points$centroid_x_ventral <- ifelse(test = reference.points$level != "level68",
                                              yes = reference.points$centroid_x_EPv,
                                              no = reference.points$centroid_x_CLAc)
reference.points$centroid_y_ventral <- ifelse(test = reference.points$level != "level68",
                                              yes = reference.points$centroid_y_EPv,
                                              no = reference.points$centroid_y_CLAc)

# Set reference points x coordinates to the right side
reference.points$centroid_x_center <- ifelse(test = reference.points$zone == "zone1",
                                             yes = reference.points$centroid_x_center,
                                             no = -reference.points$centroid_x_center)
reference.points$centroid_x_dorsal <- ifelse(test = reference.points$zone == "zone1",
                                             yes = reference.points$centroid_x_dorsal,
                                             no = -reference.points$centroid_x_dorsal)
reference.points$centroid_x_ventral <- ifelse(test = reference.points$zone == "zone1",
                                              yes = reference.points$centroid_x_ventral,
                                              no = -reference.points$centroid_x_ventral)

# Center coordinates to reference point in image
rnascope.data.l <- mapply(
  rnascope.data = split(x = rnascope.data,
                        f = rnascope.data$image),
  reference.data = split(x = reference.points,
                         f = reference.points$image),
  FUN = function(rnascope.data,
                 reference.data) {
    # Center coordinates to reference point in image
    rnascope.data$centroid_x <- rnascope.data$centroid_x - reference.data$centroid_x_center
    rnascope.data$centroid_y <- rnascope.data$centroid_y - reference.data$centroid_y_center
    
    # Return data
    return(rnascope.data)
  },
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE)

# Get list of angles between query and reference images
rotation.angle.l <- unlist(x = mapply(
  query.coords = split(x = reference.points,
                       f = reference.points$level),
  reference.coords = split(x = reference.points[reference.points$is_reference_image,],
                           f = reference.points[reference.points$is_reference_image,]$level),
  FUN = function(query.coords,
                 reference.coords) {
    # Compute angles between query and reference images
    unlist(x = lapply(
      X = split(x = query.coords,
                f = query.coords$image),
      FUN = ComputeRotationAngle,
      reference.coords = reference.coords,
      x1.name = "centroid_x_dorsal",
      x2.name = "centroid_x_ventral",
      y1.name = "centroid_y_dorsal",
      y2.name = "centroid_y_ventral"))
  },
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE))

# Rotate coordinates using rotation angle values
rnascope.data.l <- mapply(
  rnascope.data = rnascope.data.l,
  rotation.angle = rotation.angle.l[names(x = rnascope.data.l)],
  FUN = function(rnascope.data,
                 rotation.angle) {
    # Rotate coordinates
    rotated.coords <- RotateCoordinates(xy = rnascope.data[, c("centroid_x", "centroid_y")],
                                        angle = rotation.angle)
    
    # Rename columns of rotated coordinates
    colnames(x = rotated.coords) <- c("centroid_x_rot", "centroid_y_rot")
    
    # Bind rotated coordinates to image data
    rnascope.data <- cbind(rnascope.data,
                           rotated.coords)
    
    # Add rotation angle to image data
    rnascope.data$rotation_angle <- rotation.angle
    
    # Return data
    return(rnascope.data)
  },
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE)

# Bind data by rows
rnascope.data <- dplyr::bind_rows(rnascope.data.l,
                                  .id = NULL)

# Remove images to exclude (i.e. not in correspondence with other set of images)
rnascope.data <- rnascope.data[!(rnascope.data$image %in% images.to.exclude),]

# Indicate whether cell is from reference image
rnascope.data$is_reference_image <- rnascope.data$image %in% reference.images

# Define image type
rnascope.data$image_type <- ifelse(test = rnascope.data$is_reference_image,
                                   yes = "reference",
                                   no = "query")

# Bin x and y centroids in 20 micrometer bins
rnascope.data$centroid_x_disc <- discretise(x = rnascope.data$centroid_x_rot,
                                            by = 20)
rnascope.data$centroid_y_disc <- discretise(x = rnascope.data$centroid_y_rot,
                                            by = 20)

# Get axis attributes
axis.attributes <- GetAxisAttributes(x = c(rnascope.data$centroid_x, rnascope.data$centroid_x_rot), # no need to get one for clusters of interest as the change is minor
                                     y = c(rnascope.data$centroid_y, rnascope.data$centroid_y_rot),
                                     bin.size = 20)

# Define alpha values for color mixing of imputed clusters
alpha.levels <- c(1, seq(from = 0.85,
                         to = 0.5,
                         by = -0.05))
names(x = alpha.levels) <- c(0:(length(x = alpha.levels) - 2), # see code of binned map for more details
                             paste0(length(x = alpha.levels) - 1, "+"))

# Create a sham clustering column for background cell plotting
rnascope.data$sham_clustering_results <- "background"

# Define specific parameters for each plot type
figure.attributes.l <- list("plot" = list("raw_maps" = list(map.type = "raw",
                                                            centroid.x.column = "centroid_x",
                                                            centroid.y.column = "centroid_y",
                                                            clustering.column = "image_type",
                                                            cluster.colors = c("reference" = "#6C7B8B"),
                                                            cluster.breaks = c("reference" = "reference"),
                                                            legend.title = "image type (# cells)",
                                                            path = "./figures/CLA_map/raw_maps"),
                                          "rotated_maps" = list(map.type = "rotated",
                                                                centroid.x.column = "centroid_x_rot",
                                                                centroid.y.column = "centroid_y_rot",
                                                                clustering.column = "image_type",
                                                                cluster.colors = c("reference" = "#6C7B8B"),
                                                                cluster.breaks = c("reference" = "reference"),
                                                                legend.title = "image type (# cells)",
                                                                path = "./figures/CLA_map/rotated_maps"),
                                          "clustered_maps" = list(map.type = "clustered",
                                                                  centroid.x.column = "centroid_x_rot",
                                                                  centroid.y.column = "centroid_y_rot",
                                                                  clustering.column = "clustering_results",
                                                                  cluster.colors = c("CLA" = "#4B0082",
                                                                                     "Shell" = "#FF4D00",
                                                                                     "Ctgf+" = "#66CDAA",
                                                                                     "Rprm+" = "#1E90FF",
                                                                                     "Syt17+" = "#D0D2D3",
                                                                                     "Vglut2+" = "#6D6E70",
                                                                                     "Syt17+;Vglut2+" = "#A6A8AB",
                                                                                     "noise" = "#E7DECC",
                                                                                     "not clustered" = "#F0E5D3"),
                                                                  cluster.breaks = c("CLA" = "CLA",
                                                                                     "Shell" = "Shell",
                                                                                     "Ctgf+" = "Ctgf+",
                                                                                     "Rprm+" = "Rprm+",
                                                                                     "Syt17+" = "Syt17+",
                                                                                     "Vglut2+" = "Vglut2+",
                                                                                     "Syt17+;Vglut2+" = "Syt17+;Vglut2+",
                                                                                     "noise" = "noise",
                                                                                     "not clustered" = "not clustered"),
                                                                  cluster.labels = c("CLA" = "CLA",
                                                                                     "Shell" = "Shell",
                                                                                     "Ctgf+" = "*Ctgf*<sup>+</sup>",
                                                                                     "Rprm+" = "*Rprm*<sup>+</sup>",
                                                                                     "Syt17+" = "*Syt17*<sup>+</sup>",
                                                                                     "Vglut2+" = "*Vglut2*<sup>+</sup>",
                                                                                     "Syt17+;Vglut2+" = "*Syt17*<sup>+</sup>;*Vglut2*<sup>+</sup>",
                                                                                     "noise" = "noise",
                                                                                     "not clustered" = "not clustered"),
                                                                  legend.title = "neuronal clusters (# cells)",
                                                                  alpha.levels = alpha.levels,
                                                                  path = "./figures/CLA_map/clustered_maps"),
                                          "2d_density_maps" = list(map.type = "2d density",
                                                                   centroid.x.column = "centroid_x_rot",
                                                                   centroid.y.column = "centroid_y_rot",
                                                                   clustering.column = "sham_clustering_results",
                                                                   background.clusters = c("CLA", "Shell", "Ctgf+", "Rprm+", "Syt17+", "Vglut2+", "Syt17+;Vglut2+"),
                                                                   cluster.colors = c("background" = "#D3D3D3"),
                                                                   cluster.breaks = c("background" = "background"),
                                                                   cluster.labels = c("background" = "background"),
                                                                   density.column = "clustering_results",
                                                                   density.clusters = c("CLA", "Shell", "Ctgf+", "Rprm+"),
                                                                   density.colors = c("CLA" = "#4B0082",
                                                                                      "Shell" = "#FF4D00",
                                                                                      "Ctgf+" = "#66CDAA",
                                                                                      "Rprm+" = "#1E90FF"),
                                                                   density.breaks = c("CLA" = "CLA",
                                                                                      "Shell" = "Shell",
                                                                                      "Ctgf+" = "Ctgf+",
                                                                                      "Rprm+" = "Rprm+"),
                                                                   denssity.labels = c("CLA" = "CLA",
                                                                                       "Shell" = "Shell",
                                                                                       "Ctgf+" = "*Ctgf*<sup>+</sup>",
                                                                                       "Rprm+" = "*Rprm*<sup>+</sup>"),
                                                                   legend.title = "neuronal clusters",
                                                                   path = "./figures/CLA_map/2d_density_maps"),
                                          "hdr_maps" = list(map.type = "hdr",
                                                            hdr.probs = setNames(object = seq(from = 0.05,
                                                                                              to = 0.95,
                                                                                              by = 0.05),
                                                                                 nm = paste0("prob_",
                                                                                             format(x = seq(from = 0.05,
                                                                                                            to = 0.95,
                                                                                                            by = 0.05),
                                                                                                    nsmall = 2))),
                                                            centroid.x.column = "centroid_x_rot",
                                                            centroid.y.column = "centroid_y_rot",
                                                            clustering.column = "sham_clustering_results",
                                                            background.clusters = c("CLA", "Shell", "Ctgf+", "Rprm+", "Syt17+", "Vglut2+", "Syt17+;Vglut2+"),
                                                            cluster.colors = c("background" = "#D3D3D3"),
                                                            cluster.breaks = c("background" = "background"),
                                                            cluster.labels = c("background" = "background"),
                                                            density.column = "clustering_results",
                                                            density.clusters = c("CLA", "Shell", "Ctgf+", "Rprm+"),
                                                            density.colors = c("CLA" = "#4B0082",
                                                                               "Shell" = "#FF4D00",
                                                                               "Ctgf+" = "#66CDAA",
                                                                               "Rprm+" = "#1E90FF"),
                                                            density.breaks = c("CLA" = "CLA",
                                                                               "Shell" = "Shell",
                                                                               "Ctgf+" = "Ctgf+",
                                                                               "Rprm+" = "Rprm+"),
                                                            denssity.labels = c("CLA" = "CLA",
                                                                                "Shell" = "Shell",
                                                                                "Ctgf+" = "*Ctgf*<sup>+</sup>",
                                                                                "Rprm+" = "*Rprm*<sup>+</sup>"),
                                                            legend.title = "neuronal clusters",
                                                            path = "./figures/CLA_map/hdr_maps")),
                            "axis" = axis.attributes)

####################################################
# Figure: Plot violin plots with final cluster labels
# Plot and save violin plots
vln.plots.l <- mapply(
  expression.data = all.data,
  probe.set = names(x = all.data),
  channels = list("C2.Nr4a2_C3.Syt17_C4.Vglut2" = c("Nr4a2", "Syt17", "Vglut2"),
                  "C2.Nr4a2_C3.Rprm_C4.Ctgf" = c("Nr4a2", "Ctgf", "Rprm")),
  FUN = function(expression.data,
                 probe.set,
                 channels,
                 cluster.colors,
                 cluster.breaks,
                 cluster.labels) {
    # Remove cells that were not clustered
    expression.data <- expression.data[expression.data$clustering_results != "not clustered",]
    
    # Re-define cluster colors, breaks and labels (in case is needed)
    cluster.breaks <- cluster.breaks[cluster.breaks %in% expression.data$clustering_results]
    cluster.colors <- cluster.colors[cluster.breaks]
    cluster.labels <- cluster.labels[cluster.breaks]
    
    # Define y scale transformation
    y.scale.trans <- c("identity", "log10")
    names(x = y.scale.trans) <- y.scale.trans
    
    # Plot violin plots for each y scale transformation
    vln.plots.l <- lapply(
      X = y.scale.trans,
      FUN = function(y.scale.trans,
                     expression.data,
                     probe.set,
                     channels,
                     cluster.colors,
                     cluster.breaks,
                     cluster.labels) {
        # Plot violin plots
        vln.plots <- StackViolinPlots(
          data.use = t(x = expression.data[,channels]),
          genes.use = channels,
          cluster.data = expression.data[,"clustering_results", drop = FALSE],
          cellnames.as.rownames.in.cluster.data = TRUE,
          cluster.ident.name = "clustering_results",
          cluster.colors = cluster.colors,
          cluster.breaks = cluster.breaks,
          cluster.labels = cluster.labels,
          is.log.transformed = FALSE,
          log.scale = "log",
          pseudocount.use = 1,
          y.scale.trans = y.scale.trans,
          point.size = 2,
          alpha.use = 1,
          vln.border.colour = "black",
          vln.border.stroke = 0,
          plot.title.size = 5,
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
        
        # Adjust y axis labels in patchwork
        vln.plots$patches$plots[[1]] <- vln.plots$patches$plots[[1]] +
          theme(axis.text.y = element_markdown(size = 5,
                                               family = "Arial",
                                               face = "plain"))
        
        # Save plot in pdf format
        FixSizeAndSave(plot = vln.plots,
                       filename = paste("figures/CLA_map/violin_plots/e20230417_e20230425_CLA_map",
                                        probe.set,
                                        "violin_plots",
                                        y.scale.trans,
                                        "scale.pdf",
                                        sep = "_"),
                       is.ggassemble = TRUE,
                       panel.width = 0.5,
                       panel.height = 0.325 * length(x = cluster.breaks),
                       margin = 2,
                       unit.use = "cm",
                       use.ggsave = TRUE,
                       useDingbats = FALSE)
        
        # Return plot
        return(vln.plots)
      },
      expression.data = expression.data,
      probe.set = probe.set,
      channels = channels,
      cluster.colors = cluster.colors,
      cluster.breaks = cluster.breaks,
      cluster.labels = cluster.labels)
    
    # Return plots
    return(vln.plots.l)
  },
  MoreArgs = list(cluster.colors = figure.attributes.l$plot$clustered_maps$cluster.colors,
                  cluster.breaks = figure.attributes.l$plot$clustered_maps$cluster.breaks,
                  cluster.labels = figure.attributes.l$plot$clustered_maps$cluster.labels),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE)

####################################################
# Figure: Plot spatial maps before and after rotation with image and cluster labels
# Plot spatial maps
spatial.map.plots <- lapply(
  X = split(x = rnascope.data,
            f = rnascope.data$level),
  FUN = function(plot.data,
                 plot.attributes.l,
                 axis.attributes,
                 uniform.scale,
                 save.plot) {
    # Plot various spatial maps
    spatial.map.plots <- lapply(
      X = plot.attributes.l,
      FUN = function(plot.attributes,
                     plot.data,
                     axis.attributes,
                     uniform.scale,
                     save.plot) {
        # Define slice level
        slice.level <- unique(x = plot.data$level)
        
        # Define number of images
        n.images <- length(x = unique(x = plot.data$image))
        
        # Preprocess data
        if (plot.attributes$map.type %in% c("raw", "rotated")) {
          # Make queries unique
          image.types <- plot.data[!duplicated(x = plot.data$image),]$image_type
          names(x = image.types) <- plot.data[!duplicated(x = plot.data$image),]$image
          image.types <- uniqtag::make_unique(xs = image.types, sep = " ")
          
          # Add unique queries to data
          plot.data$image_type <- image.types[plot.data$image]
          
          # Number of query data
          n.query <- length(x = image.types) - 1
          
          # Define query colors
          query.colors <- colorRampPalette(colors = c("#9FB6CD", "#C6E2FF"))(n.query)
          names(x = query.colors) <- paste("query",
                                           1:n.query,
                                           sep = " ")
          
          # Define query breaks
          query.breaks <- names(x = query.colors)
          names(x = query.breaks) <- query.breaks
          
          # Define cluster colors
          plot.attributes$cluster.colors <- c(plot.attributes$cluster.colors,
                                              query.colors)
          
          # Define cluster breaks
          plot.attributes$cluster.breaks <- c(plot.attributes$cluster.breaks,
                                              query.breaks)
          
          # Define cluster labels
          plot.attributes$cluster.labels <- plot.attributes$cluster.breaks
        } else if (plot.attributes$map.type %in% c("2d density", "hdr")) {
          # Subset plot data for background clusters
          plot.data <- plot.data[plot.data[[plot.attributes$density.column]] %in% plot.attributes$background.clusters,]
          
          # Create a new dataset for density maps
          density.data <- plot.data[plot.data[[plot.attributes$density.column]] %in% plot.attributes$density.clusters,]
          
          # Split density dataset by density clusters
          density.data.l <- split(x = density.data,
                                  f = density.data[[plot.attributes$density.column]]) # only for 2d density plots
          
          # Order data by order of clusters
          density.data.l <- density.data.l[plot.attributes$density.clusters] # only for 2d density plots
          
          # Define order of overlay for hdr plots
          density.data[[plot.attributes$density.column]] <- factor(x = density.data[[plot.attributes$density.column]],
                                                                   levels = rev(x = plot.attributes$density.clusters)) # only for hdr plots
        }
        
        # Remove uniform scale if clustered or 2d density maps for level 44 are plotted
        # if (slice.level == "level44" && plot.attributes$map.type %in% c("clustered", "2d density")) {
        #   uniform.scale <- FALSE
        # }
        
        if (!uniform.scale) {
          # Get axis specifics for subsetted data
          axis.attributes <- GetAxisAttributes(x = plot.data[[plot.attributes$centroid.x.column]],
                                               y = plot.data[[plot.attributes$centroid.y.column]],
                                               bin.size = 20)
        }
        
        # Plot spatial map
        spatial.map.plot <- Plot2DEmbedding(
          embedding.data = plot.data[,c(plot.attributes$centroid.x.column, plot.attributes$centroid.y.column)],
          cluster.data = plot.data[,plot.attributes$clustering.column, drop = FALSE],
          cellnames.as.rownames.in.embedding.data = TRUE,
          cellnames.as.rownames.in.cluster.data = TRUE,
          embedding.x.name = plot.attributes$centroid.x.column,
          embedding.y.name = plot.attributes$centroid.y.column,
          embedding.x.label = "lateral",
          embedding.y.label = "dorsal",
          cluster.ident.name = plot.attributes$clustering.column,
          cluster.colors = plot.attributes$cluster.colors,
          cluster.breaks = plot.attributes$cluster.breaks,
          cluster.labels = plot.attributes$cluster.labels,
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
          arrow.type = "closed",
          axis.lineend = "butt",
          axis.gap = 0,
          range.scale = 0.05
        )
        
        # Adjust axes to standardize across all plots
        spatial.map.plot <- spatial.map.plot +
          scale_x_continuous(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                             breaks = c(axis.attributes$min.x,
                                        axis.attributes$min.x + 0.1 * axis.attributes$range.max)) + # equal arrow size for x and y axes after plot resizing
          scale_y_reverse(expand = c(0, 0), # equal spacing for x and y axes after plot resizing
                          breaks = c(axis.attributes$max.y - 0.1 * axis.attributes$range.max, # equal arrow size for x and y axes after plot resizing
                                     axis.attributes$max.y)) +
          expand_limits(x = c(axis.attributes$min.x - 0.025 * axis.attributes$range.max,
                              axis.attributes$max.x),
                        y = c(axis.attributes$min.y,
                              axis.attributes$max.y + 0.025 * axis.attributes$range.max)) +
          labs(title = paste(plot.attributes$map.type,
                             " map\n(",
                             gsub(pattern = "level",
                                  replacement = "level ",
                                  x = slice.level),
                             "; n = ",
                             n.images,
                             " images)",
                             sep = ""),
               colour = plot.attributes$legend.title) +
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
        
        # Plot clustered data twice - with and without noisy and not clustered cells
        if (plot.attributes$map.type %in% c("raw", "rotated")) {
          # Define data filtering type
          data.filtering.type <- "unfiltered"
          
          # Reorder data in plot
          spatial.map.plot$data <- spatial.map.plot$data[order(x = spatial.map.plot$data$image_type, decreasing = FALSE),]
          
          # Prepare plot for saving
          spatial.map.plot.l <- list("unfiltered" = spatial.map.plot)
        } else if (plot.attributes$map.type %in% "clustered") {
          # Define data filtering type
          data.filtering.type <- c("unfiltered", "filtered")
          
          # Reorder data in plot
          spatial.map.plot$data[[plot.attributes$clustering.column]] <- factor(x = spatial.map.plot$data[[plot.attributes$clustering.column]],
                                                                               levels = plot.attributes$cluster.breaks)
          spatial.map.plot$data <- spatial.map.plot$data[order(x = spatial.map.plot$data[[plot.attributes$clustering.column]], decreasing = TRUE),]
          
          # Prepare plot for saving
          spatial.map.plot.l <- list("unfiltered" = spatial.map.plot,
                                     "filtered" = {
                                       spatial.map.plot$data <- spatial.map.plot$data[!spatial.map.plot$data[[plot.attributes$clustering.column]] %in% c("noise", "not clustered"),]
                                       spatial.map.plot
                                     })
        } else if (plot.attributes$map.type == "2d density") {
          # Define data filtering type
          data.filtering.type <- "filtered"
          
          # New Plot2DEmbedding() creates an alpha scale that interfere with alpha scale of stat_density_2d_filled
          spatial.map.plot$layers[[1]]$aes_params$alpha <- 1 # remove existing alpha scale before producing graphs
          spatial.map.plot$scales$scales[[2]] <- NULL # remove existing alpha scale before producing graphs
          
          # Plot 2d density on spatial map
          spatial.map.plot.l <- lapply(X = density.data.l,
                                       FUN = function(density.data,
                                                      spatial.map.plot,
                                                      plot.attributes) {
                                         spatial.map.plot <- ggrastr::rasterise(input = spatial.map.plot,
                                                                                dpi = 300) +
                                           guides(colour = "none") +
                                           ggnewscale::new_scale_fill() +
                                           ggnewscale::new_scale_colour() + # or add plot.attributes$cluster.colors etc. to scale_colour_manual
                                           stat_density_2d_filled(data = density.data,
                                                                  mapping = aes(alpha = after_stat(x = level),
                                                                                colour = .data[[plot.attributes$density.column]],
                                                                                fill = .data[[plot.attributes$density.column]]),
                                                                  contour = TRUE,
                                                                  contour_var = "ndensity",
                                                                  n = 100,
                                                                  h = NULL, # bandwidth estimated using MASS::bandwidth.nrd()
                                                                  adjust = c(1, 1),
                                                                  breaks = seq(from = 0.1,
                                                                               to = 1,
                                                                               by = 0.1)) +
                                           scale_colour_manual(values = plot.attributes$density.colors,
                                                               breaks = plot.attributes$density.breaks,
                                                               labels = plot.attributes$denssity.labels) +
                                           scale_fill_manual(values = plot.attributes$density.colors,
                                                             breaks = plot.attributes$density.breaks,
                                                             labels = plot.attributes$denssity.labels) +
                                           guides(colour = "none",
                                                  alpha = "none") +
                                           labs(title = NULL,
                                                colour = NULL,
                                                fill = NULL) +
                                           theme(plot.margin = margin(t = 0,
                                                                      r = 0,
                                                                      b = 0,
                                                                      l = 0,
                                                                      unit = "mm"))
                                       },
                                       spatial.map.plot = spatial.map.plot,
                                       plot.attributes = plot.attributes)
          
          # Combine plots using patchwork
          spatial.map.plot.l <- list("filtered" = patchwork::wrap_plots(spatial.map.plot.l,
                                                                        nrow = 1,
                                                                        guides = "collect") +
                                       patchwork::plot_annotation(title = paste(plot.attributes$map.type,
                                                                                " maps\n(",
                                                                                gsub(pattern = "level",
                                                                                     replacement = "level ",
                                                                                     x = slice.level),
                                                                                "; n = ",
                                                                                n.images,
                                                                                " images)",
                                                                                sep = ""),
                                                                  theme = theme(plot.title = element_text(size = 6,
                                                                                                          family = "Arial",
                                                                                                          colour = "black",
                                                                                                          hjust = 0.5,
                                                                                                          vjust = 0,
                                                                                                          margin = margin(t = 0,
                                                                                                                          r = 0,
                                                                                                                          b = 1,
                                                                                                                          l = 0,
                                                                                                                          unit = "mm")),
                                                                                legend.position = "bottom",
                                                                                legend.box = "horizontal",
                                                                                legend.direction = "horizontal")))
        } else if (plot.attributes$map.type == "hdr") {
          # Define data filtering type
          data.filtering.type <- "filtered"
          
          # Create scale bar data
          scale.bar <- list(scale = data.frame(x = 0.99 * axis.attributes$max.x - 500, # scale of 500 micrometer
                                               xend = 0.99 * axis.attributes$max.x,
                                               y = 0.99 * axis.attributes$max.y,
                                               yend = 0.99 * axis.attributes$max.y),
                            label = data.frame(x = 0.99 * axis.attributes$max.x - 250, # centered at 250
                                               y = 0.95 * axis.attributes$max.y,
                                               label = "500 &mu;m"))
          
          # Plot hdr on spatial map
          spatial.map.plot.l <- lapply(X = plot.attributes$hdr.probs,
                                       FUN = function(hdr.prob,
                                                      density.data,
                                                      spatial.map.plot,
                                                      plot.attributes,
                                                      scale.bar,
                                                      slice.level,
                                                      n.images) {
                                         spatial.map.plot <- ggrastr::rasterise(input = spatial.map.plot,
                                                                                dpi = 300) +
                                           geom_segment(data = scale.bar$scale,
                                                        mapping = aes(x = x,
                                                                      xend = xend,
                                                                      y = y,
                                                                      yend = yend),
                                                        colour = "black",
                                                        size = 0.25,
                                                        inherit.aes = FALSE) +
                                           # ggtext::geom_richtext(data = scale.bar$label,
                                           #                       mapping = aes(x = x,
                                           #                                     y = y,
                                           #                                     label = label),
                                           #                       size = 5 * (1/72 * 25.4),
                                           #                       hjust = 0.5,
                                           #                       vjust = 0.5,
                                           #                       fill = NA,
                                           #                       label.color = NA,
                                           #                       inherit.aes = FALSE) +
                                           guides(colour = "none") +
                                           ggnewscale::new_scale_fill() +
                                           ggnewscale::new_scale_colour() + # or add plot.attributes$cluster.colors etc. to scale_colour_manual using
                                           ggdensity::stat_hdr(data = density.data,
                                                               mapping = aes(fill = .data[[plot.attributes$density.column]]),
                                                               method = "kde",
                                                               probs = hdr.prob,
                                                               n = 100,
                                                               linewidth = 0,
                                                               alpha = 0.75) +
                                           scale_colour_manual(values = plot.attributes$density.colors,
                                                               breaks = plot.attributes$density.breaks,
                                                               labels = plot.attributes$denssity.labels) +
                                           scale_fill_manual(values = plot.attributes$density.colors,
                                                             breaks = plot.attributes$density.breaks,
                                                             labels = plot.attributes$denssity.labels) +
                                           guides(colour = "none",
                                                  alpha = "none",
                                                  fill = guide_legend(override.aes = list(alpha = 1))) +
                                           labs(title = paste(plot.attributes$map.type,
                                                              " maps (prob = ",
                                                              format(x = hdr.prob,
                                                                     nsmall = 2),
                                                              ")\n(",
                                                              gsub(pattern = "level",
                                                                   replacement = "level ",
                                                                   x = slice.level),
                                                              "; n = ",
                                                              n.images,
                                                              " images)",
                                                              sep = ""),
                                                colour = NULL,
                                                fill = NULL)
                                       },
                                       spatial.map.plot = spatial.map.plot,
                                       density.data = density.data,
                                       plot.attributes = plot.attributes,
                                       scale.bar = scale.bar,
                                       slice.level = slice.level,
                                       n.images = n.images)
        }
        
        # Define file name for list of plot(s)
        file.name.l <- file.path(plot.attributes$path,
                                 paste("e20230417_e20230425_CLA_map_",
                                       data.filtering.type,
                                       "_data_",
                                       slice.level,
                                       "_",
                                       n.images,
                                       "_images_",
                                       gsub(pattern = "\\s",
                                            replacement = "_",
                                            x = plot.attributes$map.type),
                                       {
                                         if (!any(names(x = spatial.map.plot.l) %in% c("unfiltered", "filtered"))) {
                                           paste0("_",
                                                  names(x = spatial.map.plot.l))
                                         }
                                       },
                                       "_spatial_maps_",
                                       {
                                         if (uniform.scale) {
                                           "uniform_scale"
                                         } else {
                                           "nonuniform_scale"
                                         }
                                       },
                                       ".pdf",
                                       sep = ""))
        
        # Save plot in pdf format
        if (save.plot) {
          mapply(spatial.map.plot = spatial.map.plot.l,
                 file.name = file.name.l,
                 FUN = function(spatial.map.plot,
                                file.name,
                                range.ratio) {
                   # Save plot in pdf format
                   FixSizeAndSave(plot = spatial.map.plot,
                                  filename = file.name,
                                  is.ggassemble = "patchwork" %in% class(x = spatial.map.plot),
                                  # panel.width = min(7, 7 * range.ratio),
                                  # panel.height = min(7, 7 / range.ratio),
                                  panel.width = 4.2, # figure limit for width
                                  panel.height = 4.2 / range.ratio,
                                  unit.use = "cm",
                                  margin = 0,
                                  use.ggsave = TRUE,
                                  useDingbats = FALSE)
                   
                   # Return "DONE"
                   return("DONE")
                 },
                 MoreArgs = list(range.ratio = axis.attributes$range.ratio),
                 SIMPLIFY = FALSE,
                 USE.NAMES = TRUE)
        }
        
        # Return plots
        return(spatial.map.plot.l)
      },
      plot.data = plot.data,
      axis.attributes = axis.attributes,
      uniform.scale = uniform.scale,
      save.plot = save.plot)
    
    # Return plots
    return(spatial.map.plots)
  },
  plot.attributes.l = figure.attributes.l$plot,
  axis.attributes = figure.attributes.l$axis,
  uniform.scale = TRUE,
  save.plot = TRUE)

####################################################
# Figure: Plot binned spatial maps       ####################### Find a way to add axis on the data
# Plot binned claustro-insular map of broad clusters
binned.spatial.maps <- lapply(
  X = split(x = rnascope.data,
            f = rnascope.data$level),
  FUN = function(plot.data,
                 clusters.of.interest,
                 cluster.colors,
                 cluster.breaks,
                 cluster.labels,
                 alpha.levels,
                 axis.attributes,
                 impute.clusters,
                 uniform.scale,
                 save.plot,
                 path) {
    # Define slice level
    slice.level <- unique(x = plot.data$level)
    
    # Define number of images
    n.images <- length(x = unique(x = plot.data$image))
    
    if (!uniform.scale) {
      # Bin x and y centroids in 20 micrometer bins
      plot.data$centroid_x_disc <- discretise(x = plot.data$centroid_x_rot,
                                              by = 20) # include.lowest = TRUE, so the first level is not the same
      plot.data$centroid_y_disc <- discretise(x = plot.data$centroid_y_rot,
                                              by = 20) # include.lowest = TRUE, so the first level is not the same
      
      # Get axis specifics for subsetted data
      axis.attributes <- GetAxisAttributes(x = plot.data[["centroid_x_rot"]],
                                           y = plot.data[["centroid_y_rot"]],
                                           bin.size = 20)
    }
    
    # Set data as data.table and subset data for clusters of interest
    plot.data <- data.table::data.table(plot.data[plot.data$clustering_results %in% clusters.of.interest,])
    
    # Define proportion columns (used for cluster identity identification)
    prop.cols <- paste("prop",
                       gsub(pattern = "\\+|;",
                            replacement = "",
                            x = tolower(x = clusters.of.interest)),
                       sep = "_")
    
    # Define proportions to compute
    props.to.compute <- paste0("sum(clustering_results == '",
                               clusters.of.interest,
                               "') / length(x = clustering_results)")
    names(x = props.to.compute) <- prop.cols
    
    # Compute proportions of each cluster per bin
    # plot.data <- plot.data[, list(prop_cla = sum(clustering_results == "CLA") / length(x = clustering_results),
    #                               prop_shell = sum(clustering_results == "Shell") / length(x = clustering_results),
    #                               prop_ctgf = sum(clustering_results == "Ctgf+") / length(x = clustering_results),
    #                               prop_rprm = sum(clustering_results == "Rprm+") / length(x = clustering_results),
    #                               prop_syt17 = sum(clustering_results == "Syt17+") / length(x = clustering_results),
    #                               prop_vglut2 = sum(clustering_results == "Vglut2+") / length(x = clustering_results),
    #                               prop_syt17vglut2 = sum(clustering_results == "Syt17+;Vglut2+") / length(x = clustering_results)),
    #                        by = c("centroid_x_disc", "centroid_y_disc", "level")]
    plot.data <- plot.data[, lapply(X = props.to.compute,
                                    FUN = function(prop.to.compute) {
                                      eval(expr = parse(text = prop.to.compute))
                                    }),
                           by = c("centroid_x_disc", "centroid_y_disc", "level")]
    
    # Identify cluster identity of each bin
    plot.data <- plot.data[,
                           c("bin_cluster_id",
                             "bin_cluster_label",
                             "bin_cluster_color",
                             "bin_order") := GetBinAttributes(row.data = .SD,
                                                              clusters.of.interest = clusters.of.interest,
                                                              cluster.labels = cluster.labels,
                                                              cluster.colors = cluster.colors,
                                                              n.pad = length(x = clusters.of.interest)),
                           by = 1:nrow(x = plot.data),
                           .SDcols = prop.cols]
    
    # Get bin attributes
    bin.attributes <- plot.data[!is.na(x = plot.data$bin_cluster_id) & !duplicated(x = plot.data$bin_cluster_id),
                                c("bin_cluster_id", "bin_cluster_label", "bin_cluster_color", "bin_order")]
    
    # Define bin cluster breaks, labels and colors
    bin.cluster.breaks <- bin.attributes$bin_cluster_id
    names(x = bin.cluster.breaks) <- bin.cluster.breaks
    bin.cluster.labels <- bin.attributes$bin_cluster_label
    names(x = bin.cluster.labels) <- bin.cluster.breaks
    bin.cluster.colors <- bin.attributes$bin_cluster_color
    names(x = bin.cluster.colors) <- bin.cluster.breaks
    
    # Define column used for fill aesthetics
    cluster.cols <- c("pre_imputation" = "bin_cluster_id")
    
    if (impute.clusters) {
      # Get x and y bin levels from corresponding data (bins in plot data are discontinuous and can't be used as is)
      x.bin.levels.in.data <- which(x = levels(x = plot.data$centroid_x_disc) %in% plot.data$centroid_x_disc)
      y.bin.levels.in.data <- which(x = levels(x = plot.data$centroid_y_disc) %in% plot.data$centroid_y_disc)
      x.bin.levels <- levels(x = plot.data$centroid_x_disc)[min(x.bin.levels.in.data) : max(x.bin.levels.in.data)]
      y.bin.levels <- levels(x = plot.data$centroid_y_disc)[min(y.bin.levels.in.data) : max(y.bin.levels.in.data)]
      
      # Expand bin coordinates for imputation (use bin levels from data to reduce computation time when uniform scale)
      centroid.data <- cbind(expand.grid(centroid_x_disc = x.bin.levels,
                                         centroid_y_disc = y.bin.levels),
                             expand.grid(centroid_x_val = seq_along(along.with = x.bin.levels),
                                         centroid_y_val = seq_along(along.with = y.bin.levels)))
      
      # Add unique ID for centroid x and y combination
      centroid.data$bin_id <- paste0("bin_", seq_len(length.out = nrow(x = centroid.data)))
      
      # Merge expanded data to plot.data
      plot.data <- merge(x = plot.data,
                         y = centroid.data,
                         by = c("centroid_x_disc", "centroid_y_disc"),
                         all = TRUE,
                         sort = FALSE)
      
      # Compute Euclidean distance between coordinates (coordinates immediately in vicinity to each other in a matrix have distances < 2)
      dist.mat <- fields::rdist(x1 = centroid.data[,c("centroid_x_val", "centroid_y_val")],
                                x2 = NULL)
      
      # Rename column and row names of distance matrix
      colnames(x = dist.mat) <- centroid.data$bin_id
      rownames(x = dist.mat) <- centroid.data$bin_id
      
      # Prepare input data for imputation of bin cluster identities
      bin.clusters <- plot.data$bin_cluster_id
      names(x = bin.clusters) <- plot.data$bin_id
      bin.clusters <- bin.clusters[colnames(x = dist.mat)]
      
      # Keep track of the iteration at which bin cluster was imputed
      imputation.iteration.l <- list()
      
      # Impute bin cluster identities
      na.bins <- sum(is.na(x = bin.clusters))
      i <- 0
      repeat {
        i <- i + 1
        imputed.clusters.l <- ImputeNearestNeighbors(dist_mat = dist.mat,
                                                     clusters = bin.clusters)
        imputed.clusters.l <- unlist(x = imputed.clusters.l,
                                     recursive = FALSE)
        imputed.clusters <- unlist(x = lapply(X = imputed.clusters.l,
                                              FUN = "[[",
                                              "output_bin_cluster"))
        is.imputed <- unlist(x = lapply(X = imputed.clusters.l,
                                        FUN = "[[",
                                        "is_imputed"))
        imputation.iteration.l[[i]] <- ifelse(test = is.imputed,
                                              yes = i,
                                              no = 0)
        if (sum(is.na(x = imputed.clusters)) < na.bins) {
          bin.clusters <- imputed.clusters
          na.bins <- sum(is.na(x = bin.clusters))
        } else {
          break
        }
      }
      
      # Simplify imputation tracker list to array
      imputation.iteration <- simplify2array(x = imputation.iteration.l)
      imputation.iteration <- rowSums(x = imputation.iteration)
      
      # Add imputed bin cluster identities to the data
      plot.data$imputed_bin_cluster_id <- imputed.clusters[plot.data$bin_id]
      
      # Add bin cluster labels and colors to the data
      plot.data$imputed_bin_cluster_label <- bin.cluster.labels[plot.data$imputed_bin_cluster_id]
      plot.data$imputed_bin_cluster_color <- bin.cluster.colors[plot.data$imputed_bin_cluster_id]
      
      # Add imputation iteration to the data
      plot.data$imputation_iteration <- imputation.iteration[plot.data$bin_id]
      
      # Add alpha values for color mixing of imputed clusters to the data
      n.alpha.levels <- length(x = alpha.levels) - 1 # first level for non imputed bins (alpha = 1)
      plot.data$imputation_iteration_cat <- ifelse(test = plot.data$imputation_iteration >= n.alpha.levels,
                                                   yes = paste0(n.alpha.levels, "+"),
                                                   no = plot.data$imputation_iteration) # categorical transformation of imputation iteration
      plot.data$imputed_bin_cluster_color_alpha <- alpha.levels[as.character(x = plot.data$imputation_iteration_cat)]
      
      # Add imputation iteration to bin cluster identity and label
      plot.data$imputed_bin_cluster_id_iter <- ifelse(test = plot.data$imputation_iteration == 0,
                                                      yes = plot.data$imputed_bin_cluster_id,
                                                      no = paste0(plot.data$imputed_bin_cluster_id,
                                                                  " iter",
                                                                  plot.data$imputation_iteration_cat))
      plot.data$imputed_bin_cluster_label_iter <- ifelse(test = plot.data$imputation_iteration == 0,
                                                         yes = plot.data$imputed_bin_cluster_label,
                                                         no = paste0(plot.data$imputed_bin_cluster_label,
                                                                     " iter",
                                                                     plot.data$imputation_iteration_cat))
      
      # Mix bin cluster colors with white if cluster is imputed (decrease alpha values with increasing imputation iteration)
      plot.data <- plot.data[,
                             "imputed_bin_cluster_color_iter" := MixBinColors(row.data = .SD),
                             by = 1:nrow(x = plot.data),
                             .SDcols = c("imputed_bin_cluster_color", "imputed_bin_cluster_color_alpha")]
      
      # Get bin attributes
      bin.attributes <- plot.data[!is.na(x = plot.data$imputed_bin_cluster_id_iter) & !duplicated(x = plot.data$imputed_bin_cluster_id_iter),
                                  c("imputed_bin_cluster_id_iter", "imputed_bin_cluster_label_iter", "imputed_bin_cluster_color_iter")]
      
      # Order bin attributes by bin_order
      # bin.attributes <- bin.attributes[stringi::stri_order(str = bin.attributes$bin_order, numeric = TRUE),]
      
      # Define bin cluster colors (used to show only main clusters in the legend)
      bin.cluster.colors <- bin.attributes$imputed_bin_cluster_color_iter
      names(x = bin.cluster.colors) <- bin.attributes$imputed_bin_cluster_id_iter
      
      # Define column used for fill aesthetics
      cluster.cols <- c(cluster.cols,
                        "post_imputation" = "imputed_bin_cluster_id",
                        "post_imputation_alpha" = "imputed_bin_cluster_id_iter")
    }
    
    # Define file name for list of plot(s)
    file.name.l <- file.path(path,
                             paste("e20230417_e20230425_CLA_map_filtered_data_",
                                   slice.level,
                                   "_",
                                   n.images,
                                   "_images_binned_spatial_maps_",
                                   names(x = cluster.cols),
                                   {
                                     if (uniform.scale) {
                                       "_uniform_scale"
                                     } else {
                                       "_nonuniform_scale"
                                     }
                                   },
                                   ".pdf",
                                   sep = ""))
    
    # Define plot title for list of plot(s)
    plot.title <- paste("claustro-insular map\n(",
                        gsub(pattern = "level",
                             replacement = "level ",
                             x = slice.level),
                        "; n = ",
                        n.images,
                        " images)",
                        sep = "")
    
    # Plot binned claustro-insular map of broad clusters
    binned.spatial.maps <- mapply(
      cluster.col = cluster.cols,
      file.name = file.name.l,
      FUN = function(cluster.col,
                     file.name,
                     plot.data,
                     plot.title,
                     axis.attributes,
                     bin.cluster.colors,
                     cluster.labels,
                     cluster.breaks,
                     clusters.of.interest,
                     save.plot) {
        # Subset bin cluster colors (in case is needed)
        bin.cluster.colors <- bin.cluster.colors[na.omit(object = unique(x = plot.data[[cluster.col]]))]
        
        # Plot claustro-insular map of broad clusters
        binned.spatial.map <- ggplot(data = plot.data,
                                     mapping = aes_string(x = "centroid_x_disc",
                                                          y = "centroid_y_disc",
                                                          fill = cluster.col)) +
          geom_tile() +
          scale_x_discrete(limits = axis.attributes$x.bin.levels) +
          scale_y_discrete(limits = rev(x = axis.attributes$y.bin.levels)) +
          scale_fill_manual(values = bin.cluster.colors, # should be named vectors for this to work
                            labels = cluster.labels[names(x = cluster.labels) %in% clusters.of.interest],
                            breaks = cluster.breaks[names(x = cluster.breaks) %in% clusters.of.interest],
                            # values = bin.attributes$bin_cluster_color,
                            # labels = bin.attributes$bin_cluster_label,
                            # breaks = bin.attributes$bin_cluster_id,
                            na.value = "white") +
          coord_cartesian(clip = "off") +
          labs(title = plot.title,
               fill = "bin clusters") +
          theme_classic() +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                legend.text = element_markdown(size = 6,
                                               family = "Arial",
                                               colour = "black",
                                               margin = margin(t = 0,
                                                               r = 0,
                                                               b = 0,
                                                               l = 0,
                                                               unit = "mm")),
                legend.title = element_text(size = 6,
                                            family = "Arial",
                                            colour = "black"),
                axis.ticks = element_blank(),
                axis.line = element_blank(),
                plot.title = element_text(size = 6,
                                          family = "Arial",
                                          colour = "black",
                                          hjust = 0.5,
                                          vjust = 0,
                                          margin = margin(t = 0,
                                                          r = 0,
                                                          b = 1,
                                                          l = 0,
                                                          unit = "mm")),
                legend.margin = margin(t = 0,
                                       r = 0,
                                       b = 0,
                                       l = 0,
                                       unit = "mm"),
                legend.key.size = unit(x = 3,
                                       units = "mm"),
                strip.text = element_blank(),
                strip.background = element_blank(),
                panel.background = element_blank(),
                plot.background = element_blank(),
                legend.background = element_blank())
        
        if (save.plot) {
          # Save plot in pdf format
          FixSizeAndSave(plot = binned.spatial.map,
                         filename = file.name,
                         is.ggassemble = FALSE,
                         # panel.width = min(7, 7 * axis.attributes$range.ratio),
                         # panel.height = min(7, 7 / axis.attributes$range.ratio),
                         panel.width = 4.2, # figure limit for width
                         panel.height = 4.2 / axis.attributes$range.ratio,
                         unit.use = "cm",
                         margin = 0,
                         use.ggsave = TRUE,
                         useDingbats = FALSE)
        }
        
        # Return plot
        return(binned.spatial.map)
      },
      MoreArgs = list(plot.data = plot.data,
                      plot.title = plot.title,
                      axis.attributes = axis.attributes,
                      bin.cluster.colors = bin.cluster.colors,
                      cluster.labels = cluster.labels,
                      cluster.breaks = cluster.breaks,
                      clusters.of.interest = clusters.of.interest,
                      save.plot = save.plot),
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE)
    
    # Return plots
    return(binned.spatial.maps)
  },
  clusters.of.interest = clusters.of.interest,
  cluster.colors = figure.attributes.l$plot$clustered_maps$cluster.colors,
  cluster.breaks = figure.attributes.l$plot$clustered_maps$cluster.breaks,
  cluster.labels = figure.attributes.l$plot$clustered_maps$cluster.labels,
  alpha.levels = figure.attributes.l$plot$clustered_maps$alpha.levels,
  axis.attributes = figure.attributes.l$axis,
  impute.clusters = TRUE,
  uniform.scale = TRUE,
  save.plot = TRUE,
  path = "./figures/CLA_map/binned_maps")
