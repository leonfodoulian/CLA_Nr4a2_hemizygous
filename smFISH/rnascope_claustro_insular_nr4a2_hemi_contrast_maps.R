# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/")

# Set stringsAsFactors to FALSE
options(stringsAsFactors = FALSE)

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/rcpp_ImputeNearestNeighbors.R")
source(file = "/Users/leonfodoulian/scData/spatial_registration_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/rnascope_claustro_insular_maps_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/discretise.R")
source(file = "/Users/leonfodoulian/scData/CellTypeProportions.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(patchwork)
require(ggplot2)
require(extrafont)
require(ggtext)
require(ggrastr)
require(uniqtag)

# Create directory where to save files
dirs <- c("./figures/CLA_Nr4a2_hemi_map",
          "./figures/CLA_Nr4a2_hemi_map/violin_plots",
          "./figures/CLA_Nr4a2_hemi_map/raw_maps",
          "./figures/CLA_Nr4a2_hemi_map/rotated_maps",
          "./figures/CLA_Nr4a2_hemi_map/clustered_maps",
          "./figures/CLA_Nr4a2_hemi_map/contrast_maps",
          "./figures/CLA_Nr4a2_hemi_map/proportion_plots")
lapply(X = dirs,
       FUN = dir.create,
       recursive = TRUE)

# List RDS files
rds.files <- c("results/e20230605_QuPath/C2.Nr4a2_C3.Ntm_C4.Syt17/all_cells/e20230605_QuPath_C2.Nr4a2_C3.Ntm_C4.Syt17_log_transformed_hclust_wardD2_clustered_data.rds",
               "results/e20230605_QuPath/C2.Nr4a2_C3.Ryr2_C4.Syt17/all_cells/e20230605_QuPath_C2.Nr4a2_C3.Ryr2_C4.Syt17_log_transformed_hclust_wardD2_clustered_data.rds")
names(x = rds.files) <- gsub(pattern = ".*QuPath_|_log_transformed.*",
                             replacement = "",
                             x = rds.files)

# Define reference points file
# rp.file <- "reference_points/e20230605_reference_points_20240812.xlsx"
rp.file <- "reference_points/e20230605_reference_points_20240813.xlsx"

# Define cluster names
cluster.names.l <- list("C2.Nr4a2_C3.Ntm_C4.Syt17" = c("1" = "Shell",
                                                       "2" = "noise",
                                                       "3" = "Ntm+",
                                                       "4" = "CLA",
                                                       "5" = "Syt17+",
                                                       "not clustered" = "not clustered"),
                        "C2.Nr4a2_C3.Ryr2_C4.Syt17" = c("1" = "Ryr2+",
                                                        "2" = "noise",
                                                        "3" = "Syt17+",
                                                        "4" = "CLA",
                                                        "5" = "Shell",
                                                        "not clustered" = "not clustered"))

# Define reference images
# reference.images <- c("Slide3_zone2_m88672_Nr4a2.wt.wt_level44_C2.Nr4a2_C3.Ryr2_C4.Syt17.tif",
#                       "Slide5_zone1_m88671_Nr4a2.wt.wt_level44_C2.Nr4a2_C3.Ntm_C4.Syt17.tif") # per probe set
reference.images <- "Slide3_zone2_m88672_Nr4a2.wt.wt_level44_C2.Nr4a2_C3.Ryr2_C4.Syt17.tif" # all together 

# Images to exclude
images.to.exclude <- NA

# Define clusters of interest
clusters.of.interest <- c("CLA", "Shell") # can only be of length 2

# Load RNAscope data
rnascope.data.l <- lapply(X = rds.files,
                          readRDS)

# Preprocess RNAscope data
rnascope.data.l <- mapply(
  rnascope.data = rnascope.data.l,
  cluster.names = cluster.names.l,
  cluster.col = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = "k5_clusters",
                     "C2.Nr4a2_C3.Ryr2_C4.Syt17" = "k5_clusters"),
  cell.suffix = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = "_1",
                     "C2.Nr4a2_C3.Ryr2_C4.Syt17" = "_2"),
  cluster.to.remove = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = NULL,
                           "C2.Nr4a2_C3.Ryr2_C4.Syt17" = NULL),
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
    rnascope.data$centroid_x <- ifelse(test = rnascope.data$zone %in% c("zone1", "zone3"),
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
reference.points$centroid_x_dorsal <- reference.points$centroid_x_EPd
reference.points$centroid_y_dorsal <- reference.points$centroid_y_EPd
reference.points$centroid_x_ventral <- reference.points$centroid_x_EPv
reference.points$centroid_y_ventral <- reference.points$centroid_y_EPv

# Set reference points x coordinates to the right side
reference.points$centroid_x_center <- ifelse(test = reference.points$zone %in% c("zone1", "zone3"),
                                             yes = reference.points$centroid_x_center,
                                             no = -reference.points$centroid_x_center)
reference.points$centroid_x_dorsal <- ifelse(test = reference.points$zone %in% c("zone1", "zone3"),
                                             yes = reference.points$centroid_x_dorsal,
                                             no = -reference.points$centroid_x_dorsal)
reference.points$centroid_x_ventral <- ifelse(test = reference.points$zone %in% c("zone1", "zone3"),
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
                       f = reference.points$probe_set),
  # reference.coords = split(x = reference.points[reference.points$is_reference_image,],
  #                          f = reference.points[reference.points$is_reference_image,]$probe_set), # per probe set
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
  MoreArgs = list(reference.coords = reference.points[reference.points$is_reference_image,]), # all together
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

# Create a sham clustering column for background cell plotting
rnascope.data$sham_clustering_results <- "background"

# Add abbreviated genotypes to data
rnascope.data$genotype_abbr <- ifelse(test = rnascope.data$genotype == "Nr4a2.del.wt",
                                      yes = "del",
                                      no = "wt")

# Define specific parameters for each plot type
figure.attributes.l <- list("plot" = list("raw_maps" = list(map.type = "raw",
                                                            centroid.x.column = "centroid_x",
                                                            centroid.y.column = "centroid_y",
                                                            clustering.column = "image_type",
                                                            cluster.colors = c("reference" = "#6C7B8B"),
                                                            cluster.breaks = c("reference" = "reference"),
                                                            genotype.labels = c("Nr4a2.wt.wt" = "*Nr4a2<sup>wt/wt</sup>*",
                                                                                "Nr4a2.del.wt" = "*Nr4a2<sup>del/wt</sup>*"),
                                                            legend.title = "image type (# cells)",
                                                            path = "./figures/CLA_Nr4a2_hemi_map/raw_maps"),
                                          "rotated_maps" = list(map.type = "rotated",
                                                                centroid.x.column = "centroid_x_rot",
                                                                centroid.y.column = "centroid_y_rot",
                                                                clustering.column = "image_type",
                                                                cluster.colors = c("reference" = "#6C7B8B"),
                                                                cluster.breaks = c("reference" = "reference"),
                                                                genotype.labels = c("Nr4a2.wt.wt" = "*Nr4a2<sup>wt/wt</sup>*",
                                                                                    "Nr4a2.del.wt" = "*Nr4a2<sup>del/wt</sup>*"),
                                                                legend.title = "image type (# cells)",
                                                                path = "./figures/CLA_Nr4a2_hemi_map/rotated_maps"),
                                          "clustered_maps" = list(map.type = "clustered",
                                                                  centroid.x.column = "centroid_x_rot",
                                                                  centroid.y.column = "centroid_y_rot",
                                                                  clustering.column = "clustering_results",
                                                                  cluster.colors = c("CLA" = "#4B0082",
                                                                                     "Shell" = "#FF4D00",
                                                                                     "Syt17+" = "#929497",
                                                                                     "Ntm+" = "#BBBDBF",
                                                                                     "Ryr2+" = "#BBBDBF",
                                                                                     "noise" = "#E7DECC",
                                                                                     "not clustered" = "#F0E5D3"),
                                                                  cluster.breaks = c("CLA" = "CLA",
                                                                                     "Shell" = "Shell",
                                                                                     "Syt17+" = "Syt17+",
                                                                                     "Ntm+" = "Ntm+",
                                                                                     "Ryr2+" = "Ryr2+",
                                                                                     "noise" = "noise",
                                                                                     "not clustered" = "not clustered"),
                                                                  cluster.labels = c("CLA" = "CLA",
                                                                                     "Shell" = "Shell",
                                                                                     "Syt17+" = "*Syt17*<sup>+</sup>",
                                                                                     "Ntm+" = "*Ntm*<sup>+</sup>",
                                                                                     "Ryr2+" = "*Ryr2*<sup>+</sup>",
                                                                                     "noise" = "noise",
                                                                                     "not clustered" = "not clustered"),
                                                                  genotype.labels = c("Nr4a2.wt.wt" = "*Nr4a2<sup>wt/wt</sup>*",
                                                                                      "Nr4a2.del.wt" = "*Nr4a2<sup>del/wt</sup>*"),
                                                                  legend.title = "neuronal clusters (# cells)",
                                                                  path = "./figures/CLA_Nr4a2_hemi_map/clustered_maps")),
                            "axis" = axis.attributes)

####################################################
# Figure: Plot violin plots with final cluster labels
# Plot and save violin plots
vln.plots.l <- mapply(
  expression.data = all.data,
  probe.set = names(x = all.data),
  channels = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = c("Nr4a2", "Syt17", "Ntm"),
                  "C2.Nr4a2_C3.Ryr2_C4.Syt17" = c("Nr4a2", "Syt17", "Ryr2")),
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
                       filename = paste("figures/CLA_Nr4a2_hemi_map/violin_plots/e20230605_CLA_Nr4a2_hemi_map",
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
            # f = rnascope.data[,c("probe_set", "genotype")], # per probe set and genotype
            # sep = "_"),
            f = rnascope.data$genotype), # per genotype
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
        # Define probe set
        # probe.set <- unique(x = plot.data$probe_set) # per probe set and genotype
        
        # Define genotype
        genotype <- unique(x = plot.data$genotype)
        
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
          if ("reference" %in% plot.data$image_type) {
            n.query <- length(x = image.types) - 1
          } else {
            n.query <- length(x = image.types)
            plot.attributes$cluster.colors <- NULL
            plot.attributes$cluster.breaks <- NULL
          }
          
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
        } else if (plot.attributes$map.type %in% "clustered") {
          # Re-define cluster colors, breaks and labels (in case is needed)
          plot.attributes$cluster.breaks <- plot.attributes$cluster.breaks[plot.attributes$cluster.breaks %in% plot.data$clustering_results]
          plot.attributes$cluster.colors <- plot.attributes$cluster.colors[plot.attributes$cluster.breaks]
          plot.attributes$cluster.labels <- plot.attributes$cluster.labels[plot.attributes$cluster.breaks]
        }
        
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
          labs(title = paste(plot.attributes$genotype.labels[genotype],
                             " ",
                             plot.attributes$map.type,
                             " map<br>(",
                             gsub(pattern = "level",
                                  replacement = "level ",
                                  x = slice.level),
                             "; n = ",
                             n.images,
                             " images)",
                             sep = ""),
               colour = plot.attributes$legend.title) +
          theme(plot.title = element_markdown(size = 6,
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
          data.filtering.type <- c("unfiltered", "filtered", "filtered")
          
          # Reorder data in plot
          spatial.map.plot$data[[plot.attributes$clustering.column]] <- factor(x = spatial.map.plot$data[[plot.attributes$clustering.column]],
                                                                               levels = plot.attributes$cluster.breaks)
          spatial.map.plot$data <- spatial.map.plot$data[order(x = spatial.map.plot$data[[plot.attributes$clustering.column]], decreasing = TRUE),]
          
          # Prepare plot for saving
          spatial.map.plot.l <- list("unfiltered" = spatial.map.plot,
                                     "filtered" = {
                                       spatial.map.plot$data <- spatial.map.plot$data[!spatial.map.plot$data[[plot.attributes$clustering.column]] %in% c("noise", "not clustered"),]
                                       spatial.map.plot
                                     },
                                     "filtered_rasterised" = {
                                       spatial.map.plot$data <- spatial.map.plot$data[!spatial.map.plot$data[[plot.attributes$clustering.column]] %in% c("noise", "not clustered"),]
                                       ggrastr::rasterise(input = spatial.map.plot,
                                                          dpi = 1000)
                                     })
        }
        
        # Define file name for list of plot(s)
        file.name.l <- file.path(plot.attributes$path,
                                 paste("e20230605",
                                       # "_",
                                       # probe.set, # per probe set and genotype
                                       "_CLA_Nr4a2_hemi_map_",
                                       genotype,
                                       "_",
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
                                       {
                                         ifelse(test = grepl(pattern = "rasterised",
                                                             x = names(x = spatial.map.plot.l)),
                                                yes = "_rasterised_1000dpi.pdf",
                                                no = ".pdf")
                                       },
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
# Figure: Plot CLA/Shell contrast maps
# Plot CLA/Shell contrast maps
contrast.maps <- lapply(
  X = split(x = rnascope.data,
            # f = rnascope.data[,c("probe_set", "genotype")], # per probe set and genotype
            # sep = "_"),
            f = rnascope.data$genotype), # per genotype
  FUN = function(plot.data,
                 clusters.of.interest,
                 cluster.colors,
                 genotype.labels,
                 axis.attributes,
                 uniform.scale,
                 save.plot,
                 path) {
    # Define probe set
    # probe.set <- unique(x = plot.data$probe_set) # per probe set and genotype
    
    # Define genotype
    genotype <- unique(x = plot.data$genotype)
    
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
    plot.data <- plot.data[, lapply(X = props.to.compute,
                                    FUN = function(prop.to.compute) {
                                      eval(expr = parse(text = prop.to.compute))
                                    }),
                           by = c("centroid_x_disc", "centroid_y_disc", "genotype", "level")]
    
    # Compute ClA/Shell contrast per bin
    plot.data$contrast <- plot.data[[prop.cols[1]]] - plot.data[[prop.cols[2]]]
    
    # Define file name
    file.name <- file.path(path,
                           paste("e20230605",
                                 # "_",
                                 # probe.set, # per probe set and genotype
                                 "_CLA_Nr4a2_hemi_map_",
                                 genotype,
                                 "_filtered_data_",
                                 slice.level,
                                 "_",
                                 n.images,
                                 "_images_CLA_shell_contrast_maps",
                                 {
                                   if (uniform.scale) {
                                     "_uniform_scale"
                                   } else {
                                     "_nonuniform_scale"
                                   }
                                 },
                                 ".pdf",
                                 sep = ""))
    
    # Define plot title
    plot.title <- paste(genotype.labels[genotype],
                        " ",
                        "claustro-insular map<br>(",
                        gsub(pattern = "level",
                             replacement = "level ",
                             x = slice.level),
                        "; n = ",
                        n.images,
                        " images)",
                        sep = "")
    
    # Plot CLA/Shell contrast map
    contrast.map <- ggplot(data = plot.data,
                           mapping = aes_string(x = "centroid_x_disc",
                                                y = "centroid_y_disc",
                                                fill = "contrast")) +
      geom_tile() +
      scale_x_discrete(limits = axis.attributes$x.bin.levels) +
      scale_y_discrete(limits = rev(x = axis.attributes$y.bin.levels)) +
      scale_fill_gradient2(low = cluster.colors[clusters.of.interest[2]],
                           mid = "white",
                           high = cluster.colors[clusters.of.interest[1]],
                           midpoint = 0,
                           na.value = "white",
                           limits = c(-1,1),
                           breaks = c(-1,1),
                           labels = c("Shell",
                                      "CLA"),
                           guide = guide_colourbar(title = "CLA/Shell contrast",
                                                   title.position = "top",
                                                   title.hjust = 0.5,
                                                   title.vjust = 0.5,
                                                   title.theme = element_text(angle = 0,
                                                                              margin = ggplot2::margin(t = 0,
                                                                                                       r = 0,
                                                                                                       b = 1,
                                                                                                       l = 0,
                                                                                                       unit = "mm"),
                                                                              colour = "black",
                                                                              size = 6,
                                                                              family = "Arial"),
                                                   label.hjust = 0.5,
                                                   label.vjust = 0.5,
                                                   label.position = "bottom",
                                                   label.theme = element_text(angle = 0, 
                                                                              margin = ggplot2::margin(t = 0,
                                                                                                       r = 0,
                                                                                                       b = 0,
                                                                                                       l = 0,
                                                                                                       unit = "cm"),
                                                                              colour = "black",
                                                                              size = 6,
                                                                              family = "Arial"),
                                                   barwidth = 3,
                                                   barheight = 0.5,
                                                   direction = "horizontal",
                                                   position = "bottom",
                                                   ticks.linewidth = NA)) +
      coord_cartesian(clip = "off") +
      labs(title = plot.title,
           fill = "CLA/Shell contrast") +
      theme_classic() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            plot.title = element_markdown(size = 6,
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
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank())
    
    if (save.plot) {
      # Save plot in pdf format
      FixSizeAndSave(plot = contrast.map,
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
    return(contrast.map)
  },
  clusters.of.interest = clusters.of.interest,
  cluster.colors = figure.attributes.l$plot$clustered_maps$cluster.colors,
  genotype.labels = figure.attributes.l$plot$clustered_maps$genotype.labels,
  axis.attributes = figure.attributes.l$axis,
  uniform.scale = TRUE,
  save.plot = TRUE,
  path = "./figures/CLA_Nr4a2_hemi_map/contrast_maps")

####################################################
# Figure: # Calculate and plot cell type proportions
# Plot cell type proportions
cell.type.prop.l <- CellTypeProportions(
  x = NULL,
  clusters = factor(x = rnascope.data[rnascope.data$clustering_results %in% c("CLA", "Shell", "Syt17+"),]$clustering_results,
                    levels = c("CLA", "Shell", "Syt17+")),
  sample = factor(x = rnascope.data[rnascope.data$clustering_results %in% c("CLA", "Shell", "Syt17+"),]$image),
  group = factor(x = rnascope.data[rnascope.data$clustering_results %in% c("CLA", "Shell", "Syt17+"),]$genotype_abbr,
                 levels = c("del", "wt")),
  trend = FALSE,
  robust = TRUE,
  transform = "logit",
  return.plot = TRUE,
  markdown = TRUE,
  reference.group = "wt",
  pseudocount = 0,
  cluster.colors = figure.attributes.l$plot$clustered_maps$cluster.colors[1:3],
  bar.border.stroke = 0.5,
  label.size = 15,
  label.position = "top",
  point.size = 1.25,
  sd.line.size = 0.5,
  y.axis.title = "log<sub>2</sub> ratio<br>*Nr4a2<sup>del/wt</sup>* &frasl; *Nr4a2<sup>wt/wt</sup>*",
  font.family = "Arial",
  font.face = "plain",
  axis.text.size = 6,
  axis.title.size = 7,
  axis.ticks.length = 0.5,
  axis.line.size = rel(x = 0.5)
)

# Save plot in pdf format
FixSizeAndSave(plot = cell.type.prop.l$prop.plot +
                 scale_x_discrete(breaks = figure.attributes.l$plot$clustered_maps$cluster.breaks[1:3],
                                  label = figure.attributes.l$plot$clustered_maps$cluster.labels[1:3]),
               filename = "./figures/CLA_Nr4a2_hemi_map/proportion_plots/e20230605_CLA_Nr4a2_hemi_map_cell_type_proportion_ratios.pdf",
               is.ggassemble = FALSE,
               panel.width = 0.45 * 3,
               panel.height = 4.5,
               unit.use = "cm",
               margin = 0,
               use.ggsave = FALSE,
               useDingbats = FALSE)

# Plot cell type proportions for all genotypes
cell.type.prop.l$sample.prop$group <- factor(x = cell.type.prop.l$sample.prop$group,
                                             levels = c("wt", "del"))
prop.plot <- ggplot(data = cell.type.prop.l$sample.prop,
                    mapping = aes(x = clusters,
                                  y = proportion,
                                  color = group)) +
  stat_summary(geom = "col",
               size = 0.5,
               fill = NA,
               width = 0.75,
               position = position_dodge(width = 0.9,
                                         preserve = "total"),
               fun = "mean",
               show.legend = FALSE) +
  geom_point(size = 1.25,
             shape = 19,
             stroke = 0,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0,
                                             dodge.width = 0.9,
                                             seed = 47)) +
  scale_x_discrete(limits = figure.attributes.l$plot$clustered_maps$cluster.breaks[1:3],
                   labels = figure.attributes.l$plot$clustered_maps$cluster.labels[1:3]) +
  scale_y_continuous(limits = c(0, 0.1 * ceiling(x = max(cell.type.prop.l$sample.prop$proportion) / 0.1)),
                     breaks = seq(from = 0,
                                  to = 0.1 * ceiling(x = max(cell.type.prop.l$sample.prop$proportion) / 0.1),
                                  by = 0.1),
                     expand = c(0,0)) +
  scale_color_manual(values = c("wt" = "black",
                                "del" = "grey"),
                     labels = c("wt" = "*Nr4a2<sup>wt/wt</sup>*",
                                "del" = "*Nr4a2<sup>del/wt</sup>*")) +
  guides(colour = ggplot2::guide_legend(override.aes = list(size = 2,
                                                            alpha = 1,
                                                            stroke = 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "cell types",
       y = "proportion of cells",
       colour = "genotype") +
  theme_classic() +
  theme(axis.text.x = ggtext::element_markdown(size = 6,
                                               angle = 45,
                                               hjust = 1,
                                               vjust = 1,
                                               family = "Arial",
                                               face = "plain",
                                               colour = figure.attributes.l$plot$clustered_maps$cluster.colors[1:3]),
        axis.text.y = element_text(size = 6,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black"),
        axis.title = element_text(size = 7,
                                  family = "Arial",
                                  face = "plain",
                                  colour = "black"),
        legend.text = ggtext::element_markdown(size = 7,
                                               family = "Arial",
                                               face = "plain",
                                               colour = "black",
                                               margin = margin(t = 0,
                                                               r = 0,
                                                               b = 0,
                                                               l = 0,
                                                               unit = "mm")),
        legend.title = element_text(size = 7,
                                    family = "Arial",
                                    face = "plain",
                                    colour = "black"),
        axis.line = element_line(size = rel(x = 0.5),
                                 colour = "black",
                                 lineend = "square"),
        axis.ticks = element_line(size = rel(x = 0.5),
                                  colour = "black",
                                  lineend = "square"),
        axis.ticks.length = ggplot2::unit(x = 0.5,
                                          units = "mm"),
        plot.title = element_blank(),
        legend.position = "right",
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

# Save plot in pdf format
FixSizeAndSave(plot = prop.plot,
               filename = "./figures/CLA_Nr4a2_hemi_map/proportion_plots/e20230605_CLA_Nr4a2_hemi_map_cell_type_proportions.pdf",
               is.ggassemble = FALSE,
               panel.width = 0.45 * 3 * 2,
               panel.height = 4.5,
               unit.use = "cm",
               margin = 0,
               use.ggsave = FALSE,
               useDingbats = FALSE)
