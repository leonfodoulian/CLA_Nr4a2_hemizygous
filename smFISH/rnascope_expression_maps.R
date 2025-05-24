# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/")

# Set stringsAsFactors to FALSE
options(stringsAsFactors = FALSE)

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/rnascope_claustro_insular_maps_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/discretise.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(patchwork)
require(ggplot2)
require(extrafont)

# Create directory where to save files
dir <- "./figures/expression_maps"
dir.create(path = dir,
           recursive = TRUE)

# Get list of RNAscope files
rds.files <- grep(
  pattern = "all_cells",
  x = list.files(
    path = "./results",
    pattern = ".rds$",
    full.names = TRUE,
    recursive = TRUE
  ),
  value = TRUE
)

# Define specific parameters for each file
file.attributes.l <- lapply(
  X = rds.files,
  FUN = function(rds.file) {
    return(
      list(
        rds.file = rds.file,
        probe.set = gsub(pattern = ".*_QuPath/|/all_cells.*",
                         replacement = "",
                         x = rds.file),
        out.dir = gsub(pattern = ".*all_cells/|_log_transformed.*",
                       replacement = "",
                       x = rds.file),
        channels = strsplit(x = gsub(pattern = ".*_QuPath/|/all_cells.*",
                                     replacement = "",
                                     x = rds.file),
                            split = "(_)?C[0-9]+\\.",
                            perl = TRUE)[[1]][2:4]
      )
    )
  }
)
names(x = file.attributes.l) <- lapply(X = file.attributes.l, FUN = "[[", "probe.set")

# Load and preprocess RNAscope data
rnascope.data.l <- lapply(
  X = file.attributes.l,
  FUN = function(file.attributes) {
    # Load RNAscope data
    rnascope.data <- readRDS(file = file.attributes$rds.file)$full.data
    
    # Get image name from cell names
    rnascope.data$image <- gsub(
      pattern = "_cell.*",
      replacement = "",
      x = rnascope.data$cell
    )
    
    # Create a sham cluster column
    rnascope.data$sham_cluster <- "sham"
    
    # Return data
    return(rnascope.data)
  }
)

# Get axis attributes
axis.attributes.l <- lapply(
  X = rnascope.data.l,
  FUN = function(rnascope.data) {
    GetAxisAttributes(x = rnascope.data$centroid_x,
                      y = rnascope.data$centroid_y,
                      bin.size = 10)
  }
)

# Get maximum ranges of x and y coordinates
max.range.x <- max(
  unlist(
    x = lapply(
      X = axis.attributes.l,
      FUN = "[[",
      "range.x"
    )
  )
)
max.range.y <- max(
  unlist(
    x = lapply(
      X = axis.attributes.l,
      FUN = "[[",
      "range.y"
    )
  )
)

# Expand x and y coordinates for each data
axis.attributes.l <- lapply(
  X = axis.attributes.l,
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
  max.range.y = max.range.y
)

# Plot and save expression maps
expression.map.l <- mapply(
  rnascope.data = rnascope.data.l,
  file.attributes = file.attributes.l,
  axis.attributes = axis.attributes.l,
  FUN = function(rnascope.data,
                 file.attributes,
                 axis.attributes,
                 path) {
    # Create new directories
    dir.create(path = file.path(path, file.attributes$out.dir))
    
    # Plot expression map per image
    expression.map.l <- lapply(
      X = split(
        x = rnascope.data,
        f = rnascope.data$image
      ),
      FUN = function(section.data,
                     file.attributes,
                     axis.attributes,
                     path) {
        # Plot expression map per gene
        expression.map.l <- lapply(
          X = setNames(
            object = file.attributes$channels,
            nm = file.attributes$channels
          ),
          FUN = function(gene,
                         section.data,
                         file.attributes,
                         axis.attributes,
                         path) {
            # Keep only cells expressing gene of interest
            cells.to.keep <- section.data[[gene]] > 0
            section.data <- section.data[cells.to.keep,]
            
            # Skip if no cells are found
            if (nrow(x = section.data) == 0) {
              return(NULL)
            }
            
            # Scale data to 1
            section.data[[gene]] <- (section.data[[gene]] - min(section.data[[gene]])) / (max(section.data[[gene]]) - min(section.data[[gene]]))
            
            # Set label of x coordinate
            if (unique(x = section.data$zone) %in% c("zone1", "zone3")) {
              embedding.x.label <- "lateral"
            } else {
              embedding.x.label <- "medial"
            }
            if (file.attributes$out.dir == "e20230412_QuPath_C2.Lxn_C3.Oprk1_C4.Fosl2" && unique(x = section.data$image) == "Slide1_zone3_B6_wt.wt_level44") {
              embedding.x.label <- "medial"
            }
            
            # Plot expression map
            expression.map <- Plot2DEmbedding(
              embedding.data = section.data[,c("centroid_x", "centroid_y")],
              cluster.data = section.data[, "sham_cluster", drop = FALSE],
              alpha.data = section.data[, gene, drop = FALSE],
              cellnames.as.rownames.in.embedding.data = TRUE,
              cellnames.as.rownames.in.cluster.data = TRUE,
              cellnames.as.rownames.in.alpha.data = TRUE,
              embedding.x.name = "centroid_x",
              embedding.y.name = "centroid_y",
              embedding.x.label = embedding.x.label,
              embedding.y.label = "dorsal",
              cluster.ident.name = "sham_cluster",
              alpha.name = gene,
              cluster.colors = setNames(object = "black",
                                        nm = "sham"),
              cluster.breaks = "sham",
              cluster.labels = setNames(object = "sham",
                                        nm = "sham"),
              add.border.to.points = TRUE,
              point.size = 0.75,
              border.stroke = NA,
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
            ) +
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
              theme(aspect.ratio = 1/axis.attributes$range.ratio) + 
              AddGradientLegend(gradient.name = gene,
                                color = "black",
                                legend.title = gene,
                                legend.labels = c("min", "max"),
                                legend.title.face = "italic",
                                font.family = "Arial",
                                legend.title.size = 7,
                                legend.text.size = 7,
                                legend.bar.width = 3,
                                legend.bar.height = 0.5)
            
            # Remove clipping of data
            expression.map$coordinates$clip <- "off"
            
            # Reorder data in plot by level of expression
            expression.map$data <- expression.map$data[order(expression.map$data[[gene]], decreasing = FALSE),]
            
            # Save spatial map
            FixSizeAndSave(plot = expression.map,
                           filename = file.path(path,
                                                file.attributes$out.dir,
                                                paste0(file.attributes$out.dir,
                                                       "_",
                                                       unique(x = section.data$image),
                                                       "_",
                                                       gene,
                                                       "_expression_map",
                                                       ".pdf")),
                           is.ggassemble = FALSE,
                           panel.width = 4.2, # figure limit for width
                           panel.height = 4.2 / axis.attributes$range.ratio,
                           unit.use = "cm",
                           margin = 0,
                           use.ggsave = TRUE,
                           useDingbats = FALSE)
            
            # Return plot
            return(expression.map)
          },
          section.data = section.data,
          file.attributes = file.attributes,
          axis.attributes = axis.attributes,
          path = path
        )
        
        # Return list of plots
        return(expression.map.l)
      },
      file.attributes = file.attributes,
      axis.attributes = axis.attributes,
      path = path)
    
    # Return list of plots
    return(expression.map.l)
  },
  MoreArgs = list(path = dir),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)
