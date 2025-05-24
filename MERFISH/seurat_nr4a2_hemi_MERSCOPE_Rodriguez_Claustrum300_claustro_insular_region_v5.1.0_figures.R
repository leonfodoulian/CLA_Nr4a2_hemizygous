# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/sourcefunction.R")
source(file = "/Users/leonfodoulian/scData/discretise.R")
source(file = "/Users/leonfodoulian/scData/ecdfxy.R")
source(file = "/Users/leonfodoulian/scData/CellTypeProportions.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Source specific functions
GetAxisAttributes <- sourcefunction(
  file = "/Users/leonfodoulian/scData/rnascope_claustro_insular_maps_helper_functions.R",
  fun = "GetAxisAttributes"
)

# Required packages
require(Seurat)
require(extrafont)
require(ggplot2)
require(ggrastr)

# Create new directories
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp/spatial_maps")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp/spatial_maps/predicted_clusters")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp/spatial_maps/merscope_clusters")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp/expression_maps")

# Define rds file name prefix
file.prefix <- "MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region_seurat_v5.1.0_sctransform_azimuth_v0.5.0_neurons_filtered"

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Define CLA side to map coordinate correspondence
cla.side.coordinate <- c("CLA.1" = "medial",
                         "CLA.2" = "lateral")

# Load Seurat object
obj <- readRDS(file = file.path("results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Add abbreviated genotypes to Seurat object meta data
obj@meta.data$genotype_abbr <- ifelse(test = obj@meta.data$genotype == "Nr4a2.del.wt",
                                      yes = "del",
                                      no = "wt")

# Load pre-processed MERSCOPE data and extract cell meta data
cell.metadata <- dplyr::bind_rows(
  lapply(X = readRDS(file = "results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region_preprocessed_data.rds"),
         FUN = "[[",
         "cell.metadata"),
  .id = NULL
)

# Subset cell meta data for cells located in the claustro-insular region
cell.metadata <- cell.metadata[cell.metadata$clausto_insular_region,]

# Get list of predicted cell type labels from Seurat object
predicted.cell.type <- setNames(object = obj@meta.data$predicted.cell_type_abbr,
                                nm = rownames(x = obj@meta.data))

# Add predicted cell type labels to cell meta data
cell.metadata$predicted.cell_type_abbr <- unname(obj = predicted.cell.type[rownames(x = cell.metadata)])

# Replace NA cell type labels with 'filtered cells'
cell.metadata$predicted.cell_type_abbr[is.na(x = cell.metadata$predicted.cell_type_abbr)] <- "filtered cells"

# Get counts data
counts.data <- LayerData(object = obj,
                         layer = "counts",
                         assay = "MERSCOPE")

# Transform counts data to matrix
counts.data <- as.matrix(x = counts.data)

# Define genes of interest for expression maps
genes.of.interest <- c("Ccn2", "Col6a1", "Drd1", "Gnb4", "Nr4a2", "Rorb", "Rprm", "Sdk2", "Sst", "Syt17")
names(x = genes.of.interest) <- genes.of.interest

# Define colors for genes
genes.colors <- c("#66CDAA", "#FF4D00", "#7DAE39", "#4B0082", "#4B0082", "#E2CC3B", "#1E90FF", "#1E90FF", "#278081", "#FF4D00")
names(x = genes.colors) <- genes.of.interest

# Get meta data from Seurat object
meta.data <- obj@meta.data

# Create a sham cluster column
meta.data$sham_cluster <- "sham"

# Add expression levels of genes of interest to meta data
meta.data <- cbind(meta.data, t(x = counts.data[genes.of.interest,]))

# Prepare meta data for spatial maps
meta.data.l <- split(x = meta.data,
                     f = meta.data[,c("experiment", "mouse", "genotype", "level", "cla_side")],
                     sep = "_")
empty.data <- lapply(X = meta.data.l, FUN = nrow) == 0
meta.data.l <- meta.data.l[!empty.data]

# Define list of UMAP plots to produce
umap.params.l <- list(
  "all_cells" = NULL,
  "nr4a2_wt" = list(genotype = c("C57BL6", "Nr4a2.wt.wt"),
                    title = "*Nr4a2<sup>wt/wt</sup>*"),
  "nr4a2_hemi" = list(genotype = "Nr4a2.del.wt",
                      title = "*Nr4a2<sup>del/wt</sup>*")
)

# Define parameters for ECDF plots
prediction.params.l <- list(
  "broad_cell_type_score" = list(prediction.col = "predicted.broad_cell_type_abbr",
                                 score.col = "predicted.broad_cell_type_abbr.score",
                                 x.title = "prediction score",
                                 add.sec.axis = TRUE,
                                 cluster.breaks = cluster.attributes$broad.clusters$cluster.breaks,
                                 cluster.colors = cluster.attributes$broad.clusters$cluster.colors,
                                 file.suffix = "_ecdf_plot_broad_cell_type_prediction_scores.pdf"),
  "cell_type_score" = list(prediction.col = "predicted.cell_type_abbr",
                           score.col = "predicted.cell_type_abbr.score",
                           x.title = "prediction score",
                           add.sec.axis = TRUE,
                           cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
                           cluster.colors = cluster.attributes$subclusters$cluster.colors,
                           file.suffix = "_ecdf_plot_cell_type_prediction_scores.pdf"),
  "mapping_score" = list(prediction.col = "predicted.cell_type_abbr",
                         score.col = "mapping.score",
                         x.title = "mapping score",
                         add.sec.axis = FALSE,
                         cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
                         cluster.colors = cluster.attributes$subclusters$cluster.colors,
                         file.suffix = "_ecdf_plot_cell_type_mapping_scores.pdf")
)

# Define parameters for spatial maps
spatial.params.l <- list(
  "predicted_clusters" = list(splitting.col = "predicted.cell_type_abbr",
                              spatial.map.type = "predicted_clusters"),
  "merscope_clusters" = list(splitting.col = "SCT_snn_res.0.8",
                             spatial.map.type = "merscope_clusters")
)

# Plot clusters on UMAP
umap.plot.l <- mapply(
  umap.params = umap.params.l,
  umap.to.plot = names(x = umap.params.l),
  FUN = function(umap.params,
                 umap.to.plot,
                 obj,
                 cluster.attributes,
                 save.plot,
                 file.name) {
    # Prepare data
    if (!is.null(x = umap.params)) {
      # Subset Seurat object for genotype
      obj <- subset(x = obj,
                    subset = genotype %in% umap.params$genotype)
      
      # Rename file name to include genotype information
      file.name <- gsub(pattern = ".pdf$",
                        replacement = paste0("_",
                                             umap.to.plot,
                                             ".pdf"),
                        x = file.name)
    }
    
    # Plot clusters on UMAP
    umap.plot.l <- lapply(
      X = c(TRUE, FALSE, NA),
      FUN = function(show.prediction.score,
                     umap.params,
                     obj,
                     cluster.attributes,
                     save.plot,
                     file.name) {
        # Set alpha to 1 if predictions scores should not be displayed
        if (!isTRUE(x = show.prediction.score) && !is.na(x = show.prediction.score)) {
          obj@meta.data[["predicted.cell_type_abbr.score"]] <- 1
          file.name <- gsub(pattern = ".pdf$",
                            replacement = "_cell_types.pdf",
                            x = file.name)
        } else if (isTRUE(x = show.prediction.score)) {
          obj@meta.data[["predicted.cell_type_abbr"]] <- "sham"
          cluster.attributes$subclusters$cluster.colors <- c("sham" = "black")
          cluster.attributes$subclusters$cluster.breaks <- "sham"
          file.name <- gsub(pattern = ".pdf$",
                            replacement = "_prediction_scores.pdf",
                            x = file.name)
        }
        
        # Plot clusters on UMAP
        umap.plot <- Plot2DEmbedding(
          embedding.data = obj@reductions$umap_harmony_projected@cell.embeddings,
          cluster.data = obj@meta.data[, "predicted.cell_type_abbr", drop = FALSE],
          alpha.data = obj@meta.data[, "predicted.cell_type_abbr.score", drop = FALSE],
          cellnames.as.rownames.in.embedding.data = TRUE,
          cellnames.as.rownames.in.cluster.data = TRUE,
          cellnames.as.rownames.in.alpha.data = TRUE,
          embedding.x.name = "umapharmonyprojected_1",
          embedding.y.name = "umapharmonyprojected_2",
          embedding.x.label = "UMAP 1",
          embedding.y.label = "UMAP 2",
          cluster.ident.name = "predicted.cell_type_abbr",
          alpha.name = "predicted.cell_type_abbr.score",
          cluster.colors = cluster.attributes$subclusters$cluster.colors,
          cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
          cluster.labels = setNames(object = cluster.attributes$subclusters$cluster.breaks,
                                    nm = cluster.attributes$subclusters$cluster.breaks),
          add.border.to.points = TRUE,
          point.size = 0.5,
          border.stroke = NA,
          border.colour = "transparent",
          alpha.use = 1,
          legend.title.face = "plain",
          font.family = "Arial",
          import.font = FALSE,
          fix.aspect.ratio = FALSE,
          legend.point.size = 2,
          axis.title.size = 5,
          legend.title.size = 7,
          legend.text.size = 7,
          legend.text.space = unit(x = 3, units = "mm"),
          legend.ncol = 1,
          axis.line.size = rel(0.5),
          arrow.size = rel(2),
          arrow.type = "closed",
          axis.lineend = "butt",
          axis.gap = 0,
          range.scale = 0.08
        )
        
        # Plot prediction score legend
        if (isTRUE(x = show.prediction.score) | is.na(x = show.prediction.score)) {
          umap.plot <- umap.plot +
            AddGradientLegend(gradient.name = "predicted.cell_type_abbr.score",
                              color = "black",
                              legend.title = "prediction score",
                              legend.labels = c("min", "max"),
                              legend.title.face = "plain",
                              font.family = "Arial",
                              legend.title.size = 7,
                              legend.text.size = 7,
                              legend.bar.width = 3,
                              legend.bar.height = 0.5)
        }
        
        if (!is.null(x = umap.params)) {
          # Add title to UMAP plot
          umap.plot <- umap.plot +
            labs(title = umap.params$title) +
            theme(plot.title = ggtext::element_markdown(size = 7,
                                                        family = "Arial",
                                                        colour = "black",
                                                        hjust = 0.5,
                                                        vjust = 0,
                                                        margin = margin(t = 0,
                                                                        r = 0,
                                                                        b = 0,
                                                                        l = 0,
                                                                        unit = "mm"))) 
        }
        
        # Save plot in pdf format
        if (save.plot) {
          FixSizeAndSave(plot = umap.plot,
                         filename = file.name,
                         is.ggassemble = FALSE,
                         panel.width = 5.5,
                         panel.height = 5.5,
                         unit.use = "cm",
                         margin = 1,
                         use.ggsave = TRUE,
                         useDingbats = FALSE)
          
          # Rasterise UMAP plot
          umap.plot.raster <- ggrastr::rasterise(input = umap.plot,
                                                 dpi = 1000)
          
          # Save plot in pdf format
          FixSizeAndSave(plot = umap.plot.raster,
                         filename = gsub(pattern = ".pdf$",
                                         replacement = "_rasterised_1000dpi.pdf",
                                         x = file.name),
                         is.ggassemble = FALSE,
                         panel.width = 5.5,
                         panel.height = 5.5,
                         unit.use = "cm",
                         margin = 1,
                         use.ggsave = TRUE,
                         useDingbats = FALSE)
        }
        
        # Return plot
        return(umap.plot)
      },
      umap.params = umap.params,
      obj = obj,
      cluster.attributes,
      save.plot = save.plot,
      file.name = file.name)
    
    # Return list of plots
    return(umap.plot.l)
  },
  MoreArgs = list(obj = obj,
                  cluster.attributes = cluster.attributes,
                  save.plot = TRUE,
                  file.name = file.path("results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp",
                                        paste0(file.prefix, "_cluster_umap_plot.pdf"))),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# Plot violin plots of cell metrics (volume, # gene, # counts) per cluster
metrics.data <- t(x = cell.metadata[,c("volume", "ngene", "ncount")])
rownames(x = metrics.data) <- c("volume", "# of genes", "# of transcripts")
metrics.vln.plots <- StackViolinPlots(
  data.use = metrics.data,
  genes.use = rownames(x = metrics.data),
  cluster.data = cell.metadata[, "predicted.cell_type_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "predicted.cell_type_abbr",
  cluster.colors = c(cluster.attributes$subclusters$cluster.colors, "filtered cells" = "#D3D3D3"),
  cluster.breaks = c(cluster.attributes$subclusters$cluster.breaks, "filtered cells"),
  cluster.labels = setNames(object = c(cluster.attributes$subclusters$cluster.breaks, "filtered cells"),
                            nm = c(cluster.attributes$subclusters$cluster.breaks, "filtered cells")),
  is.log.transformed = FALSE,
  log.scale = "log",
  pseudocount.use = 1,
  y.scale.trans = "identity",
  y.min = "0",
  point.size = 2,
  alpha.use = 1,
  vln.border.colour = NA,
  vln.border.stroke = 0.25,
  plot.title.size = 7,
  plot.title.angle = 45,
  plot.title.face = "plain",
  font.family = "Arial",
  import.font = FALSE,
  hjust.use = 0,
  vjust.use = 0,
  x.axis.title = "0",
  axis.text.x.size = 6,
  axis.text.y.size = 7,
  axis.text.face = "plain",
  axis.ticks.length = 0.5,
  axis.line.size = rel(0.5),
  round.to.ceiling = TRUE,
  verbose = TRUE)

# Save plot in pdf format
FixSizeAndSave(plot = metrics.vln.plots,
               filename = file.path("results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_violin_plots_cell_metrics.pdf")),
               is.ggassemble = TRUE,
               panel.width = 1,
               panel.height = 0.325 * (length(x = cluster.attributes$subclusters$cluster.breaks) + 1),
               margin = 2.25,
               unit.use = "cm",
               use.ggsave = TRUE,
               useDingbats = FALSE)

# Plot ECDF per cluster for various prediction/mapping scores
ecdf.plots <- lapply(
  X = prediction.params.l,
  FUN = function(prediction.params,
                 meta.data,
                 save.plot,
                 path,
                 file.prefix) {
    # Get parameters
    prediction.col <- prediction.params$prediction.col
    score.col <- prediction.params$score.col
    x.title <- prediction.params$x.title
    add.sec.axis <- prediction.params$add.sec.axis
    cluster.breaks <- prediction.params$cluster.breaks
    cluster.colors <- prediction.params$cluster.colors
    file.suffix <- prediction.params$file.suffix
    
    # Split meta data by predicted labels
    meta.data.l <- split(x = meta.data,
                         f = meta.data[[prediction.col]])
    
    # Sort list by cluster breaks
    meta.data.l <- meta.data.l[cluster.breaks]
    
    # Get ECDF x and y values for plotting
    ecdf.data <- dplyr::bind_rows(
      lapply(
        X = meta.data.l,
        FUN = function(score.data,
                       score.col) {
          as.data.frame(x = ecdfxy(x = score.data[[score.col]],
                                   min.x = 0))
        },
        score.col = score.col
      ),
      .id = prediction.col)
    
    # Get max x values from ecdf data
    ecdf.data.max <- ecdf.data[ecdf.data$y == 1,]
    
    # Set predicted labels to factor  
    ecdf.data[[prediction.col]] <- factor(x = ecdf.data[[prediction.col]],
                                          levels = rev(x = cluster.breaks))
    
    # Plot ECDF
    ecdf.plot <- ggplot(data = ecdf.data,
                        mapping = aes_string(x = "x",
                                             y = "y",
                                             colour = prediction.col)) +
      geom_step(size = 0.5,
                show.legend = FALSE) +
      scale_x_continuous(limits = c(0, 1 + 10^-10),
                         breaks = seq(from = 0,
                                      to = 1,
                                      by = 0.25),
                         sec.axis = sec_axis(transform = ~ .,
                                             breaks = ecdf.data.max$x,
                                             labels = ecdf.data.max[[prediction.col]]),
                         expand = c(0,0)) +
      scale_y_continuous(limits = c(0,1),
                         breaks = seq(from = 0,
                                      to = 1,
                                      by = 0.25),
                         expand = c(0,0)) +
      scale_color_manual(values = cluster.colors,
                         breaks = cluster.breaks) +
      coord_cartesian(clip = "off") +
      labs(x = x.title,
           y = "Cumulative distribution") +
      theme_classic() +
      theme(axis.text = element_text(size = 6,
                                     family = "Arial",
                                     face = "plain",
                                     colour = "black"),
            axis.text.x.top = element_blank(),
            axis.title = element_text(size = 7,
                                      family = "Arial",
                                      face = "plain",
                                      colour = "black"),
            legend.text = element_blank(),
            legend.title = element_blank(),
            axis.line = element_line(size = rel(x = 0.5),
                                     colour = "black",
                                     lineend = "square"),
            axis.line.x.top = element_blank(),
            axis.ticks = element_line(size = rel(x = 0.5),
                                      colour = "black",
                                      lineend = "square"),
            axis.ticks.x.top = element_blank(),
            axis.ticks.length = ggplot2::unit(x = 0.5,
                                              units = "mm"),
            plot.title = element_blank(),
            legend.position = "none",
            legend.margin = element_blank(),
            legend.key.size = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank())
    
    # Modify secondary axis
    if (add.sec.axis) {
      ecdf.plot <- ecdf.plot +
        theme(axis.text.x.top = element_text(size = 6,
                                             family = "Arial",
                                             face = "plain",
                                             colour = cluster.colors),
              axis.ticks.x.top = element_line(size = 0.5,
                                              colour = cluster.colors,
                                              lineend = "butt"))
    }
    
    # Save plot in pdf format
    if (save.plot) {
      FixSizeAndSave(plot = ecdf.plot,
                     filename = file.path(path,
                                          paste0(file.prefix,
                                                 file.suffix)),
                     is.ggassemble = FALSE,
                     panel.width = 5.5,
                     panel.height = 5.5,
                     unit.use = "cm",
                     margin = 1,
                     use.ggsave = TRUE,
                     useDingbats = FALSE)
    }
    
    # Return plot
    return(ecdf.plot)
  },
  meta.data = meta.data,
  save.plot = TRUE,
  path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp",
  file.prefix = file.prefix
)

# # Plot violin plots of cluster markers
# markers.vln.plots <- StackViolinPlots(
#   data.use = t(x = meta.data[,"volume", drop = FALSE]),
#   genes.use = "volume",
#   cluster.data = meta.data[, "sham_cluster", drop = FALSE],
#   cellnames.as.rownames.in.cluster.data = TRUE,
#   cluster.ident.name = "sham_cluster",
#   # cluster.colors = cluster.attributes$subclusters$cluster.colors,
#   # cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
#   # cluster.labels = setNames(object = cluster.attributes$subclusters$cluster.breaks,
#   #                           nm = cluster.attributes$subclusters$cluster.breaks),
#   is.log.transformed = FALSE,
#   log.scale = "identity",
#   pseudocount.use = 1,
#   y.scale.trans = "log10",
#   y.min = "0",
#   point.size = 2,
#   alpha.use = 1,
#   vln.border.colour = NA,
#   vln.border.stroke = 0.1,
#   plot.title.size = 7,
#   plot.title.angle = 45,
#   plot.title.face = "italic",
#   font.family = "Arial",
#   import.font = FALSE,
#   hjust.use = 0,
#   vjust.use = 0,
#   x.axis.title = "max norm. UMI ",
#   axis.text.x.size = 6,
#   axis.text.y.size = 6,
#   axis.text.face = "plain",
#   axis.ticks.length = 0.5,
#   axis.line.size = rel(0.5),
#   round.to.ceiling = TRUE,
#   verbose = TRUE) +
#   geom_hline(yintercept = c(500, 3500))
# 
# # Save plot in pdf format
# FixSizeAndSave(plot = markers.vln.plots,
#                filename = file.path("results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/figures/temp",
#                                     paste0(file.prefix, "_violin_plots_cluster_marker_genes.pdf")),
#                is.ggassemble = TRUE,
#                panel.width = 0.5,
#                panel.height = 6.175,
#                margin = 2.25,
#                unit.use = "cm",
#                use.ggsave = TRUE,
#                useDingbats = FALSE)

# Plot and save spatial maps
spatial.map.l <- lapply(
  X = spatial.params.l,
  FUN = function(spatial.params,
                 meta.data.l,
                 cluster.attributes,
                 cla.side.coordinate,
                 path,
                 file.prefix) {
    # Plot and save spatial maps
    spatial.map.l <- mapply(
      section.data = meta.data.l,
      section.name = names(x = meta.data.l),
      FUN = function(section.data,
                     section.name,
                     spatial.params,
                     cluster.attributes,
                     cla.side.coordinate,
                     path,
                     file.prefix) {
        # Create new directories
        dir.create(path = file.path(path, spatial.params$spatial.map.type, section.name))
        
        # Get axis specifics for subsetted data
        axis.attributes <- GetAxisAttributes(x = section.data$centroid_x_rot,
                                             y = section.data$centroid_y_rot,
                                             bin.size = 20)
        
        # Plot spatial map per cluster
        spatial.map.l <- lapply(
          X = split(x = section.data,
                    f = section.data[[spatial.params$splitting.col]]),
          FUN = function(cluster.data,
                         spatial.params,
                         cluster.attributes,
                         cla.side.coordinate,
                         axis.attributes,
                         section.name,
                         path,
                         file.prefix) {
            # Skip if no cells are found
            if (nrow(x = cluster.data) == 0) {
              return(NULL)
            }
            
            # Get name of cluster displayed in plot
            cluster <- gsub(pattern = " |\\/",
                            replacement = "_",
                            x = unique(x = cluster.data[[spatial.params$splitting.col]])) 
            cluster <- if (spatial.params$spatial.map.type == "predicted_clusters") {
              paste0(cluster, "_cluster")
            } else {
              paste0("cluster_", cluster)
            }
            
            # Re-define cluster colors and labels
            cluster.breaks <- cluster.attributes$cluster.breaks
            cluster.breaks <- cluster.breaks[cluster.breaks %in% unique(x = cluster.data$predicted.cell_type_abbr)]
            cluster.colors <- cluster.attributes$cluster.colors[cluster.breaks]
            cluster.labels <- setNames(object = cluster.breaks,
                                       nm = cluster.breaks)
            
            # Define gradient color
            gradient.color <- if (spatial.params$spatial.map.type == "predicted_clusters") {
              cluster.colors
            } else {
              "black"
            }
            
            # Plot spatial map
            spatial.map <- Plot2DEmbedding(
              embedding.data = cluster.data[,c("centroid_x_rot", "centroid_y_rot")],
              cluster.data = cluster.data[, "predicted.cell_type_abbr", drop = FALSE],
              alpha.data = cluster.data[, "predicted.cell_type_abbr.score", drop = FALSE],
              cellnames.as.rownames.in.embedding.data = TRUE,
              cellnames.as.rownames.in.cluster.data = TRUE,
              cellnames.as.rownames.in.alpha.data = TRUE,
              embedding.x.name = "centroid_x_rot",
              embedding.y.name = "centroid_y_rot",
              embedding.x.label = cla.side.coordinate[unique(x = cluster.data$cla_side)],
              embedding.y.label = "dorsal",
              cluster.ident.name = "predicted.cell_type_abbr",
              alpha.name = "predicted.cell_type_abbr.score",
              cluster.colors = cluster.colors,
              cluster.breaks = cluster.breaks,
              cluster.labels = cluster.labels,
              add.border.to.points = TRUE,
              point.size = 1.5,
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
              AddGradientLegend(gradient.name = "predicted.cell_type_abbr.score",
                                color = gradient.color,
                                legend.title = "prediction score",
                                legend.labels = c("min", "max"),
                                legend.title.face = "plain",
                                font.family = "Arial",
                                legend.title.size = 7,
                                legend.text.size = 7,
                                legend.bar.width = 3,
                                legend.bar.height = 0.5)
            
            # Remove clipping of data
            spatial.map$coordinates$clip <- "off"
            
            # Reorder data in plot by prediction score
            spatial.map$data <- spatial.map$data[order(spatial.map$data$predicted.cell_type_abbr.score, decreasing = FALSE),]
            
            # Save spatial map
            FixSizeAndSave(plot = spatial.map,
                           filename = file.path(path,
                                                spatial.params$spatial.map.type,
                                                section.name,
                                                paste0(file.prefix,
                                                       "_",
                                                       section.name,
                                                       "_",
                                                       cluster,
                                                       "_spatial_map",
                                                       ".pdf")),
                           is.ggassemble = FALSE,
                           panel.width = 4.2, # figure limit for width
                           panel.height = 4.2 / axis.attributes$range.ratio,
                           unit.use = "cm",
                           margin = 0,
                           use.ggsave = TRUE,
                           useDingbats = FALSE)
            
            # Return plot
            return(spatial.map)
          },
          spatial.params = spatial.params,
          cluster.attributes = cluster.attributes,
          cla.side.coordinate = cla.side.coordinate,
          axis.attributes = axis.attributes,
          section.name = section.name,
          path = path,
          file.prefix = file.prefix
        )
        
        # Return list of plots
        return(spatial.map.l)
      },
      MoreArgs = list(spatial.params = spatial.params,
                      cluster.attributes = cluster.attributes,
                      cla.side.coordinate = cla.side.coordinate,
                      path = path,
                      file.prefix = file.prefix),
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE
    )
    
    # Return list of plots
    return(spatial.map.l)
  },
  meta.data.l = meta.data.l,
  cluster.attributes = cluster.attributes$subclusters,
  cla.side.coordinate = cla.side.coordinate,
  path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp/spatial_maps",
  file.prefix = file.prefix
)

# Plot and save expression maps
expression.map.l <- mapply(
  section.data = meta.data.l,
  section.name = names(x = meta.data.l),
  FUN = function(section.data,
                 section.name,
                 genes.of.interest,
                 genes.colors,
                 cla.side.coordinate,
                 path,
                 file.prefix) {
    # Create new directories
    dir.create(path = file.path(path, section.name))
    
    # Get axis specifics for subsetted data
    axis.attributes <- GetAxisAttributes(x = section.data$centroid_x_rot,
                                         y = section.data$centroid_y_rot,
                                         bin.size = 20)
    
    # Plot expression map
    expression.map.l <- lapply(
      X = genes.of.interest,
      FUN = function(gene,
                     section.data,
                     genes.colors,
                     cla.side.coordinate,
                     axis.attributes,
                     section.name,
                     path,
                     file.prefix) {
        # Keep only cells expressing gene of interest at 5 counts
        cells.to.keep <- section.data[[gene]] >= 5
        section.data <- section.data[cells.to.keep,]
        
        # Skip if no cells are found
        if (nrow(x = section.data) == 0) {
          return(NULL)
        }
        
        # Normalize data by volume and scale data to 1
        section.data[[gene]] <- log2(x = section.data[[gene]] / section.data$volume)
        section.data[[gene]] <- (section.data[[gene]] - min(section.data[[gene]])) / (max(section.data[[gene]]) - min(section.data[[gene]]))
        
        # Get color for gene of interest
        gene.color <- genes.colors[gene]
        names(x = gene.color) <- "sham"
        
        # Plot expression map
        expression.map <- Plot2DEmbedding(
          embedding.data = section.data[,c("centroid_x_rot", "centroid_y_rot")],
          cluster.data = section.data[, "sham_cluster", drop = FALSE],
          alpha.data = section.data[, gene, drop = FALSE],
          cellnames.as.rownames.in.embedding.data = TRUE,
          cellnames.as.rownames.in.cluster.data = TRUE,
          cellnames.as.rownames.in.alpha.data = TRUE,
          embedding.x.name = "centroid_x_rot",
          embedding.y.name = "centroid_y_rot",
          embedding.x.label = cla.side.coordinate[unique(x = section.data$cla_side)],
          embedding.y.label = "dorsal",
          cluster.ident.name = "sham_cluster",
          alpha.name = gene,
          cluster.colors = gene.color,
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
                            color = gene.color,
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
                                            section.name,
                                            paste0(file.prefix,
                                                   "_",
                                                   section.name,
                                                   "_",
                                                   gene,
                                                   "_expression_map",
                                                   ".pdf")),
                       is.ggassemble = FALSE,
                       panel.width = 1.5, # figure limit for width
                       panel.height = 1.5 / axis.attributes$range.ratio,
                       unit.use = "cm",
                       margin = 0,
                       use.ggsave = TRUE,
                       useDingbats = FALSE)
        
        # Return plot
        return(expression.map)
      },
      section.data = section.data,
      genes.colors = genes.colors,
      cla.side.coordinate = cla.side.coordinate,
      axis.attributes = axis.attributes,
      section.name = section.name,
      path = path,
      file.prefix = file.prefix
    )
    
    # Return list of plots
    return(expression.map.l)
  },
  MoreArgs = list(genes.of.interest = genes.of.interest,
                  genes.colors = genes.colors,
                  cla.side.coordinate = cla.side.coordinate,
                  path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp/expression_maps",
                  file.prefix = file.prefix),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# Calculate and plot cell type proportions
cell.type.prop.l <- CellTypeProportions(
  x = NULL,
  clusters = factor(x = meta.data$predicted.broad_cell_type_abbr,
                    levels = cluster.attributes$broad.clusters$cluster.breaks),
  sample = factor(x = paste(meta.data$mouse,
                            meta.data$cla_side,
                            sep = "_")),
  group = factor(x = meta.data$genotype_abbr,
                 levels = c("del", "wt")),
  trend = FALSE,
  robust = TRUE,
  transform = "logit",
  return.plot = TRUE,
  markdown = FALSE,
  reference.group = "wt",
  pseudocount = 0.5,
  cluster.colors = cluster.attributes$broad.clusters$cluster.colors,
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
FixSizeAndSave(plot = cell.type.prop.l$prop.plot,
               filename = file.path("results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_cell_type_proportion_ratios.pdf")),
               is.ggassemble = FALSE,
               panel.width = 0.45 * length(x = cluster.attributes$subclusters$cluster.breaks),
               panel.height = 4.5,
               unit.use = "cm",
               margin = 0,
               use.ggsave = FALSE,
               useDingbats = FALSE)
