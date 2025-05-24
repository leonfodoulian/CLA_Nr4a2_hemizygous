# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(extrafont)
require(ggplot2)
require(ggrastr)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_unfiltered/figures/temp")

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_unfiltered

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-20_neurons_unfiltered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_unfiltered",
                                paste0(file.prefix, ".rds")))

# Join layers of the RNA assay
if (!"RNA_full" %in% names(x = obj@assays)) {
  obj[["RNA_full"]] <- JoinLayers(object = obj[["RNA"]])
}

# Prepare data for violin plots
obj.wt <- subset(x = obj,
                 subset = genotype == "Nr4a2(WT/WT)")

# Define list of UMAP plots to produce
umap.params.l <- list(
  "all_cells" = NULL,
  "nr4a2_wt" = list(genotype = "Nr4a2(WT/WT)",
                    title = "*Nr4a2<sup>wt/wt</sup>*"),
  "nr4a2_hemi" = list(genotype = "Nr4a2(SA-IRES-Dre/WT)",
                      title = "*Nr4a2<sup>del/wt</sup>*")
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
                    subset = genotype == umap.params$genotype)
      
      # Rename file name to include genotype information
      file.name <- gsub(pattern = ".pdf$",
                        replacement = paste0("_",
                                             umap.to.plot,
                                             ".pdf"),
                        x = file.name)
    }
    
    # Plot clusters on UMAP
    umap.plot <- Plot2DEmbedding(
      embedding.data = obj@reductions$umap_harmony@cell.embeddings,
      cluster.data = obj@meta.data[, "broad_cell_type_abbr", drop = FALSE],
      cellnames.as.rownames.in.embedding.data = TRUE,
      cellnames.as.rownames.in.cluster.data = TRUE,
      embedding.x.name = "UMAP_1",
      embedding.y.name = "UMAP_2",
      embedding.x.label = "UMAP 1",
      embedding.y.label = "UMAP 2",
      cluster.ident.name = "broad_cell_type_abbr",
      cluster.colors = cluster.attributes$broad.clusters$cluster.colors,
      cluster.breaks = cluster.attributes$broad.clusters$cluster.breaks,
      cluster.labels = setNames(object = cluster.attributes$broad.clusters$cluster.breaks,
                                nm = cluster.attributes$broad.clusters$cluster.breaks),
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
    
    # Add labels or title to UMAP plot depending if data was subsetted
    if (is.null(x = umap.params)) {
      # Add labels to UMAP plot
      umap.plot <- umap.plot +
        # ggtext::geom_richtext(
        geom_text(
          data = GetLabelCentroid(label.data = umap.plot$data,
                                  label.name = "broad_cell_type_abbr",
                                  embedding.x.name = "UMAP_1",
                                  embedding.y.name =  "UMAP_2"),
          mapping = aes(x = UMAP_1,
                        y = UMAP_2,
                        label = broad_cell_type_abbr),
          colour = "black",
          fill = NA,
          alpha = 1,
          size = 7 * (1/72 * 25.4),
          family = "Arial",
          show.legend = FALSE,
          inherit.aes = FALSE,
          label.color = NA,
          label.padding = ggplot2::unit(x = c(0.5, 0.5, -0.5, 0.5),
                                        units = "pt")
        )
    } else {
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
  MoreArgs = list(obj = obj,
                  cluster.attributes = cluster.attributes,
                  save.plot = TRUE,
                  file.name = file.path("results/smim32_GTFv102/neurons_unfiltered/figures/temp",
                                        paste0(file.prefix, "_cluster_umap_plot.pdf"))),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# Plot violin plots of cell metrics (# gene, # UMI, proportion of mitochondrial gene counts) per cluster
metrics.data <- t(x = obj@meta.data[,c("nFeature_RNA", "nCount_RNA", "percent_mt")])
rownames(x = metrics.data) <- c("# of genes", "# of UMI", "% of mt- gene counts")
metrics.vln.plots <- StackViolinPlots(
  data.use = metrics.data,
  genes.use = rownames(x = metrics.data),
  cluster.data = obj@meta.data[, "broad_cell_type_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "broad_cell_type_abbr",
  cluster.colors = cluster.attributes$broad.clusters$cluster.colors,
  cluster.breaks = cluster.attributes$broad.clusters$cluster.breaks,
  cluster.labels = setNames(object = cluster.attributes$broad.clusters$cluster.breaks,
                            nm = cluster.attributes$broad.clusters$cluster.breaks),
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
               filename = file.path("results/smim32_GTFv102/neurons_unfiltered/figures/temp",
                                    paste0(file.prefix, "_violin_plots_cell_metrics.pdf")),
               is.ggassemble = TRUE,
               panel.width = 1,
               panel.height = 0.325 * length(x = cluster.attributes$broad.clusters$cluster.breaks),
               margin = 2.25,
               unit.use = "cm",
               use.ggsave = TRUE,
               useDingbats = FALSE)

# Plot violin plots of cluster markers
markers.vln.plots <- StackViolinPlots(
  data.use = LayerData(object = obj.wt,
                       layer = "data",
                       assay = "RNA_full"),
  genes.use = c("Slc17a7", "Slc32a1", # neuronal
                "Nr4a2", "Gnb4", # CLA
                "Sdk2", # shell and L6
                "Syt17", "Col6a1", # shell
                "Rprm", # L6a
                "Ccn2", # L6b
                "Tafa1", # L2/3-L5b-pir
                "Drd1", # dMSN
                "Drd2", # iMSN
                "Tshz1", # sMSN
                "C1qb" # micro
  ),
  cluster.data = obj.wt@meta.data[, "broad_cell_type_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "broad_cell_type_abbr",
  cluster.colors = cluster.attributes$broad.clusters$cluster.colors,
  cluster.breaks = cluster.attributes$broad.clusters$cluster.breaks,
  cluster.labels = setNames(object = cluster.attributes$broad.clusters$cluster.breaks,
                            nm = cluster.attributes$broad.clusters$cluster.breaks),
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  y.scale.trans = "log10",
  y.min = "0",
  point.size = 2,
  alpha.use = 1,
  vln.border.colour = NA,
  vln.border.stroke = 0.25,
  plot.title.size = 7,
  plot.title.angle = 45,
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  hjust.use = 0,
  vjust.use = 0,
  x.axis.title = "max norm. UMI ",
  axis.text.x.size = 6,
  axis.text.y.size = 7,
  axis.text.face = "plain",
  axis.ticks.length = 0.5,
  axis.line.size = rel(0.5),
  round.to.ceiling = TRUE,
  verbose = TRUE)

# Save plot in pdf format
FixSizeAndSave(plot = markers.vln.plots,
               filename = file.path("results/smim32_GTFv102/neurons_unfiltered/figures/temp",
                                    paste0(file.prefix, "_violin_plots_cluster_marker_genes.pdf")),
               is.ggassemble = TRUE,
               panel.width = 0.5,
               panel.height = 0.325 * length(x = cluster.attributes$broad.clusters$cluster.breaks),
               margin = 2.25,
               unit.use = "cm",
               use.ggsave = TRUE,
               useDingbats = FALSE)
