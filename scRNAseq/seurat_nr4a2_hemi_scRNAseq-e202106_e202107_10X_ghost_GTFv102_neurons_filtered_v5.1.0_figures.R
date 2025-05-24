# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/ExpandClusterAttributes.R")
source(file = "/Users/leonfodoulian/scData/BlendColorAlpha.R")
source(file = "/Users/leonfodoulian/scData/CountDEG.R")
source(file = "/Users/leonfodoulian/scData/CellTypeProportions.R")
source(file = "/Users/leonfodoulian/scData/Plot2DEmbedding.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(speckle)
require(extrafont)
require(ggplot2)
require(ggrastr)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp")

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Expand cluster attributes with genotype labels
hemi.cluster.attributes <- ExpandClusterAttributes(
  cluster.attributes = cluster.attributes$subclusters,
  labels = c("wt", "del"),
  alpha = c(1, 0.4),
  prefix = TRUE
)

# Define abbreviations of genotypes
genotype.abbreviations <- setNames(object = c("wt", "del"),
                                   nm = c("Nr4a2(WT/WT)", "Nr4a2(SA-IRES-Dre/WT)"))

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Load list of modulated genes
hemi.genes <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                       paste0(file.prefix, "_hemi_modulated_genes.rds")))

# Get list of observed modulated genes
hemi.observed <- lapply(
  X = hemi.genes,
  FUN = "[[",
  "observed"
)

# Count number of upregulated and downregulated genes
n.hemi.observed <- CountDEG(
  deg.list = hemi.observed,
  fc.col = "avg_log2FC",
  id.col = "cell_type_abbr"
)

# Melt data for plotting
n.hemi.observed <- reshape2::melt(data = n.hemi.observed,
                                  id.vars = "cell_type_abbr",
                                  variable.name = "modulation",
                                  value.name = "ngene")

# Order data to plot in order of cluster breaks
n.hemi.observed$cell_type_abbr <- factor(x = n.hemi.observed$cell_type_abbr,
                                         levels = cluster.attributes$subclusters$cluster.breaks)
n.hemi.observed <- n.hemi.observed[order(n.hemi.observed$cell_type_abbr, decreasing = TRUE),]
n.hemi.observed$cell_type_abbr <- as.character(x = n.hemi.observed$cell_type_abbr)

# Join layers of the RNA assay
if (!"RNA_full" %in% names(x = obj@assays)) {
  obj[["RNA_full"]] <- JoinLayers(object = obj[["RNA"]])
}

# Prepare data for violin plots
obj.wt <- subset(x = obj,
                 subset = genotype == "Nr4a2(WT/WT)")

# Add abbreviated genotypes to Seurat object meta data
obj@meta.data$genotype_abbr <- genotype.abbreviations[obj@meta.data$genotype]

# Create a column combining genotype and cell type identities of cells 
obj@meta.data$cell_type_genotype_abbr <- paste(obj@meta.data$cell_type_abbr,
                                               obj@meta.data$genotype_abbr)

# Define list of UMAP plots to produce
umap.params.l <- list(
  "all_cells" = NULL,
  "nr4a2_wt" = list(genotype = "Nr4a2(WT/WT)",
                    title = "*Nr4a2<sup>wt/wt</sup>*"),
  "nr4a2_hemi" = list(genotype = "Nr4a2(SA-IRES-Dre/WT)",
                      title = "*Nr4a2<sup>del/wt</sup>*")
)

# Define parameters for modulated genes violin and fraction plots
hemi.params.l <- list(
  "downregulated_subset" = list(genes = c("Nr4a2",
                                          "Cdh13",
                                          "Nr2f2",
                                          "Rxfp1",
                                          "Scn1b"),
                                subset = c("CLA", "shell", "L6a"),
                                subset.col = "broad_cell_type_abbr",
                                cluster.ident.name = "cell_type_genotype_abbr",
                                file.suffix = "downregulated_subset.pdf"),
  "upregulated_subset" = list(genes = c("Syt17",
                                        "Ntm",
                                        "Ryr2"),
                              subset = c("CLA", "shell", "L6a"),
                              subset.col = "broad_cell_type_abbr",
                              cluster.ident.name = "cell_type_genotype_abbr",
                              file.suffix = "upregulated_subset.pdf")
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
      cluster.data = obj@meta.data[, "cell_type_abbr", drop = FALSE],
      cellnames.as.rownames.in.embedding.data = TRUE,
      cellnames.as.rownames.in.cluster.data = TRUE,
      embedding.x.name = "UMAP_1",
      embedding.y.name = "UMAP_2",
      embedding.x.label = "UMAP 1",
      embedding.y.label = "UMAP 2",
      cluster.ident.name = "cell_type_abbr",
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
    
    # Add labels or title to UMAP plot depending if data was subsetted
    if (is.null(x = umap.params)) {
      # Add labels to UMAP plot
      umap.plot <- umap.plot +
        # ggtext::geom_richtext(
        geom_text(
          data = GetLabelCentroid(label.data = umap.plot$data,
                                  label.name = "cell_type_abbr",
                                  embedding.x.name = "UMAP_1",
                                  embedding.y.name =  "UMAP_2"),
          mapping = aes(x = UMAP_1,
                        y = UMAP_2,
                        label = cell_type_abbr),
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
                  file.name = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
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
  cluster.data = obj@meta.data[, "cell_type_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "cell_type_abbr",
  cluster.colors = cluster.attributes$subclusters$cluster.colors,
  cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
  cluster.labels = setNames(object = cluster.attributes$subclusters$cluster.breaks,
                            nm = cluster.attributes$subclusters$cluster.breaks),
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
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_violin_plots_cell_metrics.pdf")),
               is.ggassemble = TRUE,
               panel.width = 1,
               panel.height = 0.325 * length(x = cluster.attributes$subclusters$cluster.breaks),
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
                "Nr4a2", "Gnb4", "Oprk1", "Lxn", "Smim32", "Slc17a6", # CLA
                "Fosl2", "Egr2", # CLA 2 and shell 4
                "Nfib", "Sdk2", # shell and L6
                "Syt17", "Col6a1", "Nnat", "Npsr1", "Npy", # shell
                "Rprm", "Crym", # L6a
                "Ccn2", # L6b
                "Rorb", # L2/3
                "Tshz2", # L5b
                "Tafa1", # pir
                "Sst", # IN
                "Scn4b", # MSN
                "Drd1", # dMSN
                "Drd2", # iMSN
                "Tshz1" # sMSN
  ),
  cluster.data = obj.wt@meta.data[, "cell_type_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "cell_type_abbr",
  cluster.colors = cluster.attributes$subclusters$cluster.colors,
  cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
  cluster.labels = setNames(object = cluster.attributes$subclusters$cluster.breaks,
                            nm = cluster.attributes$subclusters$cluster.breaks),
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
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_violin_plots_cluster_marker_genes.pdf")),
               is.ggassemble = TRUE,
               panel.width = 0.5,
               panel.height = 0.325 * length(x = cluster.attributes$subclusters$cluster.breaks),
               margin = 2.25,
               unit.use = "cm",
               use.ggsave = TRUE,
               useDingbats = FALSE)

# Plot violin plots of Macaque cluster markers (doi: 10.1016/j.cell.2023.06.009)
markers.vln.plots <- StackViolinPlots(
  data.use = LayerData(object = obj.wt,
                       layer = "data",
                       assay = "RNA_full"),
  genes.use = c("Gnb4", "Nr4a2", # C1,3,5,8,9
                "Synpr", # C1-9
                "Nnat", # C2,4,6,7
                "Trpc3", "Cdc14a", # C2,4
                "Stk32a", # C6
                "Kank4", # C7
                "Smyd1", # C1,3,5,8,9
                "Ephx4", # C9
                "Slc17a8", # C1,9 and human
                "Sntb1", # C5,8
                "Postn", # C3,5,8
                "Tnnt1", "Adra1b", "Dcstamp", "Tfap2d" # human
  ),
  cluster.data = obj.wt@meta.data[, "cell_type_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "cell_type_abbr",
  cluster.colors = cluster.attributes$subclusters$cluster.colors,
  cluster.breaks = cluster.attributes$subclusters$cluster.breaks,
  cluster.labels = setNames(object = cluster.attributes$subclusters$cluster.breaks,
                            nm = cluster.attributes$subclusters$cluster.breaks),
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
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_violin_plots_primate_claustrum_marker_genes.pdf")),
               is.ggassemble = TRUE,
               panel.width = 0.5,
               panel.height = 0.325 * length(x = cluster.attributes$subclusters$cluster.breaks),
               margin = 2.25,
               unit.use = "cm",
               use.ggsave = TRUE,
               useDingbats = FALSE)

# Plot violin plots of modulated genes
hemi.vln.plots.l <- lapply(
  X = hemi.params.l,
  FUN = function(hemi.params,
                 obj,
                 cluster.attributes,
                 save.plot,
                 file.prefix) {
    # Get cluster attributes
    cluster.breaks <- cluster.attributes$cluster.breaks
    cluster.colors <- cluster.attributes$cluster.colors
    cluster.alpha <- cluster.attributes$cluster.alpha
    
    if (!is.null(x = hemi.params$subset)) {
      # Define cells to keep
      cells.to.keep <- rownames(x = obj@meta.data[obj@meta.data[[hemi.params$subset.col]] %in% hemi.params$subset,])
      
      # Subset Seurat object
      obj <- subset(x = obj,
                    cells = cells.to.keep)
    }
    
    # Check which clusters are present in Seurat object
    clusters.in.object <- unique(x = obj@meta.data[[hemi.params$cluster.ident.name]])
    if (any(!cluster.breaks %in% clusters.in.object)) {
      # Subset genotype x cluster breaks, colors and alpha
      cluster.breaks <- cluster.breaks[cluster.breaks %in% clusters.in.object]
      cluster.colors <- cluster.colors[cluster.breaks]
      cluster.alpha <- cluster.alpha[cluster.breaks]
    }
    
    # Plot violin plots of modulated genes
    hemi.vln.plots <- StackViolinPlots(
      data.use = LayerData(object = obj,
                           layer = "data",
                           assay = "RNA_full"),
      genes.use = hemi.params$genes,
      cluster.data = obj@meta.data[, hemi.params$cluster.ident.name, drop = FALSE],
      cellnames.as.rownames.in.cluster.data = TRUE,
      cluster.ident.name = hemi.params$cluster.ident.name,
      cluster.colors = cluster.colors,
      cluster.breaks = cluster.breaks,
      cluster.labels = setNames(object = cluster.breaks,
                                nm = cluster.breaks),
      is.log.transformed = TRUE,
      log.scale = "log",
      pseudocount.use = 1,
      y.scale.trans = "log10",
      y.min = "0",
      point.size = 2,
      alpha.use = cluster.alpha,
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
    if (save.plot) {
      FixSizeAndSave(plot = hemi.vln.plots,
                     filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                          paste0(file.prefix, "_violin_plots_hemi_", hemi.params$file.suffix)),
                     is.ggassemble = TRUE,
                     panel.width = 0.5,
                     panel.height = 0.325 * length(x = cluster.breaks),
                     margin = 2.25,
                     unit.use = "cm",
                     use.ggsave = TRUE,
                     useDingbats = FALSE)
    }
    
    # Return plot
    return(hemi.vln.plots)
  },
  obj = obj,
  cluster.attributes = hemi.cluster.attributes,
  save.plot = TRUE,
  file.prefix = file.prefix
)

# Plot number of upregulated and downregulated genes per cluster
n.observed.plot <- ggplot(data = n.hemi.observed,
                          mapping = aes(x = modulation,
                                        y = ngene,
                                        colour = cell_type_abbr)) +
  geom_point(size = 2,
             shape = 19,
             stroke = 0,
             position = position_jitter(width = 0.25,
                                        height = 0,
                                        seed = 47)) +
  scale_x_discrete(breaks = c("downregulated", "upregulated")) +
  scale_y_continuous(limits = c(0, 50 * ceiling(x = max(n.hemi.observed$ngene) / 50)),
                     breaks = seq(from = 0,
                                  to = 50 * ceiling(x = max(n.hemi.observed$ngene) / 50),
                                  by = 50),
                     expand = c(0, 0)) +
  scale_colour_manual(values = cluster.attributes$subclusters$cluster.colors,
                      breaks = cluster.attributes$subclusters$cluster.breaks) +
  guides(colour = ggplot2::guide_legend(override.aes = list(size = 2,
                                                            alpha = 1,
                                                            stroke = 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "modulation",
       y = "# of significantly modulated genes\n (adjusted p value < 0.05)",
       colour = "cell types") +
  theme_classic() +
  theme(axis.text = element_text(size = 6,
                                 family = "Arial",
                                 face = "plain",
                                 colour = "black"),
        axis.title = element_text(size = 7,
                                  family = "Arial",
                                  face = "plain",
                                  colour = "black"),
        legend.text = element_text(size = 7,
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
FixSizeAndSave(plot = n.observed.plot,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_modulated_genes_per_cluster.pdf")),
               is.ggassemble = FALSE,
               panel.width = 4,
               panel.height = 4,
               unit.use = "cm",
               margin = 2.5,
               use.ggsave = TRUE,
               useDingbats = FALSE)

# Calculate and plot cell type proportions
cell.type.prop.l <- CellTypeProportions(
  x = NULL,
  clusters = factor(x = obj@meta.data$cell_type_abbr,
                    levels = cluster.attributes$subclusters$cluster.breaks),
  sample = factor(x = obj@meta.data$mouse),
  group = factor(x = obj@meta.data$genotype_abbr,
                 levels = c("del", "wt")),
  trend = FALSE,
  robust = TRUE,
  transform = "logit",
  return.plot = TRUE,
  markdown = FALSE,
  reference.group = "wt",
  pseudocount = 0.5,
  cluster.colors = cluster.attributes$subclusters$cluster.colors,
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
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_cell_type_proportion_ratios.pdf")),
               is.ggassemble = FALSE,
               panel.width = 0.45 * length(x = cluster.attributes$subclusters$cluster.breaks),
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
  scale_x_discrete(limits = cluster.attributes$subclusters$cluster.breaks) +
  scale_y_continuous(limits = c(0, 0.05 * ceiling(x = max(cell.type.prop.l$sample.prop$proportion) / 0.05)),
                     breaks = seq(from = 0,
                                  to = 0.05 * ceiling(x = max(cell.type.prop.l$sample.prop$proportion) / 0.05),
                                  by = 0.05),
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
  theme(axis.text.x = element_text(size = 6,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   family = "Arial",
                                   face = "plain",
                                   colour = cluster.attributes$subclusters$cluster.colors),
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
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_cell_type_proportions.pdf")),
               is.ggassemble = FALSE,
               panel.width = 0.45 * length(x = cluster.attributes$subclusters$cluster.breaks) * 2,
               panel.height = 4.5,
               unit.use = "cm",
               margin = 0,
               use.ggsave = FALSE,
               useDingbats = FALSE)
