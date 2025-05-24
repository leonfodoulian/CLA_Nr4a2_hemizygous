# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/ExpandClusterAttributes.R")
source(file = "/Users/leonfodoulian/scData/BlendColorAlpha.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(extrafont)
require(ggplot2)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp")

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Expand cluster attributes with genotype labels
hemi.cluster.attributes <- ExpandClusterAttributes(
  cluster.attributes = cluster.attributes$broad.clusters,
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

# Subset Seurat object for CLA cell types
obj <- subset(x = obj,
              subset = broad_cell_type_abbr == "CLA")

# Add abbreviated genotypes to Seurat object meta data
obj@meta.data$genotype_abbr <- genotype.abbreviations[obj@meta.data$genotype]

# Create a column combining genotype and broad cell type identities of cells 
obj@meta.data$broad_cell_type_genotype_abbr <- paste(obj@meta.data$broad_cell_type_abbr,
                                                     obj@meta.data$genotype_abbr)

# Define parameters for channels/receptors violin plots
vln.plot.params.l <- list(
  "main_figure" = list(
    genes = c(
      "Kcnv1", "Kcnc1", "Trpc4", "Scn1b", # downregulated channel
      "Npffr1", # downregulated receptor
      "Kctd16", # downregulated modulator of GABAbR
      "Gnao1", "Gnb4", "Gng2", # downregulated G protein subunit
      "Dgkg", # downregulated second messenger modulator
      "Kcnf1", "Kcnj3", "Kcnj4", "Scn8a", "Cckbr", # upregulated channel
      "Cnr1", "Htr1f", "Htr2c", "Itpr1", "Gabbr2", "Gabra3", # upregulated receptor
      "Gabrb1", "Gpr26", "Gria3", "Grid2", "Grin2b", # upregulated receptor
      "Ryr2", # upregulated channel
      "Gnal", "Rgs4", # upregulated G protein subunit
      "Dgkb", "Dgkz" # upregulated second messenger modulator
    ),
    file.suffix = "main_figure.pdf"),
  "supp_figure" = list(
    genes = c(
      "Scn1a", "Scn2a", "Scn3a", "Scn4a", "Scn5a", "Scn7a",
      "Scn8a", "Scn9a", "Scn10a", "Scn11a", "Scn1b", "Kcna1",
      "Kcna2", "Kcna3", "Kcna4", "Kcna5", "Kcna6", "Kcna7",
      "Kcna10", "Kcnc1", "Kcnc2", "Kcnc3", "Kcnc4", "Kcnd1",
      "Kcnd2", "Kcnd3", "Kcnj1", "Kcnj2", "Kcnj3", "Kcnj4",
      "Kcnj5", "Kcnj6", "Kcnj8", "Kcnj9", "Kcnj10", "Kcnj11",
      "Kcnj12", "Kcnj13", "Kcnj14", "Kcnj15", "Kcnj16", "Kcnma1",
      "Kcnmb1", "Kcnmb2", "Kcnmb3", "Kcnmb4", "Kcnn1", "Kcnn2",
      "Kcnn3", "Kcnn4", "Kcnq1", "Kcnq2", "Cacna1s", "Cacna1c",
      "Cacna1d", "Cacna1f", "Cacna1a", "Cacna1b", "Cacna1e", "Cacna1g",
      "Cacna1h", "Cacna1i"
    ),
    file.suffix = "supp_figure.pdf")
)

# Plot violin plots of channels/receptors
vln.plots.l <- lapply(
  X = vln.plot.params.l,
  FUN = function(vln.plot.params,
                 obj,
                 cluster.attributes,
                 save.plot,
                 file.prefix) {
    # Plot violin plots of channels/receptors
    vln.plots <- StackViolinPlots(
      data.use = LayerData(object = obj,
                           layer = "data",
                           assay = "RNA_full"),
      genes.use = vln.plot.params$genes,
      cluster.data = obj@meta.data[, "broad_cell_type_genotype_abbr", drop = FALSE],
      cellnames.as.rownames.in.cluster.data = TRUE,
      cluster.ident.name = "broad_cell_type_genotype_abbr",
      cluster.colors = cluster.attributes$cluster.colors[1:2],
      cluster.breaks = cluster.attributes$cluster.breaks[1:2],
      cluster.labels = setNames(object = cluster.attributes$cluster.breaks[1:2],
                                nm = cluster.attributes$cluster.breaks[1:2]),
      is.log.transformed = TRUE,
      log.scale = "log",
      pseudocount.use = 1,
      y.scale.trans = "log10",
      y.min = "0",
      point.size = 2,
      alpha.use = cluster.attributes$cluster.alpha[1:2],
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
      FixSizeAndSave(plot = vln.plots,
                     filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                          paste0(file.prefix, "_violin_plots_ephys_", vln.plot.params$file.suffix)),
                     is.ggassemble = TRUE,
                     panel.width = 0.5,
                     panel.height = 0.325 * 2,
                     margin = 2.25,
                     unit.use = "cm",
                     use.ggsave = TRUE,
                     useDingbats = FALSE)
    }
    
    # Return plot
    return(vln.plots)
  },
  obj = obj,
  cluster.attributes = hemi.cluster.attributes,
  save.plot = TRUE,
  file.prefix = file.prefix
)
