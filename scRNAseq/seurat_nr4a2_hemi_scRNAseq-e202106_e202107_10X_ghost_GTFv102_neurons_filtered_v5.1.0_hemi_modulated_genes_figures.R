# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/ExpandClusterAttributes.R")
source(file = "/Users/leonfodoulian/scData/BlendColorAlpha.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/StackFractionPlots.R")
source(file = "/Users/leonfodoulian/scData/GetFractionOfExpressingCells.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(extrafont)
require(ggplot2)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp/hemi_modulated_genes",
           recursive = TRUE)

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

# Join layers of the RNA assay
if (!"RNA_full" %in% names(x = obj@assays)) {
  obj[["RNA_full"]] <- JoinLayers(object = obj[["RNA"]])
}

# Add abbreviated genotypes to Seurat object meta data
obj@meta.data$genotype_abbr <- genotype.abbreviations[obj@meta.data$genotype]

# Create a column combining genotype and cell type identities of cells 
obj@meta.data$cell_type_genotype_abbr <- paste(obj@meta.data$cell_type_abbr,
                                               obj@meta.data$genotype_abbr)

# Plot violin and fraction plots of modulated genes
hemi.plots.l <- mapply(
  cluster.hemi.genes = hemi.genes,
  cell.type = gsub(pattern = "\\s|\\/",
                   replacement = "",
                   x = names(x = hemi.genes)),
  FUN = function(cluster.hemi.genes,
                 cell.type,
                 obj,
                 cluster.attributes,
                 save.plot,
                 path,
                 file.prefix) {
    # Get list of observed modulated genes
    cluster.hemi.genes <- cluster.hemi.genes$observed
    
    # Stop function if no modulated genes are found
    if (length(x = cluster.hemi.genes) == 0) {
      return(NULL)
    }
    
    # Add modulation to the data
    cluster.hemi.genes$modulation <- ifelse(test = cluster.hemi.genes$avg_log2FC > 0,
                                            yes = "upregulated",
                                            no = "downregulated")
    
    # Plot violin and fraction plots of modulated genes for each modulation
    hemi.plots.l <- lapply(
      X = split(x = cluster.hemi.genes,
                f = cluster.hemi.genes$modulation),
      FUN = function(modulation.genes,
                     obj,
                     cell.type,
                     cluster.attributes,
                     save.plot,
                     path,
                     file.prefix) {
        # Get modulation type
        modulation <- unique(x = modulation.genes$modulation)
        
        # Create new directories
        file.dir <- file.path(path,
                              cell.type,
                              modulation)
        dir.create(path = file.dir,
                   recursive = TRUE)
        
        # Create chunks of 35 genes per modulation and per graph
        modulation.genes$chunk <- ceiling(x = seq_len(length.out = nrow(x = modulation.genes)) / 35)
        
        # Get maximum number of chunks
        max.chunk <- max(modulation.genes$chunk)
        
        # Rename chunks in data
        modulation.genes$chunk <- paste0(modulation.genes$chunk,
                                         "_of_",
                                         max.chunk)
        
        # Plot violin and fraction plots of modulated genes for each chunk
        hemi.plots.l <- lapply(
          X = split(x = modulation.genes,
                    f = modulation.genes$chunk),
          FUN = function(chunk.modulation.genes,
                         obj,
                         cell.type,
                         modulation,
                         cluster.attributes,
                         save.plot,
                         file.dir,
                         file.prefix) {
            # Define file suffix
            file.suffix <- paste(cell.type,
                                 modulation,
                                 unique(x = chunk.modulation.genes$chunk),
                                 sep = "_")
            
            # Plot violin plots of modulated genes
            hemi.vln.plots <- StackViolinPlots(
              data.use = LayerData(object = obj,
                                   layer = "data",
                                   assay = "RNA_full"),
              genes.use = rownames(x = chunk.modulation.genes),
              cluster.data = obj@meta.data[, "cell_type_genotype_abbr", drop = FALSE],
              cellnames.as.rownames.in.cluster.data = TRUE,
              cluster.ident.name = "cell_type_genotype_abbr",
              cluster.colors = cluster.attributes$cluster.colors,
              cluster.breaks = cluster.attributes$cluster.breaks,
              cluster.labels = setNames(object = cluster.attributes$cluster.breaks,
                                        nm = cluster.attributes$cluster.breaks),
              is.log.transformed = TRUE,
              log.scale = "log",
              pseudocount.use = 1,
              y.scale.trans = "log10",
              y.min = "0",
              point.size = 2,
              alpha.use = cluster.attributes$cluster.alpha,
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
                             filename = file.path(file.dir,
                                                  paste0(file.prefix, "_violin_plots_hemi_", file.suffix, ".pdf")),
                             is.ggassemble = TRUE,
                             panel.width = 0.5,
                             panel.height = 0.325 * length(x = cluster.attributes$cluster.breaks),
                             margin = 2.25,
                             unit.use = "cm",
                             use.ggsave = TRUE,
                             useDingbats = FALSE)
            }
            
            # Plot fraction plots of modulated genes
            hemi.fraction.plots <- StackFractionPlots(
              data.use = LayerData(object = obj,
                                   layer = "counts",
                                   assay = "RNA_full"),
              genes.use = rownames(x = chunk.modulation.genes),
              cluster.data = obj@meta.data[, "cell_type_genotype_abbr", drop = FALSE],
              cellnames.as.rownames.in.cluster.data = TRUE,
              cluster.ident.name = "cell_type_genotype_abbr",
              cluster.colors = cluster.attributes$cluster.colors,
              cluster.breaks = cluster.attributes$cluster.breaks,
              cluster.labels = setNames(object = cluster.attributes$cluster.breaks,
                                        nm = cluster.attributes$cluster.breaks),
              is.log.transformed = FALSE,
              log.scale = "log",
              pseudocount.use = 1,
              min.exp = 1,
              thresh.plot = 15,
              alpha.use = cluster.attributes$cluster.alpha,
              bar.border.colour = NA,
              bar.border.stroke = 0.25,
              plot.title.size = 7,
              plot.title.angle = 45,
              plot.title.face = "italic",
              font.family = "Arial",
              import.font = FALSE,
              hjust.use = 0,
              vjust.use = 0,
              x.axis.title = "Fraction of \nexpressing cells (%)",
              axis.text.x.size = 6,
              axis.text.y.size = 7,
              axis.text.face = "plain",
              axis.ticks.length = 0.5,
              axis.line.size = rel(0.5),
              verbose = TRUE)
            
            # Save plot in pdf format
            if (save.plot) {
              FixSizeAndSave(plot = hemi.fraction.plots,
                             filename = file.path(file.dir,
                                                  paste0(file.prefix, "_fraction_plots_hemi_", file.suffix, ".pdf")),
                             is.ggassemble = TRUE,
                             panel.width = 0.5,
                             panel.height = 0.325 * length(x = cluster.attributes$cluster.breaks),
                             margin = 2.25,
                             unit.use = "cm",
                             use.ggsave = TRUE,
                             useDingbats = FALSE)
            }
            
            # Return list of plot
            return(list(hemi.vln.plots = hemi.vln.plots,
                        hemi.fraction.plots = hemi.fraction.plots))
          },
          obj = obj,
          cell.type = cell.type,
          modulation = modulation,
          cluster.attributes = cluster.attributes,
          save.plot = save.plot,
          file.dir = file.dir,
          file.prefix = file.prefix)
        
        # Return list of plot
        return(hemi.plots.l)
      },
      obj = obj,
      cell.type = cell.type,
      cluster.attributes = cluster.attributes,
      save.plot = save.plot,
      path = path,
      file.prefix = file.prefix)
    
    # Return list of plot
    return(hemi.plots.l)
  },
  MoreArgs = list(obj = obj,
                  cluster.attributes = hemi.cluster.attributes,
                  save.plot = TRUE,
                  path = "./results/smim32_GTFv102/neurons_filtered/figures/temp/hemi_modulated_genes",
                  file.prefix = file.prefix),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# Save list of modulated genes in xlsx format
hemi.genes.observed <- lapply(
  X = hemi.genes,
  FUN = function(modulated.genes) {
    modulated.genes <- modulated.genes$observed
    modulated.genes$gene <- rownames(x = modulated.genes)
    rownames(x = modulated.genes) <- NULL
    modulated.genes <- modulated.genes[c("gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]
    colnames(x = modulated.genes) <- c("gene", "log2_fold_change", "p_value", "adjusted_p_value", "percentage_del", "percentage_wt")
    return(modulated.genes)
  })
hemi.genes.observed <- hemi.genes.observed[lapply(X = hemi.genes.observed, FUN = nrow) > 0]
names(x = hemi.genes.observed) <- gsub(pattern = "\\/",
                                       replacement = "-",
                                       x = paste0("modulated genes in ",
                                                  names(x = hemi.genes.observed)))
writexl::write_xlsx(x = hemi.genes.observed,
                    path = file.path("results/smim32_GTFv102/neurons_filtered",
                                     paste0(file.prefix, "_hemi_modulated_genes.xlsx")))
