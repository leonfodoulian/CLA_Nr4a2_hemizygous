# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_mice_data.R")
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/RunSeuratV5.R")
source(file = "/Users/leonfodoulian/scData/RKneeLocator.R")
source(file = "/Users/leonfodoulian/scData/ClusteringPairVI.R")
source(file = "/Users/leonfodoulian/scData/FilterUnusualCells.R")

# Required packages
require(Seurat)

# Create new directories
dir.create(path = "./results/smim32_GTFv102")
dir.create(path = "./results/smim32_GTFv102/neurons_filtered")
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures")
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/azimuth_Claustrum300_reference")

# Load list of genes for sPCA calculation
spca.features <- read.table(file = "MERFISH_data/MERFISH_Rodriguez_Claustrum300_genes.txt",
                            header = TRUE)[[1]]

# Load Seurat object
obj.ufneurons <- readRDS(file = "results/smim32_GTFv102/neurons_unfiltered/scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-20_neurons_unfiltered.rds")

# Subset Seurat object for neuronal cells
obj.neurons <- subset(x = obj.ufneurons,
                      subset = broad_cell_type_abbr %in% c("CLA", "shell", "L6a", "L6b", "L2/3-L5b-pir", "IN", "dMSN", "iMSN", "sMSN"))

# Get names of counts layers from Seurat object
layers <- grep(pattern = "counts",
               x = names(x = obj.neurons@assays$RNA@layers),
               value = TRUE)
names(x = layers) <- gsub(pattern = ".*\\.",
                          replacement = "",
                          x = layers)

# Get counts data from Seurat object
counts.mat.l <- lapply(X = layers,
                       FUN = function(layer,
                                      obj) {
                         LayerData(object = obj,
                                   layer = layer,
                                   assay = "RNA")
                       },
                       obj = obj.neurons)

# Remove loaded Seurat objects
rm(list = c("obj.ufneurons", "obj.neurons"))

# Filter out unusual cells
unusual.cells <- FilterUnusualCells(mat = Reduce(f = cbind,
                                                 x = counts.mat.l),
                                    is.expr = 0,
                                    mito.genes.pattern = "mt-",
                                    nGene.thresh = 1000,
                                    percent.mito.thresh = 0.1,
                                    loess.span = 0.5,
                                    loess.degree = 2,
                                    pval.thresh = 0.05,
                                    use.median = TRUE,
                                    verbose = TRUE)

# Filter count matrices
counts.mat.l <- lapply(X = counts.mat.l,
                       FUN = function(counts.mat,
                                      outlier.cells) {
                         # Remove outlier cells from the data
                         counts.mat <- counts.mat[,!colnames(x = counts.mat) %in% outlier.cells]
                         
                         # Return data
                         return(counts.mat)
                       },
                       outlier.cells = unusual.cells$outlier.cells)

# Run Seurat version 5 pipeline with Harmony integration
obj.merged <- RunSeuratV5(
  counts.mat.l = counts.mat.l,
  sample.data = mice.data,
  sample.col = "mouse",
  project = "neurons_filtered",
  resolution = seq(from = 0.1,
                   to = 2,
                   by = 0.1),
  cluster.attributes = cluster.attributes$neurons_filtered,
  azimuth = TRUE,
  spca.features = spca.features,
  file.prefix = "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger",
  path = "results/smim32_GTFv102"
)

# Match full and abbreviated names for broad and sub-clusters
broad.cell.type.id <- match(
  x = cluster.attributes$neurons_filtered$broad.clusters$cluster.breaks,
  table = cluster.attributes$neurons_filtered$broad.clusters$cluster.abbreviations
)
cell.type.id <- match(
  x = cluster.attributes$neurons_filtered$subclusters$cluster.breaks,
  table = cluster.attributes$neurons_filtered$subclusters$cluster.abbreviations
)

# Prepare color map for Azimuth data
colormap <- list(
  "broad_cell_type" = setNames(
    object = cluster.attributes$neurons_filtered$broad.clusters$cluster.colors,
    nm = cluster.attributes$neurons_filtered$broad.clusters$cluster.names[broad.cell.type.id]
  ),
  "broad_cell_type_abbr" = cluster.attributes$neurons_filtered$broad.clusters$cluster.colors,
  "cell_type" = setNames(
    object = cluster.attributes$neurons_filtered$subclusters$cluster.colors,
    nm = cluster.attributes$neurons_filtered$subclusters$cluster.names[cell.type.id]
  ),
  "cell_type_abbr" = cluster.attributes$neurons_filtered$subclusters$cluster.colors
)

# Create a Seurat object compatible with Azimuth
obj.ref <- Azimuth::AzimuthReference(
  object = obj.merged,
  refUMAP = "umap_harmony",
  refDR = "spca",
  refAssay = "SCT_full",
  dims = 1:obj.merged@misc$spca.kneed.res$knee,
  k.param = 31,
  plotref = "umap_harmony",
  plot.metadata = NULL,
  ori.index = NULL,
  colormap = colormap,
  assays = "SCT_full",
  metadata = intersect(x = c("broad_cell_type", "broad_cell_type_abbr", "cell_type", "cell_type_abbr"),
                       y = colnames(x = obj.merged@meta.data)),
  reference.version = "0.0.0",
  verbose = TRUE
)

# Save Azimuth references and neighbors index
Azimuth::SaveAzimuthReference(object = obj.ref,
                              folder = "results/smim32_GTFv102/neurons_filtered/azimuth_Claustrum300_reference/")
