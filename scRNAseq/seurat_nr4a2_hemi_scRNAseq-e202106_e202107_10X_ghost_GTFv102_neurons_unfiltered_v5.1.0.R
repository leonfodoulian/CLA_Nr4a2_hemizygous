# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_mice_data.R")
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/RunSeuratV5.R")
source(file = "/Users/leonfodoulian/scData/RKneeLocator.R")
source(file = "/Users/leonfodoulian/scData/ClusteringPairVI.R")

# Required packages
require(Seurat)

# Create new directories
dir.create(path = "./results/smim32_GTFv102")
dir.create(path = "./results/smim32_GTFv102/neurons_unfiltered")
dir.create(path = "./results/smim32_GTFv102/neurons_unfiltered/figures")

# Load Seurat object
obj.all.cells <- readRDS(file = "results/smim32_GTFv102/all_cells/scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_all_cells.rds")

# Subset Seurat object for neuronal cells
obj.neurons <- subset(x = obj.all.cells,
                      subset = broad_cell_type_abbr %in% c("CLA", "shell-like", "L6", "IN", "MSN"))

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
rm(list = c("obj.all.cells", "obj.neurons"))

# Run Seurat version 5 pipeline with Harmony integration
obj.merged <- RunSeuratV5(
  counts.mat.l = counts.mat.l,
  sample.data = mice.data,
  sample.col = "mouse",
  project = "neurons_unfiltered",
  resolution = seq(from = 0.1,
                   to = 1,
                   by = 0.1),
  cluster.attributes = cluster.attributes$neurons_unfiltered,
  azimuth = FALSE,
  spca.features = NULL,
  file.prefix = "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger",
  path = "results/smim32_GTFv102"
)
