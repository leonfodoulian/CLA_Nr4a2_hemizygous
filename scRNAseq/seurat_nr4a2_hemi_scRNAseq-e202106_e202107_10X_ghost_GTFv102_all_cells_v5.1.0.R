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
dir.create(path = "./results/smim32_GTFv102/all_cells")
dir.create(path = "./results/smim32_GTFv102/all_cells/figures")

# Get 10X filtered output count file directory names
path.l <- c("scRNAseq-e20210604_e20210611/cellranger_count_ghost_GTFv102",
            "scRNAseq-e20210719_e20210720/cellranger_count_ghost_GTFv102")
dirs <- unlist(x = lapply(X = path.l,
                          FUN = function(path) {
                            samples <- list.dirs(path = file.path(".", path),
                                                 full.names = FALSE,
                                                 recursive = FALSE)
                            dirs <- file.path(getwd(), path, samples, "outs/filtered_feature_bc_matrix")
                            names(x = dirs) <- samples
                            return(dirs)
                          }))

# Load data using Read10X
counts.mat.l <- mapply(dir = dirs,
                       sample.name = names(x = dirs),
                       FUN = function(dir,
                                      sample.name) {
                         # Load data
                         counts.mat <- Read10X(data.dir = dir,
                                               gene.column = 2,
                                               cell.column = 1,
                                               unique.features = TRUE,
                                               strip.suffix = FALSE)
                         
                         # Add sample name in cell names
                         colnames(x = counts.mat) <- paste("s", sample.name, "_",
                                                           colnames(x = counts.mat),
                                                           sep = "")
                         
                         # Rename Ghost as Smim32
                         rownames(x = counts.mat)[rownames(x = counts.mat) == "Ghost"] <- "Smim32"
                         
                         # Return data
                         return(counts.mat)
                       },
                       SIMPLIFY = FALSE,
                       USE.NAMES = TRUE)

# Run Seurat version 5 pipeline with Harmony integration
obj.merged <- RunSeuratV5(
  counts.mat.l = counts.mat.l,
  sample.data = mice.data,
  sample.col = "mouse",
  project = "all_cells",
  resolution = seq(from = 0.1,
                   to = 1,
                   by = 0.1),
  cluster.attributes = cluster.attributes$all_cells,
  azimuth = FALSE,
  spca.features = NULL,
  file.prefix = "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger",
  path = "results/smim32_GTFv102"
)
