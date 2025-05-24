# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/FindMarkersSCTV5.R")

# Required packages
require(Seurat)

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Prepare object to run differential expression on SCT assay with multiple models
obj.deg <- PrepSCTFindMarkers(object = obj,
                              assay = "SCT",
                              verbose = TRUE)

# Set SCT as default assay
DefaultAssay(object = obj.deg) <- "SCT"

# Keep only the SCT assay in the Seurat object
obj.deg <- DietSeurat(
  object = obj.deg,
  layers = NULL,
  features = NULL,
  assays = "SCT",
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)

# Get names of samples
samples <- unique(x = obj@meta.data$mouse)

# Define parameters for differential expression analysis
deg.params.l <- data.frame(
  permute.ident = rep(x = c(FALSE,
                            TRUE),
                      times = c(length(x = samples) + 1,
                                500)),
  remove.sample = rep(x = c(FALSE,
                            TRUE,
                            FALSE),
                      times = c(1,
                                length(x = samples),
                                500)),
  sample.to.remove = c(NA,
                       samples,
                       rep(x = NA,
                           times = 500)),
  sample.col = "mouse",
  row.names = c("observed",
                paste0("leave_one_out_",
                       samples),
                paste0("permuted_",
                       1:500))
)

# Compute down and upregulated genes in hemizygous mice for each cluster
hemi.markers.l <- lapply(
  X = setNames(object = cluster.attributes$subclusters$cluster.breaks,
               nm = cluster.attributes$subclusters$cluster.breaks),
  FUN = function(cell.type,
                 obj,
                 deg.params.l) {
    # Print DEG type to keep track of progress
    cat("Computing DEG for cluster",
        cell.type,
        "\n",
        sep = " ")
    
    # Subset Seurat object for a given cell type
    obj.ct <- subset(x = obj,
                     subset = cell_type_abbr == cell.type)
    
    # Compute differentially expressed genes between the two genotypes
    hemi.markers <- lapply(
      X = deg.params.l,
      FUN = FindMarkersSCTV5,
      obj = obj.ct,
      ident.1 = "Nr4a2(SA-IRES-Dre/WT)",
      ident.2 = "Nr4a2(WT/WT)",
      deg.ident = "genotype",
      slot = "data",
      fc.slot = "data",
      logfc.threshold = 0.25,
      test.use = "wilcox",
      min.pct = 0.1,
      min.diff.pct = -Inf,
      verbose = TRUE,
      only.pos = FALSE,
      max.cells.per.ident = Inf,
      random.seed = 1,
      latent.vars = NULL,
      min.cells.feature = 3,
      min.cells.group = 3,
      densify = FALSE,
      pseudocount.use = 1,
      norm.method = NULL,
      mean.fxn = NULL,
      fc.name = NULL,
      base = 2,
      recorrect_umi = FALSE
    )
    
    # Return list of data
    return(list("observed" = hemi.markers[["observed"]],
                "permuted" = hemi.markers[grepl(pattern = "permuted", x = names(x = hemi.markers))],
                "leave_one_out" = hemi.markers[grepl(pattern = "leave_one_out", x = names(x = hemi.markers))]))
  },
  obj = obj.deg,
  deg.params.l = split(x = deg.params.l,
                       f = rownames(x = deg.params.l))
)

# Save list of markers
saveRDS(object = hemi.markers.l,
        file = file.path("results/smim32_GTFv102/neurons_filtered",
                         paste0(file.prefix, "_hemi_modulated_genes.rds")))
