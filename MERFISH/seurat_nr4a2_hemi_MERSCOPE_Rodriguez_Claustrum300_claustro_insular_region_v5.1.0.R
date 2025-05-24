# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/spatial_registration_helper_functions.R")
source(file = "/Users/leonfodoulian/scData/PreprocessMERFISHSeuratV5.R")
source(file = "/Users/leonfodoulian/scData/RunSeuratV5.R")
source(file = "/Users/leonfodoulian/scData/RKneeLocator.R")
source(file = "/Users/leonfodoulian/scData/ClusteringPairVI.R")

# Required packages
require(Seurat)

# Set options for future globals maxSize and Azimuth ndims
options(future.globals.maxSize = 1e16,
        Azimuth.map.ndims = 14)

# Create new directories
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/all_cells")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered")
dir.create(path = "./results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/figures")

# Load mice data
mice.data <- as.data.frame(
  x = readxl::read_xlsx(path = "MERFISH_data/MERSCOPE/MERFISH_mice.xlsx",
                        na = "NaN")
)
rownames(x = mice.data) <- mice.data$mouse

# Load rotation angle data
rotation.angles <- as.data.frame(
  x = readxl::read_xlsx(path = "MERFISH_data/MERSCOPE/MERFISH_images_rotations_angles.xlsx",
                        na = "NaN")
)

# List directories of each experiment
dirs <- list.dirs(path = "./MERFISH_data/MERSCOPE",
                  full.names = TRUE,
                  recursive = FALSE)
names(x = dirs) <- gsub(pattern = ".*\\/",
                        replacement = "",
                        x = dirs)

# Preprocess MERSCOPE data
merscope.data.l <- mapply(
  dir = dirs,
  experiment = names(x = dirs),
  FUN = function(dir,
                 experiment,
                 sample.data,
                 rotation.angles) {
    # List files to load
    cell.metadata.file <- list.files(path = dir,
                                     pattern = "cell.metadata",
                                     full.names = TRUE) # cell meta data
    counts.mat.file <- list.files(path = dir,
                                  pattern = "cell.by.gene",
                                  full.names = TRUE) # count data
    cell.id.files <- list.files(path = dir,
                                pattern = "cell.IDs",
                                full.names = TRUE) # cell ID
    
    # Get attributes from experiment
    file.attributes <- unlist(x = strsplit(x = experiment, split = "_"))
    names(x = file.attributes) <- c("experiment", "mouse", "genotype", "level")[1:length(x = file.attributes)]
    if (!"level" %in% names(x = file.attributes)) {
      file.attributes["level"] <- "44"
    }
    
    # Load cell meta data
    cell.metadata <- read.csv(file = cell.metadata.file,
                              header = TRUE,
                              check.names = FALSE,
                              colClasses = c("EntityID" = "character"))
    
    # Subset cell meta data for columns of interest
    cell.metadata <- cell.metadata[,c("EntityID", "volume", "center_x", "center_y")]
    
    # Add sample name to cell names of meta data
    cell.metadata$cell <- paste(file.attributes[["mouse"]],
                                cell.metadata$EntityID,
                                sep = "_")
    
    # Add cell names to rownames of meta data
    rownames(x = cell.metadata) <- cell.metadata$cell
    
    # Add file attributes to the meta data
    cell.metadata[names(x = file.attributes)] <- as.list(x = file.attributes)
    
    # Add sample data to meta data
    sample.data <- sample.data[file.attributes[["mouse"]],,drop = FALSE]
    sample.data <- sample.data[setdiff(x = colnames(x = sample.data),
                                     y = colnames(x = cell.metadata))]
    cell.metadata[names(x = sample.data)] <- sample.data
    
    # Get rotation angle
    angle <- rotation.angles$`rotation angle`[rotation.angles$experiment == experiment]
    
    # Rotate x and y coordinates
    rotated.coords <- RotateCoordinates(xy = cell.metadata[,c("center_x", "center_y")],
                                        angle = angle * pi / 180)
    colnames(x = rotated.coords) <- c("centroid_x_rot", "centroid_y_rot")
    
    # Add rotated coordinates to meta data
    cell.metadata <- cbind(cell.metadata,
                           rotated.coords)
    
    # Select cells from the clausto-insular region
    cla.ins.cells.l <- lapply(
      X = cell.id.files,
      FUN = function(cell.id.file) {
        # Load cells to keep
        cla.ins.cells <- read.csv(file = cell.id.file,
                                  header = TRUE,
                                  check.names = FALSE,
                                  colClasses = "character")[[1]]
        
        # Get brain hemisphere information
        cla.side <- gsub(pattern = ".*\\/|\\.cell.IDs.*",
                         replacement = "",
                         x = cell.id.file)
        
        # Return data
        return(data.frame(cla_ins_cells = cla.ins.cells,
                          cla_side = cla.side))
      }
    )
    
    # Indicate in meta data whether cell is located in the clausto-insular region
    for (i in 1:length(x = cla.ins.cells.l)) {
      cell.metadata[cell.metadata$EntityID %in% cla.ins.cells.l[[i]][["cla_ins_cells"]],
                    "clausto_insular_region"] <- TRUE
      cell.metadata[cell.metadata$EntityID %in% cla.ins.cells.l[[i]][["cla_ins_cells"]],
                    "cla_side"] <- unique(x = cla.ins.cells.l[[i]][["cla_side"]])
    }
    cell.metadata[is.na(x = cell.metadata$clausto_insular_region), "clausto_insular_region"] <- FALSE
    cell.metadata[is.na(x = cell.metadata$cla_side), "cla_side"] <- "none"
    
    # Load counts data
    counts.mat <- read.csv(file = counts.mat.file,
                           header = TRUE,
                           check.names = FALSE,
                           colClasses = c("cell" = "character",
                                          "EntityID" = "character"))
    colnames(x = counts.mat)[1] <- "cell"
    
    # Add sample name to cell names of counts data
    counts.mat$cell <- paste(file.attributes[["mouse"]],
                             counts.mat$cell,
                             sep = "_")
    
    # Add cell names to rownames of counts data
    rownames(x = counts.mat) <- counts.mat$cell
    counts.mat$cell <- NULL # remove cell column
    
    # Transpose counts data
    counts.mat <- t(x = counts.mat)
    
    # Rename Ghost as Smim32
    rownames(x = counts.mat)[rownames(x = counts.mat) == "Ghost"] <- "Smim32"
    
    # Remove Blank rows
    counts.mat <- counts.mat[!grepl(pattern = "Blank", x = rownames(x = counts.mat)),]
    
    # Add number of genes and counts (excluding Blank values) to meta data
    cell.metadata$ncount <- colSums(x = counts.mat)
    cell.metadata$ngene <- colSums(x = counts.mat > 0)
    
    # Indicate whether cell has volume within norms for neurons
    cell.metadata$volume_to_keep <- cell.metadata$volume >= 500 & cell.metadata$volume <= 3500
    
    # Indicate whether cell has number of counts within norms
    cell.metadata$ncount_to_keep <- cell.metadata$ncount >= 10 & cell.metadata$ncount <= 3000
    
    # Keep cells that meet the following criteria
    cells.to.keep <- cell.metadata$clausto_insular_region & # cell in clausto-insular region
      cell.metadata$volume_to_keep & # volume within norm defined above
      cell.metadata$ncount_to_keep # number of counts within norm defined above
    
    # Preprocess MERFISH data
    obj <- PreprocessMERFISHSeuratV5(
      counts.mat = counts.mat[,cells.to.keep],
      cell.metadata = cell.metadata[cells.to.keep,],
      sample.name = file.attributes[["mouse"]],
      assay.name = "MERSCOPE"
    )
    
    # Return list of data
    return(list(merscope.data = list(counts.mat = counts.mat,
                                     cell.metadata = cell.metadata,
                                     cells.to.keep = cells.to.keep),
                seurat.object = obj))
  },
  MoreArgs = list(sample.data = mice.data,
                  rotation.angles = rotation.angles),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# Save pre-processed MERSCOPE data
merscope.data <- lapply(
  X = merscope.data.l,
  FUN = "[[",
  "merscope.data"
)
saveRDS(object = merscope.data,
        file = "results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region_preprocessed_data.rds")

# Get list of Seurat objects
obj.l <- lapply(
  X = merscope.data.l,
  FUN = "[[",
  "seurat.object"
)

# Merge Seurat objects
obj.merged <- merge(
  x = obj.l[[1]],
  y = obj.l[-1],
  add.cell.ids = NULL,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = TRUE
)

# Set all genes as variable features
obj.merged@assays$SCT@var.features <- rownames(x = obj.merged@assays$SCT@scale.data)

# Run Seurat version 5 pipeline with Harmony integration
obj.merged <- RunSeuratV5(
  obj.merged = obj.merged,
  sample.col = "mouse",
  project = "all_cells",
  resolution = seq(from = 0.1,
                   to = 1,
                   by = 0.1),
  azimuth = FALSE,
  spca.features = NULL,
  file.prefix = "MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region",
  path = "results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region"
)

# Subset Seurat object for neuronal cells
obj.merged <- subset(x = obj.merged,
                     subset = SCT_snn_res.0.1 %in% c("1", # shell-like
                                                     "2", # MSN
                                                     "5", # L6
                                                     "6", # IN
                                                     "8", # CLA
                                                     "11", # L5a
                                                     "12")) # L5b

# Get names of counts layers from Seurat object
layers <- grep(pattern = "counts",
               x = names(x = obj.merged@assays$MERSCOPE@layers),
               value = TRUE)
names(x = layers) <- gsub(pattern = ".*\\counts.",
                          replacement = "",
                          x = layers)

# Get counts data from Seurat object
counts.mat.l <- lapply(X = layers,
                       FUN = function(layer,
                                      obj) {
                         LayerData(object = obj,
                                   layer = layer,
                                   assay = "MERSCOPE")
                       },
                       obj = obj.merged)

# Preprocess MERFISH data
obj.l <- lapply(
  X = counts.mat.l,
  FUN = function(counts.mat,
                 cell.metadata) {
    # Define cells to keep
    cells.to.keep <- colnames(x = counts.mat)
    
    # Define sample name
    sample.name <- unique(
      x = gsub(
        pattern = "_.*",
        replacement = "",
        x = cells.to.keep
      )
    )
    
    # Preprocess MERFISH data
    obj <- PreprocessMERFISHSeuratV5(
      counts.mat = counts.mat[,cells.to.keep],
      cell.metadata = cell.metadata[cells.to.keep,],
      sample.name = sample.name,
      assay.name = "MERSCOPE"
    )
    
    # Return object
    return(obj)
  },
  cell.metadata = dplyr::bind_rows(
    lapply(
      X = merscope.data,
      FUN = "[[",
      "cell.metadata"),
    .id = NULL
  )
)

# Merge Seurat objects
obj.merged <- merge(
  x = obj.l[[1]],
  y = obj.l[-1],
  add.cell.ids = NULL,
  collapse = FALSE,
  merge.data = TRUE,
  merge.dr = TRUE
)

# Set all genes as variable features
obj.merged@assays$SCT@var.features <- rownames(x = obj.merged@assays$SCT@scale.data)

# Run Seurat version 5 pipeline with Harmony integration
obj.merged <- RunSeuratV5(
  obj.merged = obj.merged,
  sample.col = "mouse",
  project = "neurons_filtered",
  resolution = seq(from = 0.1,
                   to = 1,
                   by = 0.1),
  azimuth = FALSE,
  spca.features = NULL,
  file.prefix = "MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region",
  path = "results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region"
)

# Join layers of the MERSCOPE assay
obj.merged[["MERSCOPE"]] <- JoinLayers(object = obj.merged[["MERSCOPE"]])

# Run sctransform
obj.merged <- SCTransform(
  object = obj.merged,
  assay = "MERSCOPE",
  new.assay.name = "SCT_full",
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = ncol(x = obj.merged), # use all cells for parameter estimation
  residual.features = NULL, # don't input residual genes as clipping is an issue here
  variable.features.n = nrow(x = obj.merged), # set variable genes to all genes (until issue with clipping is resolved)
  variable.features.rv.th = 1.3,
  vars.to.regress = "volume",
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(-10, 10), # tailored for MERFISH data
  vst.flavor = "v2",
  conserve.memory = FALSE,
  return.only.var.genes = FALSE, # return all genes
  seed.use = 47,
  verbose = TRUE,
  method = "glmGamPoi_offset",
  exclude_poisson = TRUE,
  # latent_var = c("log_umi", "umi_per_gene"),
  n_genes = nrow(x = obj.merged), # use all genes for parameter estimation
  min_cells = 0
)

# Run azimuth
obj.merged <- Azimuth::RunAzimuth(
  query = obj.merged,
  reference = "results/smim32_GTFv102/neurons_filtered/azimuth_Claustrum300_reference/",
  query.modality = "SCT_full",
  annotation.levels = NULL,
  umap.name = "umap_harmony_projected",
  do.adt = FALSE,
  verbose = TRUE,
  assay = "SCT_full",
  k.weight = 50,
  n.trees = 20,
  mapping.score.k = 100)

# Save Seurat object
saveRDS(object = obj.merged,
        file = "results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered/MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region_seurat_v5.1.0_sctransform_azimuth_v0.5.0_neurons_filtered.rds")
