# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/MERSCOPE_Rodriguez_Claustrum300_CLA_EP_boundary.R")
source(file = "/Users/leonfodoulian/scData/FindMarkersSCTV5.R")

# Required packages
require(Seurat)

# Define rds file name prefixes
scrnaseq.file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"
merscope.file.prefix <- "MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region_seurat_v5.1.0_sctransform_azimuth_v0.5.0_neurons_filtered"

# Load Seurat objects
scrnaseq.obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                         paste0(scrnaseq.file.prefix, ".rds")))
merscope.obj <- readRDS(file = file.path("results/MERSCOPE_Rodriguez_Claustrum300/claustro_insular_region/neurons_filtered",
                                         paste0(merscope.file.prefix, ".rds")))

# Prepare scRNaseq object to run differential expression on SCT assay with multiple models
scrnaseq.obj <- PrepSCTFindMarkers(object = scrnaseq.obj,
                                   assay = "SCT",
                                   verbose = TRUE)

# Get index of CLA projection neurons from the scRNaseq dataset
cla.cell.id <- which(x = scrnaseq.obj@meta.data$broad_cell_type_abbr == "CLA")

# Prepare MERSCOPE meta data
meta.data <- merscope.obj@meta.data

# Subset MERSCOPE meta data for CLA neurons with a prediction score of 0.7
cells.to.keep <- meta.data$predicted.broad_cell_type_abbr == "CLA" &
  meta.data$predicted.broad_cell_type_abbr.score >= 0.7
meta.data <- meta.data[cells.to.keep,]

# Merge MERSCOPE meta data to CLA/EP boundary data
meta.data <- merge(
  x = meta.data,
  y = cla.ep.boundary,
  by = c("mouse", "level", "cla_side")
)

# Define where CLA neurons are located within the CLA/EP complex
meta.data$region <- ifelse(
  test = meta.data$centroid_y_rot > meta.data$boundary,
  yes = gsub(pattern = ".*\\/",
             replacement = "",
             x = meta.data$boundary_type),
  no = gsub(pattern = "\\/.*",
            replacement = "",
            x = meta.data$boundary_type)
)

# Compute proportion of CLA 1 and CLA 2 cells in each region of the CLA/EP complex
cla.cell.type.prop <- data.table::data.table(meta.data)[
  , list(cla_1 = sum(predicted.cell_type_abbr == "CLA 1") / length(x = predicted.cell_type_abbr),
         cla_2 = sum(predicted.cell_type_abbr == "CLA 2") / length(x = predicted.cell_type_abbr)),
  by = c("experiment", "genotype", "mouse", "level", "cla_side", "region")]

# Prepare nearest neighbor data
nn.data <- data.frame(merscope.obj@neighbors$query_ref.nn@nn.idx)
colnames(x = nn.data) <- paste0("nn", 1:ncol(x = nn.data))

# Bind MERSCOPE meta data and nearest neighbor data
meta.nn.data <- cbind(
  meta.data[,c("cell", "experiment", "genotype", "mouse", "level", "cla_side", "boundary_type", "region")],
  nn.data[meta.data$cell,]
)

# Define reference region for each level
reference.regions <- c(
  "44" = "CLA",
  "53" = "CLA",
  "68" = "CLA-like L6"
)

# Calculate the expected proportion of CLA neurons within the CLA/EP complex
expected.props <- lapply(
  X = split(x = meta.nn.data,
            f = meta.nn.data$level),
  FUN = function(level.data,
                 reference.regions) {
    # Define reference region for proportion calculations
    reference.region <- reference.regions[unique(x = level.data$level)]
    
    # Calculate the expected proportion of CLA neurons within the CLA/EP complex
    expected.prop <- mean(data.table::data.table(level.data)[
      ,list(prop = sum(region == reference.region) / length(x = region)),
      by = c("experiment", "genotype", "mouse", "level", "cla_side")]$prop)
    
    # Return data
    return(expected.prop)
  },
  reference.regions = reference.regions
)

# Melt MERSCOPE meta data
meta.nn.data <- reshape2::melt(
  data = meta.nn.data,
  id.vars = c("cell", "experiment", "genotype", "mouse", "level", "cla_side", "boundary_type", "region"),
  variable.name = "nn_position",
  value.name = "nn_id"
)

# Subset melted MERSCOPE meta data for CLA projection neurons from the scRNAseq dataset
meta.nn.data <- meta.nn.data[meta.nn.data$nn_id %in% cla.cell.id,]

# Add cell names of CLA projection neurons from the scRNAseq dataset
meta.nn.data$nn_name <- rownames(x = scrnaseq.obj@meta.data)[as.numeric(x = meta.nn.data$nn_id)]

# Determine where each CLA projection neuron from the scRNAseq dataset is located in the CLA/EP complex by level
region.attribute.data.l <- lapply(
  X = split(x = meta.nn.data,
            f = meta.nn.data$level),
  FUN = function(level.data,
                 expected.props,
                 reference.regions) {
    # Define reference region
    reference.region <- reference.regions[unique(x = level.data$level)]
    
    # Define query region
    query.region <- setdiff(
      x = unique(x = level.data$region),
      y = reference.region
    )
    
    # Get expected proportion of reference neurons
    expected.prop <- expected.props[[unique(x = level.data$level)]]
    
    # Calculate the observed proportion of CLA neurons among neighbors
    prop.data <- dplyr::bind_rows(
      lapply(
        X = split(x = level.data,
                  f = level.data$nn_name),
        FUN = function(cell.data,
                       reference.region,
                       query.region) {
          return(data.frame(reference = sum(cell.data$region == reference.region),
                            query = sum(cell.data$region == query.region),
                            prop = sum(cell.data$region == reference.region) / length(x = cell.data$region),
                            row.names = NULL))
        },
        reference.region = reference.region,
        query.region = query.region
      ),
      .id = "nn_name"
    )
    
    # Compute p-values in inequalities between observed and expected proportions of CLA neurons
    prop.data$p_value <- unlist(
      x = lapply(
        X = split(x = prop.data,
                  f = 1:nrow(x = prop.data)),
        FUN = function(cell.data,
                       expected.prop) {
          return(prop.test(x = as.matrix(x = cell.data[c("reference","query")]),
                           p = expected.prop)$p.value)
        },
        expected.prop = expected.prop
      )
    )
    
    # Adjust p-values for multiple comparisons
    prop.data$adjusted_p_value <- p.adjust(
      p = prop.data$p_value,
      method = "holm"
    )
    
    # Determine where each CLA projection neuron from the scRNAseq dataset is located in the CLA/EP complex
    prop.data$region <- ifelse(
      # test = (prop.data$prop < expected.prop & prop.data$adjusted_p_value < 0.05) | prop.data$prop == 0,
      test = (prop.data$prop < expected.prop & prop.data$adjusted_p_value < 0.05),
      yes = query.region,
      no = NA
    )
    prop.data$region <- ifelse(
      # test = (prop.data$prop > expected.prop & prop.data$adjusted_p_value < 0.05) | prop.data$prop == 1,
      test = (prop.data$prop > expected.prop & prop.data$adjusted_p_value < 0.05),
      yes = reference.region,
      no = prop.data$region
    )
    
    # Return data
    return(prop.data)
  },
  expected.props = expected.props,
  reference.regions = reference.regions
)

# Bind proportion data by row
region.attribute.data <- dplyr::bind_rows(
  region.attribute.data.l,
  .id = "level"
)

# Disregard cells whose region attributes could not be determined
region.attribute.data <- region.attribute.data[!is.na(x = region.attribute.data$region) & region.attribute.data$region != "CLA-EP",]

# Disregard cells with a mixed region attribute
region.attributes <- unlist(
  x = lapply(
    X = split(x = region.attribute.data,
              f = region.attribute.data$nn_name),
    FUN = function(cell.data) {
      unique.region <- unique(x = cell.data$region)
      if (length(x = unique.region) == 1) {
        return(unique.region)
      } else {
        return(NA)
      }
    }
  )
)
region.attributes <- region.attributes[!is.na(x = region.attributes)]

# Add region attribute of each cell to the scRNAseq Seurat object
scrnaseq.obj@meta.data$region <- unname(obj = region.attributes[rownames(x = scrnaseq.obj@meta.data)])

# Define CLA/EP complex region pairs for markers calculation
region.pairs <- data.frame(
  ident.1 = c("CLA", "CLA", "CLA-like L6"),
  ident.2 = c("EP", "CLA-like L6", "EP")
)
region.pairs$pair <- paste(
  region.pairs$ident.1,
  region.pairs$ident.2,
  sep = " vs ")

# Compute pairwise marker genes for CLA, EP and CLA-like L6 neurons
markers.l <- lapply(
  X = split(x = region.pairs,
            f = region.pairs$pair),
  FUN = function(pair.data,
                 scrnaseq.obj) {
    # Subset Seurat object for region pair keeping only WT cells
    obj <- subset(x = scrnaseq.obj,
                  subset = region %in% c(pair.data$ident.1, pair.data$ident.2) & genotype %in% "Nr4a2(WT/WT)")
    
    # Skip calculation if any region has less than 10 cells
    if (any(table(obj$region) < 10)) {
      # Return NULL
      return(NULL)
    } else {
      # Compute pairwise marker genes
      markers <- FindMarkersSCTV5(
        obj = obj,
        ident.1 = pair.data$ident.1,
        ident.2 = pair.data$ident.2,
        deg.ident = "region",
        deg.cond = list(),
        slot = "data",
        fc.slot = "data",
        logfc.threshold = 0.25,
        test.use = "wilcox",
        # test.use = "negbinom",
        min.pct = 0.1,
        min.diff.pct = -Inf,
        verbose = TRUE,
        only.pos = FALSE,
        max.cells.per.ident = Inf,
        random.seed = 1,
        latent.vars = NULL,
        # latent.vars = c("mouse", "gender", "age"),
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
      
      # Return data
      return(markers)
    }
  },
  scrnaseq.obj = scrnaseq.obj
)

# Save list of data
saveRDS(object = list(meta.data = meta.data,
                      cla.cell.type.prop = cla.cell.type.prop,
                      region.attribute.data = region.attribute.data,
                      region.attributes = region.attributes,
                      markers = markers.l),
        file = file.path("results/smim32_GTFv102/neurons_filtered",
                         paste0(scrnaseq.file.prefix, "_cla_ep_analysis.rds")))
