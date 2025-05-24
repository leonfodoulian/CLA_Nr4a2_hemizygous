PreprocessMERFISHSeuratV5 <- function(
    counts.mat,
    cell.metadata,
    sample.name,
    assay.name = "MERFISH"
){
  # Create a Seurat object
  obj <- CreateSeuratObject(
    counts = counts.mat,
    assay = assay.name,
    names.field = 1,
    names.delim = "_",
    meta.data = cell.metadata,
    project = sample.name,
    min.cells = 0,
    min.features = 0
  )
  
  # Run sctransform
  obj <- SCTransform(
    object = obj,
    assay = assay.name,
    new.assay.name = "SCT",
    reference.SCT.model = NULL,
    do.correct.umi = TRUE,
    ncells = ncol(x = obj), # use all cells for parameter estimation
    residual.features = NULL, # don't input residual genes as clipping is an issue here
    variable.features.n = nrow(x = obj), # set variable genes to all genes (until issue with clipping is resolved)
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
    n_genes = nrow(x = obj), # use all genes for parameter estimation
    min_cells = 0
  )
  
  # Rename SCTModel.list with sample name
  names(x = obj@assays$SCT@SCTModel.list) <- sample.name
  
  # Return object
  return(obj)
}
