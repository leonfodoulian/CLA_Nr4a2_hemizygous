# Generic function to run Seurat v5
RunSeuratV5 <- function(
    counts.mat.l,
    obj.merged,
    sample.data,
    sample.col,
    project,
    resolution,
    cluster.attributes,
    azimuth = FALSE,
    spca.features = NULL,
    file.prefix,
    path
){
  if (!missing(x = counts.mat.l)) {
    # Identify genes expressed in at least 3 cells in any dataset
    genes.to.keep <- Reduce(
      f = union,
      x = lapply(
        X = counts.mat.l,
        FUN = function(counts.mat) {
          # Identify genes expressed in at least 3 cells
          genes.to.keep <- rowSums(x = counts.mat > 0) >= 3
          
          # Return list of genes
          return(rownames(x = counts.mat)[genes.to.keep])
        }
      )
    )
    
    # Filter count matrices
    counts.mat.l <- lapply(
      X = counts.mat.l,
      FUN = function(counts.mat,
                     genes.to.keep) {
        # Remove all genes not expressed in at least 3 cells in any dataset
        counts.mat <- counts.mat[genes.to.keep,]
        
        # Return data
        return(counts.mat)
      },
      genes.to.keep = genes.to.keep
    )
    
    # Get maximum number of genes
    ngenes <- Reduce(
      f = max,
      x = lapply(
        X = counts.mat.l,
        FUN = nrow
      )
    )
    
    # Get maximum number of cells
    ncells <- Reduce(
      f = max,
      x = lapply(
        X = counts.mat.l,
        FUN = ncol
      )
    )
    
    # Identify variable genes using M3Drop's Depth-Adjusted Negative Binomial model
    future::plan(strategy = future::multisession)
    var.genes.l <- future.apply::future_lapply(
      X = counts.mat.l,
      FUN = function(counts.mat) {
        # Remove all genes not expressed in at least 3 cells
        genes.to.keep <- rowSums(x = counts.mat > 0) >= 3
        counts.mat <- counts.mat[genes.to.keep,]
        
        # Fit DANB to the raw counts matrix
        danb.fit <- M3Drop::NBumiFitModel(counts = counts.mat)
        
        # Check fit quality of DANB
        fit.plot <- R.devices::capturePlot(
          expr = {
            old_par <- par(mfrow = c(1,2))
            stats.danb <- M3Drop::NBumiCheckFitFS(
              counts = counts.mat,
              fit = danb.fit,
              suppress.plot = FALSE
            )
            par(old_par)
          }
        )
        
        # Perform dropout-based DANB feature selection
        feature.selection.plot <- R.devices::capturePlot(
          expr = {
            nb.drop.fs <- M3Drop::NBumiFeatureSelectionCombinedDrop(
              fit = danb.fit,
              ntop = NULL,
              method = "fdr",
              qval.thresh = 0.05,
              suppress.plot = FALSE
            )
          }
        )
        
        # Print "Done" to keep track of progress
        message("Done")
        
        # Return list of data
        return(list(danb.fit = danb.fit,
                    stats.danb = stats.danb,
                    nb.drop.fs = nb.drop.fs,
                    plots = list(fit.plot = fit.plot,
                                 feature.selection.plot = feature.selection.plot)))
      },
      future.seed = TRUE,
      future.scheduling = FALSE
    )
    future::plan(strategy = future::sequential)
    
    # Save variable genes in RDS format
    saveRDS(object = var.genes.l,
            file = file.path(path,
                             project,
                             paste0(file.prefix,
                                    "_M3Drop_v3.10.6_danb_variable_features_",
                                    project,
                                    ".rds")))
    
    # Get union of variable genes in all datasets
    var.genes <- Reduce(
      f = union,
      x = lapply(
        X = var.genes.l,
        FUN = function(var.genes) {
          var.genes <- as.character(x = var.genes$nb.drop.fs[, "Gene"])
          return(var.genes)
        }
      )
    )
    
    # Create Seurat objects and run sctransform
    obj.l <- mapply(
      counts.mat = counts.mat.l,
      sample.name = names(x = counts.mat.l),
      FUN = function(counts.mat,
                     sample.name,
                     ncells,
                     ngenes,
                     sample.data) {
        message("Processing sample ", sample.name)
        
        # Create a Seurat object
        obj <- CreateSeuratObject(
          counts = counts.mat,
          assay = "RNA",
          names.field = 1,
          names.delim = "_",
          meta.data = NULL,
          project = sample.name,
          min.cells = 0,
          min.features = 0
        )
        
        # Calculate the number of UMI per detected gene
        obj@meta.data$nUMI_per_gene <- obj@meta.data$nCount_RNA / obj@meta.data$nFeature_RNA
        
        # Calculate the percentage of mitochondrial gene counts
        obj <- PercentageFeatureSet(
          object = obj,
          pattern = "^mt-",
          features = NULL,
          col.name = "percent_mt",
          assay = NULL
        )
        
        # Add sample data to meta.data
        sample.data <- sample.data[sample.name,,drop = FALSE]
        sample.data <- as.data.frame(x = sample.data[obj@meta.data$orig.ident,])
        rownames(x = sample.data) <- NULL
        for (col.name in colnames(x = sample.data)) {
          obj <- AddMetaData(
            object = obj,
            metadata = sample.data[,col.name],
            col.name = col.name
          )
        }
        
        # Normalize data by library size
        obj <- NormalizeData(
          object = obj,
          assay = "RNA",
          normalization.method = "LogNormalize",
          scale.factor = 10000,
          margin = 1,
          verbose = TRUE
        )
        
        # Run sctransform
        obj <- SCTransform(
          object = obj,
          assay = "RNA",
          new.assay.name = "SCT",
          reference.SCT.model = NULL,
          do.correct.umi = TRUE,
          ncells = ncells, # use all cells for parameter estimation
          residual.features = NULL, # don't input residual genes as clipping is an issue here
          variable.features.n = ngenes, # set variable genes to all genes (until issue with clipping is resolved)
          variable.features.rv.th = 1.3,
          vars.to.regress = c("nUMI_per_gene", "percent_mt"),
          do.scale = FALSE,
          do.center = TRUE,
          clip.range = c(-sqrt(x = ncells/30),
                         sqrt(x = ncells/30)), # same clipping range for all datasets
          vst.flavor = "v2",
          conserve.memory = FALSE,
          return.only.var.genes = FALSE, # return all genes
          seed.use = 47,
          verbose = TRUE,
          method = "glmGamPoi_offset",
          exclude_poisson = TRUE,
          # latent_var = c("log_umi", "umi_per_gene"),
          n_genes = ngenes, # use all genes for parameter estimation
          min_cells = 3
        )
        
        # Rename SCTModel.list with sample name
        names(x = obj@assays$SCT@SCTModel.list) <- sample.name
        
        # Print "Done" to keep track of progress
        message("Done\n")
        
        # Return object
        return(obj)
      },
      MoreArgs = list(ncells = ncells,
                      ngenes = ngenes,
                      sample.data = sample.data),
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE
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
    
    # Subset sctransform and set variable genes to those selected with M3Drop's DANB
    obj.merged@assays$SCT@scale.data <- obj.merged@assays$SCT@scale.data[rownames(x = obj.merged@assays$SCT@scale.data) %in% var.genes,]
    obj.merged@assays$SCT@var.features <- rownames(x = obj.merged@assays$SCT@scale.data)
  } else if (missing(x = counts.mat.l) && missing(x = obj.merged)) {
    stop("Either a list of count matrices or a merged Seurat object should be provided to the function")
  }
  
  # Run Principal Component Analysis
  obj.merged <- RunPCA(
    object = obj.merged,
    assay = "SCT",
    features = NULL,
    npcs = 100,
    rev.pca = FALSE,
    weight.by.var = TRUE,
    verbose = TRUE,
    ndims.print = 1:5,
    nfeatures.print = 30,
    reduction.name = "pca",
    reduction.key = "PC_",
    seed.use = 47
  )
  
  # Detect knee in PC sd distribution using kn.KneeLocator()
  pca.kneed.res <- RKneeLocator(
    x = seq_len(length.out = length(x = obj.merged@reductions$pca@stdev)),
    y = obj.merged@reductions$pca@stdev,
    S = 1,
    curve = "convex",
    direction = "decreasing",
    online = TRUE
  )
  
  # Add kneed results to Seurat object
  obj.merged@misc$pca.kneed.res <- pca.kneed.res
  
  # Define a vector of embedding dimensions for downstream analysis
  embedding.dims <- 1:obj.merged@misc$pca.kneed.res$knee
  
  # Subset PCA embeddings to dimensions of interest
  obj.merged@reductions$pca@cell.embeddings <- obj.merged@reductions$pca@cell.embeddings[,embedding.dims]
  obj.merged@reductions$pca@feature.loadings <- obj.merged@reductions$pca@feature.loadings[,embedding.dims]
  
  # Integrate data using Harmony
  obj.merged@misc$harmony.convergence.plot <- R.devices::capturePlot(
    expr = {
      set.seed(seed = 47) # for reproducibility
      obj.merged <- harmony::RunHarmony(
        object = obj.merged,
        group.by.vars = sample.col,
        reduction.use = "pca",
        dims.use = embedding.dims,
        reduction.save = "harmony",
        project.dim = TRUE,
        theta = NULL,
        sigma = 0.1,
        lambda = NULL,
        nclust = NULL,
        max_iter = 50,
        early_stop = TRUE,
        ncores = 1,
        plot_convergence = TRUE,
        verbose = TRUE,
        .options = harmony::harmony_options(
          alpha = 0.2,
          tau = 0,
          block.size = 0.05,
          max.iter.cluster = 200,
          epsilon.cluster = 1e-05,
          epsilon.harmony = 1e-05
        )
      )
    }
  )
  
  # Construct a Shared Nearest Neighbor (SNN) Graph
  obj.merged <- FindNeighbors(
    object = obj.merged,
    reduction = "harmony",
    dims = embedding.dims,
    assay = NULL,
    features = NULL,
    k.param = 15,
    return.neighbor = FALSE,
    compute.SNN = TRUE,
    prune.SNN = 1/15,
    nn.method = "annoy",
    n.trees = 50,
    annoy.metric = "euclidean",
    nn.eps = 0,
    verbose = TRUE,
    do.plot = FALSE,
    graph.name = NULL,
    l2.norm = FALSE,
    cache.index = FALSE
  )
  
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
  obj.merged <- FindClusters(
    object = obj.merged,
    graph.name = NULL,
    cluster.name = NULL,
    modularity.fxn = NULL,
    initial.membership = NULL,
    node.sizes = NULL,
    resolution = resolution,
    method = "igraph",
    algorithm = 4, # Leiden algorithm
    n.start = NULL,
    n.iter = -1L, # Leidenalg parameter: runs until optimal clustering is reached
    random.seed = 0,
    group.singletons = TRUE,
    temp.file.location = NULL,
    edge.file.name = NULL,
    verbose = TRUE
  )
  
  # Run UMAP
  obj.merged <- RunUMAP(
    object = obj.merged,
    dims = embedding.dims,
    reduction = "harmony",
    features = NULL,
    graph = NULL,
    assay = NULL,
    nn.name = NULL,
    slot = NULL,
    umap.method = "uwot",
    reduction.model = NULL,
    return.model = TRUE,
    n.neighbors = 15L,
    n.components = 2L,
    metric = "cosine",
    n.epochs = NULL,
    learning.rate = 1,
    min.dist = 0.1,
    spread = 1,
    set.op.mix.ratio = 1,
    local.connectivity = 1L,
    repulsion.strength = 1,
    negative.sample.rate = 5L,
    a = NULL,
    b = NULL,
    uwot.sgd = FALSE,
    seed.use = 47L,
    metric.kwds = NULL,
    angular.rp.forest = FALSE,
    densmap = FALSE,
    dens.lambda = 2,
    dens.frac = 0.3,
    dens.var.shift = 0.1,
    verbose = TRUE,
    reduction.name = "umap_harmony",
    reduction.key = "UMAP_"
  )
  
  # Calculate VI from clustering pairs
  clustering.resolutions <- sort(x = grep(pattern = "SCT_snn_res.",
                                          x = colnames(x = obj.merged@meta.data),
                                          value = TRUE))
  pairs.data <- data.frame(comm1 = clustering.resolutions[-length(x = clustering.resolutions)],
                           comm2 = clustering.resolutions[-1],
                           stringsAsFactors = FALSE)
  obj.merged@misc$clustering.pair.vi <- dplyr::bind_rows(
    lapply(
      X = split(x = pairs.data,
                f = 1:nrow(x = pairs.data)),
      FUN = ClusteringPairVI,
      meta.data = obj.merged@meta.data
    )
  )
  
  # Add cluster attributes to the meta data
  if (!missing(x = cluster.attributes)) {
    clustering.result <- cluster.attributes$clustering.result
    obj.merged@meta.data$seurat_clusters <- obj.merged@meta.data[[clustering.result]]
    if ("broad.clusters" %in% names(x = cluster.attributes)) {
      cluster.names <- cluster.attributes$broad.clusters$cluster.names
      cluster.abbreviations <- cluster.attributes$broad.clusters$cluster.abbreviations
      obj.merged@meta.data$broad_cell_type <- cluster.names[obj.merged@meta.data[[clustering.result]]]
      obj.merged@meta.data$broad_cell_type_abbr <- cluster.abbreviations[obj.merged@meta.data[[clustering.result]]]
    }
    if ("subclusters" %in% names(x = cluster.attributes)) {
      cluster.names <- cluster.attributes$subclusters$cluster.names
      cluster.abbreviations <- cluster.attributes$subclusters$cluster.abbreviations
      obj.merged@meta.data$cell_type <- cluster.names[obj.merged@meta.data[[clustering.result]]]
      obj.merged@meta.data$cell_type_abbr <- cluster.abbreviations[obj.merged@meta.data[[clustering.result]]]
    }
  }
  
  if (azimuth) {
    # Join layers of the RNA assay
    obj.merged[["RNA_full"]] <- JoinLayers(object = obj.merged[["RNA"]])
    
    # Run sctransform
    obj.merged <- SCTransform(
      object = obj.merged,
      assay = "RNA_full",
      new.assay.name = "SCT_full",
      reference.SCT.model = NULL,
      do.correct.umi = TRUE,
      ncells = ncol(x = obj.merged), # use all cells for parameter estimation
      residual.features = NULL, # don't input residual genes as clipping is an issue here
      variable.features.n = nrow(x = obj.merged), # set variable genes to all genes (until issue with clipping is resolved)
      variable.features.rv.th = 1.3,
      vars.to.regress = c("nUMI_per_gene", "percent_mt"),
      do.scale = FALSE,
      do.center = TRUE,
      clip.range = c(-sqrt(x = ncol(x = obj.merged)/30),
                     sqrt(x = ncol(x = obj.merged)/30)),
      vst.flavor = "v2",
      conserve.memory = FALSE,
      return.only.var.genes = FALSE, # return all genes
      seed.use = 47,
      verbose = TRUE,
      method = "glmGamPoi_offset",
      exclude_poisson = TRUE,
      # latent_var = c("log_umi", "umi_per_gene"),
      n_genes = nrow(x = obj.merged), # use all genes for parameter estimation
      min_cells = 3
    )
    
    # Run supervised PCA
    obj.merged <- RunSPCA(
      object = obj.merged,
      assay = "SCT_full",
      features = intersect(x = spca.features,
                           y = rownames(x = obj.merged@assays$SCT_full@scale.data)),
      npcs = 100,
      reduction.name = "spca",
      reduction.key = "SPC_",
      graph = "SCT_snn",
      verbose = TRUE,
      seed.use = 47
    )
    
    # Detect knee in sPC sd distribution using kn.KneeLocator()
    spca.kneed.res <- RKneeLocator(
      x = seq_len(length.out = length(x = obj.merged@reductions$spca@stdev)),
      y = obj.merged@reductions$spca@stdev,
      S = 1,
      curve = "convex",
      direction = "decreasing",
      online = TRUE
    )
    
    # Add kneed results to Seurat object
    obj.merged@misc$spca.kneed.res <- spca.kneed.res
    
    # Subset sPCA embeddings to dimensions of interest
    obj.merged@reductions$spca@cell.embeddings <- obj.merged@reductions$spca@cell.embeddings[,1:obj.merged@misc$spca.kneed.res$knee]
    obj.merged@reductions$spca@feature.loadings <- obj.merged@reductions$spca@feature.loadings[,1:obj.merged@misc$spca.kneed.res$knee]
    
    # Create a Seurat object compatible with Azimuth
    # obj.ref <- Azimuth::AzimuthReference(
    #   object = obj.merged,
    #   refUMAP = "umap_harmony",
    #   refDR = "spca",
    #   refAssay = "SCT_full",
    #   dims = 1:obj.merged@misc$spca.kneed.res$knee,
    #   k.param = 31,
    #   plotref = "umap_harmony",
    #   plot.metadata = NULL,
    #   ori.index = NULL,
    #   colormap = NULL,
    #   assays = "SCT_full",
    #   metadata = intersect(x = c("broad_cell_type", "broad_cell_type_abbr", "cell_type", "cell_type_abbr"),
    #                        y = colnames(x = obj.merged@meta.data)),
    #   reference.version = "0.0.0",
    #   verbose = TRUE
    # )
    
    # Save Azimuth references and neighbors index
    # Azimuth::SaveAzimuthReference(object = obj.ref,
    #                               folder = file.path(path, project, "azimuth_reference/"))
  }
  
  # Save Seurat object
  saveRDS(object = obj.merged,
          file = file.path(path,
                           project,
                           paste0(file.prefix,
                                  {
                                    if (!missing(x = counts.mat.l)) {
                                      "_M3Drop_v3.10.6_danb"
                                    }
                                  },
                                  "_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-",
                                  obj.merged@misc$pca.kneed.res$knee,
                                  "_",
                                  project,
                                  ".rds")))
  
  # Print "Seurat object saved"
  message("Seurat object saved")
  
  # Return object
  return(obj.merged)
}
