# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Required packages
require(Seurat)
require(reticulate)

# Python packages
bi <- reticulate::import_builtins(convert = TRUE)
dask <- reticulate::import(module = "dask", convert = TRUE)
pandas <- reticulate::import(module = "pandas", convert = TRUE)
arboreto <- reticulate::import(module = "arboreto", convert = FALSE)
ctxcore.rnkdb <- reticulate::import(module = "ctxcore.rnkdb", convert = FALSE)
pyscenic <- reticulate::import(module = "pyscenic", convert = FALSE)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/pyscenic")

# List files to load for pySCENIC
tf.names.file <- "/Users/leonfodoulian/scData/scenic/cistarget/tf_lists/allTFs_mm.txt"
tf.ranking.files <- list.files(path = "/Users/leonfodoulian/scData/scenic/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/",
                               pattern = "genes_vs_motifs.rankings.feather$",
                               full.names = TRUE)
names(x = tf.ranking.files) <- gsub(pattern = ".*/|\\.feather",
                                    replacement = "",
                                    x = tf.ranking.files)
motif.annotations.file <- "/Users/leonfodoulian/scData/scenic/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load files for pySCENIC
tf.names <- arboreto$utils$load_tf_names(path = tf.names.file)
tf.ranking.l <- mapply(fname = tf.ranking.files,
                       name = names(x = tf.ranking.files),
                       FUN = function(fname,
                                      name) {
                         return(ctxcore.rnkdb$FeatherRankingDatabase(fname = fname,
                                                                     name = name))
                       },
                       SIMPLIFY = FALSE,
                       USE.NAMES = TRUE)

# Get names of genes for which we have TF motif rankings
genes.with.motifs <- Reduce(
  f = union,
  x = lapply(
    X = tf.ranking.l,
    FUN = function(tf.ranking) {
      return(bi$list(tf.ranking$genes))
    }
  )
)

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Subset Seurat object for CLA, shell and L6a cell types
obj <- subset(x = obj,
              subset = broad_cell_type_abbr %in% c("CLA", "shell", "L6a"))

# Load list of modulated genes in CLA clusters
hemi.genes <- Reduce(
  f = union,
  x = lapply(
    X = readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                 paste0(file.prefix, "_hemi_modulated_genes.rds")))[c("CLA 1", "CLA 2")],
    FUN = function(modulated.genes) {
      return(rownames(x = modulated.genes$observed))
    }
  )
)

# Prepare counts data for pySCENIC
counts.data <- obj@assays$SCT_full@data
genes.to.keep <- intersect(
  # x = rownames(x = counts.data)[rowSums(x = counts.data > 0) >= 3], # for all genes
  x = hemi.genes, # for modulated genes only
  y = genes.with.motifs
)
# genes.to.keep <- union(
#   x = genes.to.keep,
#   y = bi$list(tf.names) # union of all TF genes and modulated genes with motifs
# ) 
counts.data <- counts.data[rownames(x = counts.data) %in% genes.to.keep,] # subset for genes to keep
counts.data.df <- as.data.frame(x = t(x = counts.data)) # transform counts to data frame

# Subset TF genes for genes in counts data
tf.names <- bi$list(tf.names)
tf.names <- tf.names[tf.names %in% rownames(x = counts.data)]

# Remove Sema4a from list of TF genes
tf.names <- tf.names[tf.names != "Sema4a"]

# Run GENIE3 for inference of co-expression modules (gene regulatory network)
client <- dask$distributed$Client(address = dask$distributed$LocalCluster(n_workers = 1L,
                                                                          threads_per_worker = 10L))
adjacencies <- arboreto$algo$genie3(
  expression_data = t(x = counts.data),
  gene_names = rownames(x = counts.data),
  tf_names = tf.names,
  client_or_address = client,
  limit = NULL,
  seed = 47L, # for reproducibility
  verbose = TRUE
)
client$cluster$close()

# Save adjacencies data
adjacencies$to_csv(file.path("results/smim32_GTFv102/neurons_filtered/pyscenic/",
                             paste0(file.prefix, "_pyscenic_genie3_adjacencies.csv")))

# Create modules from a dataframe containing weighted adjacencies between a TF and its target genes
modules <- pyscenic$utils$modules_from_adjacencies(
  adjacencies = adjacencies, # pd.DataFrame
  ex_mtx = counts.data.df, # pd.DataFrame
  thresholds = c(0.90, 0.95),
  top_n_targets = c(25L, 50L),
  top_n_regulators = c(5L, 10L, 50L),
  min_genes = 20L,
  absolute_thresholds = FALSE, # more robust
  rho_dichotomize = TRUE,
  keep_only_activating = FALSE, # get also repressed modules
  rho_threshold = 0.03,
  rho_mask_dropouts = FALSE # use all cells
)

# Calculate a list of enriched motifs and the corresponding target genes for all modules
client <- dask$distributed$Client(address = dask$distributed$LocalCluster(n_workers = 6L,
                                                                          threads_per_worker = 1L))
prune2df <- pyscenic$prune$prune2df(
  rnkdbs = c(tf.ranking.l$mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings,
             tf.ranking.l$mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings),
  modules = modules,
  motif_annotations_fname = motif.annotations.file,
  rank_threshold = 500L,
  auc_threshold = 0.01,
  nes_threshold = 2.75,
  motif_similarity_fdr = 0.001,
  orthologuous_identity_threshold = 0,
  weighted_recovery = FALSE,
  client_or_address = client,
  num_workers = NULL,
  module_chunksize = 100,
  filter_for_annotation = TRUE
)
client$cluster$close()

# Save enriched motifs and corresponding target genes data
prune2df$to_csv(file.path("results/smim32_GTFv102/neurons_filtered/pyscenic/",
                          paste0(file.prefix, "_pyscenic_genie3_prune2df.csv")))

# Create regulons from table of enriched motifs
regulons <- pyscenic$prune$df2regulons(df = prune2df)

# Create a list of regulons compatible with R
regulons.l <- lapply(
  X = bi$list(regulons),
  FUN = function(regulon) {
    return(list(name = regulon$name,
                transcription_factor = regulon$transcription_factor,
                context = bi$list(regulon$context),
                gene_weights = setNames(object = unlist(x = regulon$weights),
                                        nm = unlist(x = regulon$genes)),
                score = regulon$score))
  }
)

# Calculate enrichment of gene signatures for single cells
aucell <- pyscenic$aucell$aucell(
  exp_mtx = counts.data.df,
  signatures = regulons,
  auc_threshold = 0.01,
  noweights = FALSE,
  normalize = FALSE,
  seed = 47L, # for reproducibility
  num_workers = 7L
)

# Save regulons and aucell data
saveRDS(object = list(regulons = regulons.l,
                      aucell = pandas$DataFrame(data = aucell)),
        file = file.path("results/smim32_GTFv102/neurons_filtered/",
                         paste0(file.prefix, "_pyscenic_genie3.rds")))
