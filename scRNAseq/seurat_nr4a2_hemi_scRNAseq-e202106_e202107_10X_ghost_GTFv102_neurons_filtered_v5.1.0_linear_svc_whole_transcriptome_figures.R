# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(reticulate)
require(extrafont)
require(ggplot2)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp")

# Python packages
svm <- reticulate::import(module = "sklearn.svm")
ms <- reticulate::import(module = "sklearn.model_selection")

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Subset Seurat object for CLA, shell and L6a cell types
obj <- subset(x = obj,
              subset = broad_cell_type_abbr %in% c("CLA", "shell", "L6a"))

# Load predicted labels data using linear SVC
svm.data.l <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                       paste0(file.prefix, "_linear_svc_predicted_labels.rds")))

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

# Get names of wt cells
wt.cells <- rownames(x = obj@meta.data[obj@meta.data$genotype == "Nr4a2(WT/WT)",])

# Prepare counts and labels data for SVM
counts.data <- t(x = obj@assays$SCT_full@scale.data)
broad.cell.type <- setNames(object = obj@meta.data[,"broad_cell_type_abbr"],
                            nm = rownames(x = obj@meta.data))

# Identify genes expressed in at least 3 cells
genes.to.keep <- intersect(
  x = rownames(x = obj@assays$SCT@counts)[rowSums(x = obj@assays$SCT@counts > 0) >= 3],
  y = rownames(x = obj@assays$SCT_full@scale.data)
)

# Set linear SVC model
svm.model <- svm$LinearSVC(
  penalty = svm.data.l$best.parameters$penalty, # "l2"
  loss = "squared_hinge",
  dual = svm.data.l$best.parameters$dual, # "auto"
  tol = 1e-5,
  C = svm.data.l$best.parameters$C, # 9e-04
  multi_class = "ovr",
  fit_intercept = TRUE,
  intercept_scaling = 1,
  class_weight = NULL,
  verbose = 0L,
  random_state = NULL,
  max_iter = 50000L
)

# Fit SVC model with all wt cells
svm.fit <- svm.model$fit(X = counts.data[wt.cells, genes.to.keep],
                         y = broad.cell.type[wt.cells])

# Compute gene weights
gene.weights <- setNames(object = colSums(x = abs(x = svm.fit$coef_)),
                         nm = genes.to.keep)

# Get gene weights threshold
gene.weights.thresh <- median(x = gene.weights) + 3 * mad(x = gene.weights)
percentile <- 1 - sum(gene.weights >= gene.weights.thresh) / length(x = gene.weights)

# Prepare data frame for statistics and plotting
gene.weights <- data.frame(gene_weight = gene.weights)
gene.weights$gene <- rownames(x = gene.weights)
gene.weights$is_modulated <- gene.weights$gene %in% hemi.genes
gene.weights$is_high_weighted <- gene.weights$gene_weight >= gene.weights.thresh

# Check whether modulated genes are enriched among high-weighted genes
prop.res <- prop.test(
  x = table(gene.weights[gene.weights$is_high_weighted,]$is_modulated)[c("TRUE", "FALSE")],
  p = sum(gene.weights$is_modulated) / nrow(x = gene.weights)
)

# Compute statistics on the data using Wilcoxon's rank sum test
wilcox.res <- wilcox.test(
  formula = gene_weight ~ is_modulated,
  data = gene.weights
)

# Plot gene weights across modulated and non-modulated genes
gene.weights.plot <- ggplot(data = gene.weights,
                          mapping = aes(x = is_modulated,
                                        y = gene_weight,
                                        colour = is_modulated)) +
  geom_boxplot(outlier.shape = 19,
               outlier.size = 2,
               outlier.stroke = 0,
               show.legend = FALSE) +
  scale_x_discrete(limits = c("TRUE", "FALSE"),
                   breaks = c("TRUE", "FALSE"),
                   labels = c("TRUE" = "modulated",
                              "FALSE" = "unmodulated")) +
  scale_y_continuous(limits = c(0, 0.01 * ceiling(x = max(gene.weights$gene_weight) / 0.01)),
                     breaks = seq(from = 0,
                                  to = 0.01 * ceiling(x = max(gene.weights$gene_weight) / 0.01),
                                  by = 0.01),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("TRUE" = "#333333",
                                 "FALSE" = "#999999"),
                      breaks = c("TRUE", "FALSE")) +
  coord_cartesian(clip = "off") +
  labs(x = "genes",
       y = "SVC gene weights") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,
                                 angle = 45,
                                 hjust = 1,
                                 vjust = 1,
                                 family = "Arial",
                                 face = "plain",
                                 colour = "black"),
        axis.text.y = element_text(size = 6,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black"),
        axis.title = element_text(size = 7,
                                  family = "Arial",
                                  face = "plain",
                                  colour = "black"),
        legend.text = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(size = rel(x = 0.5),
                                 colour = "black",
                                 lineend = "square"),
        axis.ticks = element_line(size = rel(x = 0.5),
                                  colour = "black",
                                  lineend = "square"),
        axis.ticks.length = ggplot2::unit(x = 0.5,
                                          units = "mm"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.margin = element_blank(),
        legend.key.size = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank())

# Save plot in pdf format
FixSizeAndSave(plot = gene.weights.plot,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_linear_svc_whole_transcriptome_gene_weights.pdf")),
               is.ggassemble = FALSE,
               panel.width = 1,
               panel.height = 4.2,
               unit.use = "cm",
               margin = 1,
               use.ggsave = TRUE,
               useDingbats = FALSE)
