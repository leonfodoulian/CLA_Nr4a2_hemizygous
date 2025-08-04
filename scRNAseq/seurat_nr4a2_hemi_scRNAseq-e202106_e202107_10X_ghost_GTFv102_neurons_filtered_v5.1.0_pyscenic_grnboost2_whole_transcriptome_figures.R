# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/SCENICRegulonDataFrame.R")
source(file = "/Users/leonfodoulian/scData/ExpandClusterAttributes.R")
source(file = "/Users/leonfodoulian/scData/BlendColorAlpha.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(extrafont)
require(ggplot2)
require(ggtext)

# Python packages
pandas <- reticulate::import(module = "pandas", convert = TRUE)
series <- reticulate::import(module = "pandas", convert = FALSE)$Series
pyscenic <- reticulate::import(module = "pyscenic", convert = FALSE)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp")

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Expand cluster attributes with genotype labels
hemi.cluster.attributes <- lapply(
  X = cluster.attributes[1:2],
  FUN = ExpandClusterAttributes,
  labels = c("wt", "del"),
  alpha = c(1, 0.4),
  prefix = TRUE
)

# Define abbreviations of genotypes
genotype.abbreviations <- setNames(object = c("wt", "del"),
                                   nm = c("Nr4a2(WT/WT)", "Nr4a2(SA-IRES-Dre/WT)"))

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Subset Seurat object for CLA, shell and L6a cell types
obj <- subset(x = obj,
              subset = broad_cell_type_abbr %in% c("CLA", "shell", "L6a"))

# Add abbreviated genotypes to Seurat object meta data
obj@meta.data$genotype_abbr <- genotype.abbreviations[obj@meta.data$genotype]

# Create a column combining genotype and cell type identities of cells 
obj@meta.data$cell_type_genotype_abbr <- paste(obj@meta.data$cell_type_abbr,
                                               obj@meta.data$genotype_abbr)

# Load list of modulated genes
hemi.genes <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                       paste0(file.prefix, "_hemi_modulated_genes.rds")))[["CLA 1"]]$observed
hemi.genes$gene <- rownames(x = hemi.genes)

# Load pySCENIC results
pyscenic.res <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                         paste0(file.prefix, "_pyscenic_grnboost2_whole_transcriptome.rds")))

# Get cell type identities of cells 
cell.type <- setNames(object = obj@meta.data$broad_cell_type_abbr,
                      nm = rownames(x = obj@meta.data))

# Reorder aucell rows with names of cells in Seurat object
aucell <- pyscenic.res$aucell[names(x = cell.type),]

# Calculate regulon specificity score for each genotype
rss.l <- lapply(
  X = split(x = obj@meta.data[, "genotype", drop = FALSE],
            f = obj@meta.data$genotype),
  FUN = function(genotype.data,
                 aucell,
                 cell.type) {
    # Calculate RSS
    rss <- pyscenic$rss$regulon_specificity_scores(
      auc_mtx = aucell[rownames(x = genotype.data),],
      cell_type_series = series(data = cell.type[rownames(x = genotype.data)])
    )
    
    # Return data
    return(pandas$DataFrame(data = rss))
  },
  aucell = aucell,
  cell.type = cell.type
)

# Calculate RSS contrast (wt / del)
wt.rss <- rss.l[["Nr4a2(WT/WT)"]]
hemi.rss <- rss.l[["Nr4a2(SA-IRES-Dre/WT)"]]
hemi.rss <- hemi.rss[rownames(x = wt.rss), colnames(x = wt.rss)]
# rss.contrast <- (wt.rss - hemi.rss) / (wt.rss + hemi.rss)
rss.contrast <- (wt.rss - hemi.rss)

# Melt data for plotting
rss.contrast.m <- setNames(
  object = reshape2::melt(data = as.matrix(x = rss.contrast)),
  nm = c("cell_type", "regulon", "rss_contrast")
) # rss contrast
rss.m <- setNames(
  object = rbind(
    cbind(reshape2::melt(data = as.matrix(x = wt.rss)),
          "genotype" = "wt"),
    cbind(reshape2::melt(data = as.matrix(x = hemi.rss)),
          "genotype" = "del")
  ),
  nm = c("cell_type", "regulon", "rss", "genotype")
) # rss

# Remove regulons with NA rss or rss contrast
rss.contrast.m <- rss.contrast.m[!is.na(x = rss.contrast.m$rss_contrast),] # rss contrast
rss.m <- rss.m[!is.na(x = rss.m$rss),] # rss

# Get ordered regulons for CLA cell type
regulon.levels <- as.character(
  x = subset(
    x = rss.contrast.m[order(rss.contrast.m$rss_contrast, decreasing = TRUE),],
    subset = cell_type == "CLA"
  )$regulon
)

# Convert regulon to factor using order in CLA cell type
rss.contrast.m$regulon <- factor(x = rss.contrast.m$regulon,
                                 levels = regulon.levels) # rss contrast
rss.m$regulon <- factor(x = rss.m$regulon,
                        levels = regulon.levels) # rss

# Convert cell type to factor
rss.contrast.m$cell_type <- factor(x = rss.contrast.m$cell_type,
                                   levels = c("CLA", "shell", "L6a")) # rss contrast
rss.m$cell_type <- factor(x = rss.m$cell_type,
                          levels = c("CLA", "shell", "L6a")) # rss

# Order data by cell type
rss.contrast.m <- rss.contrast.m[order(rss.contrast.m$cell_type, decreasing = TRUE),] # rss contrast
rss.m <- rss.m[order(rss.m$cell_type, decreasing = TRUE),] # rss

# Create a column combining genotype and cell type identities
rss.m$cell_type_genotype <- paste(rss.m$cell_type,
                                  rss.m$genotype)

# Get maximum contrast value
max.contrast <- max(abs(x = rss.contrast.m$rss_contrast))
max.contrast.round <- 0.1 * ceiling(x = max.contrast / 0.1)

# Get maximum rss value
max.rss <- 0.1 * ceiling(x = max(rss.m$rss) / 0.1)

# Generate alpha through min-max normalization of RSS
rss.contrast.m$alpha <- (rss.contrast.m$rss_contrast + max.contrast) / (2 * max.contrast)
rss.contrast.m$alpha <- (rss.contrast.m$alpha * (1 - 0.4)) + 0.4 # scale to [0.4, 1]

# Get ordered TFs
tf.levels <- unique(
  x = gsub(
    pattern = "\\(.*",
    replacement = "",
    x = regulon.levels
  )
)

# Prepare regulon data frame
regulons <- dplyr::bind_rows(
  lapply(
    X = pyscenic.res$regulons,
    FUN = SCENICRegulonDataFrame,
    genes = hemi.genes$gene
  ),
  .id = NULL
)

# Bind regulon data to modulated genes data
regulons <- merge(
  x = regulons,
  y = hemi.genes,
  by = "gene",
  all = TRUE
)

# Subset regulons for modulated genes
regulons <- regulons[regulons$gene %in% hemi.genes$gene,]

# Convert TFs to factor
regulons$transcription_factor <- factor(x = regulons$transcription_factor,
                                        levels = tf.levels)

# Add modulation to the data and convert to factor
regulons$modulation <- factor(
  x = ifelse(test = regulons$avg_log2FC < 0,
             yes = "downregulated",
             no = "upregulated"),
  levels = c("downregulated",
             "upregulated")
)

# Order of genes by fold change 
regulons <- regulons[order(regulons$avg_log2FC, decreasing = FALSE),]
regulons$gene <- factor(x = regulons$gene,
                        levels = unique(x = regulons$gene))

# Order regulons by TFs
regulons <- regulons[order(regulons$transcription_factor, regulons$regulon_context, decreasing = FALSE),]
regulons$regulon_name <- factor(x = regulons$regulon_name,
                                levels = rev(x = unique(x = regulons$regulon_name)))

# Define regulon breaks and labels
regulon.labels <- setNames(
  object = ifelse(test = grepl(pattern = "\\+",
                               x = rev(x = levels(x = regulons$regulon_name))),
                  yes = "act.",
                  no = "rep."),
  nm = rev(x = levels(x = regulons$regulon_name))
)

# Define parameters for regulon plots
regulon.params.l <- list(
  "enriched_regulons" = list(
    genes.to.highlight = c(
      "Car3", "Col11a1", "Cux1", "Gfra1", "Gnb4", # downregulated
      "Gng2", "Nr2f2", "Nr4a2", "Satb1", "Synpr", # downregulated
      "Cnr1", "Cpeb1", "Nfib", "Nos1ap", "Ryr2" # upregulated
    ),
    tf.to.keep = tf.levels[1:3], # top 3 CLA regulons enriched in WT cells
    panel.width = 13,
    panel.height = 3,
    file.suffix = "CLA_wt_enriched_regulons.pdf"
  ),
  "nr4a2_regulon" = list(
    genes.to.highlight = as.character(x = regulons$gene[regulons$transcription_factor == "Nr4a2" & !is.na(x = regulons$weight)]),
    tf.to.keep = "Nr4a2",
    panel.width = 13,
    panel.height = 0.5,
    file.suffix = "nr4a2_regulon.pdf"
  )
)

# Plot RSS contrast for each cell type
rss.contrast.plots.l <- lapply(
  X = split(x = rss.contrast.m,
            f = rss.contrast.m$cell_type),
  FUN = function(plot.data,
                 max.contrast,
                 cluster.blended.colors,
                 save.plot,
                 file.prefix) {
    # Set factors to character
    plot.data$cell_type <- as.character(x = plot.data$cell_type)
    plot.data$regulon <- as.character(x = plot.data$regulon)
    
    # Order data by RSS contrast
    plot.data <- plot.data[order(plot.data$rss_contrast, decreasing = TRUE),]
    plot.data$regulon <- factor(x = plot.data$regulon,
                                levels = plot.data$regulon)
    
    # Get name of cell type
    cell.type <- unique(x = plot.data$cell_type)
    
    # Plot RSS contrast
    rss.contrast.plot <- ggplot(data = plot.data,
                                mapping = aes(x = regulon,
                                              y = rss_contrast,
                                              colour = rss_contrast)) +
      geom_point(size = 2,
                 shape = 19,
                 stroke = 0,
                 show.legend = FALSE) +
      scale_x_discrete(breaks = levels(x = plot.data$regulon),
                       label = setNames(object = gsub(pattern = "\\(.*",
                                                      replacement = "",
                                                      x = levels(x = plot.data$regulon)),
                                        nm = levels(x = plot.data$regulon))) +
      scale_y_continuous(limits = c(-max.contrast,
                                    max.contrast),
                         breaks = c(-max.contrast,
                                    0,
                                    max.contrast),
                         expand = c(0, 0)) +
      scale_colour_gradient2(low = cluster.blended.colors[paste0(cell.type, " del")],
                             mid = BlendColorAlpha(color = cluster.blended.colors[paste0(cell.type, " wt")],
                                                   alpha = 0.7),
                             high = cluster.blended.colors[paste0(cell.type, " wt")],
                             midpoint = 0,
                             limits = c(-max.contrast,
                                        max.contrast)) +
      coord_cartesian(clip = "off") +
      labs(x = "regulon",
           y = "rss contrast<br>*Nr4a2<sup>wt/wt</sup>* &minus; *Nr4a2<sup>del/wt</sup>*",
           title = cell.type,
           colour = "cell types") +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6,
                                       angle = 90,
                                       hjust = 1,
                                       vjust = 0.5,
                                       family = "Arial",
                                       face = "plain",
                                       colour = ifelse(test = grepl(pattern = "\\+",
                                                                    x = levels(x = plot.data$regulon)),
                                                       yes = "#00008B",
                                                       no = "#CD0000")),
            axis.text.y = element_text(size = 6,
                                       family = "Arial",
                                       face = "plain",
                                       colour = "black"),
            axis.title = ggtext::element_markdown(size = 7,
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
            plot.title = element_text(size = 7,
                                      family = "Arial",
                                      colour = cluster.blended.colors[paste0(cell.type, " wt")],
                                      hjust = 0.5,
                                      vjust = 0,
                                      margin = margin(t = 0,
                                                      r = 0,
                                                      b = 0,
                                                      l = 0,
                                                      unit = "mm")),
            legend.position = "none",
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank())
    
    # Save plot in pdf format
    if (save.plot) {
      FixSizeAndSave(plot = rss.contrast.plot,
                     filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                          paste0(file.prefix, "_pyscenic_grnboost2_whole_transcriptome_rss_contrast_", cell.type, ".pdf")),
                     is.ggassemble = FALSE,
                     panel.width = 6.5,
                     panel.height = 4.2,
                     unit.use = "cm",
                     margin = 2.25,
                     use.ggsave = FALSE,
                     useDingbats = FALSE)
    }
    
    # Return plot
    return(rss.contrast.plot)
  },
  max.contrast = max.contrast.round,
  cluster.blended.colors = hemi.cluster.attributes$broad.clusters$cluster.blended.colors,
  save.plot = TRUE,
  file.prefix = file.prefix
)

# Plot RSS contrast for all data combined
rss.contrast.plot <- ggplot(data = rss.contrast.m,
                            mapping = aes(x = regulon,
                                          y = rss_contrast,
                                          colour = cell_type,
                                          alpha = alpha)) +
  geom_point(size = 2,
             shape = 19,
             stroke = 0,
             show.legend = TRUE) +
  scale_x_discrete(breaks = levels(x = rss.contrast.m$regulon),
                   label = setNames(object = gsub(pattern = "\\(.*",
                                                  replacement = "",
                                                  x = levels(x = rss.contrast.m$regulon)),
                                    nm = levels(x = rss.contrast.m$regulon))) +
  scale_y_continuous(limits = c(-max.contrast.round,
                                max.contrast.round),
                     breaks = c(-max.contrast.round,
                                0,
                                max.contrast.round),
                     expand = c(0, 0)) +
  scale_colour_manual(values = cluster.attributes$broad.clusters$cluster.colors[1:3],
                      breaks = cluster.attributes$broad.clusters$cluster.breaks[1:3]) +
  scale_alpha(limits = c(0,1),
              range = c(0.005,1),
              guide = "none") +
  guides(colour = ggplot2::guide_legend(override.aes = list(size = 2,
                                                            alpha = 1,
                                                            stroke = 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "regulon",
       y = "rss contrast<br>*Nr4a2<sup>wt/wt</sup>* &minus; *Nr4a2<sup>del/wt</sup>*",
       colour = "cell types") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   family = "Arial",
                                   face = "plain",
                                   colour = ifelse(test = grepl(pattern = "\\+",
                                                                x = levels(x = rss.contrast.m$regulon)),
                                                   yes = "#00008B",
                                                   no = "#CD0000")),
        axis.text.y = element_text(size = 6,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black"),
        axis.title = ggtext::element_markdown(size = 7,
                                              family = "Arial",
                                              face = "plain",
                                              colour = "black"),
        legend.text = element_text(size = 7,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black",
                                   margin = margin(t = 0,
                                                   r = 0,
                                                   b = 0,
                                                   l = 0,
                                                   unit = "mm")),
        legend.title = element_text(size = 7,
                                    family = "Arial",
                                    face = "plain",
                                    colour = "black"),
        axis.line = element_line(size = rel(x = 0.5),
                                 colour = "black",
                                 lineend = "square"),
        axis.ticks = element_line(size = rel(x = 0.5),
                                  colour = "black",
                                  lineend = "square"),
        axis.ticks.length = ggplot2::unit(x = 0.5,
                                          units = "mm"),
        plot.title = element_blank(),
        legend.position = "right",
        legend.margin = margin(t = 0,
                               r = 0,
                               b = 0,
                               l = 0,
                               unit = "mm"),
        legend.key.size = unit(x = 3,
                               units = "mm"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank())

# Save plot in pdf format
FixSizeAndSave(plot = rss.contrast.plot,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_pyscenic_grnboost2_whole_transcriptome_rss_contrast.pdf")),
               is.ggassemble = FALSE,
               panel.width = 6.5,
               panel.height = 4.2,
               unit.use = "cm",
               margin = 2.25,
               use.ggsave = FALSE,
               useDingbats = FALSE)

# Plot RSS for all data combined
rss.plot <- ggplot(data = rss.m,
                   mapping = aes(x = regulon,
                                 y = rss,
                                 colour = cell_type_genotype)) +
  geom_jitter(size = 2,
              shape = 19,
              stroke = 0,
              position = position_jitter(width = 0.2,
                                         height = 0,
                                         seed = 47),
              show.legend = TRUE) +
  scale_x_discrete(breaks = levels(x = rss.m$regulon),
                   label = setNames(object = gsub(pattern = "\\(.*",
                                                  replacement = "",
                                                  x = levels(x = rss.m$regulon)),
                                    nm = levels(x = rss.m$regulon))) +
  scale_y_continuous(limits = c(0,
                                max.rss),
                     breaks = seq(from = 0,
                                  to = max.rss,
                                  by = 0.2),
                     expand = c(0, 0)) +
  scale_colour_manual(values = hemi.cluster.attributes$broad.clusters$cluster.blended.colors[1:6],
                      breaks = hemi.cluster.attributes$broad.clusters$cluster.breaks[1:6]) +
  guides(colour = ggplot2::guide_legend(override.aes = list(size = 2,
                                                            alpha = 1,
                                                            stroke = 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "regulon",
       y = "regulon specificity score (rss)",
       colour = "cell types") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   family = "Arial",
                                   face = "plain",
                                   colour = ifelse(test = grepl(pattern = "\\+",
                                                                x = levels(x = rss.m$regulon)),
                                                   yes = "#00008B",
                                                   no = "#CD0000")),
        axis.text.y = element_text(size = 6,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black"),
        axis.title = ggtext::element_markdown(size = 7,
                                              family = "Arial",
                                              face = "plain",
                                              colour = "black"),
        legend.text = element_text(size = 7,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black",
                                   margin = margin(t = 0,
                                                   r = 0,
                                                   b = 0,
                                                   l = 0,
                                                   unit = "mm")),
        legend.title = element_text(size = 7,
                                    family = "Arial",
                                    face = "plain",
                                    colour = "black"),
        axis.line = element_line(size = rel(x = 0.5),
                                 colour = "black",
                                 lineend = "square"),
        axis.ticks = element_line(size = rel(x = 0.5),
                                  colour = "black",
                                  lineend = "square"),
        axis.ticks.length = ggplot2::unit(x = 0.5,
                                          units = "mm"),
        plot.title = element_blank(),
        legend.position = "right",
        legend.margin = margin(t = 0,
                               r = 0,
                               b = 0,
                               l = 0,
                               unit = "mm"),
        legend.key.size = unit(x = 3,
                               units = "mm"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank())

# Save plot in pdf format
FixSizeAndSave(plot = rss.plot,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_pyscenic_grnboost2_whole_transcriptome_rss.pdf")),
               is.ggassemble = FALSE,
               panel.width = 6.5,
               panel.height = 4.2,
               unit.use = "cm",
               margin = 2.25,
               use.ggsave = FALSE,
               useDingbats = FALSE)

# Plot regulons and their target genes
tf.target.plot.l <- lapply(
  X = regulon.params.l,
  FUN = function(regulon.params,
                 regulons,
                 regulon.labels,
                 file.prefix) {
    # Subset regulon data for TFs of interest
    regulons <- regulons[regulons$transcription_factor %in% regulon.params$tf.to.keep,]
    
    # Drop levels from regulon data
    regulons$regulon_name <- droplevels(x = regulons$regulon_name)
    regulons$transcription_factor <- droplevels(x = regulons$transcription_factor)
    
    # Subset regulon breaks and labels for those found in data
    regulon.labels <- regulon.labels[rev(x = levels(x = regulons$regulon_name))]
    
    # Get list of genes to highlight on regulon plot
    genes.to.highlight <- regulon.params$genes.to.highlight
    
    # Plot regulons and their target genes
    tf.target.plot <- ggplot(data = regulons,
                             mapping = aes(x = gene,
                                           y = regulon_name,
                                           fill = regulon_context)) +
      geom_tile(colour = NA) +
      scale_y_discrete(breaks = names(x = regulon.labels),
                       labels = regulon.labels) +
      scale_fill_manual(values = c("activating" = "#00008B",
                                   "repressing" = "#CD0000"),
                        breaks = c("activating",
                                   "repressing"),
                        na.value = "white") +
      guides(x.sec = ggh4x::guide_axis_manual(breaks = genes.to.highlight, 
                                              labels = genes.to.highlight,
                                              label_family = "Arial",
                                              label_face = "italic",
                                              label_colour = "black",
                                              label_size = 7,
                                              angle = 90)) +
      coord_cartesian(clip = "off") +
      labs(x = "modulated genes",
           y = "transcription factor") +
      facet_grid(rows = transcription_factor ~ modulation,
                 scales = "free",
                 space = "free",
                 switch = "both") +
      theme_classic() +
      # theme(axis.text.x = element_blank(),
      #       axis.text.y = element_text(size = 6,
      #                                  family = "Arial",
      #                                  face = "plain",
      #                                  colour = "black"),
      theme(axis.text = element_blank(),
            axis.text.x.top = element_text(),
            axis.title = element_text(size = 7,
                                      family = "Arial",
                                      face = "plain",
                                      colour = "black"),
            legend.text = element_text(size = 7,
                                       family = "Arial",
                                       face = "plain",
                                       colour = "black",
                                       margin = margin(t = 1,
                                                       r = 0,
                                                       b = 1,
                                                       l = 1,
                                                       unit = "mm")),
            legend.title = element_blank(),
            # axis.line = element_line(size = rel(x = 0.5),
            #                          colour = "black",
            #                          lineend = "square"),
            # axis.ticks = element_line(size = rel(x = 0.5),
            #                           colour = "black",
            #                           lineend = "square"),
            # axis.ticks.x = element_blank(),
            # axis.ticks.length = ggplot2::unit(x = 0.5,
            #                                   units = "mm"),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.x.top = element_line(size = rel(x = 0.5),
                                            colour = "black",
                                            lineend = "square"),
            axis.ticks.length.x.top = ggplot2::unit(x = 0.5,
                                                    units = "mm"),
            plot.title = element_blank(),
            legend.position = "right",
            legend.margin = margin(t = 0,
                                   r = 0,
                                   b = 0,
                                   l = 0,
                                   unit = "mm"),
            legend.key.size = unit(x = 1,
                                   units = "mm"),
            strip.text.x.bottom = element_text(size = 7,
                                               family = "Arial",
                                               face = "plain",
                                               colour = "black",
                                               margin = margin(t = 0,
                                                               r = 0,
                                                               b = 0,
                                                               l = 0,
                                                               unit = "mm")),
            strip.text.y.left = element_text(size = 7,
                                             angle = 0,
                                             hjust = 1,
                                             vjust = 0.5,
                                             family = "Arial",
                                             face = "italic",
                                             colour = "black",
                                             margin = margin(t = 0,
                                                             r = 0,
                                                             b = 0,
                                                             l = 0,
                                                             unit = "mm")),
            strip.placement = "outside",
            strip.clip = "off",
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black",
                                        fill = NA,
                                        linewidth = rel(x = 0.5),
                                        linetype = "solid"),
            panel.background = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank())
    
    # Save plot in pdf format
    FixSizeAndSave(plot = patchwork::wrap_plots(tf.target.plot),
                   filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                        paste0(file.prefix, "_pyscenic_grnboost2_tf_target_", regulon.params$file.suffix)),
                   is.ggassemble = TRUE,
                   panel.width = regulon.params$panel.width,
                   panel.height = regulon.params$panel.height,
                   unit.use = "cm",
                   margin = 2.5,
                   use.ggsave = TRUE,
                   useDingbats = FALSE)
    
    # Return plot
    return(tf.target.plot)
  },
  regulons = regulons,
  regulon.labels = regulon.labels,
  file.prefix = file.prefix
)

# Plot violin plots of modulated TFs
hemi.vln.plots <- StackViolinPlots(
  data.use = LayerData(object = obj,
                       layer = "data",
                       assay = "RNA_full"),
  genes.use = regulon.params.l$enriched_regulons$tf.to.keep,
  cluster.data = obj@meta.data[, "cell_type_genotype_abbr", drop = FALSE],
  cellnames.as.rownames.in.cluster.data = TRUE,
  cluster.ident.name = "cell_type_genotype_abbr",
  cluster.colors = hemi.cluster.attributes$subclusters$cluster.colors[1:18],
  cluster.breaks = hemi.cluster.attributes$subclusters$cluster.breaks[1:18],
  cluster.labels = setNames(object = hemi.cluster.attributes$subclusters$cluster.breaks[1:18],
                            nm = hemi.cluster.attributes$subclusters$cluster.breaks[1:18]),
  is.log.transformed = TRUE,
  log.scale = "log",
  pseudocount.use = 1,
  y.scale.trans = "log10",
  y.min = "0",
  point.size = 2,
  alpha.use = hemi.cluster.attributes$subclusters$cluster.alpha[1:18],
  vln.border.colour = NA,
  vln.border.stroke = 0.25,
  plot.title.size = 7,
  plot.title.angle = 45,
  plot.title.face = "italic",
  font.family = "Arial",
  import.font = FALSE,
  hjust.use = 0,
  vjust.use = 0,
  x.axis.title = "max norm. UMI ",
  axis.text.x.size = 6,
  axis.text.y.size = 7,
  axis.text.face = "plain",
  axis.ticks.length = 0.5,
  axis.line.size = rel(0.5),
  round.to.ceiling = TRUE,
  verbose = TRUE)

# Save plot in pdf format
FixSizeAndSave(plot = hemi.vln.plots,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                    paste0(file.prefix, "_pyscenic_grnboost2_violin_plots_hemi_modulated_tf.pdf")),
               is.ggassemble = TRUE,
               panel.width = 0.5,
               panel.height = 0.325 * 18,
               margin = 2.25,
               unit.use = "cm",
               use.ggsave = TRUE,
               useDingbats = FALSE)
