# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/CellTypeProportions.R")
source(file = "/Users/leonfodoulian/scData/StackViolinPlots.R")
source(file = "/Users/leonfodoulian/scData/StackFractionPlots.R")
source(file = "/Users/leonfodoulian/scData/GetFractionOfExpressingCells.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(Seurat)
require(extrafont)
require(ggplot2)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp/cla_ep_analysis",
           recursive = TRUE)

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Load list of data from CLA/EP complex analysis
cla.ep.analysis <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                            paste0(file.prefix, "_cla_ep_analysis.rds")))

# Join layers of the RNA assay
if (!"RNA_full" %in% names(x = obj@assays)) {
  obj[["RNA_full"]] <- JoinLayers(object = obj[["RNA"]])
}

# Prepare data for violin plots
obj.wt <- subset(x = obj,
                 subset = genotype == "Nr4a2(WT/WT)")

# Add region attribute of each cell to the scRNAseq Seurat object
obj.wt@meta.data$region <- unname(obj = cla.ep.analysis$region.attributes[rownames(x = obj.wt@meta.data)])

# Plot violin and fraction plots of region-specific genes
cla.ep.plots.l <- mapply(
  cla.ep.genes = cla.ep.analysis$markers,
  region.pair = names(x = cla.ep.analysis$markers),
  FUN = function(cla.ep.genes,
                 region.pair,
                 obj,
                 cluster.colors,
                 cluster.labels,
                 regions.to.keep,
                 save.plot,
                 path,
                 file.prefix) {
    # Stop function if no modulated genes are found
    if (length(x = cla.ep.genes) == 0) {
      return(NULL)
    }
    
    # Add region attribute to the data
    cla.ep.genes$region <- ifelse(test = cla.ep.genes$avg_log2FC > 0,
                                  yes = gsub(pattern = " vs.*",
                                             replacement = "",
                                             x = region.pair),
                                  no = gsub(pattern = ".*vs ",
                                            replacement = "",
                                            x = region.pair))
    
    # Plot violin and fraction plots of region-specific genes for each region
    cla.ep.plots.l <- lapply(
      X = split(x = cla.ep.genes,
                f = cla.ep.genes$region),
      FUN = function(region.genes,
                     obj,
                     region.pair,
                     cluster.colors,
                     cluster.labels,
                     regions.to.keep,
                     save.plot,
                     path,
                     file.prefix) {
        # Subset Seurat object for regions to keep
        obj <- subset(x = obj,
                      subset = region %in% regions.to.keep)
        
        # Get query region attribute
        query.region <- unique(x = region.genes$region)
        
        # Get reference region attribute
        reference.region <- setdiff(x = unlist(x = strsplit(x = region.pair,
                                                            split = " vs ")),
                                    y = query.region)
        
        # Define file suffix
        file.suffix <- gsub(pattern = "\\-|\\s",
                            replacement = "_",
                            x = paste(query.region,
                                      "specific_genes_vs",
                                      reference.region,
                                      sep = "_"))
        
        # Plot violin plots of region-specific genes
        vln.plots <- StackViolinPlots(
          data.use = LayerData(object = obj,
                               layer = "data",
                               assay = "RNA_full"),
          genes.use = rownames(x = region.genes),
          cluster.data = obj@meta.data[, "region", drop = FALSE],
          cellnames.as.rownames.in.cluster.data = TRUE,
          cluster.ident.name = "region",
          cluster.colors = cluster.colors[regions.to.keep],
          cluster.breaks = names(x = cluster.labels[regions.to.keep]),
          cluster.labels = cluster.labels[regions.to.keep],
          is.log.transformed = TRUE,
          log.scale = "log",
          pseudocount.use = 1,
          y.scale.trans = "log10",
          y.min = "0",
          point.size = 2,
          alpha.use = 1,
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
        if (save.plot) {
          FixSizeAndSave(plot = vln.plots,
                         filename = file.path(path,
                                              paste0(file.prefix, "_cla_ep_analysis_violin_plots_", file.suffix, ".pdf")),
                         is.ggassemble = TRUE,
                         panel.width = 0.5,
                         panel.height = 0.325 * length(x = regions.to.keep),
                         margin = 2.5,
                         unit.use = "cm",
                         use.ggsave = TRUE,
                         useDingbats = FALSE)
        }
        
        # Plot fraction plots of region-specific genes
        fraction.plots <- StackFractionPlots(
          data.use = LayerData(object = obj,
                               layer = "counts",
                               assay = "RNA_full"),
          genes.use = rownames(x = region.genes),
          cluster.data = obj@meta.data[, "region", drop = FALSE],
          cellnames.as.rownames.in.cluster.data = TRUE,
          cluster.ident.name = "region",
          cluster.colors = cluster.colors[regions.to.keep],
          cluster.breaks = names(x = cluster.labels[regions.to.keep]),
          cluster.labels = cluster.labels[regions.to.keep],
          is.log.transformed = FALSE,
          log.scale = "log",
          pseudocount.use = 1,
          min.exp = 1,
          thresh.plot = 15,
          alpha.use = 1,
          bar.border.colour = NA,
          bar.border.stroke = 0.25,
          plot.title.size = 7,
          plot.title.angle = 45,
          plot.title.face = "italic",
          font.family = "Arial",
          import.font = FALSE,
          hjust.use = 0,
          vjust.use = 0,
          x.axis.title = "Fraction of \nexpressing cells (%)",
          axis.text.x.size = 6,
          axis.text.y.size = 7,
          axis.text.face = "plain",
          axis.ticks.length = 0.5,
          axis.line.size = rel(0.5),
          verbose = TRUE)
        
        # Save plot in pdf format
        if (save.plot) {
          FixSizeAndSave(plot = fraction.plots,
                         filename = file.path(path,
                                              paste0(file.prefix, "_cla_ep_analysis_fraction_plots_", file.suffix, ".pdf")),
                         is.ggassemble = TRUE,
                         panel.width = 0.5,
                         panel.height = 0.325 * length(x = regions.to.keep),
                         margin = 2.5,
                         unit.use = "cm",
                         use.ggsave = TRUE,
                         useDingbats = FALSE)
        }
        
        # Return list of plot
        return(list(vln.plots = vln.plots,
                    fraction.plots = fraction.plots))
      },
      obj = obj,
      region.pair = region.pair,
      cluster.colors = cluster.colors,
      cluster.labels = cluster.labels,
      regions.to.keep = regions.to.keep,
      save.plot = save.plot,
      path = path,
      file.prefix = file.prefix)
    
    # Return list of plot
    return(cla.ep.plots.l)
  },
  MoreArgs = list(obj = obj.wt,
                  cluster.colors = c("CLA" = "#EB008D",
                                     "EP" = "#00A69C",
                                     "CLA-like L6" = "#C9B3DA"),
                  cluster.labels = c("CLA" = "CLA",
                                     "EP" = "EP",
                                     "CLA-like L6" = "CLA-like L6"),
                  regions.to.keep = c("CLA", "EP"),
                  save.plot = TRUE,
                  path = "./results/smim32_GTFv102/neurons_filtered/figures/temp/cla_ep_analysis",
                  file.prefix = file.prefix),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)

# Calculate and plot cell type proportions
meta.data <- cla.ep.analysis$meta.data
meta.data$sample <- paste(meta.data$experiment,
                          meta.data$genotype,
                          meta.data$mouse,
                          meta.data$level,
                          meta.data$cla_side,
                          meta.data$region,
                          sep = "_")
meta.data <- meta.data[meta.data$genotype != "Nr4a2.del.wt" & meta.data$level == "44",]
cell.type.prop.l <- CellTypeProportions(
  x = NULL,
  clusters = factor(x = meta.data$predicted.cell_type_abbr,
                    levels = c("CLA 1", "CLA 2")),
  sample = factor(x = meta.data$sample),
  group = factor(x = meta.data$region,
                 levels = c("EP", "CLA")),
  trend = FALSE,
  robust = FALSE,
  transform = "logit",
  return.plot = TRUE,
  markdown = FALSE,
  reference.group = "CLA",
  pseudocount = 0,
  cluster.colors = cluster.attributes$subclusters$cluster.colors[c("CLA 1", "CLA 2")],
  bar.border.stroke = 0.5,
  label.size = 15,
  label.position = "top",
  point.size = 1.25,
  sd.line.size = 0.5,
  y.axis.title = "log<sub>2</sub> ratio",
  font.family = "Arial",
  font.face = "plain",
  axis.text.size = 6,
  axis.title.size = 7,
  axis.ticks.length = 0.5,
  axis.line.size = rel(x = 0.5)
)

# Prepare CLA 1 and CLA 2 proportion data for plotting
cla.cell.type.prop <- cla.ep.analysis$cla.cell.type.prop
cla.cell.type.prop <- cla.cell.type.prop[cla.cell.type.prop$genotype != "Nr4a2.del.wt" & cla.cell.type.prop$level == "44",]
cla.cell.type.prop <- reshape2::melt(
  data = cla.cell.type.prop,
  id.vars = c("experiment", "genotype", "mouse", "level", "cla_side", "region"),
  variable.name = "cla_cell_type",
  value.name = "proportion"
)
cla.cell.type.prop$cla_cell_type <- toupper(
  x = gsub(
    pattern = "_",
    replacement = " ",
    x = cla.cell.type.prop$cla_cell_type 
  )
)

# Plot proportion of CLA 1 and CLA 2 cells in each region of the CLA/EP complex
prop.plot.r <- ggplot(data = cla.cell.type.prop,
                      mapping = aes(x = region,
                                    y = proportion,
                                    color = cla_cell_type)) +
  stat_summary(geom = "col",
               size = 0.5,
               fill = NA,
               width = 0.75,
               position = position_dodge(width = 0.9,
                                         preserve = "total"),
               fun = "mean",
               show.legend = FALSE) +
  stat_summary(geom = "segment",
               size = 0.5,
               position = position_dodge(width = 0.9,
                                         preserve = "total"),
               fun.data = function(x) {
                 return(c(y = mean(x = x) - sd(x = x),
                          yend = mean(x = x) + sd(x = x)))
               },
               show.legend = FALSE) +
  geom_point(size = 1.25,
             shape = 19,
             stroke = 0,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0,
                                             dodge.width = 0.9,
                                             seed = 53)) +
  scale_x_discrete(limits = c("CLA", "EP")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(from = 0,
                                  to = 1,
                                  by = 0.1),
                     expand = c(0,0)) +
  scale_color_manual(values = cluster.attributes$subclusters$cluster.colors[c("CLA 1", "CLA 2")]) +
  guides(colour = ggplot2::guide_legend(override.aes = list(size = 2,
                                                            alpha = 1,
                                                            stroke = 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "localization",
       y = "proportion of cells",
       colour = "cell types") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   family = "Arial",
                                   face = "plain",
                                   colour = c("CLA" = "#EB008D",
                                              "EP" = "#00A69C")),
        axis.text.y = element_text(size = 6,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black"),
        axis.title = element_text(size = 7,
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
FixSizeAndSave(plot = prop.plot.r,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp/cla_ep_analysis",
                                    paste0(file.prefix, "_CLA_cell_type_proportions_by_localization.pdf")),
               is.ggassemble = FALSE,
               panel.width = 0.45 * 2 * 2,
               panel.height = 4.5,
               unit.use = "cm",
               margin = 0,
               use.ggsave = FALSE,
               useDingbats = FALSE)

# Plot proportion of CLA/EP complex regions distribution across CLA 1 and CLA 2 cell types
prop.plot.c <- ggplot(data = cla.cell.type.prop,
                      mapping = aes(x = cla_cell_type,
                                    y = proportion,
                                    color = region)) +
  stat_summary(geom = "col",
               size = 0.5,
               fill = NA,
               width = 0.75,
               position = position_dodge(width = 0.9,
                                         preserve = "total"),
               fun = "mean",
               show.legend = FALSE) +
  stat_summary(geom = "segment",
               size = 0.5,
               position = position_dodge(width = 0.9,
                                         preserve = "total"),
               fun.data = function(x) {
                 return(c(y = mean(x = x) - sd(x = x),
                          yend = mean(x = x) + sd(x = x)))
               },
               show.legend = FALSE) +
  geom_point(size = 1.25,
             shape = 19,
             stroke = 0,
             position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0,
                                             dodge.width = 0.9,
                                             seed = 53)) +
  scale_x_discrete(limits = c("CLA 1", "CLA 2")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(from = 0,
                                  to = 1,
                                  by = 0.1),
                     expand = c(0,0)) +
  scale_color_manual(values = c("CLA" = "#EB008D",
                                "EP" = "#00A69C")) +
  guides(colour = ggplot2::guide_legend(override.aes = list(size = 2,
                                                            alpha = 1,
                                                            stroke = 0))) +
  coord_cartesian(clip = "off") +
  labs(x = "cell types",
       y = "proportion of cells",
       colour = "localization") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 6,
                                   angle = 45,
                                   hjust = 1,
                                   vjust = 1,
                                   family = "Arial",
                                   face = "plain",
                                   colour = cluster.attributes$subclusters$cluster.colors[c("CLA 1", "CLA 2")]),
        axis.text.y = element_text(size = 6,
                                   family = "Arial",
                                   face = "plain",
                                   colour = "black"),
        axis.title = element_text(size = 7,
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
FixSizeAndSave(plot = prop.plot.c,
               filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp/cla_ep_analysis",
                                    paste0(file.prefix, "_localization_proportions_by_CLA_cell_type.pdf")),
               is.ggassemble = FALSE,
               panel.width = 0.45 * 2 * 2,
               panel.height = 4.5,
               unit.use = "cm",
               margin = 0,
               use.ggsave = FALSE,
               useDingbats = FALSE)
