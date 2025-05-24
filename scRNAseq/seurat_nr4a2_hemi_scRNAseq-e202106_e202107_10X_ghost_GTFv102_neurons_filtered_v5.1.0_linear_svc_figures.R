# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/scRNAseq-e202106_e202107_10X_cluster_attributes.R")
source(file = "/Users/leonfodoulian/scData/ExpandClusterAttributes.R")
source(file = "/Users/leonfodoulian/scData/BlendColorAlpha.R")
source(file = "/Users/leonfodoulian/scData/FixSizeAndSave.R")

# Required packages
require(glmmTMB)
require(extrafont)
require(ggplot2)
require(ggtext)
require(ggh4x)

# Create new directories
dir.create(path = "./results/smim32_GTFv102/neurons_filtered/figures/temp")

# Define cluster attributes to be used
cluster.attributes <- cluster.attributes$neurons_filtered

# Expand cluster attributes with genotype labels
hemi.cluster.attributes <- ExpandClusterAttributes(
  cluster.attributes = cluster.attributes$broad.clusters,
  labels = c("wt", "del"),
  alpha = c(1, 0.4),
  prefix = TRUE
)

# Define abbreviations of genotypes
genotype.abbreviations <- setNames(object = c("wt", "del"),
                                   nm = c("Nr4a2(WT/WT)", "Nr4a2(SA-IRES-Dre/WT)"))

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load predicted labels data using linear SVC
svm.data.l <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                       paste0(file.prefix, "_linear_svc_predicted_labels.rds")))

# Compute proportions of predictions by actual label for observed and permuted data
accuracy.data.l <- lapply(
  X = svm.data.l$predictions,
  FUN = function(predicted.labels) {
    # Compute proportions of predictions by actual label for each genotype
    accuracy.data <- dplyr::bind_rows(
      lapply(
        X = split(x = predicted.labels,
                  f = predicted.labels[, c("iteration", "genotype")],
                  sep = ":"),
        FUN = function(iteration.data) {
          # Compute confusion matrix of actual and predicted labels
          con.mat <- table(actual_label = iteration.data$actual,
                           predicted_label = iteration.data$predicted)
          
          # Compute proportions of predictions by actual label
          con.mat <- con.mat / rowSums(x = con.mat)
          
          # Transform table to data frame
          con.mat <- as.data.frame(x = con.mat,
                                   responseName = "proportion")
          
          # Return data
          return(con.mat)
        }
      ),
      .id = "group"
    )
    
    # Add iteration and genotype information to the data
    accuracy.data[, c("iteration", "genotype")] <- setNames(
      object = data.frame(stringr::str_split(string = accuracy.data$group,
                                             pattern = ":",
                                             simplify = TRUE)),
      nm = c("iteration", "genotype"))
    
    # Return data
    return(accuracy.data)
  }
)

# Compute statistics on the data using a binomial regression with logit family
lrt.data.l <- lapply(
  X = accuracy.data.l,
  FUN = function(accuracy.data) {
    # Compute statistics on the data using a binomial regression with logit family
    models.l <- lapply(
      X = split(x = accuracy.data,
                f = accuracy.data[, c("actual_label", "predicted_label")],
                sep = "_"),
      FUN = function(pair.data) {
        # Set WT genotype as the reference
        pair.data$genotype <- factor(x = pair.data$genotype,
                                     levels = c("Nr4a2(WT/WT)", "Nr4a2(SA-IRES-Dre/WT)"))
        
        # Run glmmTMB
        model <- glmmTMB::glmmTMB(formula = proportion ~ genotype,
                                  data = pair.data,
                                  family = binomial(link = "logit"))
        
        # Return model
        return(model)
      }
    )
    
    # Compute p values using a LRT
    lrt.l <- lapply(X = models.l,
                    FUN = drop1,
                    test = "Chisq")
    
    # Transform to data frame
    lrt.data <- dplyr::bind_rows(
      lapply(
        X = lrt.l,
        FUN = function(lrt) {
          return(data.frame(lrt["genotype",],
                            row.names = NULL,
                            check.names = FALSE))
        }
      ),
      .id = "group"
    )
    
    # Add actual and predicted labels information to the data
    lrt.data[, c("actual_label", "predicted_label")] <- setNames(
      object = data.frame(stringr::str_split(string = lrt.data$group,
                                             pattern = "_",
                                             simplify = TRUE)),
      nm = c("actual_label", "predicted_label"))
    
    # Adjust p values using the holm correction method
    lrt.data$adjusted_p_value <- p.adjust(p = lrt.data[, "Pr(>Chi)"],
                                          method = "holm")
    
    # Return list of data
    return(list(models = models.l,
                lrt = lrt.l,
                lrt.data = lrt.data))
  }
)

# Plot proportions of correctly classified cells per actual-predicted labels pair and genotype
accuracy.plot.l <- mapply(
  plot.data = accuracy.data.l,
  label.data = lapply(X = lrt.data.l,
                      FUN = "[[",
                      "lrt.data"),
  file.suffix = names(x = accuracy.data.l),
  FUN = function(plot.data,
                 label.data,
                 file.suffix,
                 genotype.abbreviations,
                 cluster.attributes,
                 hemi.cluster.attributes,
                 nested,
                 save.plot,
                 file.prefix) {
    # Add group column to plot data
    plot.data$group <- paste(plot.data$actual_label,
                             plot.data$predicted_label,
                             sep = "_")
    
    # Add abbreviated genotypes to plot data
    plot.data$genotype_abbr <- genotype.abbreviations[plot.data$genotype]
    
    # Create a column combining genotype and cell type identities of cells 
    plot.data$cell_type_genotype_abbr <- paste(plot.data$actual_label,
                                               plot.data$genotype_abbr)
    
    # Transform actual label, predicted label and genotype to factor
    plot.data$actual_label <- factor(x = plot.data$actual_label,
                                     levels = cluster.attributes$broad.clusters$cluster.breaks[1:3])
    plot.data$predicted_label <- factor(x = plot.data$predicted_label,
                                        levels = cluster.attributes$broad.clusters$cluster.breaks[1:3])
    plot.data$genotype <- factor(x = plot.data$genotype,
                                 levels = c("Nr4a2(WT/WT)",
                                            "Nr4a2(SA-IRES-Dre/WT)"))
    
    # Order plot data by actual label, predicted label and genotype
    plot.data <- plot.data[order(plot.data$actual_label, plot.data$predicted_label, plot.data$genotype, decreasing = FALSE),]
    
    # Get unique group values from plot data
    group <- setNames(object = plot.data$group[!duplicated(x = plot.data$group)],
                      nm = plot.data$predicted_label[!duplicated(x = plot.data$group)])
    
    # Transform group and cell type-abbreviated genotype pairs to factor
    plot.data$group <- factor(x = plot.data$group,
                              levels = group)
    plot.data$cell_type_genotype_abbr <- factor(x = plot.data$cell_type_genotype_abbr,
                                                levels = unique(x = plot.data$cell_type_genotype_abbr))
    
    # Add * label for significant differences in classification
    label.data$label <- ifelse(test = label.data$adjusted_p_value < 0.05,
                               yes = "*",
                               no = "")
    
    # Order label data to match order of unique group values
    label.data <- label.data[match(x = group, table = label.data$group),]
    
    # Plot proportions of correctly classified cells per actual-predicted labels pair and genotype
    if (isTRUE(x = nested)) {
      # Nested plot
      accuracy.plot <- ggplot(data = plot.data,
                              mapping = aes(x = interaction(predicted_label,
                                                            actual_label), # nested
                                            y = proportion,
                                            colour = cell_type_genotype_abbr,
                                            fill = cell_type_genotype_abbr)) +
        geom_point(data = plot.data[!duplicated(x = plot.data$genotype),],
                   mapping = aes(x = interaction(predicted_label,
                                                 actual_label), # nested
                                 y = proportion,
                                 alpha = genotype),
                   shape = NA, # trick to add point legend without plotting the data
                   show.legend = TRUE,
                   inherit.aes = FALSE) +
        geom_vline(xintercept = c(3.5, 6.5), # separation between nested plots
                   colour = "#D3D3D3",
                   linetype = "dashed",
                   linewidth = 0.25) + # nested
        scale_x_discrete(guide = "axis_nested") + # nested
        guides(x.sec = ggh4x::guide_axis_manual(labels = label.data$label,
                                                breaks = interaction(label.data$predicted_label,
                                                                     label.data$actual_label))) + # nested
        labs(x = "actual-predicted cell type pair") # nested
    } else {
      # Facetted plot
      accuracy.plot <- ggplot(data = plot.data,
                              mapping = aes(x = group, # facet
                                            y = proportion,
                                            colour = cell_type_genotype_abbr,
                                            fill = cell_type_genotype_abbr)) +
        geom_point(data = plot.data[!duplicated(x = plot.data$genotype),],
                   mapping = aes(x = group, # facet
                                 y = proportion,
                                 alpha = genotype),
                   shape = NA, # trick to add point legend without plotting the data
                   show.legend = TRUE,
                   inherit.aes = FALSE) +
        scale_x_discrete(breaks = group,
                         labels = names(x = group)) + # facet
        guides(x.sec = ggh4x::guide_axis_manual(labels = label.data$label,
                                                breaks = label.data$group)) + # facet
        labs(x = "predicted cell type", # facet
             title = "actual cell type") + # facet
        ggh4x::facet_wrap2(
          facets = ~ actual_label,
          nrow = 1,
          scales = "free_x",
          strip = ggh4x::strip_themed(
            background_x = element_blank(),
            text_x = ggh4x::elem_list_text(size = rep(x = 7,
                                                      times = 3),
                                           vjust = rep(x = 7,
                                                       times = 3),
                                           family = rep(x = "Arial",
                                                        times = 3),
                                           face = rep(x = "plain",
                                                      times = 3),
                                           colour = cluster.attributes$broad.clusters$cluster.colors[1:3])
          )
        ) # facet
    }
    
    # Final plot
    accuracy.plot <- accuracy.plot +
      geom_violin(position = position_dodge(width = 0.85),
                  scale = "width",
                  adjust = 1,
                  size = 1,
                  width = 0.5,
                  show.legend = FALSE) + # violin plot
      scale_y_continuous(limits = c(0,1),
                         breaks = seq(from = 0,
                                      to = 1,
                                      by = 0.25),
                         labels = seq(from = 0,
                                      to = 100,
                                      by = 25),
                         expand = c(0,0)) +
      scale_colour_manual(values = hemi.cluster.attributes$cluster.blended.colors[1:6]) +
      scale_fill_manual(values = hemi.cluster.attributes$cluster.blended.colors[1:6]) +
      scale_alpha_manual(values = c("Nr4a2(WT/WT)" = 1,
                                    "Nr4a2(SA-IRES-Dre/WT)" = 0.4),
                         breaks = c("Nr4a2(WT/WT)",
                                    "Nr4a2(SA-IRES-Dre/WT)"),
                         labels = c("*Nr4a2<sup>wt/wt</sup>*",
                                    "*Nr4a2<sup>del/wt</sup>*")) +
      guides(colour = "none",
             fill = "none",
             alpha = guide_legend(override.aes = list(colour = "black",
                                                      size = 2,
                                                      shape = 19,
                                                      stroke = 0))) +
      coord_cartesian(clip = "off") +
      labs(y = "percentage of classified cells",
           alpha = "genotype") +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6,
                                       family = "Arial",
                                       face = "plain",
                                       colour = cluster.attributes$broad.clusters$cluster.colors[1:3]),
            axis.text.x.top = element_text(size = 15,
                                           vjust = ifelse(test = nested, yes = 0, no = -7),
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
            legend.text = ggtext::element_markdown(size = 7,
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
            axis.line.x.top = element_blank(),
            axis.ticks = element_line(size = rel(x = 0.5),
                                      colour = "black",
                                      lineend = "square"),
            axis.ticks.x.top = element_blank(),
            axis.ticks.length = ggplot2::unit(x = 0.5,
                                              units = "mm"),
            legend.position = "right",
            legend.margin = margin(t = 0,
                                   r = 0,
                                   b = 0,
                                   l = 0,
                                   unit = "mm"),
            legend.key.size = unit(x = 3,
                                   units = "mm"),
            strip.background = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank())
    
    if (isTRUE(x = nested)) {
      accuracy.plot <- accuracy.plot +
        theme(ggh4x.axis.nestline = element_line(linewidth = rel(x = 1),
                                                 colour = cluster.attributes$broad.clusters$cluster.colors[1:3]), # nested
              plot.title = element_blank(), # nested
              strip.text = element_blank()) # nested
    } else {
      accuracy.plot <- accuracy.plot +
        theme(plot.title = ggtext::element_markdown(size = 7,
                                                    family = "Arial",
                                                    colour = "black",
                                                    hjust = 0.5,
                                                    vjust = 0,
                                                    margin = margin(t = 0,
                                                                    r = 0,
                                                                    b = 0,
                                                                    l = 0,
                                                                    unit = "mm")), # facet
              panel.spacing = unit(x = 2, units = "mm")) # facet
    }
    
    # Save plot in pdf format
    if (save.plot) {
      FixSizeAndSave(plot = patchwork::wrap_plots(accuracy.plot),
                     filename = file.path("results/smim32_GTFv102/neurons_filtered/figures/temp",
                                          paste0(file.prefix, "_linear_svc_plot_", file.suffix, ".pdf")),
                     is.ggassemble = TRUE,
                     panel.width = 6.6,
                     panel.height = 4.2,
                     margin = 0,
                     unit.use = "cm",
                     use.ggsave = TRUE,
                     useDingbats = FALSE)
    }
    
    # Return plot
    return(accuracy.plot)
  },
  MoreArgs = list(genotype.abbreviations = genotype.abbreviations,
                  cluster.attributes = cluster.attributes,
                  hemi.cluster.attributes = hemi.cluster.attributes,
                  nested = TRUE,
                  save.plot = TRUE,
                  file.prefix = file.prefix),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)
