# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Load required packages
require(reshape2)
require(data.table)
require(glmmTMB)

# Load data
modulated.genes <- readxl::read_xlsx(path = "files_to_submit/Figure_4GH/e20221114-e20230605_Nr4a2_modulated_genes_raw_data.xlsx",
                                     na = "NA",
                                     col_types = rep(x = c("guess", "numeric"),
                                                     times = c(10, 9)))

# Keep cells having > 0 expression of either Oprk1 or Nr4a2 or both
cells.to.keep <- (modulated.genes$Oprk1 > 0 & !is.na(x = modulated.genes$Oprk1)) | (modulated.genes$Nr4a2 > 0 & !is.na(x = modulated.genes$Nr4a2))
modulated.genes <- modulated.genes[cells.to.keep,]

# Subset data for columns of interest
tested.genes <- c("Nr2f2_cytoplasmic", "Scn1b", "Ntm", "Syt17", "Ryr2", "Cdh13", "Rxfp1")
grouping.columns <- c("Image", "experiment", "genotype", "group")
modulated.genes <- modulated.genes[,c(grouping.columns, tested.genes)]

# Melt data
modulated.genes <- reshape2::melt(data = modulated.genes,
                                  id.vars = grouping.columns,
                                  variable.name = "gene",
                                  value.name = "expression_level")
modulated.genes$gene <- as.character(x = modulated.genes$gene)

# Remove rows with NA expression levels
modulated.genes <- modulated.genes[!is.na(x = modulated.genes$expression_level),]

# Round expression levels to convert data to integers (for model building)
modulated.genes$expression_level <- round(x = modulated.genes$expression_level)

# Set negative expression levels to 0
modulated.genes$expression_level <- ifelse(test = modulated.genes$expression_level < 0,
                                           yes = 0,
                                           no = modulated.genes$expression_level)

# Set WT genotype as the reference
modulated.genes$genotype <- factor(x = modulated.genes$genotype,
                                   levels = c("Nr4a2.wt.wt", "Nr4a2.del.wt"))

# Compute statistics on the data
lrt.data.l <- lapply(
  X = list("modulated.genes" = NA),
  FUN = function(void,
                 modulated.genes) {
    # Run glmmTMB
    models.l <- lapply(
      X = split(x = modulated.genes,
                f = modulated.genes$gene),
      FUN = function(gene.data) {
        # Run glmmTMB
        model <- glmmTMB::glmmTMB(formula = expression_level ~ genotype + (1 | Image),
                                  data = gene.data,
                                  family = glmmTMB::nbinom2(link = "log"))
        
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
      .id = "gene"
    )
    
    # Adjust p values using the holm correction method
    lrt.data$adjusted_p_value <- p.adjust(p = lrt.data[, "Pr(>Chi)"],
                                          method = "holm")
    
    # Return list of data
    return(list(model = models.l,
                lrt = lrt.l,
                lrt.data = lrt.data))
  },
  modulated.genes = modulated.genes
)

# Compute mean expression levels per image
modulated.genes.m <- data.table(modulated.genes)
modulated.genes.m <- modulated.genes.m[, list(mean_expression_level = mean(x = expression_level)),
                                       by = c(grouping.columns, "gene")]

# Compute ratio of expression level with the mean of WT animals
modulated.genes.m.l <- split(x = modulated.genes.m,
                             f = modulated.genes.m[,c("group", "gene")],
                             sep = "_")
empty.data <- lapply(X = modulated.genes.m.l, FUN = nrow) == 0
modulated.genes.m.l <- modulated.genes.m.l[!empty.data]
modulated.genes.m.l <- lapply(X = modulated.genes.m.l,
                              FUN = function(gene.data) {
                                # Compute mean of WT images
                                gene.data$mean_expression_level_wt <- mean(x = gene.data$mean_expression_level[gene.data$genotype == "Nr4a2.wt.wt"])
                                # Compute the ratio for each image relative to the mean of WT images
                                gene.data$ratio <- gene.data$mean_expression_level / gene.data$mean_expression_level_wt
                                # Transform ratio to log2
                                gene.data$log2_ratio <- log2(x = gene.data$ratio)
                                # Return data
                                return(gene.data)
                              })
modulated.genes.m <- dplyr::bind_rows(modulated.genes.m.l,
                                      .id = NULL)

# Save mean and ratio data
write.csv(x = as.data.frame(x = modulated.genes.m),
          file = "files_to_submit/Figure_4GH/e20221114-e20230605_Nr4a2_modulated_genes_mean_data_20240911.csv",
          quote = FALSE,
          row.names = FALSE)
