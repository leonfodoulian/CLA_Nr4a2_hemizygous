# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Load required packages
require(glmmTMB)

# List files
files <- list.files(path = "files_to_submit/Figure_4IJ/",
                    pattern = "clustered_data.csv",
                    full.names = TRUE)
names(x = files) <- gsub(pattern = ".*QuPath_|_clustered.*",
                         replacement = "",
                         x = files)

# Compute statistics on the data
lrt.data.l <- mapply(
  file = files,
  cluster.names = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = c("1" = "Shell",
                                                      "2" = "noise",
                                                      "3" = "Ntm+",
                                                      "4" = "CLA",
                                                      "5" = "Syt17+"),
                       "C2.Nr4a2_C3.Ryr2_C4.Syt17" = c("1" = "Ryr2+",
                                                       "2" = "noise",
                                                       "3" = "Syt17+",
                                                       "4" = "CLA",
                                                       "5" = "Shell")),
  count.col = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = "Ntm",
                   "C2.Nr4a2_C3.Ryr2_C4.Syt17" = "Ryr2"),
  FUN = function(file,
                 cluster.names,
                 count.col,
                 cluster.col,
                 cluster.to.remove) {
    # Load data
    rnascope.data <- read.csv(file = file,
                              header = TRUE,
                              sep = ",")
    
    # Remove not clustered cells
    rnascope.data <- rnascope.data[rnascope.data$is_clustered, ]
    
    # Define image name for each cell
    rnascope.data$image <- paste(rnascope.data$slide,
                                 rnascope.data$zone,
                                 rnascope.data$mouse,
                                 rnascope.data$genotype,
                                 rnascope.data$level,
                                 paste0(rnascope.data$probe_set,
                                        ".tif"),
                                 sep = "_")
    
    # Rename clusters
    rnascope.data$clustering_results <- cluster.names[rnascope.data[[cluster.col]]]
    
    # Remove clusters from data if necessary
    rnascope.data <- rnascope.data[!(rnascope.data$clustering_results %in% cluster.to.remove),]
    
    # Set WT genotype as the reference
    rnascope.data$genotype <- factor(x = rnascope.data$genotype,
                                     levels = c("Nr4a2.wt.wt", "Nr4a2.del.wt"))
    
    # Round counts to convert data to integers (for model building)
    rnascope.data$expression_level <- round(x = rnascope.data[[count.col]])
    
    # Run glmmTMB
    models.l <- lapply(
      X = split(x = rnascope.data,
                f = rnascope.data$clustering_results),
      FUN = function(cluster.data,
                     count.col) {
        # Run glmmTMB
        model <- glmmTMB::glmmTMB(formula = expression_level ~ genotype + (1 | image),
                                  data = cluster.data,
                                  family = glmmTMB::nbinom2(link = "log"))
        
        # Return model
        return(model)
      },
      count.col = count.col
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
      .id = "cluster"
    )
    
    # Adjust p values using the holm correction method
    lrt.data$adjusted_p_value <- p.adjust(p = lrt.data[, "Pr(>Chi)"],
                                          method = "holm")
    
    # Return list of data
    return(list(model = models.l,
                lrt = lrt.l,
                lrt.data = lrt.data))
  },
  MoreArgs = list(cluster.col = "k5_clusters",
                  cluster.to.remove = "noise"),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)
