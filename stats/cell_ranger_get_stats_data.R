setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")
path.l <- c("scRNAseq-e20210604_e20210611/cellranger_count_ghost_GTFv102",
            "scRNAseq-e20210719_e20210720/cellranger_count_ghost_GTFv102")

dirs <- unlist(x = lapply(X = path.l,
                          FUN = function(path) {
                            samples <- list.dirs(path = file.path(".", path),
                                                 full.names = FALSE,
                                                 recursive = FALSE)
                            dirs <- file.path(getwd(), path, samples, "outs/filtered_feature_bc_matrix")
                            names(x = dirs) <- samples
                            return(dirs)
                          }))

files <- list.files(path = path.l,
                    recursive = TRUE,
                    full.names = TRUE)
files <- grep(pattern = "outs/metrics_summary.csv", x = files, value = TRUE)
names(x = files) <- files
x <- dplyr::bind_rows(lapply(X = files,
                             FUN = read.csv,
                             colClasses = "character"),
                      .id = "sample")
round(mean(as.numeric(gsub(pattern = ",", replacement = "", x = x$Estimated.Number.of.Cells))))
round(mean(as.numeric(gsub(pattern = ",", replacement = "", x = x$Mean.Reads.per.Cell))))
round(mean(as.numeric(gsub(pattern = ",", replacement = "", x = x$Median.Genes.per.Cell))))
round(mean(as.numeric(gsub(pattern = ",", replacement = "", x = x$Median.UMI.Counts.per.Cell))))
