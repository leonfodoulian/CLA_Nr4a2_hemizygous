# seurat_nr4a2_hemi_MERSCOPE_Rodriguez_Claustrum300_claustro_insular_region_v5.1.0_figures.R

# Add real genotypes to meta data
obj@meta.data$real_genotype <- c("C57BL6" = "C57BL6",
                                 "Nr4a2.del.wt" = "Nr4a2(SA-IRES-Dre/WT)",
                                 "Nr4a2.wt.wt" = "Nr4a2(WT/WT)")[obj@meta.data$genotype]

# Add bregma to meta data
obj@meta.data$bregma <- c("44" = "+1.045",
                          "53" = "+0.145",
                          "68" = "-1.055")[obj@meta.data$level]

# Prepare source data
source.data <- list(
  "MERSCOPE" = data.frame(
    cbind(
      setNames(
        object = obj@meta.data[,c("mouse", "gender", "real_genotype", "age", "bregma", "cla_side", "centroid_x_rot", "centroid_y_rot",
                                  "volume", "ngene", "ncount", "mapping.score",
                                  "predicted.broad_cell_type_abbr", "predicted.broad_cell_type_abbr.score",
                                  "predicted.cell_type_abbr", "predicted.cell_type_abbr.score")],
        nm = c("mouse", "gender", "genotype", "age", "bregma", "cla_side", "center_x", "center_y", "volume", "n_gene", "n_count", "mapping_score",
               "predicted_cell_type", "predicted_cell_type_score", "predicted_cell_subtype", "predicted_cell_subtype_score")
      ),
      t(
        x = as.matrix(
          x = LayerData(
            object = obj,
            layer = "counts",
            assay = "MERSCOPE",
            features = c("Ccn2", "Col6a1", "Drd1", "Gnb4", "Nr4a2", "Rorb", "Rprm", "Sdk2", "Sst", "Syt17")
          )
        )
      ),
      obj@reductions$umap_harmony@cell.embeddings
    )
  )
)

# Save file in xlsx format
for (n in names(x = source.data)) {
  do::write_xlsx(source.data[[n]],
                 sheet = n,
                 file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_LF.xlsx",
                 col.names = TRUE,
                 row.names = TRUE,
                 overwrite = FALSE,
                 append = TRUE) 
}

# Save file in xlsx format
do::write_xlsx(source.data$MERSCOPE,
               sheet = names(x = source.data),
               file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_MERSCOPE_updated_LF.xlsx",
               col.names = TRUE,
               row.names = TRUE,
               overwrite = TRUE)
