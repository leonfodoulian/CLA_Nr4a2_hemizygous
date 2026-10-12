# seurat_nr4a2_hemi_scRNAseq-e202106_e202107_10X_ghost_GTFv102_all_cells_v5.1.0_figures.R

# Prepare source data
source.data <- list(
  "scRNAseq_all_cells" = data.frame(
    cbind(
      setNames(
        object = obj@meta.data[,c("mouse", "gender", "genotype", "age", "nFeature_RNA", "nCount_RNA", "percent_mt", "broad_cell_type_abbr")],
        nm = c("mouse", "gender", "genotype", "age", "n_gene", "n_umi", "percent_mt", "cell_type")
      ),
      t(
        x = as.matrix(
          x = LayerData(
            object = obj,
            layer = "data",
            assay = "RNA_full",
            features = c("C1qb", "Dcn", "Flt1", "Gja1", "Gnb4",
                         "Mobp", "Nr4a2", "Pdgfra", "Rgs5", "Rprm",
                         "Scn4b", "Slc17a7", "Slc32a1", "Syt17")
          )
        )
      ),
      obj@reductions$umap_harmony@cell.embeddings
    )
  )
)

# Save file in xlsx format
do::write_xlsx(source.data$scRNAseq_all_cells,
               sheet = names(x = source.data),
               file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_LF.xlsx",
               col.names = TRUE,
               row.names = TRUE,
               overwrite = TRUE)
