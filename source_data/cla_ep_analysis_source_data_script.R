# seurat_nr4a2_hemi_scRNAseq-e202106_e202107_10X_ghost_GTFv102_neurons_filtered_v5.1.0_cla_ep_analysis_figures.R

# Subset Seurat object for regions of interest
obj <- subset(x = obj.wt,
              subset = region %in% c("CLA", "EP"))

# Add replicate to cell type proportions data
cla.cell.type.prop$replicate <- gsub(pattern = ".*\\.",
                                     replacement = "",
                                     x = cla.cell.type.prop$cla_side)

# Prepare source data
source.data <- list(
  "Fig9b" = data.frame(
    cbind(
      setNames(
        object = obj@meta.data[,c("region"),drop = FALSE],
        nm = c("predicted_localization")
      ),
      t(
        x = as.matrix(
          x = LayerData(
            object = obj,
            layer = "data",
            assay = "RNA_full",
            features = sort(x = rownames(x = cla.ep.analysis$markers$`CLA vs EP`))
          )
        )
      )
    )
  ),
  "Fig9c" = setNames(
    object = cla.cell.type.prop[,c("mouse", "replicate", "region", "cla_cell_type", "proportion")],
    nm = c("mouse", "replicate", "localization", "predicted_cell_subtype", "proportion")
  )
)

# Save file in xlsx format
i <- 1
for (n in names(x = source.data)) {
  do::write_xlsx(source.data[[n]],
                 sheet = n,
                 file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_Fig9_LF.xlsx",
                 col.names = TRUE,
                 row.names = n == "Fig9b",
                 overwrite = i == 1,
                 append = i != 1)
  i <- i + 1
}
