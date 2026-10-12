# seurat_nr4a2_hemi_scRNAseq-e202106_e202107_10X_ghost_GTFv102_neurons_filtered_v5.1.0_figures.R
# rnascope_claustro_insular_nr4a2_hemi_contrast_maps.R
paste(
  unlist(
    x = lapply(
      X = split(x = cell.type.prop.l$cell.type.prop,
                f = 1:nrow(x = cell.type.prop.l$cell.type.prop)),
      FUN = function(row.data) {
        paste0(
          row.data[["BaselineProp.clusters"]], ": ",
          # "mean Nr4a2(wt/wt) = ", signif(x = row.data[["PropMean.del"]],
          #                                digits = 3),
          # ", mean Nr4a2(del/wt) = ", signif(x = row.data[["PropMean.wt"]],
          #                                   digits = 3),
          "ratio = ", signif(x = row.data[["PropRatio"]],
                               digits = 3),
          ", t-value = ", format(x = round(x = row.data[["Tstatistic"]],
                                         digits = 2),
                               nsmall = 2),
          ", p-value = ", signif(x = row.data[["FDR"]],
                                 digits = 3)
        )
      }
    )
  ),
  collapse = "; "
)
