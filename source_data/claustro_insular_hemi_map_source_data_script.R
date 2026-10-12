# rnascope_claustro_insular_nr4a2_hemi_contrast_maps.R

# Define abbreviations of genotypes
genotype.abbreviations <- setNames(object = c("wt", "del"),
                                   nm = c("Nr4a2(WT/WT)", "Nr4a2(SA-IRES-Dre/WT)"))

# Prepare gene expression data
gene.exp.data <- dplyr::bind_rows(
  mapply(
    probeset.data = all.data,
    cols.to.add = list("C2.Nr4a2_C3.Ntm_C4.Syt17" = "Ryr2",
                       "C2.Nr4a2_C3.Ryr2_C4.Syt17" = "Ntm"),
    FUN = function(probeset.data,
                   cols.to.add,
                   cols.to.keep) {
      for (col in cols.to.add) {
        probeset.data[[col]] <- NA
      }
      return(probeset.data[,cols.to.keep])
    },
    MoreArgs = list(cols.to.keep = c("cell", "Nr4a2", "Syt17", "Ntm", "Ryr2")),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  ),
  .id = "probe_set"
)

# Add gene expression to data
rnascope.data <- merge(
  x = rnascope.data,
  y = gene.exp.data,
  by = c("probe_set", "cell"),
  all.x = TRUE,
  sort = FALSE
)

# Add real genotypes to data
rnascope.data$real_genotype <- names(x = genotype.abbreviations[match(x = rnascope.data$genotype_abbr, table = genotype.abbreviations)])
cell.type.prop.l$sample.prop$genotype <- names(x = genotype.abbreviations[match(x = cell.type.prop.l$sample.prop$group, table = genotype.abbreviations)])

# Prepare source data
source.data <- list(
  "Fig6f" = setNames(
    object = cell.type.prop.l$sample.prop[,c("sample", "genotype", "clusters", "log2_prop_ratio")],
    nm = c("image", "genotype", "cell_type", "log2_prop_ratio")
  ),
  "Fig6g" = setNames(
    object = rnascope.data[,c("cell", "slide", "zone", "mouse", "real_genotype", "level",
                              "centroid_x_rot", "centroid_y_rot", "centroid_x_disc", "centroid_y_disc", "image_type",
                              "Nr4a2", "Syt17", "Ntm", "Ryr2", "clustering_results")],
    nm = c("cell", "slide", "zone", "mouse", "genotype", "level",
           "centroid_x", "centroid_y", "centroid_x_disc", "centroid_y_disc", "image_type",
           "Nr4a2", "Syt17", "Ntm", "Ryr2", "clustering_results")
  )
)

# Save file in xlsx format
for (n in names(x = source.data)) {
  do::write_xlsx(source.data[[n]],
                 sheet = n,
                 file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_LF.xlsx",
                 col.names = TRUE,
                 row.names = FALSE,
                 overwrite = FALSE,
                 append = TRUE) 
}
