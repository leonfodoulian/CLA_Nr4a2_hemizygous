# rnascope_claustro_insular_maps.R

# Create a vector of new cluster names
new.cluster.names <- c(
  "CLA" = "CLA",
  "Ctgf+" = "Ccn2+",
  "noise" = "noise",
  "not clustered" = "not clustered",
  "Rprm+" = "Rprm+",
  "Shell" = "shell",
  "Syt17+" = "Syt17+",
  "Syt17+;Vglut2+" = "Syt17+;Slc17a6+",
  "Vglut2+" = "Slc17a6+"
)

# Prepare gene expression data
gene.exp.data <- dplyr::bind_rows(
  mapply(
    probeset.data = all.data,
    cols.to.add = list("C2.Nr4a2_C3.Syt17_C4.Vglut2" = c("Rprm", "Ctgf"),
                       "C2.Nr4a2_C3.Rprm_C4.Ctgf" = c("Syt17", "Vglut2")),
    FUN = function(probeset.data,
                   cols.to.add,
                   cols.to.keep) {
      for (col in cols.to.add) {
        probeset.data[[col]] <- NA
      }
      return(probeset.data[,cols.to.keep])
    },
    MoreArgs = list(cols.to.keep = c("cell", "Nr4a2", "Rprm", "Ctgf", "Syt17", "Vglut2")),
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

# Add real genotype to data
rnascope.data$real_genotype <- "C57BL6"

# Rename Ctgf and Slc17a6 clusters
rnascope.data$clustering_results_new_names <- new.cluster.names[rnascope.data$clustering_results]

# Prepare source data
source.data <- list(
  "Fig2" = setNames(
    object = rnascope.data[,c("cell", "slide", "zone", "mouse", "real_genotype", "level",
                              "centroid_x_rot", "centroid_y_rot", "centroid_x_disc", "centroid_y_disc", "image_type",
                              "Nr4a2", "Rprm", "Ctgf", "Syt17", "Vglut2", "clustering_results_new_names")],
    nm = c("cell", "slide", "zone", "mouse", "genotype", "level",
           "centroid_x", "centroid_y", "centroid_x_disc", "centroid_y_disc", "image_type",
           "Nr4a2", "Rprm", "Ccn2", "Syt17", "Slc17a6", "clustering_results")
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
