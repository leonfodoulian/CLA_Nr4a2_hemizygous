# rnascope_subcluster_spatial_maps.R

# Add real genotype to data
rnascope.data.l <- lapply(
  X = rnascope.data.l,
  FUN = function(rnascope.data,
                 cols.to.keep) {
    rnascope.data$genotype <- "C57BL6"
    rownames(x = rnascope.data) <- NULL
    return(rnascope.data[,c(cols.to.keep, colnames(x = rnascope.data)[7:9])])
  },
  cols.to.keep = c("cell", "slide", "zone", "mouse", "genotype", "level", "centroid_x", "centroid_y")
)

# Assign figure panels to data
figure.panels <- c(
  "Fig3d" = "e20241217_QuPath_C2.Nr4a2_C3.Nnat_C4.Rprm_log_transformed\nSlide1_zone3_B6_wt.wt_level44\nk7_clusters_hclust_wardD2_cluster7_spatial_maps_uniform_scale.pdf",
  "Fig3g" = "e20241217_QuPath_C2.Nr4a2_C3.Nnat_C4.Rprm_log_transformed\nSlide1_zone3_B6_wt.wt_level44\nk7_clusters_hclust_wardD2_cluster5_spatial_maps_uniform_scale.pdf",
  "Fig3h" = "e20241217_QuPath_C2.Nr4a2_C3.Nnat_C4.Rprm_log_transformed\nSlide1_zone3_B6_wt.wt_level44\nk7_clusters_hclust_wardD2_cluster2_spatial_maps_uniform_scale.pdf"
)

# Prepare source data
source.data <- rnascope.data.l[figure.panels]
names(x = source.data) <- names(x = figure.panels)

# Save file in xlsx format
i <- 1
for (n in names(x = source.data)) {
  do::write_xlsx(source.data[[n]],
                 sheet = n,
                 file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_Fig3_updated_LF.xlsx",
                 col.names = TRUE,
                 row.names = FALSE,
                 overwrite = i == 1,
                 append = i != 1)
  i <- i + 1
}
