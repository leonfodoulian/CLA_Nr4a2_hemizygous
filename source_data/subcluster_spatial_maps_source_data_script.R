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
  "Fig3b" = "e20230412_QuPath_C2.Lxn_C3.Oprk1_C4.Fosl2_log_transformed\nSlide2_zone4_B6_wt.wt_level44\nk4_clusters_hclust_wardD2_cluster2_spatial_maps_uniform_scale.pdf",
  "Fig3c" = "e20230412_QuPath_C2.Lxn_C3.Oprk1_C4.Fosl2_log_transformed\nSlide2_zone4_B6_wt.wt_level44\nk4_clusters_hclust_wardD2_cluster4_spatial_maps_uniform_scale.pdf",
  "Fig3d" = "e20220929_QuPath_C2.Npsr1_C3.Rorb_C4.Nr4a2_log_transformed\nSlide3_zone1_B6_wt.wt_level44\nk8_clusters_hclust_wardD2_cluster8_spatial_maps_uniform_scale.pdf",
  "Fig3e" = "e20220929_QuPath_C2.Npy_C3.Nr4a2_C4.Gad1_log_transformed\nSlide2_zone1_B6_wt.wt_level44\nk12_clusters_hclust_wardD2_cluster12_spatial_maps_uniform_scale.pdf",
  "Fig3f" = "e20220929_QuPath_C2.Egr2_C3.Slc17a6_C4.Nnat_log_transformed\nSlide1_zone3_B6_wt.wt_level44\nk11_clusters_hclust_wardD2_cluster7_spatial_maps_uniform_scale.pdf",
  "Fig3g" = "e20220913_QuPath_C2.Nr4a2_C3.Rprm_C4.Col6a1_log_transformed\nSlide4_zone1_B6_wt.wt_level44\nk6_clusters_hclust_wardD2_cluster3_spatial_maps_uniform_scale.pdf",
  "Fig3h" = "e20220913_QuPath_C2.Nr4a2_C3.Rprm_C4.Col6a1_log_transformed\nSlide4_zone1_B6_wt.wt_level44\nk6_clusters_hclust_wardD2_cluster1_spatial_maps_uniform_scale.pdf",
  "Fig3i" = "e20220929_QuPath_C2.Npsr1_C3.Rorb_C4.Nr4a2_log_transformed\nSlide3_zone1_B6_wt.wt_level44\nk8_clusters_hclust_wardD2_cluster4_spatial_maps_uniform_scale.pdf",
  "FigSupp6b_right" = "e20220913_QuPath_C2.Nr4a2_C3.Rprm_C4.Col6a1_log_transformed\nSlide4_zone1_B6_wt.wt_level44\nk6_clusters_hclust_wardD2_cluster4_spatial_maps_uniform_scale.pdf",
  "FigSupp6b_left" = "e20220913_QuPath_C2.Nr4a2_C3.Rprm_C4.Col6a1_log_transformed\nSlide4_zone3_B6_wt.wt_level44\nk5_clusters_hclust_wardD2_cluster3_spatial_maps_uniform_scale.pdf"
)

# Prepare source data
source.data <- rnascope.data.l[figure.panels]
names(x = source.data) <- names(x = figure.panels)

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
