# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/RNAscope_data_final/results/")

# Create new directories
dir.create(path = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH")

# Define raw data to submit
raw.data.l <- list(
  "Fig2_Nr4a2_Rprm_Ctgf" = "e20230425_QuPath/C2.Nr4a2_C3.Rprm_C4.Ctgf/all_cells",
  "Fig2_Nr4a2_Syt17_Slc17a6" = "e20230417_QuPath/C2.Nr4a2_C3.Syt17_C4.Vglut2/all_cells",
  "Fig3b_Fig3c_Lxn_Oprk1_Fosl2" = "e20230412_QuPath/C2.Lxn_C3.Oprk1_C4.Fosl2/all_cells",
  "Fig3d_Fig3g_Fig3h_Nr4a2_Nnat_Rprm" = "e20241217_QuPath/C2.Nr4a2_C3.Nnat_C4.Rprm/all_cells",
  "Fig3e_Nr4a2_Npy" = "e20220929_QuPath/C2.Npy_C3.Nr4a2_C4.Gad1/all_cells",
  "Fig3f_Egr2_Slc17a6_Nnat" = "e20220929_QuPath/C2.Egr2_C3.Slc17a6_C4.Nnat/all_cells",
  "Fig3i_Npsr1_Rorb_Nr4a2" = "e20220929_QuPath/C2.Npsr1_C3.Rorb_C4.Nr4a2/all_cells",
  "Fig4hi_Fig6fg_Nr4a2_Ntm_Syt17" = "e20230605_QuPath/C2.Nr4a2_C3.Ntm_C4.Syt17/all_cells",
  "Fig4hj_Fig6fg_Nr4a2_Ryr2_Syt17" = "e20230605_QuPath/C2.Nr4a2_C3.Ryr2_C4.Syt17/all_cells"
)

# Save data
mapply(
  raw.data = raw.data.l,
  file.prefix = names(x = raw.data.l),
  FUN = function(raw.data,
                 file.prefix) {
    # Define RNAscope rds file
    rds.file <- list.files(path = raw.data,
                           pattern = ".rds",
                           full.names = TRUE)
    
    # Load RNAscope data
    rnascope.data <- readRDS(file = rds.file)$full.data
    
    # Subset RNAscope data for columns of interest
    columns.to.remove <- grepl(pattern = "_scaled|_cluster",
                               x = colnames(x = rnascope.data))
    rnascope.data <- rnascope.data[,!columns.to.remove]
    
    # Replace Vglut2 with Slc17a6
    if (any(grepl(pattern = "Vglut2", x = colnames(x = rnascope.data)))) {
      colnames(x = rnascope.data) <- gsub(pattern = "Vglut2",
                                          replacement = "Slc17a6",
                                          colnames(x = rnascope.data))
    }
    
    # Save RNAscope data to submission folder
    write.csv(x = rnascope.data,
              file = file.path(
                "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH",
                paste0(file.prefix,
                       "_raw_data.csv")
              ),
              quote = FALSE,
              row.names = FALSE)
    
    # Return "DONE"
    return("DONE")
  },
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)
