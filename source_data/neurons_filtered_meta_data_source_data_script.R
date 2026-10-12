# seurat_nr4a2_hemi_scRNAseq-e202106_e202107_10X_ghost_GTFv102_neurons_filtered_v5.1.0_figures.R
# seurat_nr4a2_hemi_scRNAseq-e202106_e202107_10X_ghost_GTFv102_neurons_filtered_v5.1.0_linear_svc_figures.R

# Add real genotypes to cell type proportions data
cell.type.prop.l$sample.prop$genotype <- names(x = genotype.abbreviations[match(x = cell.type.prop.l$sample.prop$group, table = genotype.abbreviations)])

# Prepare source data
source.data <- list(
  "scRNAseq_neurons_filtered" = data.frame(
    cbind(
      setNames(
        object = obj@meta.data[,c("mouse", "gender", "genotype", "age", "nFeature_RNA", "nCount_RNA", "percent_mt", "broad_cell_type_abbr", "cell_type_abbr")],
        nm = c("mouse", "gender", "genotype", "age", "n_gene", "n_umi", "percent_mt", "cell_type", "cell_subtype")
      ),
      t(
        x = as.matrix(
          x = LayerData(
            object = obj,
            layer = "data",
            assay = "RNA_full",
            features = c("Cacna1a", "Cacna1b", "Cacna1c", "Cacna1d", "Cacna1e",
                         "Cacna1f", "Cacna1g", "Cacna1h", "Cacna1i", "Cacna1s",
                         "Car3", "Cckbr", "Ccn2", "Cdh13", "Cnr1",
                         "Col6a1", "Cpeb1", "Crym", "Cux1", "Dgkb",
                         "Dgkg", "Dgkz", "Drd1", "Drd2", "Egr2",
                         "Fhl2", "Fosl2", "Gabbr2", "Gabra3", "Gabrb1",
                         "Gfra1", "Gnal", "Gnao1", "Gnb4", "Gng2",
                         "Gpr26", "Gria3", "Grid2", "Grin2b", "Htr1f",
                         "Htr2c", "Itpr1", "Kcna1", "Kcna10", "Kcna2",
                         "Kcna3", "Kcna4", "Kcna5", "Kcna6", "Kcnc1",
                         "Kcnc2", "Kcnc3", "Kcnc4", "Kcnd1", "Kcnd2",
                         "Kcnd3", "Kcnf1", "Kcnj10", "Kcnj11", "Kcnj12",
                         "Kcnj13", "Kcnj14", "Kcnj16", "Kcnj2", "Kcnj3",
                         "Kcnj4", "Kcnj5", "Kcnj6", "Kcnj8", "Kcnj9",
                         "Kcnma1", "Kcnmb1", "Kcnmb2", "Kcnmb3", "Kcnmb4",
                         "Kcnn1", "Kcnn2", "Kcnn3", "Kcnq1", "Kcnq2",
                         "Kcnv1", "Kctd16", "Lxn", "Nfib", "Nnat",
                         "Nos1ap", "Npffr1", "Npsr1", "Npy", "Nr2f2",
                         "Nr4a2", "Nt5dc3", "Ntm", "Oprk1", "Rgs4",
                         "Rorb", "Rprm", "Rxfp1", "Ryr2", "Scn10a",
                         "Scn1a", "Scn1b", "Scn2a", "Scn3a", "Scn4b",
                         "Scn5a", "Scn7a", "Scn8a", "Scn9a", "Sdk2",
                         "Slc17a6", "Slc17a7", "Slc30a3", "Slc32a1", "Smim32",
                         "Sst", "Synpr", "Syt17", "Tafa1", "Trpc4",
                         "Tshz1", "Tshz2")
          )
        )
      ),
      obj@reductions$umap_harmony@cell.embeddings
    )
  ),
  "Fig4d" = setNames(
    object = CountDEG(
      deg.list = hemi.observed,
      fc.col = "avg_log2FC",
      id.col = "cell_type_abbr"
    ),
    nm = c("cell_type", "n_downregulated", "n_upregulated")
  ),
  "Fig6c" = setNames(
    object = cell.type.prop.l$sample.prop[,c("sample", "genotype", "clusters", "log2_prop_ratio")],
    nm = c("mouse", "genotype", "cell_type", "log2_prop_ratio")
  ),
  "Fig7a" = accuracy.data.l$observed[stringr::str_order(x = accuracy.data.l$observed$iteration, numeric = TRUE),-1],
  "FigSupp11" = accuracy.data.l$permuted[stringr::str_order(x = accuracy.data.l$permuted$iteration, numeric = TRUE),-1]
)

# Save file in xlsx format
for (n in names(x = source.data)) {
  if (n == "scRNAseq_neurons_filtered") {
    row.names <- TRUE
  } else {
    row.names <- FALSE
  }
  do::write_xlsx(source.data[[n]],
                 sheet = n,
                 file = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/source_data/Nr4a2_hemi_paper_source_data_LF.xlsx",
                 col.names = TRUE,
                 row.names = row.names,
                 overwrite = FALSE,
                 append = TRUE) 
}
