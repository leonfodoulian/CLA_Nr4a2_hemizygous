# Find markers using the SCT assay from Seurat v5
FindMarkersSCTV5 <- function(
    obj,
    ident.1,
    ident.2,
    deg.ident,
    deg.cond = list(),
    ...
){
  # Permute variable
  if (isTRUE(x = deg.cond$permute.ident)) {
    # Set seed to NULL as Seurat::FindMarkers() sets a seed too which is problematic for the shuffling
    set.seed(seed = NULL)
    
    # Shuffle deg.ident identity of cells
    obj@meta.data[[deg.ident]] <- sample(x = obj@meta.data[[deg.ident]],
                                         size = nrow(x = obj@meta.data),
                                         replace = FALSE,
                                         prob = NULL)
  }
  
  # Remove sample from Seurat object
  if (isTRUE(x = deg.cond$remove.sample) && !is.na(x = deg.cond$sample.to.remove) && !is.na(x = deg.cond$sample.col)) {
    # Sample to remove
    sample.to.remove <- deg.cond$sample.to.remove[[1]]
    
    # Sample column name
    sample.col <- deg.cond$sample.col[[1]]
    
    # Get names of cells to keep
    cells.to.keep <- rownames(x = obj@meta.data[obj@meta.data[[sample.col]] != sample.to.remove,])
    
    # Subset Seurat object
    obj <- subset(x = obj,
                  cells = cells.to.keep)
  }
  
  # Set deg.ident as cell identities
  obj <- SetIdent(object = obj,
                  value = deg.ident)
  
  # Compute DEG markers
  deg.markers <- tryCatch(expr = FindMarkers(object = obj,
                                             ident.1 = ident.1,
                                             ident.2 = ident.2,
                                             assay = "SCT",
                                             ...),
                          error = function(cond) {
                            cat("No modulated genes found\n")
                            return(NULL)
                          })
  
  if (!is.null(x = deg.markers)) {
    # Get fold-change column name
    fc.col <- grep(pattern = "avg_",
                   x = colnames(x = deg.markers),
                   value = TRUE)
    
    # Clean list of markers and keep only genes with a significant p adjusted value
    deg.markers <- deg.markers[deg.markers$p_val_adj <= 0.05,,drop = FALSE]
    
    # Sort list of markers by decreasing fold change value
    deg.markers <- deg.markers[order(abs(x = deg.markers[[fc.col]]), decreasing = TRUE),]
  }
  
  # Return data
  return(deg.markers)
}
