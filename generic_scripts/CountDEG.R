# Count number of upregulated and downregulated genes
CountDEG <- function(
    deg.list,
    fc.col,
    id.col
){
  # Count number of modulated genes
  n.modulated.genes <- dplyr::bind_rows(
    lapply(
      X = deg.list,
      FUN = function(deg,
                     fc.col) {
        # Compute number of modulated genes
        n_downregulated <- sum(deg[[fc.col]] < 0)
        n_upregulated <- sum(deg[[fc.col]] > 0)
        
        # Return data
        return(data.frame(downregulated = n_downregulated,
                          upregulated = n_upregulated))
      },
      fc.col = fc.col
    ),
    .id = id.col
  )
  
  # Return data
  return(n.modulated.genes)
}
