# Prepare pySCENIC regulon data frame
SCENICRegulonDataFrame <- function(
    regulon,
    genes = NULL
){
  # Prepare regulon data frame 
  regulon.data <- data.frame(
    regulon_name = regulon$name,
    regulon_context = ifelse(test = grepl(pattern = "\\+",
                                          x = regulon$name),
                             yes = "activating",
                             no = "repressing"),
    transcription_factor = regulon$transcription_factor,
    gene = names(x = regulon$gene_weights),
    weight = regulon$gene_weights,
    norm_weight = regulon$gene_weights / sum(regulon$gene_weights),
    row.names = NULL
  )
  
  # Add to data genes not present in regulon list of genes
  if (!is.null(x = genes)) {
    # Get list of genes not in regulon data
    genes <- setdiff(x = genes,
                     y = regulon.data$gene)
    
    # Create data frame with remainder of genes
    if (length(x = genes) > 0) {
      regulon.data <- rbind(
        regulon.data,
        data.frame(
          regulon_name = regulon$name,
          regulon_context = NA,
          transcription_factor = regulon$transcription_factor,
          gene = genes,
          weight = NA,
          norm_weight = NA,
          row.names = NULL
        )
      )
    }
  }
  
  # Return data
  return(regulon.data)
}
