# Partition data into train and test subsets
PartitionData <- function(
    data,
    ...
){
  # Set seed to NULL
  set.seed(seed = NULL)
  
  # Create balanced partitions of the data
  partitions <- groupdata2::partition(
    data = partition.data,
    ...
  )
  
  # Prepare data to return
  partitions <- as.data.frame(x = partitions)
  partitions$partition_type <- c("1" = "test", "2" = "train")[as.character(x = partitions$.partitions)]
  partitions.l <- split(x = partitions,
                        f = partitions$partition_type)
  
  # Return a list of data
  return(partitions.l)
}
