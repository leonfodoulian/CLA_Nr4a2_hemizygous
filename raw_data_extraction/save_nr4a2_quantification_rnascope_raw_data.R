# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/smFISH_raw_data/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/rnascope_helper_functions.R")

# Create new directories
dir.create(path = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH")

# Define experiment directory on which function will be executed
dir <- "Fig4c_Fig6a_Fig6b_Nr4a2"
names(x = dir) <- dir

# Define attributes to extract from file names
attributes <- c("slide",
                "zone",
                "mouse",
                "genotype",
                "level",
                "Channel 2",
                "Channel 3",
                "Channel 4")

# List files to load
files <- list.files(path = dir,
                    pattern = ".xlsx",
                    full.names = TRUE,
                    recursive = TRUE)
names(x = files) <- files

# Keep only RNAscope data
files <- files[grepl(pattern = "Slide", x = files)]

# Load RNAscope data
rnascope.data.l <- lapply(
  X = files,
  FUN = function(file,
                 attributes) {
    # Extract file attributes from file name
    file.attributes <- GetFileAttributes(
      file = file,
      attributes = attributes,
      channel.attributes = 6:8,
      split.pattern = "_",
      prefix = gsub(pattern = "/Slide.*",
                    replacement = "",
                    x = file),
      suffix = ".xlsx",
      path = NULL,
      subdirectory = NULL
    )
    
    # Get list of channels (probe set)
    channels <- gsub(pattern = ".*\\.",
                     replacement = "",
                     x = file.attributes$channels)
    
    # Load RNAscope data
    rnascope.data <- LoadRNAscopeData(
      file = file.attributes$file,
      objects.to.keep = "PathCellObject",
      meta.data = file.attributes$meta.data,
      channels = channels
    )
    
    # Add probe set information to RNAscope data
    rnascope.data$probe_set <- paste(file.attributes$channels,
                                     collapse = "_")
    
    # Subset RNAscope data for columns of interest
    columns.to.remove <- colnames(x = rnascope.data) %in% setdiff(x = channels,
                                                                  y = "Nr4a2")
    rnascope.data <- rnascope.data[,!columns.to.remove]
    
    # Return RNAscope data
    return(rnascope.data)
  },
  attributes = attributes
)

# Bind RNAscope data by row
rnascope.data <- dplyr::bind_rows(rnascope.data.l,
                                  .id = NULL)

# Save RNAscope data to submission folder
write.csv(x = rnascope.data,
          file = file.path(
            "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH",
            paste0("Fig4c_Fig6a_Fig6b_Nr4a2",
                   "_raw_data.csv")
          ),
          quote = FALSE,
          row.names = FALSE)
