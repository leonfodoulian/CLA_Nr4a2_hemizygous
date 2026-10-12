# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/smFISH_raw_data/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/rnascope_helper_functions.R")

# Create new directories
dir.create(path = "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH")

# List experiment directories on which function will be executed
dirs <- grep(
  pattern = "Fig4g",
  x = list.dirs(path = ".",
                full.names = FALSE,
                recursive = FALSE),
  value = TRUE
)
names(x = dirs) <- dirs

# Define attributes to extract from file names
attributes <- c("slide",
                "zone",
                "mouse",
                "genotype",
                "level",
                "Channel 2",
                "Channel 3",
                "Channel 4")

# Load and preprocess RNAscope data
mapply(
  dir = dirs,
  file.prefix = names(x = dirs),
  FUN = function(dir,
                 file.prefix,
                 attributes) {
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
          suffix = ".xlsx$|_cellular.xlsx$|_nuclear.xlsx$",
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
        
        # Return RNAscope data
        return(rnascope.data)
      },
      attributes = attributes
    )
    
    if (all(grepl(pattern = "_cellular.xlsx$|_nuclear.xlsx$", x = names(x = rnascope.data.l)))) {
      # Split RNAscope data by count type
      rnascope.data.l <- split(x = rnascope.data.l,
                               f = gsub(pattern = ".*_|.xlsx$",
                                        replacement = "",
                                        x = names(x = rnascope.data.l)))
      
      # Bind RNAscope data by row
      rnascope.data.l <- lapply(X = rnascope.data.l,
                                FUN = dplyr::bind_rows,
                                .id = NULL)
      
      # Substract nuclear Nr2f2 counts from cellular Nr2f2 counts
      if(nrow(x = rnascope.data.l$cellular) != nrow(x = rnascope.data.l$nuclear)) {
        stop("cellular and nuclear data do not match")
      }
      rnascope.data.l$cellular$Nr2f2 <- rnascope.data.l$cellular$Nr2f2 - rnascope.data.l$nuclear$Nr2f2
      
      # Set negative Nr2f2 counts to 0
      negative.counts <- rnascope.data.l$cellular$Nr2f2 < 0
      if (sum(negative.counts) > 0) {
        message("negative counts in data (n = ", sum(negative.counts), ")")
        rnascope.data.l$cellular$Nr2f2[negative.counts] <- 0 
      }
      
      # Save RNAscope data to submission folder
      write.csv(x = rnascope.data.l$cellular,
                file = file.path(
                  "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH",
                  paste0(file.prefix,
                         "_raw_data.csv")
                ),
                quote = FALSE,
                row.names = FALSE)
    } else {
      # Bind RNAscope data by row
      rnascope.data <- dplyr::bind_rows(rnascope.data.l,
                                        .id = NULL)
      
      # Save RNAscope data to submission folder
      write.csv(x = rnascope.data,
                file = file.path(
                  "/Users/leonfodoulian/scData/Nr4a2_hemi_paper/raw_data/smFISH",
                  paste0(file.prefix,
                         "_raw_data.csv")
                ),
                quote = FALSE,
                row.names = FALSE)
    }
    
    # Return "DONE"
    return("DONE")
  },
  MoreArgs = list(attributes = attributes),
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE
)
