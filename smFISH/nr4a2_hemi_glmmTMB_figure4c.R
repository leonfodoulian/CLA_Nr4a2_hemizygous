# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Load required packages
require(glmmTMB)

# Load data
nr4a2.puncta <- readxl::read_xlsx(path = "files_to_submit/Figure_4C/e20230202-e20230605_nb_Nr4a2_puncta_raw_data.xlsx",
                                  na = "NaN")

# Set WT genotype as the reference
nr4a2.puncta$genotype <- factor(x = nr4a2.puncta$genotype,
                                levels = c("Nr4a2(wt/wt)", "Nr4a2(del/wt)"))

# Compute statistics on the data
lrt.data.l <- lapply(
  X = list("nr4a2.puncta" = "Nr4a2"),
  FUN = function(count.col,
                 nr4a2.puncta) {
    # Round counts to convert data to integers (for model building)
    nr4a2.puncta[[count.col]] <- round(x = nr4a2.puncta[[count.col]])
    
    # Run glmmTMB
    model <- glmmTMB::glmmTMB(formula = as.formula(object = paste0(count.col, " ~ genotype + (1 | probeset) + (1 | Image)")),
                              data = nr4a2.puncta,
                              family = glmmTMB::nbinom2(link = "log"))
    
    # Compute p values using a LRT
    lrt <- drop1(object = model,
                 test = "Chisq")
    
    # Transform to data frame
    lrt.data <- data.frame(lrt["genotype",],
                           row.names = NULL,
                           check.names = FALSE)
    
    # Return list of data
    return(list(model = model,
                lrt = lrt,
                lrt.data = lrt.data))
  },
  nr4a2.puncta = nr4a2.puncta
)
