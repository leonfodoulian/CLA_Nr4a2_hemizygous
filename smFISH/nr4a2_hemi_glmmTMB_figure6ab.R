# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Load required packages
require(glmmTMB)

# Load data
cell.numbers <- readxl::read_xlsx(path = "files_to_submit/Figure_5AB/e20230202-e20230605_nb_Nr4a2_cells_ratios.xlsx",
                                  na = "NaN")

# Set WT genotype as the reference
cell.numbers$genotype <- factor(x = cell.numbers$genotype,
                                levels = c("Nr4a2(wt/wt)", "Nr4a2(del/wt)"))

# Compute statistics on the data
lrt.data.l <- lapply(
  X = list("total.cells" = "count.total",
           "nr4a2.cells" = "count.Nr4a2_pos"),
  FUN = function(count.col,
                 cell.numbers) {
    # Run glmmTMB
    model <- glmmTMB::glmmTMB(formula = as.formula(object = paste0(count.col, " ~ genotype + (1 | probeset)")),
                              data = cell.numbers,
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
  cell.numbers = cell.numbers
)
