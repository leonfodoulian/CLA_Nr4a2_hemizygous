# Filter cells based on various criteria
FilterUnusualCells <- function(
  mat,
  reads.per.cell,
  is.expr = 0,
  mito.genes.pattern = "mt-",
  nGene.thresh = 1000,
  percent.mito.thresh = 0.1,
  loess.span = 0.5,
  loess.degree = 2,
  pval.thresh = 0.05,
  use.median = TRUE,
  verbose = TRUE
) {
  
  if (any(unique(x = apply(X = mat, MARGIN = 2, FUN = class)) == "character")) {
    stop("'mat' should be in wide format, with genes as rownames, cells as colnames, and entries as numeric")
  }
  
  if (nGene.thresh < 0 | is.numeric(x = nGene.thresh) && !isInteger(x = nGene.thresh) | is.character(x = nGene.thresh)) {
    stop("'nGene.thresh' should be a positive integer")
  }
  
  if (percent.mito.thresh < 0 | percent.mito.thresh > 1) {
    stop("'percent.mito.thresh' should be comprised between 0 and 1")
  }
  
  # Calculate metrics
  nGene <- colSums(x = mat > is.expr, na.rm = TRUE) # calculate number of genes
  nUMI <- colSums(x = mat, na.rm = TRUE) # calculate number of UMI
  mito.genes <- grep(pattern = mito.genes.pattern, x = rownames(x = mat), value = TRUE)
  percent.mito <- colSums(x = mat[mito.genes, ], na.rm = TRUE) / colSums(x = mat, na.rm = TRUE) # calculate percent.mito
  
  # Filter cells based on nGene.thresh
  nGene.selected.cells <- nGene >= nGene.thresh # select all cells having at least or more than nGene.thresh genes
  nGene.filtered.cells <- nGene.selected.cells[!nGene.selected.cells]
  if (length(x = nGene.filtered.cells) != 0) {
    message.verbose(x = paste("Cells ", PasteNames(n = names(x = nGene.filtered.cells)), 
                              " have less than ", nGene.thresh, " expressed genes: removing those cells from downstream analysis\n", sep = ""))
  } else if (length(x = nGene.filtered.cells) == 0) {
    message.verbose(x = paste("All cells have at least ", nGene.thresh, " expressed genes: not removing any cell\n", sep = ""))
  }
  
  # Filter cells based on percent.mito.thresh
  percent.mito.selected.cells <- percent.mito <= percent.mito.thresh # select all cells having at most or less than percent.mito.thresh mitochondrial abundance
  percent.mito.filtered.cells <- percent.mito.selected.cells[!percent.mito.selected.cells]
  if (length(x = percent.mito.filtered.cells) != 0) {
    message.verbose(x = paste("Cells ", PasteNames(n = names(x = percent.mito.filtered.cells)), 
                              " have more than ", percent.mito.thresh * 100, "% of mitochondrial gene abundance: removing those cells from downstream analysis\n", sep = ""))
  } else if (length(x = percent.mito.filtered.cells) == 0) {
    message.verbose(x = paste("All cells have less than ", percent.mito.thresh * 100, "% of mitochondrial gene abundance: not removing any cell\n", sep = ""))
  }
  
  # Select cells based on nGene and percent.mito
  selected.cells <- intersect(x = names(x = nGene.selected.cells[nGene.selected.cells]), 
                              y = names(x = percent.mito.selected.cells[percent.mito.selected.cells]))
  
  if (length(x = selected.cells) == 0) {
    stop("All cells are removed for downstream filtering: function cannot proceed")
  }
  
  # Compute log10 of metrics after filtering
  nGene.log10 <- log10(x = nGene[selected.cells]) # get number of genes in log10 scale
  nUMI.log10 <- log10(x = nUMI[selected.cells]) # get number of UMI in log10 scale
  percent.mito.log10 <- log10(x = percent.mito[selected.cells]) # get percent.mito in log10 scale
  
  if (!all(names(x = nGene.log10) == names(x = nUMI.log10))) {
    cat("Order of nGene.log10 and nUMI.log10 don't match: indexing nGene.log10 as in nUMI.log10", sep = "\n")
    nGene.log10 <- nGene.log10[names(x = nUMI.log10)]
  }
  
  # Fit a loess curve to remove cells having an unusual number of genes given their number of UMI
  loess.nGene.given.nUMI.log10 <- loess(formula = nGene.log10 ~ nUMI.log10, span = loess.span, degree = loess.degree)
  loess.nGene.given.nUMI.resid.log10 <- loess.nGene.given.nUMI.log10$residuals
  
  if (!missing(x = reads.per.cell)) {
    
    if (dim(x = reads.per.cell)[2] != 2 | any(!colnames(x = reads.per.cell) %in% c("Cell", "nRead"))) {
      stop("'reads.per.cell' should be a 2 column data.frame, with Cell and nRead as colnames")
    }
    
    message.verbose(x = "'reads.per.cell' is not missing: the number of reads per cell will be used for cell filtering!\n")
    
    nRead.log10 <- log10(x = reads.per.cell$nRead) # get number of reads per cell in log10 scale
    names(x = nRead.log10) <- reads.per.cell$Cell # name reads.per.cells.log10 by cell name
    nRead.log10 <- nRead.log10[names(x = nUMI.log10)] # keep only cells in mat after filtering based on nGene and percent.mito
    
    if (!all(names(x = nRead.log10) == names(x = nUMI.log10))) {
      cat("Order of nRead.log10 and nUMI.log10 don't match: indexing nRead.log10 as in nUMI.log10", sep = "\n")
      nRead.log10 <- nRead.log10[names(x = nUMI.log10)]
    }
    
    # Fit a loess curve to remove cells having an unusual number of UMI given their number of reads
    loess.nUMI.given.nRead.log10 <- loess(formula = nUMI.log10 ~ nRead.log10, span = loess.span, degree = loess.degree)
    loess.nUMI.given.nRead.resid.log10 <- loess.nUMI.given.nRead.log10$residuals
  } else {
    message.verbose(x = "'reads.per.cell' is missing:  the number of reads per cell will not be used for cell filtering!\n")
    nRead.log10 <- NULL
    loess.nUMI.given.nRead.log10 <- NULL
    loess.nUMI.given.nRead.resid.log10 <- NULL
  }
  
  # Get list of metrics
  metrics.list <- list(nGene.log10 = nGene.log10,
                       nUMI.log10 = nUMI.log10,
                       percent.mito.log10 = percent.mito.log10,
                       loess.nGene.given.nUMI.resid.log10 = loess.nGene.given.nUMI.resid.log10,
                       nRead.log10 = nRead.log10,
                       loess.nUMI.given.nRead.resid.log10 = loess.nUMI.given.nRead.resid.log10)
  
  # Get list of outlier cells
  outlier.list <- mapply(m = metrics.list, 
                         n = names(x = metrics.list),
                         pval.thresh = pval.thresh,
                         use.median = use.median,
                         FUN = DetectOutlier, 
                         SIMPLIFY = FALSE, 
                         USE.NAMES = TRUE)
  
  # Get names of outlier cells
  outlier.cell.names <- sort(x = unique(x = c(names(x = unlist(x = unname(obj = lapply(X = outlier.list, FUN = function(l) { return(l$outliers) })))),
                                              names(x = nGene.filtered.cells),
                                              names(x = percent.mito.filtered.cells))))
  
  # Get lower and upper limits for each metric
  metrics.limits <- lapply(X = outlier.list,
                           FUN = function(l) {
                             return(l$limits)
                           })
  
  # Plot metrics distribution
  metrics.plots <- mapply(m = metrics.list,
                          n = names(x = metrics.list),
                          lim = metrics.limits,
                          FUN = PlotMetrics,
                          SIMPLIFY = FALSE, 
                          USE.NAMES = TRUE)
  
  null.plots <- unlist(x = lapply(X = metrics.plots,
                                  FUN = is.null))
  
  metrics.plots <- metrics.plots[!null.plots]
  
  metrics.plot <- patchwork::wrap_plots(plot = metrics.plots, ncol = ceiling(x = sqrt(x = length(x = metrics.plots))))
  
  return(list(nGene.selected.cells = nGene.selected.cells,
              percent.mito.selected.cells = percent.mito.selected.cells,
              mito.genes = mito.genes,
              outlier.list = outlier.list,
              outlier.cells = outlier.cell.names,
              metrics = list(nGene = nGene,
                             nUMI = nUMI,
                             percent.mito = percent.mito),
              metrics.log10 = list(nGene.log10 = nGene.log10,
                                   nUMI.log10 = nUMI.log10,
                                   percent.mito.log10 = percent.mito.log10,
                                   loess.nGene.given.nUMI.log10 = loess.nGene.given.nUMI.log10,
                                   nRead.log10 = nRead.log10,
                                   loess.nUMI.given.nRead.log10 = loess.nUMI.given.nRead.log10),
              metrics.log10.plot = metrics.plot))
  
}

# Detect cells with metric 3*SD away from the mean
DetectOutlierMean <- function(m, n) {
  mean.m <- mean(x = m, na.rm = TRUE) # get mean
  sd.m <- sd(x = m, na.rm = TRUE) # get standard deviation
  upper.limit <- mean.m + 3*sd.m
  lower.limit <- mean.m - 3*sd.m
  outliers.high <- which(x = m > upper.limit) # detect outliers on the upper limit
  outliers.low <- which(x = m < lower.limit) #detect outliers on the lower limit
  
  if (length(x = outliers.high) != 0) {
    message.verbose(x = paste("Removing cell(s)", PasteNames(n = names(x = outliers.high)), "with high", getMetric(n = n)))
  }
  
  if (length(x = outliers.low) != 0) {
    message.verbose(x = paste("Removing cell(s)", PasteNames(n = names(x = outliers.low)), "with low", getMetric(n = n)))
  }
  
  if (length(x = c(outliers.high, outliers.low)) == 0) {
    message.verbose(x = "Not removing any cell")
  }
  
  cat("", sep = "\n")
  
  return(list(outliers = c(outliers.high, outliers.low),
              limits = c("lower.limit" = lower.limit,
                         "upper.limit" = upper.limit)))
}

# Detect cells with metric 3*MAD away from the median
DetectOutlierMedian <- function(m, n) {
  median.m <- median(x = m, na.rm = TRUE) # get median
  mad.m <- mad(x = m, center = median.m, na.rm = TRUE) # get median absolute deviation
  upper.limit <- median.m + 3*mad.m
  lower.limit <- median.m - 3*mad.m
  outliers.high <- which(x = m > upper.limit) # detect outliers on the upper limit
  outliers.low <- which(x = m < lower.limit) # detect outliers on the lower limit
  
  if (length(x = outliers.high) != 0) {
    message.verbose(x = paste("Removing cell(s)", PasteNames(n = names(x = outliers.high)), "with high", getMetric(n = n)))
  }
  
  if (length(x = outliers.low) != 0) {
    message.verbose(x = paste("Removing cell(s)", PasteNames(n = names(x = outliers.low)), "with low", getMetric(n = n)))
  }
  
  if (length(x = c(outliers.high, outliers.low)) == 0) {
    message.verbose(x = "Not removing any cell")
  }
  
  cat("", sep = "\n")
  
  return(list(outliers = c(outliers.high, outliers.low),
              limits = c("lower.limit" = lower.limit,
                         "upper.limit" = upper.limit)))
}

# Get appropriate metric name for message output
getMetric <- function(n) {
  if (n == "nGene.log10") { return("number of genes") }
  if (n == "nUMI.log10") { return("number of UMI") }
  if (n == "percent.mito.log10") { return("mitochondrial gene abundances") }
  if (n == "loess.nGene.given.nUMI.resid.log10") { return("number of genes given their number of UMI") }
  if (n == "nRead.log10") { return("number of reads") }
  if (n == "loess.nUMI.given.nRead.resid.log10") { return("number of UMI given their number of reads") }
  return(n)
}

# Generic function for outlier cell detection (see DetectOutlierMean and DetectOutlierMedian)
DetectOutlier <- function(m, n, pval.thresh, use.median) {
  if (use.median) {
    message.verbose(x = paste(n, ": function forced for using median(x = x) and mad(x = x, center = median(x)) for outliers detection", sep = ""))
    outliers <- DetectOutlierMedian(m = m, n = n)
  } else {
    if (shapiro.test(x = m)$p.val < pval.thresh) {
      message.verbose(x = paste("Shapiro-Wilk Normality Test p-value = ", signif(x = shapiro.test(x = m)$p.val, digits = 3), " (<", pval.thresh, "):\n", n, " is not normally distributed: using median(x = x) and mad(x = x, center = median(x)) for outliers detection", sep = ""))
      outliers <- DetectOutlierMedian(m = m, n = n) # for non normally distributed data
    } else if (shapiro.test(x = m)$p.val >= pval.thresh) {
      message.verbose(x = paste("Shapiro-Wilk Normality Test p-value = ", signif(x = shapiro.test(x = m)$p.val, digits = 3), " (>=", pval.thresh, "):\ncannot reject that ", n, " is normally distributed: using mean(x = x) and sd(x = x) for outliers detection", sep = ""))
      outliers <- DetectOutlierMean(m = m, n = n) # for normally distributed data
    }
  }
  return(outliers)
}

# Function to paste names
PasteNames <- function(n) {
  if(length(x = n) > 1) {
    m <- paste(paste(n[-length(x = n)], collapse = ", "),
               n[length(x = n)], sep = " and ")
  } else {
    m <- n
  }
  return(m)
}

# Plot distribution of metric and set vertical lines indicating the filtered values
PlotMetrics <- function(m, n, lim) {
  if (!is.null(x = m)) {
    n <- gsub(pattern = "\\.", replacement = " ", x = n) # remove points from metric name
    max.density <- ceiling(x = max(density(x = m)$y) * 2) / 2 # round to the closest 1 or 0.5 
    lim <- sort(x = lim, decreasing = FALSE) # sort lim by decreasing order
    x.lim.min <- ifelse(test = lim[1] > min(m),
                        yes = min(m),
                        no = lim[[1]]) # get lower limit of x axis values
    x.lim.max <- ifelse(test = lim[[2]] < max(m),
                        yes = max(m),
                        no = lim[[2]]) # get upper limit of x axis values
    plot <- ggplot2::ggplot(data = as.data.frame(x = m, 
                                                 nm = "metric"), 
                            mapping = ggplot2::aes_string(x = "metric",
                                                          y = "..density..")) +
      ggplot2::geom_density(fill = "dodgerblue1", alpha = 0.5) +
      ggplot2::labs(x = "",
                    y = "Density",
                    title = n) +
      ggplot2::scale_x_continuous(breaks = lim,
                                  labels = round(x = lim, digits = 2)) + 
      ggplot2::scale_y_continuous(expand = c(0,0),
                                  breaks = c(seq(from = 0, 
                                                 to = max.density, 
                                                 by = 0.5))) +
      ggplot2::expand_limits(x = c(x.lim.min, x.lim.max),
                             y = c(0, max.density)) +
      ggplot2::geom_vline(xintercept = lim,
                          color = "red",
                          size = ggplot2::rel(0.8)) 
    return(plot)
  } else {
    NULL
  }
}

# 'message' if 'verbose' is 'TRUE'
message.verbose <- function(x, verbose = TRUE) {
  if (isTRUE(x = verbose)) {
    message(x)
  }
}

# Check if value is integer (same as isTRUE)
isInteger <- function(x){
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if(test == TRUE){ return(TRUE) }
  else { return(FALSE) }
}
