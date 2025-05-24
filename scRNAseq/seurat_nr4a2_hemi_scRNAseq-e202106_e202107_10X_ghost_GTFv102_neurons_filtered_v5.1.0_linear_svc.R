# Set working directory
setwd(dir = "/Users/leonfodoulian/scData/nr4a2_hemi_scRNAseq-e202106_e202107/")

# Functions to load from disk
source(file = "/Users/leonfodoulian/scData/PartitionData.R")

# Required packages
require(Seurat)
require(reticulate)

# Python packages
svm <- reticulate::import(module = "sklearn.svm")
ms <- reticulate::import(module = "sklearn.model_selection")

# Define rds file name prefix
file.prefix <- "scRNAseq-e20210604_e20210611_scRNAseq-e20210719_e20210720_smim32_GTFv102_cellranger_M3Drop_v3.10.6_danb_seurat_v5.1.0_sctransform_harmony_v1.2.0_unsupervised_clustering_PC1-19_neurons_filtered"

# Load Seurat object
obj <- readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                paste0(file.prefix, ".rds")))

# Subset Seurat object for CLA, shell and L6a cell types
obj <- subset(x = obj,
              subset = broad_cell_type_abbr %in% c("CLA", "shell", "L6a"))

# Load list of modulated genes in CLA clusters
hemi.genes <- Reduce(
  f = union,
  x = lapply(
    X = readRDS(file = file.path("results/smim32_GTFv102/neurons_filtered",
                                 paste0(file.prefix, "_hemi_modulated_genes.rds")))[c("CLA 1", "CLA 2")],
    FUN = function(modulated.genes) {
      return(rownames(x = modulated.genes$observed))
    }
  )
)

# Get names of wt and hemizygous cells
wt.cells <- rownames(x = obj@meta.data[obj@meta.data$genotype == "Nr4a2(WT/WT)",])
hemi.cells <- rownames(x = obj@meta.data[obj@meta.data$genotype == "Nr4a2(SA-IRES-Dre/WT)",])

# Prepare counts and labels data for SVM
counts.data <- t(x = obj@assays$SCT_full@scale.data)
broad.cell.type <- setNames(object = obj@meta.data[,"broad_cell_type_abbr"],
                            nm = rownames(x = obj@meta.data))
cell.type <- setNames(object = obj@meta.data[,"cell_type_abbr"],
                      nm = rownames(x = obj@meta.data))

# Prepare data for balanced partitioning
partition.data <- data.frame(cell = wt.cells,
                             broad_cell_type = broad.cell.type[wt.cells],
                             cell_type = cell.type[wt.cells])

# Define range of regularization parameter C values for grid search
c.values <- c(as.vector(x = outer(X = 1:9,
                                  Y = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
                                  FUN = "*")),
              1)

# Define grid search parameters for best parameters estimation
grid.parameters <- list(list(C = c.values,
                             penalty = list("l1"),
                             dual = list("auto", FALSE)),
                        list(C = c.values,
                             penalty = list("l2"),
                             dual = list("auto", TRUE)))

# Set model selection fit for parameter estimation using grid search
ms.grid <- ms$GridSearchCV(
  estimator = svm$LinearSVC(
    loss = "squared_hinge",
    tol = 1e-5,
    multi_class = "ovr",
    fit_intercept = TRUE,
    intercept_scaling = 1,
    class_weight = NULL,
    verbose = 0L,
    random_state = NULL,
    max_iter = 50000L
  ),
  param_grid = grid.parameters,
  scoring = "balanced_accuracy",
  n_jobs = 7L,
  refit = TRUE,
  cv = 10L,
  verbose = 0L,
  pre_dispatch = 1L,
  error_score = "raise",
  return_train_score = FALSE
)

# Run model selection fit using grid search
ms.grid.fit <- ms.grid$fit(X = counts.data[wt.cells, hemi.genes],
                           y = broad.cell.type[wt.cells])

# Set linear SVC model
svm.model <- svm$LinearSVC(
  penalty = ms.grid.fit$best_params_$penalty, # "l2"
  loss = "squared_hinge",
  dual = ms.grid.fit$best_params_$dual, # "auto"
  tol = 1e-5,
  C = ms.grid.fit$best_params_$C, # 9e-04
  multi_class = "ovr",
  fit_intercept = TRUE,
  intercept_scaling = 1,
  class_weight = NULL,
  verbose = 0L,
  random_state = NULL,
  max_iter = 50000L
)

# Run linear SVC fit with 1000 iterations
iterations <- paste0("iteration_", 1:1000)
names(x = iterations) <- iterations
predicted.labels.l <- lapply(
  X = iterations,
  FUN = function(iteration,
                 svm.model,
                 counts.data,
                 cell.type,
                 partition.data,
                 hemi.cells) {
    # Print iteration to keep track of progress
    message("Processing ", iteration)
    
    # Partition data into 80% train, 20% test set
    partitions.l <- PartitionData(
      data = partition.data,
      p = 0.2,
      cat_col = "cell_type",
      force_equal = FALSE,
      list_out = FALSE
    )
    
    # Get list of train and test cells
    train.cells <- partitions.l$train$cell
    test.cells <- partitions.l$test$cell
    
    # Fit SVC model with train cells and predict labels
    set.seed(seed = NULL)
    predicted.labels.l <- lapply(
      X = list("observed" = unname(obj = cell.type[train.cells]),
               "permuted" = unname(obj = sample(x = cell.type[train.cells]))),
      FUN = function(y.labels,
                     svm.model,
                     counts.data,
                     train.cells,
                     test.cells,
                     hemi.cells,
                     cell.type) {
        # Fit SVC model with train cells
        svm.fit <- svm.model$fit(X = counts.data[train.cells,],
                                 y = y.labels)
        
        # Predict labels of test and hemi cells
        predicted.labels <- data.frame(
          cell = c(test.cells, hemi.cells),
          actual = cell.type[c(test.cells, hemi.cells)],
          predicted = svm.fit$predict(X = counts.data[c(test.cells, hemi.cells),]),
          genotype = rep(x = c("Nr4a2(WT/WT)",
                               "Nr4a2(SA-IRES-Dre/WT)"),
                         times = c(length(x = test.cells),
                                   length(x = hemi.cells))),
          row.names = NULL
        )
        
        # Return data
        return(predicted.labels)
      },
      svm.model = svm.model,
      counts.data = counts.data,
      train.cells = train.cells,
      test.cells = test.cells,
      hemi.cells = hemi.cells,
      cell.type = cell.type)
    
    # Return data
    return(predicted.labels.l)
  },
  svm.model = svm.model,
  counts.data = counts.data[,hemi.genes],
  cell.type = broad.cell.type,
  partition.data = partition.data,
  hemi.cells = hemi.cells
)

# Format data for downstream analysis
svm.data.l <- list(
  predictions = list(
    "observed" = dplyr::bind_rows(
      lapply(
        X = predicted.labels.l,
        FUN = "[[",
        "observed"
      ),
      .id = "iteration"
    ),
    "permuted" = dplyr::bind_rows(
      lapply(
        X = predicted.labels.l,
        FUN = "[[",
        "permuted"
      ),
      .id = "iteration"
    )
  ),
  best.parameters = ms.grid$best_params_
)

# Save data
saveRDS(object = svm.data.l,
        file = file.path("results/smim32_GTFv102/neurons_filtered",
                         paste0(file.prefix, "_linear_svc_predicted_labels.rds")))
