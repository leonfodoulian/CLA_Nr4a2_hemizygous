# nr4a2_hemi_glmmTMB_figure4c.R
# nr4a2_hemi_glmmTMB_figure4gh.R
# nr4a2_hemi_glmmTMB_figure4ij.R
# nr4a2_hemi_glmmTMB_figure5ab.R
# seurat_nr4a2_hemi_scRNAseq-e202106_e202107_10X_ghost_GTFv102_neurons_filtered_v5.1.0_linear_svc_figures.R
lapply(
  X = lrt.data.l,
  FUN = function(lrt,
                 is.svm) {
    lrt.data <- lrt$lrt.data
    if ("adjusted_p_value" %in% colnames(x = lrt.data)) {
      p.val.col <- "adjusted_p_value"
    } else {
      p.val.col <- "Pr(>Chi)"
    }
    out.data <- unlist(
      x = lapply(
        X = split(x = lrt.data,
                  f = 1:nrow(x = lrt.data)),
        FUN = function(row.data) {
          p.val <- ifelse(test = row.data[[p.val.col]] != 1,
                          yes = signif(x = row.data[[p.val.col]],
                                       digits = 3),
                          no = row.data[[p.val.col]])
          paste0(
            # "AIC = ", format(x = round(x = row.data[["AIC"]],
            #                            digits = 2),
            #                  nsmall = 2),
            # ", LRT = ", format(x = round(x = row.data[["LRT"]],
            #                              digits = 2),
            #                    nsmall = 2),
            ", LRT = ", format(x = round(x = row.data[["LRT"]],
                                         digits = 2),
                               nsmall = 2),
            ", p-value = ", p.val
          )
        }
      )
    )
    
    if (p.val.col == "adjusted_p_value") {
      estimates <- unlist(
        x = lapply(
          X = as.list(x = lrt$model),
          FUN = function(model) {
            signif(
              x = coef(
                object = summary(
                  object = model
                )
              )$cond[2, "Estimate"],
              digits = 3
            )
          }
        )
      )
    } else {
      estimates <- signif(
        x = coef(
          object = summary(
            object = lrt$model
          )
        )$cond[2, "Estimate"],
        digits = 3
      )
    }
    
    out.data <- paste0(
      "estimate = ", estimates,
      out.data
    )
    
    if (colnames(x = lrt.data)[1] != "Df") {
      if (all(lrt.data[[1]] != names(x = estimates))) {
        stop("ERROR")
      }
      out.data <- paste0(
        if (!is.svm) {
          lrt.data[[1]]
        } else {
          gsub(pattern = "_",
               replacement = "-",
               x = lrt.data[[1]])
        },
        ": ",
        out.data
      )
    }
    paste(
      out.data,
      collapse = "; "
    )
  },
  is.svm = FALSE
)
