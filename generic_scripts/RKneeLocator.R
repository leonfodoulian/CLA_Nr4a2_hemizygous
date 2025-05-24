# R function to compute KneeLocator from the kneed python package
RKneeLocator <- function(
    x,
    y,
    S,
    curve,
    direction,
    online,
    ...
){
  # import python build-in functions as bi
  bi <- reticulate::import_builtins()
  # import kneed as kn
  kn <- reticulate::import(module = "kneed")
  # from matplotlib import pyplot as plt
  plt <- reticulate::import(module = "matplotlib.pyplot")
  
  # Check if curve argument is valid
  curve <- match.arg(arg = curve,
                     choices = c("concave", "convex"),
                     several.ok = FALSE)
  
  # Check if direction argument is valid
  direction <- match.arg(arg = direction,
                         choices = c("increasing", "decreasing"),
                         several.ok = FALSE)
  
  # Check if online argument is logical
  if (!is.logical(x = online)) {
    stop('argument "online" should be logical')
  }
  
  # Detect knee using kn.KneeLocator()
  kneed.res <- kn$KneeLocator(x = x,
                              y = y,
                              S = S,
                              curve = curve,
                              direction = direction,
                              online = online,
                              ...)
  
  # Return output in R format
  return(list(all_elbows = unlist(x = bi$list(kneed.res$all_elbows)),
              all_elbows_y = unlist(x = kneed.res$all_elbows_y),
              all_knees = unlist(x = bi$list(kneed.res$all_knees)),
              all_knees_y = unlist(x = kneed.res$all_knees_y),
              all_norm_elbows = unlist(x = bi$list(kneed.res$all_norm_elbows)),
              all_norm_elbows_y = unlist(x = kneed.res$all_norm_elbows_y),
              all_norm_knees = unlist(x = bi$list(kneed.res$all_norm_knees)),
              all_norm_knees_y = unlist(x = kneed.res$all_norm_knees_y),
              curve = kneed.res$curve,
              direction = kneed.res$direction,
              Ds_y = kneed.res$Ds_y,
              elbow = kneed.res$elbow,
              elbow_y = kneed.res$elbow_y,
              find_knee = unlist(x = kneed.res$find_knee()),
              knee = kneed.res$knee,
              knee_y = kneed.res$knee_y,
              maxima_indices = kneed.res$maxima_indices,
              minima_indices = kneed.res$minima_indices,
              N = kneed.res$N,
              norm_elbow = kneed.res$norm_elbow,
              norm_elbow_y = kneed.res$norm_elbow_y,
              norm_knee = kneed.res$norm_knee,
              norm_knee_y = kneed.res$norm_knee_y,
              online = kneed.res$online,
              plot_knee = R.devices::capturePlot(
                expr = {
                  # plt$rcParams["font.family"] <- "Arial"
                  # kneed.res$plot_knee()
                  # plt$show()
                  plot(x = kneed.res$x,
                       y = kneed.res$y,
                       type = "l",
                       main = "Knee Point",
                       xlab = "x",
                       ylab = "y",
                       lty = "solid",
                       lwd = 2,
                       col = "blue")
                  abline(v = kneed.res$knee,
                         lty = "longdash",
                         lwd = 2,
                         col = "skyblue3")
                  legend(x = min(kneed.res$x),
                         y = max(kneed.res$y),
                         legend = c("data",
                                    "knee/elbow"),
                         lty = c("solid",
                                 "longdash"),
                         lwd = c(2,2),
                         bty = "n",
                         col = c("blue",
                                 "skyblue3"))
                }
              ),
              plot_knee_normalized = R.devices::capturePlot(
                expr = {
                  # plt$rcParams["font.family"] <- "Arial"
                  # kneed.res$plot_knee_normalized()
                  # plt$xlim(0, 1)
                  # plt$ylim(0, 1)
                  # plt$show()
                  plot(x = kneed.res$x_normalized,
                       y = kneed.res$y_normalized,
                       type = "l",
                       xlim = c(0,1),
                       ylim = c(0,1),
                       main = "Normalized Knee Point",
                       xlab = "x_normalized",
                       ylab = "y_normalized",
                       lty = "solid",
                       lwd = 2,
                       col = "blue")
                  lines(x = kneed.res$x_difference,
                        y = kneed.res$y_difference,
                        type = "l",
                        lty = "solid",
                        lwd = 2,
                        col = "red")
                  abline(v = kneed.res$norm_knee,
                         lty = "longdash",
                         lwd = 2,
                         col = "skyblue3")
                  legend(x = 0,
                         y = 1,
                         legend = c("normalized curve",
                                    "difference curve",
                                    "knee/elbow"),
                         lty = c("solid",
                                 "solid",
                                 "longdash"),
                         lwd = c(2,2,2),
                         bty = "n",
                         col = c("blue",
                                 "red",
                                 "skyblue3"))
                }
              ),
              polynomial_degree = kneed.res$polynomial_degree,
              S = kneed.res$S,
              Tmx = kneed.res$Tmx,
              x = kneed.res$x,
              x_difference = kneed.res$x_difference,
              x_difference_maxima = kneed.res$x_difference_maxima,
              x_difference_minima = kneed.res$x_difference_minima,
              x_normalized = kneed.res$x_normalized,
              y = kneed.res$y,
              y_difference = kneed.res$y_difference,
              y_difference_maxima = kneed.res$y_difference_maxima,
              y_difference_minima = kneed.res$y_difference_minima,
              y_normalized = kneed.res$y_normalized))
}
