# Source a specific function from a file (adapted from https://stat.ethz.ch/pipermail/r-help/2007-June/134495.html)
sourcefunction <- function(
    file,
    fun,
    ...
) {
  local({
    source(file = file,
           local = TRUE,
           ...)
    fun <- get(x = fun)
    environment(fun = fun) <- .GlobalEnv
    fun
  })
}
