# Get x and y values from ECDF
ecdfxy <- function(
    x,
    min.x = min(x)
){
  ecdf.fun <- ecdf(x = x)
  x <- sort(x = unique(x = x))
  y <- ecdf.fun(x)
  if (min.x < min(x)) {
    x <- c(min.x, x)
    y <- c(0, y)
  }
  return(list(x = x,
              y = y))
}
