# Discretise values
discretise <- function(x,
                       by = 20) {
  cut(x = x,
      breaks = seq(from = (by * floor(x = min(x) / by)),
                   to = (by * ceiling(x = max(x) / by)),
                   by = by),
      include.lowest = TRUE,
      dig.lab = 50)
}
