# Blend color with alpha
BlendColorAlpha <- function(
    color,
    alpha
){
  color.rgb <- grDevices::col2rgb(col = color) / 255
  blended.color <- grDevices::rgb(
    red = alpha * color.rgb[1,] + (1 - alpha),
    green = alpha * color.rgb[2,] + (1 - alpha),
    blue = alpha * color.rgb[3,] + (1 - alpha)
  )
  return(blended.color)
}
