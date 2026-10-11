cmyk2rgb <- function(cmyk){
  rgb(
    red   = (cmyk[, 1L] - 1.0) * (cmyk[, 4L] - 1.0),
    green = (cmyk[, 2L] - 1.0) * (cmyk[, 4L] - 1.0),
    blue  = (cmyk[, 3L] - 1.0) * (cmyk[, 4L] - 1.0)
  )
}
