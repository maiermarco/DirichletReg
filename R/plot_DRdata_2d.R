plot_DRdata_2d <- function(y, rug, main, ylim, colr, lwd, lty){
  colr <- if(is.null(colr)) "black" else colr
  y.dens <- density(y, from = 0.0, to = 1.0)
  ylims <- c(ifelse(is.null(ylim[1L]), 0.0, ylim[1L]),
             ifelse(is.null(ylim[2L]), max(y.dens$y), ylim[2L]))
  plot(y.dens, type = "l", ylim = ylims, main = main, col = colr, lwd = lwd, lty = lty)
  if(rug) rug(y)
}
