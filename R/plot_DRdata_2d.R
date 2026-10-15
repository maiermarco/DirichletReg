plot_DRdata_2d <- function(y, rug, main, ylim, colr, lwd, lty){
  colr   <- colr %||% "black"
  y.dens <- density(y, from = 0.0, to = 1.0)
  ylims  <- c(ylim[1L] %||% 0.0, ylim[2L] %||% max(y.dens$y))
  plot(y.dens, type = "l", ylim = ylims, main = main, col = colr, lwd = lwd, lty = lty)
  if(rug) rug(y)
}
