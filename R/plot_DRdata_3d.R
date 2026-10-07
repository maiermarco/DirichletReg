plot_DRdata_3d <- function(
  x,
  entropy.contours,
  colored,
  c.grid,
  ticks,
  dim.labels,
  col.scheme,
  .main             = "Ternary Plot",
  .col              = NULL,
  .pch              = 16L,
  .cex              = 1.0,
  .lwd              = 1.0,
  .lty              = 1L
){
  
  xy <- toTernary(x$Y)
  
  par(mai = rep(0.0, 4L))
  plot(NULL, axes = FALSE, ann = FALSE, xlim = c(-.0866025, 2.0 / sqrt(3.0) + .0866025), ylim = c(0.0 - 0.05, 1.0 + 0.1), asp = 1.0)
###
###
###
###
### TODO FIX ENTROPY!
  if(entropy.contours){
    if(entropy.colors){
      heatcols <- c("#D9545F", "#E06951", "#E47C42", "#E88F33", "#EAA428", "#E9B62D", "#E8C842", "#E5D961", "#E2E6BD")   # colorspace : heat_hcl(100)[round(seq(0, 100, length.out = 10), 0L)]
      for(i in 9:1) polygon(.entropy.coordinates[[i]][, 1L], .entropy.coordinates[[i]][, 2L], col = heatcols[i], border = NA)
    } else {
      for(i in 9:1) lines(.entropy.coordinates[[i]], lwd = 0.25)
    }
  }

  if(colored){
    colorz <- hsv(c(0.0, 1.0, 2.0) / 3.0, 0.8, 0.4)
    entropy.colors <- c("#D33F6A", "#D34269", "#D44468", "#D54667", "#D54866", "#D64A65", "#D74C63", "#D74E62", "#D85061", "#D95260", "#D9545F", "#DA565E", "#DB585D", "#DB5A5B", "#DC5C5A", "#DC5E59", "#DD6058", "#DD6256", "#DE6355", "#DF6554", "#DF6753", "#E06951", "#E06B50", "#E16D4F", "#E16E4D", "#E2704C", "#E2724B", "#E27449", "#E37548", "#E37747", "#E47945", "#E47B44", "#E47C42", "#E57E41", "#E58040", "#E6823E", "#E6833D", "#E6853B", "#E6873A", "#E78939", "#E78A37", "#E78C36", "#E88E35", "#E88F33", "#E89132", "#E89331", "#E89430", "#E9962E", "#E9982D", "#E99A2C", "#E99B2B", "#E99D2B", "#E99F2A", "#E9A029", "#E9A229", "#EAA428", "#EAA528", "#EAA728", "#EAA928", "#EAAA28", "#EAAC28", "#EAAE29", "#EAAF29", "#EAB12A", "#EAB32B", "#E9B42C", "#E9B62D", "#E9B82F", "#E9B930", "#E9BB32", "#E9BC33", "#E9BE35", "#E9C037", "#E8C139", "#E8C33B", "#E8C53D", "#E8C640", "#E8C842", "#E7C944", "#E7CB47", "#E7CD4A", "#E7CE4C", "#E6D04F", "#E6D152", "#E6D355", "#E5D458", "#E5D65B", "#E5D85E", "#E5D961", "#E4DB64", "#E4DC68", "#E4DE6C", "#E3DF6F", "#E3E173", "#E3E278", "#E2E37D", "#E2E582", "#E2E688", "#E2E791", "#E2E6BD") # colorspace : heat_hcl(100) # nolint
  } else {
    colorz <- rep("black", 3L)
  }

  if(c.grid || ticks){
    g.main <- (1:9) / 10
    g.aux1 <- 1.0 - g.main

    c1.grid <- cbind(toTernaryVectors(g.main, g.aux1, 0), toTernaryVectors(g.main, 0, g.aux1))
    c2.grid <- cbind(toTernaryVectors(0, g.main, g.aux1), toTernaryVectors(g.aux1, g.main, 0))
    c3.grid <- cbind(toTernaryVectors(0, g.aux1, g.main), toTernaryVectors(g.aux1, 0, g.main))

#DEBUG for(a in 0:10) abline(v = a / 10, lty = 2L)

    if(c.grid){
      segments(c1.grid[, 1L], c1.grid[, 2L], c1.grid[, 3L], c1.grid[, 4L], lwd = 0.5, lty = "36", col = colorz[1L])
      segments(c2.grid[, 1L], c2.grid[, 2L], c2.grid[, 3L], c2.grid[, 4L], lwd = 0.5, lty = "36", col = colorz[2L])
      segments(c3.grid[, 1L], c3.grid[, 2L], c3.grid[, 3L], c3.grid[, 4L], lwd = 0.5, lty = "36", col = colorz[3L])
    }
  }
  if(colored){
    if(!is.null(.col)){
      points(xy, pch = .pch, lwd = .lwd, col = .col, cex = .cex)
    } else if(col.scheme == "dims"){
      points(xy, pch = .pch, lwd = .lwd, col = rgb(x$Y), cex = .cex)
    } else if(col.scheme == "entropy"){
      ent.col <- entropy.colors[round(99.0 * (rowSums(-x$Y * log(x$Y)) / log(ncol(x$Y))), 0L) + 1L]
      points(xy, pch = .pch, lwd = .lwd, col = ent.col, cex = .cex)

      y.cuts <- seq(sqrt(3.0) / 2.0, sqrt(3.0) / 6.0, length.out = 101L)
      rect(0.95, y.cuts[-101L]  , 1.0, y.cuts[-1L], col = entropy.colors, border = NA)
      rect(0.95, y.cuts[1L]     , 1.0, y.cuts[101L])
      text(0.90, 1.0 / sqrt(3.0), srt = 90, labels = "Entropy")
    }
  } else if(!colored) {
    if(!is.null(.col)){
      points(xy, pch = .pch, lwd = .lwd, col = .col, cex = .cex)
    } else {
      points(xy, pch = .pch, lwd = .lwd, cex = .cex)
    }
  } else { stop("error! specify color = TRUE or FALSE") }

  if(ticks){
    tk.coo <- (1:9) / 10
    tk.lab <- sub("^0\\.", "\\.", tk.coo)

    segments(c3.grid[, 3L],
             c3.grid[, 4L],
             c3.grid[, 3L] + 1.0 / 60.0,
             c3.grid[, 4L], col = colorz[1L])
    segments(c1.grid[, 1L],
             c1.grid[, 2L],
             c1.grid[, 1L] - (1.0 / 120.0),
             c1.grid[, 2L] + 1.0 / (40.0 * sqrt(3.0)), col = colorz[2L])
    segments(c3.grid[, 1L],
                       0.0,
             c3.grid[, 1L] - (1.0 / 120.0),
             -1.0 / (40.0 * sqrt(3.0)), col = colorz[3L]) # sin(-60 deg)/60 and /30 for text

    text(c3.grid[, 3L] + (1.0 / 30.0), c3.grid[, 4L]                           , labels = rev(tk.lab), cex = 0.8, col = colorz[1L])
    text(c1.grid[, 1L] - (1.0 / 60.0), c1.grid[, 2L] + 1.0 / (20.0 * sqrt(3.0)), labels = rev(tk.lab), cex = 0.8, col = colorz[2L])
    text(c3.grid[, 1L] - (1.0 / 60.0), -1.0 / (20.0 * sqrt(3.0))               , labels = tk.lab     , cex = 0.8, col = colorz[3L])
  }

  # dimension labels
  text( 1.0 / sqrt(3.0)                           ,  1.0 + (1.0 / 15.0), labels = dim.labels[1L], font = 2, col = colorz[1L])
  text(-1.0 / (10.0 * sqrt(3.0))                  , -1.0 / 30.0        , labels = dim.labels[2L], font = 2, col = colorz[2L])
  text( 2.0 / sqrt(3.0) + 1.0 / (10.0 * sqrt(3.0)), -1.0 / 30.0        , labels = dim.labels[3L], font = 2, col = colorz[3L])

  ### bounding triangle and white space
  polygon(c(0.0, 2.0 / sqrt(3.0), 1.0 / sqrt(3.0), 0.0),
          c(0.0,             0.0,             1.0, 0.0))

}
