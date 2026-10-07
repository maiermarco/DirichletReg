plot_DRdata_4d <- function(
  x,
  dim.labels,
  ref.lines,
  main,
  cex,
  args.3d,
  theta,
  phi
){

  theta     <- if(is.null(theta)) 40.0 else theta
  phi       <- if(is.null(phi))   25.0 else phi
  ref.lines <- if(is.null(ref.lines)) NULL else ref.lines

  transp <- as.hexmode(round(255.0 * ifelse(is.null(args.3d$transp), 0.25, args.3d$transp), 0L))
  rgl    <- if(is.null(args.3d$rgl)) TRUE else args.3d$rgl

  xyz <- toQuaternary(x$Y)

  corners <- toQuaternary(diag(4L))
  corner.connect <- structure(c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4), dim = c(6L, 2L))

  coo.lab <- 2 * corners - toQuaternary(diag(4) * (0.9 - 0.1 / 3.0) + 0.1 / 3.0)
  lab.col <- cmyk2rgb(diag(4L) + cbind(0.0, 0.0, 0.0, c(0.2, 0.2, 0.2, 0.0)))

  # reference - axes
  ref_axes                  <- matrix(1.0 / 3.0, ncol = 4L, nrow = 4L)
  ref_axes[cbind(1:4, 1:4)] <- 0
  ref_axes_xyz              <- toQuaternary(ref_axes)


# reference points orthogonal to the planes
    .ref_pts <- list(xyz, xyz, xyz, xyz)

    .ref_pts[[1L]] <- xyz + toQuaternary(cbind(1.0 - x$Y[, 1L],
                                                     x$Y[, 1L] / 3.0,
                                                     x$Y[, 1L] / 3.0,
                                                     x$Y[, 1L] / 3.0)) - toQuaternary(matrix(rep(c(1, 0, 0, 0), nrow(x$Y)), ncol = 4L, byrow = TRUE))
    .ref_pts[[2L]] <- xyz + toQuaternary(cbind(      x$Y[, 2L] / 3.0,
                                               1.0 - x$Y[, 2L],
                                                     x$Y[, 2L] / 3.0,
                                                     x$Y[, 2L] / 3.0)) - toQuaternary(matrix(rep(c(0, 1, 0, 0), nrow(x$Y)), ncol = 4L, byrow = TRUE))
    .ref_pts[[3L]] <- xyz + toQuaternary(cbind(      x$Y[, 3L] / 3.0,
                                                     x$Y[, 3L] / 3.0,
                                               1.0 - x$Y[, 3L],
                                                     x$Y[, 3L] / 3.0)) - toQuaternary(matrix(rep(c(0, 0, 1, 0), nrow(x$Y)), ncol = 4L, byrow = TRUE))
    .ref_pts[[4L]] <- xyz + toQuaternary(cbind(      x$Y[, 4L] / 3.0,
                                                     x$Y[, 4L] / 3.0,
                                                     x$Y[, 4L] / 3.0,
                                               1.0 - x$Y[, 4L]      )) - toQuaternary(matrix(rep(c(0, 0, 0, 1), nrow(x$Y)), ncol = 4L, byrow = TRUE))


  if(rgl){

    rgl::view3d(theta = theta, phi = phi)

    rgl::segments3d(
        x = as.vector(rbind(corners[corner.connect[, 1L], 1L], corners[corner.connect[, 2L], 1L])),
        y = as.vector(rbind(corners[corner.connect[, 1L], 2L], corners[corner.connect[, 2L], 2L])),
        z = as.vector(rbind(corners[corner.connect[, 1L], 3L], corners[corner.connect[, 2L], 3L])),
        aspect = 1.0, xlim = 1.0 / sqrt(3.0) + c(-3.0, 3.0) / 4.0, ylim = 1.0 / 3.0 + c(-3.0, 3.0) / 4.0, zlim = 0.25 + c(-3.0, 3.0) / 4.0,
#       xlim = c(0.0, 1.0), ylim = c(-sqrt(3.0) / 6.0, sqrt(3.0) / 2.0), zlim = c(-sqrt(3.0) / 6.0, sqrt(3.0) / 2.0),
        line_antialias = TRUE)

    # axes - references
    rgl::segments3d(
        x = as.vector(rbind(ref_axes_xyz[, 1L], corners[, 1L])),
        y = as.vector(rbind(ref_axes_xyz[, 2L], corners[, 2L])),
        z = as.vector(rbind(ref_axes_xyz[, 3L], corners[, 3L])),
        lwd = 1.0, lty = 2L, col = rep(lab.col, each = 2L), line_antialias = TRUE)

    if(!is.null(ref.lines)){
      browser()
      for(i in ref.lines){
      rgl::segments3d(t3d(.ref_pts[[i]], VTrans)$x, t3d(.ref_pts[[i]], VTrans)$y,
               t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = paste(lab.col[i], transp, sep = "", collapse = ""))
    }}

#    segments(.ref_xy[[2L]]$x, .ref_xy[[2L]]$y, t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = "#BF00BF40")
#    segments(.ref_xy[[3L]]$x, .ref_xy[[3L]]$y, t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = "#BFBF0040")
#    segments(.ref_xy[[4L]]$x, .ref_xy[[4L]]$y, t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = "#00000040")

    rgl::points3d(xyz, cex = cex, col = cmyk2rgb(x$Y), point_antialias = TRUE)

    rgl::text3d(x = coo.lab[, 1L], y = coo.lab[, 2L], z = coo.lab[, 3L], texts = dim.labels, font = 2, col = lab.col, line_antialias = TRUE)

  } else {

    VTrans <- make.VT(theta = theta, phi = phi, d = 1e10, r = 1e10, origin = c(0.5, sqrt(3.0) / 6.0, 0.0))

    xy.corners <- as.data.frame(t3d(corners, VTrans))
    xy.coo.lab <- t3d(coo.lab, VTrans)

    par(mai = rep(0, 4L))
    plot(NULL, xlim = range(xy.coo.lab$x), ylim = range(xy.coo.lab$y), asp = 1.0, axes = FALSE, xlab = "", ylab = "")

    segments(xy.corners[corner.connect[, 1L], 1L], xy.corners[corner.connect[, 1L], 2L],
             xy.corners[corner.connect[, 2L], 1L], xy.corners[corner.connect[, 2L], 2L])

    text(xy.coo.lab, labels = dim.labels, font = 2, col = lab.col)

    ref_axes_xy <- t3d(ref_axes_xyz, VTrans)
    for(i in 1:4) segments(ref_axes_xy$x[i], ref_axes_xy$y[i], xy.corners$x[i], xy.corners$y[i], lwd = 1.0, lty = 2L, col = lab.col[i])

    if(!is.null(ref.lines)){
      for(i in ref.lines){
      segments(t3d(.ref_pts[[i]], VTrans)$x, t3d(.ref_pts[[i]], VTrans)$y,
               t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = paste(lab.col[i], transp, sep = "", collapse = ""))
    }}

#    segments(.ref_xy[[2L]]$x, .ref_xy[[2L]]$y, t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = "#BF00BF40")
#    segments(.ref_xy[[3L]]$x, .ref_xy[[3L]]$y, t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = "#BFBF0040")
#    segments(.ref_xy[[4L]]$x, .ref_xy[[4L]]$y, t3d(xyz, VTrans)$x, t3d(xyz, VTrans)$y, lwd = 0.5, col = "#00000040")

    points(t3d(xyz, VTrans), pch = 16L, cex = cex, col = cmyk2rgb(x$Y))

  } # END if(rgl)



}
