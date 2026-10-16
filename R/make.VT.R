D2R <- function(degrees) degrees * 0.0174532925199432957692369076849 # = degrees * pi / 180



t3d <- function(XYZ, ViewTrans3D) trans3d(XYZ[, 1L], XYZ[, 2L], XYZ[, 3L], ViewTrans3D)



make.VT <- function(
  theta  = 0.0,
  phi    = 15.0,
  r      = sqrt(3.0),
  d      = 1.0,
  expand = 1.0,
  origin = c(0.0, 0.0, 0.0),
  scale  = c(1.0, 1.0, 1.0)
){

  xc <- origin[1L]
  yc <- origin[2L]
  zc <- origin[3L]

  xs <- scale[1L]
  ys <- scale[2L]
  zs <- scale[3L]

  VT <- diag(4L)   # initialize

  TT <- diag(4L)   # center @ origin
  TT[4L, 1L:3L] <- c(-xc, -yc, -zc)
  VT <- VT %*% TT


  TT <- diag(4L)   # scale extents to [-1, 1]
  TT[cbind(1L:3L, 1L:3L)] <- c(1.0 / xs, 1.0 / ys, expand / zs)
  VT <- VT %*% TT


  TT <- diag(4L)   # rotate x-y plane to horizontal
  TT[cbind(c(2L, 3L, 3L, 2L), c(2L, 2L, 3L, 3L))] <- c(cos(D2R(-90.0)), -sin(D2R(-90.0)),
                                                       cos(D2R(-90.0)),  sin(D2R(-90.0)))
  VT <- VT %*% TT


  TT <- diag(4L)   # azimuthal rotation (theta)
  TT[cbind(c(1L, 3L, 3L, 1L), c(1L, 1L, 3L, 3L))] <- c(cos(D2R(-theta)),  sin(D2R(-theta)),
                                                       cos(D2R(-theta)), -sin(D2R(-theta)))
  VT <- VT %*% TT


  TT <- diag(4L)   # elevation rotation (phi)
  TT[cbind(c(2L, 3L, 3L, 2L), c(2L, 2L, 3L, 3L))] <- c(cos(D2R(phi)), -sin(D2R(phi)),
                                                       cos(D2R(phi)),  sin(D2R(phi)))
  VT <- VT %*% TT


  TT <- diag(4L)   # translate eyepoint to origin
  TT[4L, 1L:3L] <- c(0.0, 0.0, -r - d)
  VT <- VT %*% TT


  TT <- diag(4L)   # perspective
  TT[3L, 4L] <- -1.0 / d
  VT <- VT %*% TT

  return(VT)
}
