plot_DRdata_4d <- function(x,dim.labels,ref.lines,
main, cex,
args.3d, theta, phi){

  theta <- if(is.null(theta)) 40 else theta
  phi <- if(is.null(phi)) 25 else phi
  ref.lines <- if(is.null(ref.lines)) NULL else ref.lines

  transp <- as.hexmode(round(255*ifelse(is.null(args.3d$transp), .25, args.3d$transp),0))
  rgl <- if(is.null(args.3d$rgl)) TRUE else args.3d$rgl

  xyz <- coord.trafo(x$Y)

  corners <- coord.trafo(diag(4))
  corner.connect <- structure(c(1,1,1,2,2,3,2,3,4,3,4,4), .Dim=c(6,2))

  coo.lab <- 2 * corners - coord.trafo(diag(4)*(.9-.1/3)+.1/3)
  lab.col <- cmyk2rgb(diag(4)+cbind(0,0,0,c(.2,.2,.2,0)))

  
  ref_axes                 <- matrix(1/3, ncol=4, nrow=4)
  ref_axes[cbind(1:4,1:4)] <- 0
  ref_axes_xyz             <- coord.trafo(ref_axes)



    .ref_pts <- list(xyz, xyz, xyz, xyz)

    .ref_pts[[1]] <- xyz + coord.trafo(cbind(1-x$Y[,1],
                                               x$Y[,1]/3,
                                               x$Y[,1]/3,
                                               x$Y[,1]/3)) - coord.trafo(matrix(rep(c(1,0,0,0),nrow(x$Y)),ncol=4,byrow=T))
    .ref_pts[[2]] <- xyz + coord.trafo(cbind(  x$Y[,2]/3,
                                             1-x$Y[,2],
                                               x$Y[,2]/3,
                                               x$Y[,2]/3)) - coord.trafo(matrix(rep(c(0,1,0,0),nrow(x$Y)),ncol=4,byrow=T))
    .ref_pts[[3]] <- xyz + coord.trafo(cbind(  x$Y[,3]/3,
                                               x$Y[,3]/3,
                                             1-x$Y[,3],
                                               x$Y[,3]/3)) - coord.trafo(matrix(rep(c(0,0,1,0),nrow(x$Y)),ncol=4,byrow=T))
    .ref_pts[[4]] <- xyz + coord.trafo(cbind(  x$Y[,4]/3,
                                               x$Y[,4]/3,
                                               x$Y[,4]/3,
                                             1-x$Y[,4])) - coord.trafo(matrix(rep(c(0,0,0,1),nrow(x$Y)),ncol=4,byrow=T))

                                             
  if(rgl){

    view3d(theta=theta, phi=phi)
    
    segments3d(x=as.vector(rbind(corners[corner.connect[,1],1], corners[corner.connect[,2],1])),
               y=as.vector(rbind(corners[corner.connect[,1],2], corners[corner.connect[,2],2])),
               z=as.vector(rbind(corners[corner.connect[,1],3], corners[corner.connect[,2],3])),
               xlim=c(0,1), ylim=c(-sqrt(3)/6, sqrt(3)/2), zlim=c(-sqrt(3)/6, sqrt(3)/2),
               line_antialias=TRUE)
               
    
    segments3d(x=as.vector(rbind(ref_axes_xyz[,1],corners[,1])),
               y=as.vector(rbind(ref_axes_xyz[,2],corners[,2])),
               z=as.vector(rbind(ref_axes_xyz[,3],corners[,3])), lwd=1, lty=2, col=rep(lab.col,each=2), line_antialias=TRUE)

    if(!is.null(ref.lines)){
    browser()
      for(i in ref.lines){
      segments3d(t3d(.ref_pts[[i]],VTrans)$x, t3d(.ref_pts[[i]],VTrans)$y,
               t3d(xyz,VTrans)$x,t3d(xyz,VTrans)$y, lwd=.5, col=paste(lab.col[i], transp, sep="", collapse=""))
    }}





    points3d(xyz, cex=cex, col=cmyk2rgb(x$Y), point_antialias=TRUE)

    text3d(x=coo.lab[,1], y=coo.lab[,2], z=coo.lab[,3], texts=dim.labels, font=2, col=lab.col, line_antialias=TRUE)

  } else {

    VTrans <- make.VT(theta=theta, phi=phi, d=1e10, r=1e10, origin=c(.5, sqrt(3)/6, 0))

    xy.corners <- as.data.frame(t3d(corners,VTrans))
    xy.coo.lab <- t3d(coo.lab,VTrans)
              
    par(mai=rep(0,4))
    plot(NULL, xlim=range(xy.coo.lab$x), ylim=range(xy.coo.lab$y), asp=1,
         axes=F, xlab="", ylab="")

    segments(xy.corners[corner.connect[,1],1], xy.corners[corner.connect[,1],2],
             xy.corners[corner.connect[,2],1], xy.corners[corner.connect[,2],2])

    text(xy.coo.lab, labels=dim.labels, font=2, col=lab.col)

    ref_axes_xy <- t3d(ref_axes_xyz, VTrans)
    for(i in 1:4) segments(ref_axes_xy$x[i], ref_axes_xy$y[i], xy.corners$x[i], xy.corners$y[i], lwd=1, lty=2, col=lab.col[i])

    if(!is.null(ref.lines)){
      for(i in ref.lines){
      segments(t3d(.ref_pts[[i]],VTrans)$x, t3d(.ref_pts[[i]],VTrans)$y,
               t3d(xyz,VTrans)$x,t3d(xyz,VTrans)$y, lwd=.5, col=paste(lab.col[i], transp, sep="", collapse=""))
    }}





    points(t3d(xyz,VTrans), pch=16, cex=cex, col=cmyk2rgb(x$Y))

  } 



}
