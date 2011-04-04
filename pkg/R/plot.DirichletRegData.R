plot.DirichletRegData <- function(x, ...,
                                  dims,   
                                  explore,   
                                  model,   
                                  c.grid=TRUE,   
                                  ticks=TRUE,   
                                  colored=TRUE,   
                                  col.scheme="classes",   
                                  entropy.contours=FALSE,   
                                  entropy.colors=FALSE,   
                                  pt.size=.75,   
                                  dim.labels,
                                  args.3d=list(rgl=TRUE, theta=NULL, phi=NULL, ref.lines=NULL, ...)  
                                  ){

  if(missing(dims))    dims    <- NULL
  if(missing(explore)) explore <- NULL
  if(missing(model))   model   <- NULL
  
  if(class(x) != "DirichletRegData") stop("data must be prepared by 'DR.data()'")



  plot.char <-   16
  line.wdth <-    1
  point.col <- NULL
  
  add.args <- list(...)
  
  if(length(add.args) > 0){
    if(!is.null(add.args$pch)){ plot.char <- add.args$pch }
    if(!is.null(add.args$lwd)){ line.wdth <- add.args$lwd }
    if(!is.null(add.args$col)){ point.col <- add.args$col }
  }



  .marginal <- FALSE
  
  if(is.null(dims)){
    if(x$dims > 4){
      x$Y         <- x$Y[, 1:4]
      x$dims      <- 4
      x$dim.names <- x$dim.names[1:4]
      warning("data contains > 4 variables. the first four are being used. to change this, set the 'dims' argument appropriately.")
      .marginal <- TRUE
    } else {
      dims <- 1:x$dims
    }
  } else {
    if((length(dims) < 2) | (length(dims) > 4)) stop("the argument 'dims' must have 2, 3 or 4 elements")
    x$Y         <- x$Y[, dims]
    x$dims      <- length(dims)
    x$dim.names <- x$dim.names[dims]
    .marginal <- TRUE
  }
  

 
  if(missing(dim.labels)){
    dim.labels <- x$dim.names
  }


  
  if(!is.null(explore) & !is.null(model)) stop("model and explore cannot be specified at the same time.")
  
  if(!is.null(explore)){
    frmla <- explore$formula
    dta   <- explore$data

    ex.vars <- model.frame(frmla, dta)
    if( (length(ex.vars) < 1) | (length(ex.vars) > 3) ) stop("only 1 - 3 variables allowed in 'explore'.")
    if( (sum(sapply(ex.vars, is.factor)) > 2) | (sum(sapply(ex.vars, is.numeric)) > 1) ) stop("check your explore-variables")
    
    logit.Y <- as.matrix(inv.logit(x$Y))
    
    ex.model <- lm(update(logit.Y~1, frmla), dta)
    
    .to.plot <- predict(ex.model, data.frame(depth=seq(0,100,length.out=100)))
    .to.plot <- exp(.to.plot) / (1 + exp(.to.plot))
    .to.plot <- .to.plot / rowSums(.to.plot)
  } else {
    .to.plot <- NULL
  }
  
  if(!is.null(model)){
    cat("\nnot implemented yet!\n")
  }




  if(x$dims == 2){










    warning("not implemented yet")





  } else if(x$dims == 3) {











    xy <- coord.trafo(x$Y[,c(2, 3, 1)])

    par(mai=rep(0,4))
    plot(NULL,asp=1,
         xlim=c(-.0866025,1+.0866025),ylim=c(0-.05,sqrt(3)/2+.1),
         axes=F,xlab="",ylab="")
    if(entropy.contours){
      if(entropy.colors){
        
      } else {
        for(i in 9:1) lines(.entropy.coordinates[[i]], lwd=.25)
      }
    }
         
    if(colored){
      colorz <- hex(HLS(c(0,120,240),0.25,.75), fixup=TRUE)
      entropy.colors <- heat_hcl(100)
    } else {
      colorz <- rep("black", 3)
    }
    
    if(c.grid | ticks){
      g.main <- seq(.1,.9,by=.1)
      g.aux1 <- 1 - g.main
  
      y.grid <- cbind(coord.trafo(cbind(g.aux1,g.main,0)),
                      coord.trafo(cbind(0,g.main,g.aux1)))
      z.grid <- cbind(coord.trafo(cbind(g.aux1,0,g.main)),
                      coord.trafo(cbind(0,g.aux1,g.main)))
      x.grid <- cbind(coord.trafo(cbind(g.main,g.aux1,0)),
                      coord.trafo(cbind(g.main,0,g.aux1)))
  
      if(c.grid){
        segments(z.grid[,1],z.grid[,2],z.grid[,3],z.grid[,4], lwd=.5, lty=2, col=colorz[1])
        segments(x.grid[,1],x.grid[,2],x.grid[,3],x.grid[,4], lwd=.5, lty=2, col=colorz[2])
        segments(y.grid[,1],y.grid[,2],y.grid[,3],y.grid[,4], lwd=.5, lty=2, col=colorz[3])
      }
    }
  
    if(colored == TRUE){
      if(!is.null(point.col)){
        points(xy, pch=plot.char, lwd=line.wdth, col=point.col, cex=pt.size)
      } else if(col.scheme == "classes"){
        points(xy, pch=plot.char, lwd=line.wdth, col=rgb(x$Y), cex=pt.size)
      } else if(col.scheme == "entropy"){
        ent.col <- entropy.colors[round(99*(rowSums(-x$Y*log(x$Y))/log(ncol(x$Y))),0)+1]
        points(xy, pch=plot.char, lwd=line.wdth, col=ent.col,cex=pt.size)
        
        y.cuts <- seq(sqrt(3)/2,sqrt(3)/6,length.out=101)
        rect(.95,y.cuts[-101],1,y.cuts[-1],col=entropy.colors,border=NA)
        rect(.95,y.cuts[1],1,y.cuts[101])
        text(.90,1/sqrt(3),srt=90,labels="Entropy")
      }
    } else {
      if(!is.null(point.col)){
        points(xy, pch=plot.char, lwd=line.wdth, col=point.col, cex=pt.size)
      } else {
        points(xy, pch=plot.char, lwd=line.wdth, cex=pt.size)
      }
    }

    
    polygon(c(-1,-1,2,2),c(0,-1,-1,0),col="white",border=NA)
    polygon(c(-1,0,.5,.5,-1),c(0,0,sqrt(3)/2,3,3),col="white",border=NA)
    polygon(c(2,1,.5,.5,2),c(0,0,sqrt(3)/2,3,3),col="white",border=NA)
      polygon(c(0,1,.5,0),c(0,0,sqrt(3)/2,0))

    if(ticks){
      tk.coo <- seq(.1,.9,by=.1)
      tk.lab <- substr(tk.coo,2,3)

      segments(y.grid[,3],y.grid[,4],
               y.grid[,3]+1/60,y.grid[,4], col=colorz[1])
      segments(z.grid[,1],z.grid[,2],
               z.grid[,1]-(1/120),z.grid[,2]+1/(40*sqrt(3)), col=colorz[2])
      segments(tk.coo,0,
               tk.coo-(1/120),-1/(40*sqrt(3)), col=colorz[3]) 

      text(y.grid[,3]+(1/30),y.grid[,4],labels=rev(tk.lab),cex=.8, col=colorz[1])
      text(z.grid[,1]-(1/60),z.grid[,2]+1/(20*sqrt(3)),labels=rev(tk.lab),cex=.8, col=colorz[2])
      text(tk.coo-(1/60),-1/(20*sqrt(3)),labels=tk.lab,cex=.8, col=colorz[3])
    }

    if(!is.null(.to.plot)){

      .arrow.trans <- coord.trafo(.to.plot)
      .end <- nrow(.arrow.trans)
      lines(.arrow.trans[-.end,1], .arrow.trans[-.end,2], lwd=2)
      arrows(.arrow.trans[.end-1,1], .arrow.trans[.end-1,2],
             .arrow.trans[.end  ,1], .arrow.trans[.end  ,2],
             angle=20, length=.1, lwd=2)
    
    }

    
    text(.5,sqrt(3)/2+(1/15),labels=dim.labels[1],font=2, col=colorz[1])
    text(-1/(10*sqrt(3)),-1/30,labels=dim.labels[2],srt=30,font=2, col=colorz[2])
    text(1+1/(10*sqrt(3)),-1/30,labels=dim.labels[3],srt=-30,font=2, col=colorz[3])
    
  } else if(x$dims == 4){










    theta <- if(is.null(args.3d$theta)) 40 else args.3d$theta
    phi <- if(is.null(args.3d$phi)) 25 else args.3d$phi
    transp <- as.hexmode(round(255*ifelse(is.null(args.3d$transp), .25, args.3d$transp),0))
    ref.lines <- if(is.null(args.3d$ref.lines)) NULL else args.3d$ref.lines
    rgl <- if(is.null(args.3d$rgl)) TRUE else args.3d$rgl

    xyz <- coord.trafo(x$Y)

    corners <- coord.trafo(diag(4))
    corner.connect <- structure(c(1,1,1,2,2,3,2,3,4,3,4,4), .Dim=c(6,2))

    coo.lab <- 2 * corners - coord.trafo(diag(4)*(.9-.1/3)+.1/3)
    lab.col <- cmyk2rgb(diag(4)+cbind(0,0,0,c(.2,.2,.2,0)))

    
    ref_axes                 <- matrix(1/3, ncol=4, nrow=4)
    ref_axes[cbind(1:4,1:4)] <- 0
    ref_axes_xyz             <- coord.trafo(ref_axes)

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

      points3d(xyz, cex=pt.size, col=cmyk2rgb(x$Y), point_antialias=TRUE)

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
                                               
      if(!is.null(ref.lines)){
        for(i in ref.lines){
        segments(t3d(.ref_pts[[i]],VTrans)$x, t3d(.ref_pts[[i]],VTrans)$y,
                 t3d(xyz,VTrans)$x,t3d(xyz,VTrans)$y, lwd=.5, col=paste(lab.col[i], transp, sep="", collapse=""))
      }}
  
  
  
  
  
      points(t3d(xyz,VTrans), pch=16, cex=pt.size, col=cmyk2rgb(x$Y))

    } 

  }

}











