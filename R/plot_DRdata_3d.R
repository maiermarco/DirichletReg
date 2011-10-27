plot_DRdata_3d <- function(x, entropy.contours, colored, c.grid, ticks, dim.labels, col.scheme,
.main="Ternary Plot", .col=NULL, .pch=16, .cex=1, .lwd=1, .lty=1){

  xy <- coord.trafo(x$Y[,c(2, 3, 1)])
  
  par(mai=rep(0,4))
  plot(NULL,axes=F,xlab="",ylab="", xlim=c(-.0866025,1+.0866025), ylim=c(0-.05,sqrt(3)/2+.1), asp=1)
  if(entropy.contours){
    if(entropy.colors){
      ### FIX ME
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

#DEBUG for(a in 0:10) abline(v=a/10, lty=2)

    if(c.grid){
      segments(z.grid[,1],z.grid[,2],z.grid[,3],z.grid[,4], lwd=.5, lty=2, col=colorz[1])
      segments(x.grid[,1],x.grid[,2],x.grid[,3],x.grid[,4], lwd=.5, lty=2, col=colorz[2])
      segments(y.grid[,1],y.grid[,2],y.grid[,3],y.grid[,4], lwd=.5, lty=2, col=colorz[3])
    }
  }
  if(colored){
    if(!is.null(.col)){
      points(xy, pch=.pch, lwd=.lwd, col=.col, cex=.cex)
    } else if(col.scheme == "dims"){
      points(xy, pch=.pch, lwd=.lwd, col=rgb(x$Y), cex=.cex)
    } else if(col.scheme == "entropy"){
      ent.col <- entropy.colors[round(99*(rowSums(-x$Y*log(x$Y))/log(ncol(x$Y))),0)+1]
      points(xy, pch=.pch, lwd=.lwd, col=ent.col,cex=.cex)
      
      y.cuts <- seq(sqrt(3)/2,sqrt(3)/6,length.out=101)
      rect(.95,y.cuts[-101],1,y.cuts[-1],col=entropy.colors,border=NA)
      rect(.95,y.cuts[1],1,y.cuts[101])
      text(.90,1/sqrt(3),srt=90,labels="Entropy")
    }
  } else if(!colored) {
    if(!is.null(.col)){
      points(xy, pch=.pch, lwd=.lwd, col=.col, cex=.cex)
    } else {
      points(xy, pch=.pch, lwd=.lwd, cex=.cex)
    }
  } else { stop("error! specify color=TRUE or FALSE") }

  ### bounding triangle and white space
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
             tk.coo-(1/120),-1/(40*sqrt(3)), col=colorz[3]) # sin(-60°)/60 and /30 for text

    text(y.grid[,3]+(1/30),y.grid[,4],labels=rev(tk.lab),cex=.8, col=colorz[1])
    text(z.grid[,1]-(1/60),z.grid[,2]+1/(20*sqrt(3)),labels=rev(tk.lab),cex=.8, col=colorz[2])
    text(tk.coo-(1/60),-1/(20*sqrt(3)),labels=tk.lab,cex=.8, col=colorz[3])
  }

  # dimension labels
  text(.5,sqrt(3)/2+(1/15),labels=dim.labels[1],font=2, col=colorz[1])
  text(-1/(10*sqrt(3)),-1/30,labels=dim.labels[2],font=2, col=colorz[2])
  text(1+1/(10*sqrt(3)),-1/30,labels=dim.labels[3],font=2, col=colorz[3])

}
