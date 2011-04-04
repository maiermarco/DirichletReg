starplot <- function(x, persons, vars){

  resid.mat <- cbind(runif(10),runif(10),runif(10),runif(10),runif(10),runif(10),runif(10))

  if(missing(persons)) persons <- 1:nrow(resid.mat)
  if(missing(vars)) vars <- 1:ncol(resid.mat)

  resid.mat <- resid.mat[persons, vars, drop=F]

  limit <- c(min(resid.mat), max(resid.mat)) + c(-.1, .1)*sum(abs(range(resid.mat)))
  
  offs <- limit[1]
  limit <- limit[2] - offs
  
  plot.mat <- resid.mat - offs

  nvar <- ncol(resid.mat)
  nper <- nrow(resid.mat)

  unit.x <- cos(D2R(90-((1:nvar)-1)*(360/nvar)))
  unit.y <- sin(D2R(90-((1:nvar)-1)*(360/nvar)))
  
  par(mai=c(0,0,0,0))
  
  plot(NULL, xlim=range(1.25*limit*unit.x), ylim=range(1.25*limit*unit.y), asp=1)
  
  circles <- seq(0, 2*pi, length.out=200)+pi
  circles <- cbind(cos(circles), sin(circles))
  
  rnd <- -5; tixx <- unique(round(as.numeric(resid.mat),rnd))
  
  while(length(tixx) < 3){ tixx <- unique(round(as.numeric(resid.mat),rnd))
                           rnd <- rnd + 1 }
  
  for(i in tixx) lines((i - offs)*circles, lty=2)

  for(i in 1:nvar) segments(0, 0, limit*unit.x[i], limit*unit.y[i], lwd=3)
  
  for(i in 1:nper) polygon(plot.mat[i,]*unit.x, plot.mat[i,]*unit.y, lwd=2)

}
