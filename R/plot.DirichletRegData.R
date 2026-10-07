get_or_else <- function(what, or_else, where) if(any(what %in% names(where))) where[[what, exact = TRUE]] else or_else



plot.DirichletRegData <- function(
  x,
  dims,                                       # which dimensions to plot
  ticks              = TRUE,                  # plot ternary ticks?
  ref.lines          = NULL,                  # reference lines for 2d and 3d plots?
  dim.labels,
  a2d = list(
    colored          = TRUE,
    c.grid           = TRUE,                  # plot a grid?
    col.scheme       = c("dims", "entropy"),  # if colors: which scheme?
    entropy.contours = FALSE,                 # plot entropy-contour lines?
    entropy.colors   = FALSE                  # if entropy-contours: plot colored regions?
  ),
  a3d                = list(rgl = TRUE, ...), # theta and phi for the viewport
  rug                = TRUE,
  reset_par          = TRUE,
  ...
){

## ADAPT TO THE NEW DATA STRUCTURE
  full_obs <- nrow(x)
  if(anyNA(x)){
    include <- which(rowSums(is.na(x)) == 0L)
    .x <- x

    x <- as.matrix(x[include, ])
    attr(x, "Y.original") <- as.data.frame(attr(.x, "Y.original")[include, ])
    attr(x, "dims") <- attr(.x, "dims")
    attr(x, "obs") <- nrow(x)
    attr(x, "valid_obs") <- nrow(x)
    attr(x, "normalized") <- attr(.x, "normalized")
    attr(x, "transformed") <- attr(.x, "transformed")
    attr(x, "base") <- attr(.x, "base")
    class(x) <- "DirichletRegData"

    rm(.x)
  } else {
    include <- seq_len(nrow(x))
  }

  nx <- x
  x <- attributes(x)
  x$Y <- as.data.frame(unclass(nx))
  x$dim.names <- x$dimnames[[2L]]
  class(x) <- "DirichletRegData"

  if(missing(dims)) dims <- NULL

  if(!inherits(x, "DirichletRegData")) stop("data must be prepared by 'DR_data()'")

  colored <- get_or_else("colored", TRUE, a2d)
  c.grid <- get_or_else("c.grid", TRUE, a2d)
  col.scheme <- match.arg(get_or_else("col.scheme", "dims", a2d), c("dims", "entropy"))
  entropy.contours <- get_or_else("entropy.contours", FALSE, a2d)
  entropy.colors <- get_or_else("entropy.colors", FALSE, a2d)

  dotlist <- list(...)

  .main <- get_or_else("main", NULL, dotlist)
  .xlim <- get_or_else("xlim", NULL, dotlist)
  .ylim <- get_or_else("ylim", NULL, dotlist)
   .col <- get_or_else("col", NULL, dotlist)  ; if(length(.col) == full_obs) .col <- .col[include]
   .pch <- get_or_else("pch", 16L, dotlist)   ; if(length(.pch) == full_obs) .pch <- .pch[include]
   .cex <- get_or_else("cex", 1.0, dotlist)   ; if(length(.cex) == full_obs) .cex <- .cex[include]
   .lwd <- get_or_else("lwd", 1.0, dotlist)
   .lty <- get_or_else("lty", 1L, dotlist)

      theta <- get_or_else("theta", NULL, dotlist)
        phi <- get_or_else("phi", NULL, dotlist)
  ref.lines <- get_or_else("ref.lines", NULL, dotlist)

  .marginal <- FALSE

  if(is.null(dims)){
    if(x$dims > 4L){
      x$Y         <- x$Y[, seq_len(4)]
      x$dims      <- 4L
      x$dim.names <- x$dim.names[seq_len(4)]
      warning("data contains > 4 variables. the first four are being used. to change this, set the 'dims' argument appropriately.")
      .marginal <- TRUE
    } else {
      dims <- seq_len(x$dims)
    }
  } else {
    if((length(dims) < 2L) || (length(dims) > 4L)) stop("the argument 'dims' must have 2, 3 or 4 elements")
    x$Y         <- x$Y[, dims]
    x$dims      <- length(dims)
    x$dim.names <- x$dim.names[dims]
    .marginal <- TRUE
  }



  if(missing(dim.labels)){
    dim.labels <- x$dim.names
  }





  if(x$dims == 2L){
    if(is.null(.main)) .main <- "Density Plot of a Beta-Distributed Variable"
    plot_DRdata_2d(y = x$Y[, 2L], rug = rug, main = .main, ylim = .ylim, colr = .col, lwd = .lwd, lty = .lty)

  } else if(x$dims == 3L) {
    if(!all(is.null(c(.xlim, .ylim)))) warning("xlim and ylim not useable in a ternary plot. arguments ignored.")
    if(is.null(.main)) .main <- "Ternary Plot"

    if(reset_par){ # reset the current pars after plotting
      old.par <- par(no.readonly = TRUE)
      on.exit(par(old.par))
    }

    plot_DRdata_3d(x = x, entropy.contours = entropy.contours, colored = colored, c.grid = c.grid, ticks = ticks, dim.labels = dim.labels, col.scheme = col.scheme,
    .main = .main, .col = .col, .pch = .pch, .cex = .cex, .lwd = .lwd, .lty = .lty)

  } else if(x$dims == 4L){
    if( !("rgl" %in% loadedNamespaces()) && !("rgl" %in% utils::installed.packages()[, "Package"]) ){ stop('The "rgl" package could not be found. You can try installing it using:\ninstall.packages("rgl")') }
    requireNamespace("rgl")

    plot_DRdata_4d(x = x, dim.labels = dim.labels, ref.lines = ref.lines,
    main = .main, cex = .cex,
    args.3d = a3d, theta = theta, phi = phi
    )
  }

}
