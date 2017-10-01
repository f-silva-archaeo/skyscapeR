devtools::use_package("plotrix")

#' Polar plot of orientations (azimuths)
#'
#' This function creates a polar plot of azimuthal data. It is a wrapper for
#' \code{\link[plotrix]{polar.plot}}
#' @param az Array of azimuths. Values outside the [0, 360] range will be ignored.
#' @param ... Any other parameters to be passed unto \code{\link[plotrix]{polar.plot}}
#' @export
#' @seealso \code{\link[plotrix]{polar.plot}}
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' orplot(az)
orplot = function(az, ...) {
  require(plotrix)

  ind <- which(az > 360 | az < 0)
  if (length(ind) > 0) {
    az <- az[-ind]
    warning('Ignoring values outside of azimuth range [0, 360]')
    }

  n <- NROW(az)

  testlen <- c(0.01,rep(1,n))
  testpos <- c(0,az)
  label.pos <- c(0,45,90,135,180,225,270,315)
  labels <- c("N","","E","","S","","W","")
  plotrix::polar.plot(testlen,testpos,main="Orientation Plot",lwd=1.2,line.col='black', start=90, clockwise=T, labels= labels, label.pos = label.pos, show.grid.labels= F, ...)
}

#' Plot a curvigram
#'
#' This function creates a plot of a curvigram
#' @param curv Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @seealso \code{\link{curvigram}}
#' @examples
#' # Plot some azimuth data:
#' decs <- data(RugglesRSC)
#' curv <- curvigram(decs, sd)
#' plot(curv)
plot.skyscapeR.curv = function(curv, xlim=NULL, ...) {
  if (is.null(xlim)) { xlim <- c(min(curv$dec),max(curv$dec)) }
  plot.default(-100,-100, xlab='Declination', ylab='Density', xlim=xlim, ylim=c(0,max(curv$density)), axes=F, ...)
  axis(1); axis(2)
  lines(curv$dec, curv$density, lwd=1.5, col='blue')
}


#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param hor Object of \emph{skyscapeR.horizon} format.
#' @param show.az Boolean that controls whether to display azimuth values on horizontal axis.
#'  Defaults to FALSE.
#' @param max.alt Maximum altitude to display. Defaults to 45 degrees.
#' @param az0 Leftmost azimuth of plot. Defaults to 0 degrees, i.e. North at the left.
#' @param zoom Boolean that controls whether to provide a zoomed-in view of 100 degrees in
#' azimuth and 5 degrees of altitude above the horizonline. Defaults to FALSE.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @seealso \code{\link{download.HWT}}
#' @examples
#' hor <- download.HWT('HIFVTBGK')
#' plot(hor)
plot.skyscapeR.horizon <- function(hor, show.az=F, max.alt, az0 = 0, zoom=F, ...) {
  # rejiggle so plot starts at given start point
  if (az0 < -360) { az0 <- az0 + 360 }
  ind <- which(hor$az < az0)
  if (NROW(ind)>0) {
    hor$az <- c(hor$az, hor$az[ind] + 360)
    hor$alt <- c(hor$alt, hor$alt[ind])
  }

  yl <- floor(min(hor$alt, na.rm=T))-5
  if (missing(max.alt)) { yl[2] <- 45 } else { yl[2] <- max.alt }

  if (zoom == T) {
    xl <- c(az0, az0+100)
    yl <- c(-5,5)
    show.az <- T
  } else {
    xl <- c(az0, az0 + 360)
  }

  if (show.az == T) {
    xx <- "azimuth"
    ll <- seq(-360,720,by=10)
    at <- ll
  } else {
    xx <- ""
    ll <- c("N","E","S","W","N","E","S","W","N","E","S","W","N")
    at <- seq(-360,720,by=90)
  }

  plot.default(hor$az, hor$alt, type='l', xlab = xx, ylab = "", yaxs='i', xaxs='i', ylim = yl, xlim = xl, axes=F, lwd=5, ...)
  x <- c(hor$az, rev(hor$az))
  y <- c(hor$alt, rep(-20,NROW(hor$az)))
  polygon(x ,y, col=rgb(217/255,95/255,14/255,1))

  if (show.az == T) {
    axis(1, at = at, labels = ll)
  } else {
    axis(1, at = at, labels = ll)
    axis(1, at = seq(-495,720,by=90), labels = F, tcl=-0.25)
  }
  # axis(2, at = seq(0,yl[2],by=10))
  box()
}


#' Plot visible path of celestial object on horizon
#'
#' This function adds the visible path of a celestial object to a horizon plot.
#' @param orbit Object of \emph{skyscapeR.orbit} format.
#' @param col String with colour to plot path in. Defaults to red.
#' @export
#' @seealso \code{\link{plot.skyscapeR.horizon}}, \code{\link{orbit}}
#' @examples
#' hor <- download.HWT('HIFVTBGK')
#' plot(hor)
#'
#' path <- orbit(dS(-2500),hor)
#' plot(path, col='blue')
plot.skyscapeR.orbit <- function(orbit, col, ...) {
  if (missing(col)) { col <- "red" }

  # TO DO: if no horizon displayed assume flat horizon
  options(warn=2)
  test <- try (par(new=T), silent=T)
  options(warn=0)
  if (class(test)=='try-error') {
    hor <- createHor(c(0,180,360), c(0,0,0), orbit$georef[1], orbit$georef[2], 'Flatland')
    plot(hor)
  }

  ind.break <- which(abs(diff(orbit$az)) > 1)
  if (NROW(ind.break) > 0) {
    for (i in 1:(NROW(ind.break)+1)) {
      orb1 <- c()
      xx <- c()
      if (i==1) { xx <- seq(1,ind.break[i],1) }
      if (i>1 & i<=NROW(ind.break)) {
        xx <- seq(ind.break[i-1]+1, ind.break[i], 1)
      }
      if (i>NROW(ind.break)) {
        xx <- seq(ind.break[i-1]+1, NROW(orbit$az), 1)
      }
      orb1$az <- orbit$az[xx]
      orb1$alt <- orbit$alt[xx]
      lines(orb1$az, orb1$alt, col=col)
      lines(orb1$az - 360, orb1$alt, col=col)
      lines(orb1$az + 360, orb1$alt, col=col)
    }
  } else {
    lines(orbit$az, orbit$alt, col=col)
    lines(orbit$az - 360, orbit$alt, col=col)
    lines(orbit$az + 360, orbit$alt, col=col)
  }
}

