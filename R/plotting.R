devtools::use_package("plotrix")

#' Polar plot of orientations (azimuths)
#'
#' This function creates a polar plot of azimuthal data. It is a wrapper for
#' \code{\link[plotrix]{polar.plot}}
#' @param az Array of azimuths. Values outside the [0, 360] range will be ignored.
#' @param objects (Optional) A skyscapeR.object object created with \code{\link{object}}
#' for displaying the azimuths of clestial objects. Beware that this assumes a single
#' location (given by parameter loc) and a flat horizon of zero degrees.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a skyscapeR.horizon object. Only necessary for plotting potential
#' celestial targets.
#' @param ... Any other parameters to be passed unto \code{\link[plotrix]{polar.plot}}
#' @export
#' @seealso \code{\link[plotrix]{polar.plot}}, \code{\link{object}}
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' plotAz(az)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- object(c('sun','moon'), epoch=-2000, lty=c(2,3))
#' plotAz(az, tt, loc=c(35,-8))
plotAz = function(az, objects, loc, ...) {

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
  par(mar=c(4, 4, 2, 2) + 0.1)
  plotrix::polar.plot(testlen, testpos, lwd=1.2, line.col='black', start=90, clockwise=T, labels= labels, label.pos = label.pos, show.grid.labels= F, ...)
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # objects
  if (!missing(objects) & !missing(loc)) {
    rise <- c(); set <- c()
    for (i in 1:objects$n) {
      orb <- orbit(objects$decs[i], loc)
      forb <- splinefun(orb$az, orb$alt)
      rise[i] <- uniroot(forb, interval=c(0, 180))$root
      set[i] <- uniroot(forb, interval=c(180, 360))$root
    }
    tt <- c(rise, set)
    plotrix::polar.plot(c(0.01,rep(1,NROW(tt))), c(0,tt), lwd=c(0,rep(objects$lwd,2)), line.col=c(0,rep(objects$col,2)), lty=c(1,rep(objects$lty,2)), start=90, clockwise=T, add=T)
    plotrix::radial.plot.labels(c(0.01,rep(1,NROW(rise))), c(0,rise), units='polar', start=pi*90/180, clockwise=T, labels=c('', colnames(objects$decs)), cex=0.6, col=c('black',objects$col), pos=4, offset=0.5)
    plotrix::radial.plot.labels(c(0.01,rep(1,NROW(set))), c(0,set), units='polar', start=pi*90/180, clockwise=T, labels=c('', colnames(objects$decs)), cex=0.6, col=c('black',objects$col), pos=2, offset=0.5)
  }
}

#' Plot a curvigram
#'
#' This function creates a plot of a curvigram
#' @param curv Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.
#' @param objects (Optional) A skyscapeR.object object created with \code{\link{object}}
#' for displaying the dclination of clestial objects.
#' @param xlim Array of two values restricting the horizontal range of the plot.,
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @seealso \code{\link{curvigram}}, \code{\link{object}}
#' @examples
#' # Plot the curvigram of Recumbent Stone Circles:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, sd=2)
#' plot(curv, xlim=c(-40,0))
#'
#' # Redo the plot to include lunar extreme declinations:
#' LEx <- objects(c('smnLX','sMjLX'), -2000, col='red', lty=2)
#' plot(curv, LEx, xlim=c(-40,0))
plot.skyscapeR.curv = function(curv, objects, objects.label=T, xlim=NULL, ...) {
  par(mar=c(4, 4, 2, 2) + 0.1)
  if (is.null(xlim)) { xlim <- c(min(curv$dec),max(curv$dec)) }
  plot.default(-100,-100, xlab='Declination', ylab='Density', xlim=xlim, ylim=c(0,max(curv$density)), axes=F, ...)
  axis(1); axis(2)
  lines(curv$dec, curv$density, lwd=1.5, col='blue')
  box()
  mtext(paste0('skyscapeR ',packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # objects
  if (!missing(objects)) {
    for (i in 1:objects$n) {
      abline(v=objects$dec[i], col=objects$col[i], lwd=objects$lwd[i], lty=objects$lty[i])
      if (objects.label) {
        text(objects$dec[i], .95*par('usr')[4], colnames(objects$decs)[i], col=objects$col[i], pos=4, offset=0.2, cex=0.7)
      }
    }
  }
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
#' @param objects (Optional) A skyscapeR.object object created with \code{\link{object}}
#' for displaying the paths of clestial objects.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @seealso \code{\link{download.HWT}}, \code{\link{object}}
#' @examples
#' # Plot a horizon retrieved from HeyWhatsThat:
#' hor <- download.HWT('HIFVTBGK')
#' plot(hor)
#'
#' # Add the paths of the ssolstice sun in the year 1999 BC:
#' tt <- object(c('jS','dS'), -2000, 'blue')
#' plot(hor, objects=tt)
plot.skyscapeR.horizon <- function(hor, show.az=F, max.alt, az0 = 0, zoom=F, objects, ...) {
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

  par(mar=c(2,1,1,1))
  plot(-99999,-99999, xlab = xx, ylab = "", yaxs='i', xaxs='i', ylim = yl, xlim = xl, axes=F, lwd=5, ...)
  if (show.az == T) {
    axis(1, at = at, labels = ll)
  } else {
    axis(1, at = at, labels = ll)
    axis(1, at = seq(-495,720,by=90), labels = F, tcl=-0.25)
  }
  # axis(2, at = seq(0,yl[2],by=10))
  box()
  mtext(paste0('skyscapeR ',packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # objects
  if (!missing(objects)) {
    for (i in 1:objects$n) {
      orb <- orbit(objects$dec[i], hor, res=0.5)
      plot(orb, objects$col[i])
    }
  }

  # Horizon line
  line(hor$az, hor$alt)
  x <- c(hor$az, rev(hor$az))
  y <- c(hor$alt, rep(-20,NROW(hor$az)))
  polygon(x ,y, col=rgb(217/255,95/255,14/255,1))

}


#' Plot visible path of celestial object on horizon
#'
#' This function adds the visible path of a celestial object to a horizon plot.
#' @param orbit Object of \emph{skyscapeR.orbit} format.
#' @param col String with colour to plot path in. Defaults to red.
#' @seealso \code{\link{plot.skyscapeR.horizon}}, \code{\link{orbit}}
#' @examples
#' hor <- download.HWT('HIFVTBGK')
#' plot(hor)
#'
#' path <- orbit(dS(-2500),hor)
#' plot(path, col='blue')
#' @noRd
plot.skyscapeR.orbit <- function(orbit, col) {

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

