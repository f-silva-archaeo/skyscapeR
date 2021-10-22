#' Polar plot of orientations (azimuths)
#'
#' @param az (Optional) Array of azimuths. Can be omitted is \emph{obj} is given.
#' @param col (Optional) Single color or color palette to use for plotting measurements.
#' @param lwd (Optional) Line width to plot measurements. Defaults to 1.
#' @param lty (Optional) Line type to plot measurements. Defaults to 1.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the azimuths of celestial objects. Note that this assumes a single
#' location and a flat horizon of zero degrees.
#' @param show.obj.labels (Optional) Boolean to control whether to display celestial objects names.
#' Defaults to TRUE.
#' @export
#' @import utils stats plotrix
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' plotAzimuth(az)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- sky.objects('solar extremes', epoch=-2000, loc=c(35,-8), col='red')
#' plotAzimuth(az, obj=tt)
#'
#' # To display only celestial objects
#' plotAzimuth(obj=tt)
plotAzimuth = function(az, col='blue', lwd=1.5, lty=1, obj, show.obj.labels=T) {
  options(warn=-1)
  oldpar <- par('mar','mfrow')
  par(mar=c(1,1,1,1))

  if (!missing(obj)) {
    # Add celestial objects
    plotrix::polar.plot(0,0, radial.lim=c(0,5), lwd=lwd, line.col=col,
                        start=-270, clockwise = T,
                        labels = c('N', 'NE','E', 'SE', 'S', 'SW', 'W', 'NW'), label.pos = 45*seq(0,8)/180*pi,
                        show.grid.labels=F)

    for (i in 1:obj$n) {
      rise <- c(); set <- c();
      orb <- orbit(obj$decs[i], obj$loc, refraction=FALSE)
      forb <- splinefun(orb$az, orb$alt)
      rise <- uniroot(forb, interval=c(min(orb$az), 180))$root
      set <- uniroot(forb, interval=c(180, max(orb$az)))$root
      tt <- round(c(rise, set),1)

      sl <- TRUE
      for (j in 1:length(tt)) {
        plotrix::polar.plot(5, tt[j], radial.lim=c(0,5), lwd=obj$lwd[i], lty=obj$lty[i], line.col=obj$col[i],
                            start=-270, clockwise = T, add=T)

        if (show.obj.labels) {
          radial.plot.labels(5.9, tt[j], radial.lim=c(0,5), units="polar",
                             labels=colnames(obj$decs)[i], start=-270, clockwise = T, col=obj$col[i], cex=0.7)
        }
      }
    }
  }

  # Measurements
  if (!missing(az)) {
    if (class(az)=='data.frame') {
      names <- az$Name
      az <- az$True.Azimuth
    } else {
      names <- as.character(1:NROW(az))
    }

    if (length(col)==1) { col <- rep(col, length(az)) }

    plotrix::polar.plot(rep(5,length(az)), az, radial.lim=c(0,5), lwd=lwd, line.col=col,
                        start=-270, clockwise = T,
                        labels = c('N', 'NE','E', 'SE', 'S', 'SW', 'W', 'NW'), label.pos = 45*seq(0,7)/180*pi,
                        show.grid.labels=F, add=F+!missing(obj))

  }
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
  par(oldpar)
  options(warn=0)
}


#' Bar plot of orientations (declination)
#'
#' @param val Array of declination or azimuth values.
#' @param unc Single value or array of measurement uncertainty
#' @param names (Optional) Array of names of measurements in \emph{val}
#' @param unit (Optional). Either 'Declination' or 'Azimuth'. Defaults to 'Declination'.
#' @param col (Optional) Color to plot measurements in. Defaults to \emph{blue}.
#' @param shade (Optional) Boolean to control whether to shade a polygon of measurements. Defaults
#' to \emph{TRUE}
#' @param mark (Optional) Boolean to control whether to mark the declination value. Defaults to
#' \emph{TRUE}
#' @param sort (Optional) Boolean to control whether to sort the measurements by their declination
#' value. Defaults to \emph{FALSE}
#' @param xrange (Optional) Array of limits for x-axis.
#' @param yrange (Optional) Array of limits for y-axis.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination values of celestial objects.
#' @param show.obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the polar plot. Defaults to \emph{TRUE}.
#' @export
#' @import utils stats graphics
#' @seealso \code{\link{sky.objects}}
#' @examples
#' # Plot some declination data:
#' decs <- c(10, 12, -5, 4)
#' plotBars(decs, unc=5)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- sky.objects(c('solar extremes','lunar extremes'), epoch=-2000, lty=c(2,3))
#' plotBars(decs, unc=5, obj=tt)
plotBars <- function(val, unc, names, unit='Declination', col='blue', shade=TRUE, mark=TRUE, sort=FALSE, xrange, yrange, obj, show.obj.label=TRUE) {
  if (min(val)< -90 | max(val)>90 & unit=='Declination') {
    message('It appears that val includes azimuth values. Changing unit to Azimuth')
    unit <- 'Azimuth'}
  if (NROW(unc)==1) { unc <- rep(unc, NROW(val)) }
  if (sort) {
    ind <- sort(val, decreasing=TRUE, index.return=TRUE)$ix
    val <- val[ind]
    unc <- unc[ind]
    if (!missing(names)) { names <- names[ind] }
  }

  if (!missing(names)) { par(mar=c(4, 9, 2, 2) + 0.1) } else { par(mar=c(4, 2, 2, 2) + 0.1) }
  if (missing(xrange)) { xrange <- c(min(val-unc)-5, max(val+unc)+5) }
  if (missing(yrange)) { yrange <- c(0.5, NROW(val)+0.5) }
  plot.default(-100,-100, xlab=paste(unit,'(\u00b0)'), ylab='', xlim=xrange, ylim=yrange, axes=FALSE, yaxs='i')
  axis(1, at=pretty(seq(xrange[1],xrange[2])))
  axis(1, at=0, labels = 0)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 2) { axis(1, at=seq(-90,360,0.5), lwd=0.2, labels=FALSE) }
  if (scale <= 5 & scale > 1) { axis(1, at=seq(-90,360,1), lwd=0.5, labels=FALSE) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360,5), lwd=0.5, labels=FALSE) }
  if (scale > 10) { axis(1, at=seq(-90,360,10), lwd=0.5, labels=FALSE) }
  if (!missing(names)) { axis(2, at=1:NROW(names), labels=names, las=2, cex=0.7) } else { axis(2, at=1:NROW(val), las=2) }

  box()

  if (shade) { border <- NA } else { border <- col; col <- NA  }

  for (i in 1:NROW(val)) {
    xp <- c(val[i]-unc[i], val[i]+unc[i], val[i]+unc[i], val[i]-unc[i])
    yp <- c(i-0.4,i-0.4,i+0.4,i+0.4)
    polygon(xp, yp, col=MESS::col.alpha(col,0.5), border=border)
    if (mark) { lines(c(val[i],val[i]), c(i-0.4,i+0.4), col=col, lwd=2) }
  }

  # objects
  if (!missing(obj)) {
    if (unit == 'Azimuth') { message('Celestial objects can only currently be plotted when using declination values.') } else {
      for (i in 1:obj$n) {
        if (length(obj$epoch)==1) {
          abline(v=obj$decs[i], col=obj$col[i], lwd=obj$lwd[i], lty=obj$lty[i])
          if (show.obj.label) { text(obj$decs[i], .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
        } else {
          xp <- c(obj$decs[3:4,i], rev(obj$decs[3:4,i]))
          yp <- c(-1,-1,2,2)
          polygon(xp, yp, border=obj$col[i], col=MESS::col.alpha(obj$col[i],.3))
          if (show.obj.label) { text(mean(obj$decs[3:4,i]), .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
        }
      }
    }
  }
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
}




#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param x Object of \emph{skyscapeR.horizon} format.
#' @param show.az (Optional) Boolean that controls whether to display azimuth values or cardinal
#' directions. Defaults to FALSE.
#' @param xlim (Optional) Azimuth rage for plotting.
#' @param ylim (Optional) Altitude rage for plotting.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the paths of celestial objects.
#' @param refraction (Optional) Boolean switch controlling whether to take refraction into account
#' when displaying the paths of celestial objects.
#' @param col.ground (Optional) Color of the ground. Defaults to \emph{#fdae61}.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.horizon
#' @export
#' @import utils stats graphics grDevices
#' @seealso \code{\link{downloadHWT}}, \code{\link{sky.objects}}
#' @examples
#' # Plot a horizon retrieved from HeyWhatsThat:
#' hor <- downloadHWT('HIFVTBGK')
#' plot(hor)
#'
#' # Add the paths of the solstices and equinoxes sun in the year 1999 BC:
#' tt <- sky.objects('solar extremes', epoch=-2000, col='blue')
#' plot(hor, obj=tt)
plot.skyscapeR.horizon <- function(x, show.az=F, xlim, ylim, obj, refraction=F, col.ground='#fdae61', ...) {
  if (missing(xlim)) { xlim <- c(0,360) }
  if (missing(ylim)) { ylim <- c(floor(min(x$data$alt, na.rm=T))-5,45) }

  par(mar=c(2,1,1,1))
  plot(-99999,-99999, xlab = "", ylab = "", yaxs='i', xaxs='i', axes=F, lwd=5, xlim=xlim, ylim=ylim, ...)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-40,360+40,0.1), lwd=0.2, labels=F) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-40,360+40,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-40,360+40,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-40,360+40,5), lwd=0.5, labels=F) }
  if (scale < 90 & scale >= 10) { axis(1, at=seq(-40,360+40,10), lwd=0.5, labels=F) }

  if (show.az == T) {
    if (scale >= 10 & scale < 45) { scale <- 10 }
    if (scale >= 45 & scale < 90) { scale <- 45 }
    if (scale >= 90) { scale <- 90 }
    axis(1, at = seq(-90,360+90,scale), labels = seq(-90,360+90,scale), lwd=0)

  } else {

    ll <- c("N","NE","E","SE","S","SW","W","NW","N","NE","E","SE","S","SW","W","NW","N","NE","E","SE","S","SW","W","NW","N")
    axis(1, at = seq(-360,720,by=45), labels = ll, lwd=0.5)
  }

  # objects
  if (!missing(obj)) {
    ind <- sort(obj$decs[1,], decreasing=T, index.return=T)$ix
    for (i in ind) {
      if (length(obj$epoch)==1) {
        orb <- orbit(obj$decs[i], x, res=0.5, refraction=refraction)
        lines(orb$az, orb$alt, col=obj$col[i], lty=obj$lty[i], lwd=obj$lwd[i])
      } else {
        orb1 <- orbit(obj$decs[3,i], x, res=0.5)
        orb2 <- orbit(obj$decs[4,i], x, res=0.5)
        lines(orb1$az, orb1$alt, col=obj$col[i], lty=obj$lty[i], lwd=obj$lwd[i])
        lines(orb2$az, orb2$alt, col=obj$col[i], lty=obj$lty[i], lwd=obj$lwd[i])
      }
    }
  }

  # Horizon line
  line(x$data$az, x$data$alt)
  xx <- c(x$data$az, rev(x$data$az))
  yy <- c(x$data$alt, rep(-20,NROW(x$data$az)))
  polygon(xx, yy, col=col.ground, lwd=1.5)

  box()
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
}


#' Plot orientation probability distributions
#'
#' @param x A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param index (Optional) A value to indicate which distribution to plot, when only one is desired.
#' @param hdr (Optional) High density region to highlight. Defaults to 0.954
#' @param show.az (Optional) Whether to show the azimuth and transformation curve (horizon) when
#' displaying declination distributions. Defaults to TRUE.
#' @param xlim (Optional) Range of x-axis to plot.
#' @param col (Optional) Color of high density region to highlight. Defaults to blue with 0.5 alpha.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.pdf
#' @export
plot.skyscapeR.pdf <- function(x, index, hdr=0.954, show.az=T, xlim, col=MESS::col.alpha('blue',0.5), ...) {
  if (missing(index)) { index <- 1:NROW(x$data) }

  if (length(index)==1 & x$metadata$coord == 'dec') {
    # Coordinate Transformation Panels
    xx <- x$metadata$az.x$data[[index]]$x
    dens <- x$metadata$az.x$data[[index]]$y
    xl <- range(x$data[[index]]$x)
    yl <- range(xx);  yl[1] <- max(yl[1],0); yl[2] <- min(yl[2],360)

    layout(t(matrix(c(2,3,1,4), nrow=2)), widths = c(1,3), heights=c(3,1))

    # Information Panel
    par(mar=c(0,0,0,0))
    plot(-999,-999, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',axes=F, xlab='', ylab='', ...)
    hd <- hpdi(x$data[[index]], hdr)
    if (exists('name', where=x$metadata)) { text(.3, .8, labels=x$metadata$name[index], pos=4, font=2, cex=1.2)  }
    prob <- paste0(hdr*100, '% probability')
    for (j in 1:NROW(hd)) {
      prob <- paste0(prob,'\n   ',round(hd[j,1],2),'\u00b0  ::  ',round(hd[j,2],2),'\u00b0')
    }
    text(.3, .6, labels=prob, pos=4, font=1, cex=0.8)

    # Azimuth Distribution
    par(mar=c(0,4,1,0))
    plot(dens, xx, type='l', ylab='Azimuth (\u00b0)', xlab='', lwd=2, col='black', xlim=c(0, 1.2*max(dens)), ylim=yl, xaxs='i', yaxs='i', axes=F, ...); box()
    polygon(c(dens,rep(0,length(dens))), c(xx, rev(xx)), col='grey', border=NA)
    axis(2, at=pretty(seq(par('usr')[3],par('usr')[4])))
    axis(2, at=0, labels = 0)
    scale <- mean(diff(pretty(seq(par('usr')[3],par('usr')[4]))))
    if (scale <= 1) { axis(2, at=seq(-45,360+45,0.1), lwd=0.2, labels=F) }
    if (scale <= 2 & scale > 1) { axis(2, at=seq(-45,360+45,0.5), lwd=0.2, labels=F) }
    if (scale <= 5 & scale > 2) { axis(2, at=seq(-45,360+45,1), lwd=0.5, labels=F) }
    if (scale <= 20 & scale > 5) { axis(2, at=seq(-45,360+45,5), lwd=0.5, labels=F) }
    if (scale > 10) { axis(2, at=seq(-45,360+45,10), lwd=0.5, labels=F) }

    # Transformation Curve
    par(mar=c(0,0,1,1))
    plot(-999,999, axes=F, xlim=xl, ylim=yl, xaxs='i', yaxs='i', ...); box()

    ## process horizon profile
    hh <- x$metadata$horizon[[index]]$data
    alt <- approx(hh$az, hh$alt, xout=xx)$y
    alt.unc <- approx(hh$az, hh$alt.unc, xout=xx)$y
    refraction <- x$metadata$param$refraction
    atm <- x$metadata$param$atm
    temp <- x$metadata$param$temp

    dec0 <- az2dec(xx, x$metadata$horizon[[index]], alt, refraction=refraction, atm=atm, temp=temp)
    dec1 <- az2dec(xx, x$metadata$horizon[[index]], alt+alt.unc, refraction=refraction, atm=atm, temp=temp)
    dec2 <- az2dec(xx, x$metadata$horizon[[index]], alt-alt.unc, refraction=refraction, atm=atm, temp=temp)
    dec3 <- az2dec(xx, x$metadata$horizon[[index]], alt+2*alt.unc, refraction=refraction, atm=atm, temp=temp)
    dec4 <- az2dec(xx, x$metadata$horizon[[index]], alt-2*alt.unc, refraction=refraction, atm=atm, temp=temp)

    yp <- c(xx, rev(xx))
    xp <- c(dec3, rev(dec4))
    polygon(xp, yp, border=NA, col=MESS::col.alpha('blue', 0.25))
    yp <- c(xx, rev(xx))
    xp <- c(dec1, rev(dec2))
    polygon(xp, yp, border=NA, col=MESS::col.alpha('blue', 0.5))
    lines(dec0, xx, lwd=0.5, col='black')


    # Declination Distribution
    par(mar=c(4,0,0,1))
    xx <- x$data[[index]]
    plot(xx$x, xx$y, type='l', main='', xlab='Declination (\u00b0)', lwd=1.2, col='black', xlim=xl, ylim=c(0, 1.2*max(x$data[[index]]$y)), xaxs='i', yaxs='i', axes=F, ...); box()

    for (j in 1:NROW(hd)) {
      ind <- which(x$data[[index]]$x >= min(hd[j,]) & x$data[[index]]$x <= max(hd[j,]))
      xvals <- x$data[[index]]$x[ind]; xvals <- c(xvals, rev(xvals))
      yvals <- x$data[[index]]$y[ind]; yvals <- c(yvals, rep(0, length(yvals)))
      polygon(xvals, yvals, col=col, border=NA)

    }
    lines(x$data[[index]]$x, x$data[[index]]$y)

    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
    scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
    if (scale <= 1) { axis(1, at=seq(-90,90,0.1), lwd=0.2, labels=F) }
    if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
    if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
    if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
    if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0.5, labels=F) }

    box()
    mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)

  } else {
    ## Single or Multiple Distributions
    if (x$metadata$coord == 'dec') { lab <- 'Declination (\u00b0)'; xl <- c(-90,90) }
    if (x$metadata$coord == 'az') { lab <- 'Azimuth (\u00b0)'; xl <- c(0,360)}
    par(mar=c(4,1,1,1), mfrow=c(1,1))
    # if(!missing(xlim)) { xl <- xlim } else { xl <- range(hpdi(x, 1)) }
    plot(-9999, -9999, type='l', main='', xlab=lab, ylab='', lwd=2, col='black', xlim=xl, ylim=c(0,length(index)), xaxs='i', yaxs='i', axes=F, ...); box()
    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
    scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
    if (scale <= 1) { axis(1, at=seq(-90,360+50,0.1), lwd=0.2, labels=F) }
    if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,360+50,0.5), lwd=0.2, labels=F) }
    if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,360+50,1), lwd=0.5, labels=F) }
    if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360+50,5), lwd=0.5, labels=F) }
    if (scale > 10) { axis(1, at=seq(-90,360+50,10), lwd=0.5, labels=F) }

    for (i in 1:length(index)) {
      xx <- x$data[[index[i]]]
      hd <- hpdi(xx, hdr)
      for (j in 1:NROW(hd)) {
        ind <- which(xx$x >= min(hd[j,]) & xx$x <= max(hd[j,]))
        xvals <- xx$x[ind]; xvals <- c(xvals, rev(xvals))
        yvals <- xx$y[ind]; yvals <- yvals/max(yvals)*.8
        yvals <- c(yvals, rep(0, length(yvals)))
        polygon(xvals, yvals + (i-1) + 0.1, col=col, border=NA)
      }
      ind <- which(diff(xx$x) > 5*mean(diff(xx$x)))
      if (length(ind)>0) { xx[ind,2] <- NA }
      lines(xx$x, xx$y/max(xx$y, na.rm=T)*.8 + (i-1) + 0.1, lwd=1.2, col='black')
      # text(par('usr')[1] - 20, par('usr')[3]+i-0.5, labels=x$metadata$name[index[i]], xpd=T)
      text(max(hd), par('usr')[3]+i-0.5, labels=x$metadata$name[index[i]], xpd=T, pos=4)
      abline(h=i, col='grey', lty=3)
      mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
    }

  }
}


#' Plot orientation summed probability density
#'
#' @param x A \emph{skyscapeR.spd} object created with \code{\link{spd}}
#' @param xlim (Optional) Range of x-axis to plot.
#' @param ylim (Optional) Range of y-axis to plot.
#' @param title (Optional) Title to add to the plot.
#' @param col (Optional) Color of summed probability density. Defaults to blue with 0.5 alpha.
#' @param shading (Optional) Boolean to control whether to color the entire distribution. Defaults to TRUE
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.spd
#' @export
plot.skyscapeR.spd <- function(x, xlim, ylim, title=NULL, col='blue', shading=T, ...) {
  par(mar=c(5, 4, 1, 1) + 0.1)
  if (missing(xlim)) { xlim <- sort(x$data$x[c(min(which(x$data$y >= 1e-12)), max(which(x$data$y >= 1e-12)))]) }
  if (missing(ylim)) { ylim <- c(0, max(x$data$y)) }

  if (x$metadata$coord == 'dec') { label <- 'Declination (\u00b0)' }
  if (x$metadata$coord == 'az') {
    label <- 'Azimuth (\u00b0)'
    xlim[1] <- max(xlim[1], 0)
    xlim[2] <- min(xlim[2], 360)
  }

  plot(x$data$x, x$data$y, type='l', main=title, xlab=label, ylab='Density', lwd=2, col=col, xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', axes=F, ...); box()
  if (shading) {
    xp <- c(x$data$x, rev(x$data$x))
    yp <- c(x$data$y, rep(0, length(x$data$x)))
    polygon(xp, yp, col=MESS::col.alpha(col,0.5), border=NA)
  }
  if (x$metadata$coord == 'dec') {
    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
  } else if (x$metadata$coord == 'az') {
    axis(1, at=c(0, 45, 90, 135, 180, 225, 270, 315, 360))
  }
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-90,360,0.1), lwd=0.2, labels=F) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,360,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,360,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,360,10), lwd=0.5, labels=F) }
  axis(2)
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)

}


#' Plot significance test of orientations
#'
#' @param x A \emph{skyscapeR.sigTest} object created with \code{\link{randomTest}}
#' @param xlim (Optional) Range of x-axis to plot.
#' @param title (Optional) Title to add to the plot.
#' @param show.pval (Optional) Boolean to control whether to print the global p-value. Default is TRUE.
#' @param show.local (Optional) Boolean to control whether to show local regions of significance. Default is FALSE
#' @param pal (Optional) Color palette for local regions of significance. Default is brewer.pal(5, 'PRGn')
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.sigTest
#' @export
plot.skyscapeR.sigTest <- function(x, xlim, title=NULL, show.pval=T, show.local=F, pal=brewer.pal(5, 'PRGn')[c(1,5)], ...) {
  # empirical spd
  spd <- x$data$empirical
  if (missing(xlim)) { xlim <- sort(spd$x[c(min(which(spd$y >= 1e-12)), max(which(spd$y >= 1e-12)))]) }
  if (show.local) { ylim <- c(-max(spd$y)*.05,max(spd$y)) } else { ylim <- c(0,max(spd$y))}
  if (x$metadata$coord == 'dec') { label <- 'Declination (\u00b0)' }
  if (x$metadata$coord == 'az') {
    label <- 'Azimuth (\u00b0)'
    xlim[1] <- max(xlim[1], 0)
    xlim[2] <- min(xlim[2], 360)
  }
  plot(spd$x, spd$y, type='l', main=title, xlab=label, ylab='Density', lwd=2, col='blue', xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', axes=F, ...); box()
  xp <- c(spd$x, rev(spd$x))
  yp <- c(spd$y, rep(0, length(spd$x)))
  polygon(xp, yp, col=MESS::col.alpha('blue',0.5), border=NA)
  if (x$metadata$coord == 'dec') {
    axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
    axis(1, at=0, labels = 0)
  } else if (x$metadata$coord == 'az') {
    axis(1, at=c(0, 45, 90, 135, 180, 225, 270, 315, 360))
  }
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-90,360,0.1), lwd=0.2, labels=F) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,360,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,360,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,360,10), lwd=0.5, labels=F) }
  axis(2)

  # sigtest
  # lines(spd$x, sig$data$null.hyp$CE.mean, col='grey')
  xp <- c(spd$x, rev(spd$x))
  if (x$metadata$tails==1) {
    yp <- c(x$data$null.hyp$CE.upper, rep(0, length(spd$x)))
  } else {
    yp <- c(x$data$null.hyp$CE.upper, rev(x$data$null.hyp$CE.lower))
  }
  polygon(xp, yp, col=MESS::col.alpha('grey', 0.5), border=NA)

  # global p-value
  if (show.pval) {
    if (x$metadata$global.pval == 0 ) {
      pval <- paste0("global p-value < ", round(1/(x$metadata$nsims+1),4))
    } else {
      pval <- paste0('global p-value = ', x$metadata$global.pval)
    }
    text(par('usr')[2], abs(diff(par('usr')[3:4]))*.90, pos=2, pval, cex=1.2, font=2)
  }

  # regions of significance
  if (show.local) {
    abline(0,0, lwd=1)
    aux <- as.matrix(x$metadata$local.pval[,2:4])
    if (NROW(aux)>0) {
      for (i in 1:NROW(aux)) {
        if (x$metadata$local.pval[i,1] == '+') { col <- pal[2] } else { col <- pal[1] }
        if (aux[i,1] < xlim[1] & aux[i,2] > xlim[1]) { aux[i,1] <- xlim[1] }
        if (aux[i,2] > xlim[2] & aux[i,1] < xlim[2]) { aux[i,2] <- xlim[2] }
        xp <- c(aux[i,1], aux[i,1], aux[i,2], aux[i,2])
        yp <- c(ylim[1], 0, 0, ylim[1])
        polygon(xp, yp, col=MESS::col.alpha(col, .4), border=NA)
        text(mean(aux[i,1:2]), ylim[1]/2, labels=pval2stars(aux[i,3]), cex=0.9, srt=90)
      }
    }
  }

  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)

}


#' Plot stellar phase and seasonality
#'
#' This function creates a plot of stellar seasonality and phases/events.
#' @param x Object of \emph{skyscapeR.starphases} format.
#' @param ... Additional arguments to be passed to \emph{plot}.
#' @rdname plot.skyscapeR.starphases
#' @export
#' @import MESS RColorBrewer
#' @seealso \code{\link{star.phases}}
#' @examples
#' # Plot the seasonality of Aldebaran for 3999 BCE:
#' \dontrun{
#' ss <- star.phases('Aldebaran',-4000, c(35,-8, 200))
#' plot(ss)
#' }
plot.skyscapeR.starphases = function(x, ...) {

  col <- RColorBrewer::brewer.pal(4,'Accent'); col[5] <- col[4]
  code <- c('RS', 'R', 'S', 'V', 'I')
  seasons <- c("Rise and Set","Rise Only","Set Only","","")
  events <- c('acronycal rising','heliacal setting','acronycal setting','heliacal rising')
  if (x$metadata$type == 'curtailed passage' | x$metadata$type == 'dual phase') {
    seasons[4] <- 'Curtailed Passage'
  } else if (x$metadata$type == 'arising and lying hidden' | x$metadata$type == 'dual phase') {
    seasons[5] <- 'Arising and Lying Hidden'
  } else {
    seasons[4] <- 'Circumpolar'
    seasons[5] <- 'Invisible'
  }

  par(mar=c(3,1,1,1), mgp=c(3,1.5,0))
  plot(-100,-100, xlim=c(1,365), ylim=c(0,1), main=paste0(x$metadata$star$name,' at ', BC.AD(x$metadata$year)), xlab="", ylab="", axes=FALSE, ...)
  axis(1, at=c(0,91,182,273,365), labels=c('December\nSolstice','March\nEquinox','June\nSolstice','September\nEquinox','December\nSolstice'), lwd.ticks=2)
  axis(1, at=seq(0,365,30.4), labels=NA)

  for (i in 1:length(code)) {
    ind <- which(x$data$phase == code[i])
    ind.bb <- split(ind,cumsum(c(1,abs(diff(ind))>1)))

    for (j in 1:NROW(ind.bb)) {
      ind.i <- ind.bb[[j]]
      if (length(ind.i)>0) {
        x.poly <- c(ind.i[1]-.5,tail(ind.i,1)+.5,tail(ind.i,1)+.5,ind.i[1]-.5)
        y.poly <- c(0,0,1,1)
        polygon(x.poly, y.poly, col=col[i], border=NA)
        text(mean(ind.i),0.5, seasons[i], font=2, srt=90)
      }
    }

  }

  if (x$metadata$type == 'curtailed passage' | x$metadata$type == 'arising and lying hidden' | x$metadata$type == 'dual phase') {
    ind <- which(x$data$phase == code[4])
    if (length(ind) > 0) {
      lines(rep(min(ind)-1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(min(ind)-1,0.02, events[1], cex=0.8, font=1, srt=90 ,pos=2)

      lines(rep(max(ind)+1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(max(ind)+1,0.02, events[2], cex=0.8, font=1, srt=90 ,pos=4)
    }

    ind <- which(x$data$phase == code[5])
    if (length(ind) > 0) {
      lines(rep(min(ind)-1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(min(ind)-1,0.02, events[3], cex=0.8, font=1, srt=90 ,pos=4)

      lines(rep(max(ind)+1,2), c(0,1), lty=2, lwd=1.5, col='black')
      text(max(ind)+1,0.02, events[4], cex=0.8, font=1, srt=90 ,pos=4)
    }
  }
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
}


#' Prints significance test results
#'
#' This function prints the results of \code{\link{randomTest}}.
#' @param x Object of \emph{skyscapeR.sigTest} format.
#' @param ... Additional arguments to be passed to \emph{print}.
#' @export
#' @seealso \code{\link{randomTest}}
print.skyscapeR.sigTest <- function (x, ...) {
  cat("\n*** Results of Significance Test ***\n\n")
  cat(paste0(x$metadata$tails,'-tailed test at ',x$metadata$conf*100,'% confidence, based on ', x$metadata$nsims, ' simulations.\n'))

  if (x$metadata$global.pval == 0) {
    p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3))
  } else { p.value <- x$metadata$global.pval }
  cat(paste0('global p-value: ', p.value, ' (',pval2stars(p.value),')\n'))
  cat('local p-values:\n')
  for (i in 1:NROW( x$metadata$local.pval)) {
    if (x$metadata$local.pval[i,]$p.value == 0) {
      p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3))
    } else { p.value <- x$metadata$local.pval[i,]$p.value }
    cat(paste0('      ',x$metadata$local.pval[i,]$type,'  dec range [',round(x$metadata$local.pval[i,]$start,2), ', ',round(x$metadata$local.pval[i,]$end,2),'] :: p-value: ', p.value, ' (',pval2stars(p.value),')\n'))
  }
}
