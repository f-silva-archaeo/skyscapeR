## Correction for plotrix, sent off to package maintainer
#' @noRd
radial.plot.labels <- function (lengths, radial.pos = NULL, units = "radians", radial.lim = NULL,
                                start = 0, clockwise = FALSE, labels, adj = NULL, pos = NULL,
                                ...)
{
  npoints <- length(lengths)
  if (is.null(radial.pos))
    radial.pos <- seq(0, pi * (2 - 2/npoints), length.out = npoints)
  else {
    if (units == "clock24") {
      radial.pos <- pi * (450 - radial.pos * 15)/180
      start <- pi * (450 - start * 15)/180
    }
    if (units == "polar"){
      radial.pos <- pi * radial.pos/180
      start <- pi * start/180
    }
  }
  if (clockwise && units != "clock24")
    radial.pos <- -radial.pos
  if (start && units != "clock24")
    radial.pos <- radial.pos + start
  if (is.null(radial.lim))
    radial.lim <- range(lengths)
  lengths <- lengths - radial.lim[1]
  xpos <- cos(radial.pos) * lengths
  ypos <- sin(radial.pos) * lengths
  text(x = xpos, y = ypos, labels = labels, adj = adj, pos = pos,
       ...)
}




#' Polar plot of orientations (azimuths)
#'
#' This function creates a polar plot of azimuthal data using \emph{plotrix}.
#' @param az Array of azimuths.
#' @param col (Optional) Single colour or colour pallete to use for plotting measurements.
#' @param lwd (Optional) Line width to plot measurements. Defaults to 1.
#' @param lty (Optional) Line type to plot measurements. Defaults to 1.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the azimuths of celestial objects. Beware that this assumes a single
#' location (given by parameter loc) and a flat horizon of zero degrees.
#' @param show.obj.labels (Optional) Boolean to control whether to display celestial objects names.
#' Defaults to TRUE.
#' @export
#' @import utils stats plotrix
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' plotAzimuth(az)
#'
#' # To visualise this data against the common solar and lunar targets:
#' tt <- sky.objects('solar extremes', epoch=-2000, loc=c(35,-8), col='red')
#' plotAzimuth(az, obj=tt)
#'
#' # To display only celestial objects
#' plotAzimuth(az=NULL, obj=tt)
plotAzimuth = function(az, col='blue', lwd=1, lty=1, obj, show.obj.labels=T) {
  options(warn=-1)
  oldpar <- par('mar','mfrow')
  par(mar=c(1,1,1,1))

  if (class(az)=='data.frame') {
    names <- az$Name
    az <- az$True.Azimuth
  } else {
    names <- as.character(1:NROW(az))
  }
  if (length(col)==1) { col <- rep(col, length(az)) }

  if (!missing(obj)) {
    # Add celestial objects
    plotrix::polar.plot(0,0, radial.lim=c(0,5), lwd=lwd, line.col=col,
                        start=-270, clockwise = T,
                        labels = c('N', 'NE','E', 'SE', 'S', 'SW', 'W', 'NW'), label.pos = 45*seq(0,8),
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
  if (!is.null(az)) {
    plotrix::polar.plot(rep(5,length(az)), az, radial.lim=c(0,5), lwd=lwd, line.col=col,
                        start=-270, clockwise = T,
                        labels = c('N', 'NE','E', 'SE', 'S', 'SW', 'W', 'NW'), label.pos = 45*seq(0,8),
                        show.grid.labels=F, add=F+!missing(obj))

  }
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
  par(oldpar)
  options(warn=0)
}


#' Bar plot of orientations (declinations)
#'
#' This function creates a bar plot of orientation data.
#' @param val Array of declinations or azimuths.
#' @param unc Single value or array of measurement uncertainty
#' @param names (Optional) Array of names of measurements in \emph{val}
#' @param unit (Optional). Either 'Declination' or 'Azimuth'. Defaults to 'Declination'.
#' @param col (Optional) Colour to plot measurments in. Defaults to \emph{blue}.
#' @param shade (Optional) Boolean to control whether to shade a polygon of measurements. Defaults
#' to \emph{TRUE}
#' @param mark (Optional) Boolean to control whether to mark the declination value. Defaults to
#' \emph{TRUE}
#' @param sort (Optional) Boolean to control whether to sort the measurements by their declination
#' value. Defaults to \emph{FALSE}
#' @param xrange (Optional) Array of limits for x-axis.
#' @param yrange (Optional) Array of limits for y-axis.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declinations of celestial objects.
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
  plot.default(-100,-100, xlab=unit, ylab='', xlim=xrange, ylim=yrange, axes=FALSE, yaxs='i')
  axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
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




#' Plot a curvigram
#'
#' This function creates a plot of a curvigram
#' @param cc Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.
#' @param col (Optional) Colour to plot the curvigram in. Defaults to blue.
#' @param shading (Optional) Whether to shade the curvigram. Defaults to TRUE.
#' @param xrange (Optional) Array of two values restricting the horizontal range of the plot.
#' @param yrange (Optional) Array of two values restricting the vertical range of the plot.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param show.obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the polar plot. Defaults to \emph{TRUE}.
#' @export
#' @seealso \code{\link{curvigram}}, \code{\link{histogram}}, \code{\link{sky.objects}}
#' @examples
#' # Plot the curvigram of Recumbent Stone Circles:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, unc=2)
#' plot(curv)
#'
#' # Redo the plot to include lunar extreme declinations:
#' LEx <- sky.objects('lunar extremes', -2000, col='red', lty=2)
#' plot(curv, obj=LEx)
plot.skyscapeR.curv <- function(cc, col='blue', shading=T, xrange, yrange, obj, show.obj.label=T){
  par(mar=c(5, 4, 1, 1) + 0.1)
  if (missing(xrange)) { xrange <- sort(cc$data$dec[c(min(which(cc$data$density >= 1e-12)), max(which(cc$data$density >= 1e-12)))]) }
  if (missing(yrange)) { yrange <- c(0, diff(range(cc$data$density))*1.05) }
  plot.default(cc$data$dec, cc$data$density, type='l', main='', xlab='Declination', ylab='Density', lwd=2,
               col=col, xlim=xrange, ylim=yrange, xaxs='i', yaxs='i', axes=F); box()
  if (shading) {
    xp <- c(cc$data$dec, rev(cc$data$dec))
    yp <- c(cc$data$density, rep(0, length(cc$data$dec)))
    polygon(xp, yp, col=MESS::col.alpha(col,0.5), border=NA)
  }
  axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
  axis(1, at=0, labels = 0)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-45,360+45,0.1), lwd=0.2, labels=F) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0., labels=F) }
  axis(2)

  # objects
  if (!missing(obj)) {
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

  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
}



#' Plot significance test results
#'
#' This function creates a plot of the results of significance test
#' @param sig Object of \emph{skyscapeR.sigTest} format, created using \code{\link{sigTest}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param col (Optional) Colour to plot the curvigram in. Defaults to blue.
#' @param xrange (Optional) Array of two values restricting the horizontal range of the plot.
#' @param show.pval (Optional) Boolean to control whether the global p-value is written.
#' Default is True.
#' @param show.local (Optional) Boolean to control whether the local regions of significance
#' are plotted. Default is False.
#' @param legend (Optional) Show legend. Defaults to FALSE.
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' mag.az <- c(89.5, 105, 109.5)
#' data <- reduct.compass(loc, mag.az, "2016/04/02", alt=c(1,2,0))
#' data$Azimuth.Uncertainty <- 2  # adds the information on the preision of the azimuthal meaurement
#' sig <- sigTest(data)
#'
#' plot(sig, show.local=TRUE)
#' }
plot.skyscapeR.sigTest <- function(sig, obj, col='blue', xrange, show.pval=TRUE, show.local=FALSE, legend=FALSE){
  stop('Functionality remove ahead of a major update ')
}


#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param hor Object of \emph{skyscapeR.horizon} format.
#' @param show.az (Optional) Boolean that controls whether to display azimuth values or cardinal
#' directions. Defaults to FALSE.
#' @param xlim (Optional) Azimuth rage for plotting.
#' @param ylim (Optional) Altitude rage for plotting.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the paths of celestial objects.
#' @param refraction (Optional) Wheter to take refraction into account when displayingisplaying
#' the paths of celestial objects.
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
plot.skyscapeR.horizon <- function(hor, show.az=F, xlim, ylim, obj, refraction=F) {
  if (missing(xlim)) { xlim <- c(0,360) }
  if (missing(ylim)) { ylim <- c(floor(min(hor$data$alt, na.rm=T))-5,45) }

  par(mar=c(2,1,1,1))
  plot(-99999,-99999, xlab = "", ylab = "", yaxs='i', xaxs='i', axes=F, lwd=5, xlim=xlim, ylim=ylim)
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
        orb <- orbit(obj$decs[i], hor, res=0.5, refraction=refraction)
        lines(orb$az, orb$alt, col=obj$col[i], lty=obj$lty[i], lwd=obj$lwd[i])
      } else {
        orb1 <- orbit(obj$decs[3,i], hor, res=0.5)
        orb2 <- orbit(obj$decs[4,i], hor, res=0.5)
        lines(orb1$az, orb1$alt, col=obj$col[i], lty=obj$lty[i], lwd=obj$lwd[i])
        lines(orb2$az, orb2$alt, col=obj$col[i], lty=obj$lty[i], lwd=obj$lwd[i])
      }
    }
  }

  # Horizon line
  line(hor$data$az, hor$data$alt)
  x <- c(hor$data$az, rev(hor$data$az))
  y <- c(hor$data$alt, rep(-20,NROW(hor$data$az)))
  polygon(x , y, col=rgb(217/255,95/255,14/255,1), lwd=1.5)

  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
  box()
}



#' Plot stellar phase and seasonality
#'
#' This function creates a plot of stellar seasonality and phases/events.
#' @param starphase Object of \emph{skyscapeR.starphase} format.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @import MESS RColorBrewer
#' @seealso \code{\link{star.phases}}
#' @examples
#' # Plot the seasonality of Aldebaran for 3999 BCE:
#' \dontrun{
#' ss <- star.phases('Aldebaran',-4000, c(35,-8, 200))
#' plot(ss)
#' }
plot.skyscapeR.starphase = function(starphase, ...) {

  col <- RColorBrewer::brewer.pal(4,'Accent')
  seasons <- c("Rise and Set","Rise Only","Set Only","")
  events <- c()
  for (i in 1:NROW(starphase$phase)) {
    if (starphase$phase[i] == 'Curtailed Passage') {
      events <- c(events,'Acronycal\nRising','Heliacal\nSetting')
      seasons[4] <- 'Curtailed Passage'
    }
    if (starphase$phase[i] == 'Arising and Lying Hidden') {
      events <- c(events,'Acronycal\nSetting','Heliacal\nRising')
      seasons[4] <- 'Arising and Lying Hidden'
    }
  }

  par(mar=c(3,1,1,1), mgp=c(3,1.5,0))
  plot(-100,-100, xlim=c(1,365), ylim=c(0,1), main=paste0(starphase$star$name,' at ', BC.AD(starphase$year)), xlab="", ylab="", axes=FALSE, ...)
  axis(1, at=c(0,91,182,273,365), labels=c('December\nSolstice','March\nEquinox','June\nSolstice','September\nEquinox','December\nSolstice'), lwd.ticks=2)
  axis(1, at=seq(0,365,30.4), labels=NA)

  for (i in 1:4) {
    ind <- starphase$raw$seasons[[i]]
    ind.bb <- split(ind,cumsum(c(1,abs(diff(ind))>3)))

    for (j in 1:NROW(ind.bb)) {
      ind.i <- ind.bb[[j]]
      if (length(ind.i)>0) {
        x.poly <- c(ind.i[1]-.5,tail(ind.i,1)+.5,tail(ind.i,1)+.5,ind.i[1]-.5)
        y.poly <- c(0,0,1,1)
        polygon(x.poly, y.poly, col=col[i], border=NA)
        text(mean(ind.i),0.5,seasons[i], font=2, srt=90)
      }
    }

  }

  for (i in 1:(NROW(starphase$phase)*2)) {
    ind <- starphase$raw$events[[i]]
    ind.bb <- split(ind,cumsum(c(1,abs(diff(ind))>3)))

    for (j in 1:NROW(ind.bb)) {
      ind.i <- ind.bb[[j]]
      x.poly <- c(ind.i[1]-.5,tail(ind.i,1)+.5,tail(ind.i,1)+.5,ind.i[1]-.5)
      y.poly <- c(0,0,1,1)
      polygon(x.poly, y.poly, col=NA, border='black', lty=2, lwd=1.5)
      text(mean(ind.i),0.02,events[i], cex=0.8, font=1, srt=90 ,pos=4)
    }
  }
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'), side=3, at=par('usr')[2], cex=0.5, adj=1)
}


#' Prints significance test results
#'
#' This function prints the results of \code{\link{sigTest}}.
#' @param x Object of \emph{skyscapeR.sigTest} format.
#' @export
#' @seealso \code{\link{sigTest}}
print.skyscapeR.sigTest <- function (x) {
  cat("\n*** Results of Significance Test ***\n\n")
  cat(paste0('2-tailed test at ',x$metadata$conf*100,'% confidence, based on ', x$metadata$nsims, ' simulations.\n'))

  if (x$metadata$global.p.value == 0) {
    p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3))
  } else { p.value <- x$metadata$global.p.value }
  cat(paste0('global p-value: ', p.value, ' (',stars.pval(p.value),')\n'))
  cat('local p-values:\n')
  for (i in 1:NROW( x$metadata$local)) {
    if (x$metadata$local[i,]$p.value == 0) {
      p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3))
    } else { p.value <- x$metadata$local[i,]$p.value }
    cat(paste0('      ',x$metadata$local[i,]$type,'  dec range [',round(x$metadata$local[i,]$startDec,2), ', ',round(x$metadata$local[i,]$endDec,2),'] :: p-value: ', p.value, ' (',stars.pval(p.value),')\n'))
  }
}
