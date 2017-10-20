#' Polar plot of orientations (azimuths)
#'
#' This function creates a polar plot of azimuthal data. It is a wrapper for
#' \code{\link[plotrix]{polar.plot}}
#' @param az Array of azimuths or data frame with column named \emph{True.Azimuth}. Values
#' outside the [0, 360] range will be ignored.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the azimuths of celestial objects. Beware that this assumes a single
#' location (given by parameter loc) and a flat horizon of zero degrees.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. Only necessary for plotting potential
#' celestial targets.
#' @param obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the polar plot. Defaults to \emph{TRUE}.
#' @param ... Any other parameters to be passed unto \code{\link[plotrix]{polar.plot}}
#' @export
#' @import utils stats graphics
#' @seealso \code{\link[plotrix]{polar.plot}}, \code{\link{sky.objects}}
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' plotAz(az)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- sky.objects(c('sun','moon'), epoch=-2000, lty=c(2,3))
#' plotAz(az, tt, loc=c(35,-8))
plotAz = function(az, obj, loc, obj.label=T, ...) {
  if (class(az)=='data.frame') { az <- az$True.Azimuth }

  ind <- which(az > 360 | az < 0)
  if (length(ind) > 0) {
    az <- az[-ind]
    warning('Ignoring values outside of azimuth range [0, 360]')
  }

  n <- NROW(az)

  oldpar <- par()
  testlen <- c(0.01,rep(1,n))
  testpos <- c(0,az)
  label.pos <- c(0,45,90,135,180,225,270,315)
  labels <- c("N","","E","","S","","W","")
  par(mar=c(4, 4, 2, 2) + 0.1)
  plotrix::polar.plot(testlen, testpos, lwd=1.2, line.col='black', start=90, clockwise=T, labels= labels, label.pos = label.pos, show.grid.labels= F, grid.col='white', ...)
  plotrix::draw.circle(0,0,1, border='grey')
  mtext(paste0('skyscapeR ', packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # objects
  if (!missing(obj) & !missing(loc)) {

    for (i in 1:obj$n) {
      if (length(obj$epoch)==1) {
        rise <- c(); set <- c();
        orb <- orbit(obj$decs[i], loc, nutate_=F, refract_=F, aberration_=F)
        forb <- splinefun(orb$az, orb$alt)
        rise <- uniroot(forb, interval=c(0, 180))$root
        set <- uniroot(forb, interval=c(180, 360))$root
        tt <- c(rise, set)
        plotrix::polar.plot(c(0.01,rep(1,NROW(tt))), c(0,tt), lwd=c(0,rep(obj$lwd[i],2)), line.col=c(0,rep(obj$col[i],2)), lty=c(1,rep(obj$lty[i],2)), start=90, clockwise=T, add=T)
        if (obj.label) {
          plotrix::radial.plot.labels(c(0.01,1), c(0,rise), units='polar', start=pi*90/180, clockwise=T, labels=c('', colnames(obj$decs)[i]), cex=0.6, col=c('black',obj$col[i]), pos=4, offset=0.5)
          plotrix::radial.plot.labels(c(0.01,1), c(0,set), units='polar', start=pi*90/180, clockwise=T, labels=c('', colnames(obj$decs)[i]), cex=0.6, col=c('black',obj$col[i]), pos=2, offset=0.5)
        }
      } else {
        rise1 <- c(); rise2 <- c(); set1 <- c(); set2 <- c()
        orb1 <- orbit(obj$decs[3,i], loc, nutate_=F, refract_=F, aberration_=F)
        forb <- splinefun(orb1$az, orb1$alt)
        rise1 <- uniroot(forb, interval=c(0, 180))$root
        set1 <- uniroot(forb, interval=c(180, 360))$root

        orb2 <- orbit(obj$decs[4,i], loc, nutate_=F, refract_=F, aberration_=F)
        forb <- splinefun(orb2$az, orb2$alt)
        rise2 <- uniroot(forb, interval=c(0, 180))$root
        set2 <- uniroot(forb, interval=c(180, 360))$root

        rise <- seq(min(c(rise1,rise2)), max(c(rise1,rise2)), by=0.1)
        set <- seq(min(c(set1,set2)), max(c(set1,set2)), length.out=NROW(rise))
        tt <- rbind(rep(0,NROW(rise)), rise, set)
        dd <- cbind(rbind(rep(0.01,NROW(rise)),rep(0,NROW(rise)),rep(0,NROW(rise))), rbind(rep(0.01,NROW(rise)),rep(1,NROW(rise)),rep(1,NROW(rise))))

        plotrix::polar.plot(dd, polar.pos=tt, rp.type='p', poly.col=c(0,rep(MESS::col.alpha(obj$col[i],0.3),2)), line.col=c(0,rep(obj$col[i],2)), start=90, clockwise=T, add=T)
        if (obj.label) {
          plotrix::radial.plot.labels(c(0.01,1), c(0,mean(c(rise1,rise2))), units='polar', start=pi*90/180, clockwise=T, labels=c('', colnames(obj$decs)[i]), cex=0.6, col=c('black', obj$col[i]), pos=4, offset=0.5)
          plotrix::radial.plot.labels(c(0.01,1), c(0,mean(c(set1,set2))), units='polar', start=pi*90/180, clockwise=T, labels=c('', colnames(obj$decs)[i]), cex=0.6, col=c('black', obj$col[i]), pos=2, offset=0.5)
        }
      }
    }
    plotrix::polar.plot(testlen, testpos, lwd=1.2, line.col='black', start=90, clockwise=T, add = T)
  }
  options(warn=-2); par(oldpar); options(warn=0)
}

#' Plot a curvigram
#'
#' This function creates a plot of a curvigram.
#' @param curv Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the curvigram. Defaults to \emph{TRUE}.
#' @param signif (Optional) A \emph{skyscapeR.sig} object created with \code{\link{sigTest}}
#' for displaying confidence envelope around the chosen null hypothesis and overall p-value.
#' @param xlim Array of two values restricting the horizontal range of the plot.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @import utils stats graphics
#' @export
#' @seealso \code{\link{curvigram}}, \code{\link{sky.objects}}, \code{\link{sigTest}}
#' @examples
#' # Plot the curvigram of Recumbent Stone Circles:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, unc=2)
#' plotCurv(curv, xlim=c(-40,0))
#'
#' # Redo the plot to include lunar extreme declinations:
#' LEx <- sky.objects('moon', -2000, col='red', lty=2)
#' plotCurv(curv, objects=LEx, xlim=c(-40,0))
#'
#' # Add significance testing information:
#' \dontrun{
#' sig <- sigTest(curv, nh.Uniform(c(57,2)))
#' plotCurv(curv, objects=LEx, signif=sig, xlim=c(-40,0))
#' }
plotCurv = function(curv, obj, obj.label=T, signif, xlim=NULL, ...) {
  if (class(curv)!='skyscapeR.curv') { stop('No skyscapeR.curv object found.') }

  par(mar=c(4, 4, 2, 2) + 0.1)
  if (is.null(xlim)) { xlim <- c(min(curv$dec)-5, max(curv$dec)+5) }
  plot.default(-100,-100, xlab='Declination', ylab='Density', xlim=xlim, ylim=c(0,max(curv$density)), axes=F, ...)
  axis(1); axis(2)
  lines(curv$dec, curv$density, lwd=1.5, col='blue')
  box()
  mtext(paste0('skyscapeR ',packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # significance
  if (!missing(signif)) {
    x.polygon <- c(signif$null.hyp[1,], rev(signif$null.hyp[1,]))
    y.polygon <- c(signif$null.hyp[4,], rev(signif$null.hyp[2,]))
    polygon(x.polygon, y.polygon, border=NA, col=MESS::col.alpha('grey', 0.7))
    lines(signif$null.hyp[1,], signif$null.hyp[3,], col='grey3')
    abline(h=0)

    tail <- signif$type
    options(scipen=5)
    if (signif$p.value > 0) {
      text(xlim[2], max(curv$density), labels=bquote(paste('p'[.(tail)]*' = ', .(signif$p.value))), pos=2, cex=1.2)
    } else {
      text(xlim[2], max(curv$density), labels=bquote(paste('p'[.(tail)]*' < ', .(1/signif$nsims))), pos=2, cex=1.2)
    }
  }

  # objects
  if (!missing(obj)) {
    for (i in 1:obj$n) {
      if (length(obj$epoch)==1) {
        abline(v=obj$decs[i], col=obj$col[i], lwd=obj$lwd[i], lty=obj$lty[i])
        if (obj.label) { text(obj$decs[i], .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
      } else {
        xp <- c(obj$decs[3:4,i], rev(obj$decs[3:4,i]))
        yp <- c(-1,-1,2,2)
        polygon(xp, yp, border=obj$col[i], col=MESS::col.alpha(obj$col[i],.3))
        if (obj.label) { text(mean(obj$decs[3:4,i]), .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
      }
    }
  }
}



#' Plot a z-score transformed curvigram
#'
#' This function creates a plot of a z-score transformed curvigram, which is to say the
#' curvigram transformed into sigma units, based on a previously generated significance test.
#' @param signif A \emph{skyscapeR.sig} object created with \code{\link{sigTest}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the curvigram. Defaults to \emph{TRUE}.
#' @param xlim Array of two values restricting the horizontal range of the plot.
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, unc=2)
#' sig <- sigTest(curv, nh.Uniform(c(57,2)))
#'
#' plotZscore(sig)
#' }
plotZscore = function(signif, obj, obj.label=T, xlim=NULL) {
  if (class(signif)!='skyscapeR.sig') { stop('No skyscapeR.sig object found.') }

  # if (is.null(xlim)) { xlim <- c(min(signif$null.hyp.z[1,])-5, max(signif$null.hyp.z[1,])+5) }
  if (is.null(xlim)) { xlim <- signif$data.range }
  par(mar=c(4, 4, 2, 2) + 0.1)
  plot(-100,100, axes=F, xlim=xlim, ylim=c(-2, max(signif$maxima[2,])+1), xlab="Declination", ylab="")
  axis(1, at = seq(-90,90,10))
  mtext("Standard Deviations", side=2, line=2)
  axis(2, at = seq(-10,100,2), labels = seq(-10,100,2))
  lines(signif$null.hyp.z[1,], signif$null.hyp.z[3,], col='blue')
  x.polygon <- c(signif$null.hyp.z[1,], rev(signif$null.hyp.z[1,]))
  y.polygon <- c(signif$null.hyp.z[4,], rev(signif$null.hyp.z[2,]))
  polygon(x.polygon, y.polygon, border=NA, col=MESS::col.alpha('grey', 0.7))
  abline(h=0)
  box()
  mtext(paste0('skyscapeR ',packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)


  for (i in 1:NCOL(signif$maxima)) {
    points(signif$maxima[1,i], signif$maxima[2,i], pch=3)
    text(signif$maxima[1,i], signif$maxima[2,i], labels= substitute(paste(s, sigma), list(s = round(signif$maxima[2,i],2))), pos=4)
  }
  abline(h=0)

  # objects
  if (!missing(obj)) {
    for (i in 1:obj$n) {
      if (length(obj$epoch)==1) {
        abline(v=obj$decs[i], col=obj$col[i], lwd=obj$lwd[i], lty=obj$lty[i])
        if (obj.label) { text(obj$decs[i], .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
      } else {
        xp <- c(obj$decs[3:4,i], rev(obj$decs[3:4,i]))
        yp <- c(-1,-1,2,2)
        polygon(xp, yp, border=obj$col[i], col=MESS::col.alpha(obj$col[i],.3))
        if (obj.label) { text(mean(obj$decs[3:4,i]), .95*par('usr')[4], colnames(obj$decs)[i], col=obj$col[i], pos=4, offset=0.2, cex=0.7) }
      }
    }
  }

}



#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param hor Object of \emph{skyscapeR.horizon} format.
#' @param show.az Boolean that controls whether to display azimuth values on horizontal axis.
#'  Defaults to \emph{FALSE}.
#' @param max.alt Maximum altitude to display. Defaults to 45 degrees.
#' @param az0 Leftmost azimuth of plot. Defaults to 0 degrees, i.e. North at the left.
#' @param zoom Boolean that controls whether to provide a zoomed-in view of 100 degrees in
#' azimuth and 5 degrees of altitude above the horizon line. Defaults to \emph{FALSE}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the paths of celestial objects.
#' @param measure (Optional) A data.frame object with columns \emph{True.Azimuth} and
#' \emph{Altitude} such as the ones produced by \code{\link{reduct.compass}} or
#' \code{\link{reduct.theodolite}}. If no column named \emph{Altitude} is found then it will
#' plot all azimuths are zero degrees altitude.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @import utils stats graphics grDevices
#' @seealso \code{\link{download.HWT}}, \code{\link{sky.objects}}
#' @examples
#' # Plot a horizon retrieved from HeyWhatsThat:
#' hor <- download.HWT('HIFVTBGK')
#' plotHor(hor)
#'
#' # Add the paths of the solstices and equinoxes sun in the year 1999 BC:
#' tt <- sky.objects('sun', -2000, 'blue')
#' plotHor(hor, objects=tt)
plotHor <- function(hor, show.az=F, max.alt, az0 = 0, zoom=F, obj, measure, ...) {
  if (class(hor)!='skyscapeR.horizon') { stop('No skyscapeR.horizon object found.') }

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
  box()
  mtext(paste0('skyscapeR ',packageVersion('skyscapeR'),' Fabio Silva (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # objects
  if (!missing(obj)) {
    ind <- sort(obj$decs[1,], decreasing=T, index.return=T)$ix
    for (i in ind) {
      if (length(obj$epoch)==1) {
        orb <- orbit(obj$decs[i], hor, res=0.5)
        plotOrb(orb, obj$col[i])
      } else {
        orb1 <- orbit(obj$decs[3,i], hor, res=0.5)
        orb2 <- orbit(obj$decs[4,i], hor, res=0.5)
        plotOrb(orb1, obj$col[i])
        plotOrb(orb2, obj$col[i])
      }
    }
  }

  # Horizon line
  line(hor$az, hor$alt)
  x <- c(hor$az, rev(hor$az))
  y <- c(hor$alt, rep(-20,NROW(hor$az)))
  polygon(x ,y, col=rgb(217/255,95/255,14/255,1))

  # measurements
  if (!missing(measure)) {
    if (class(measure)=='data.frame') {
      if("True.Azimuth" %in% colnames(measure)) { az <- measure$True.Azimuth }
      if("Altitude" %in% colnames(measure)) { alt <- measure$Altitude } else { alt <- rep(0,NROW(az)) }
      points(az, alt)
    } else { message('Measurement data not in correct format. Check ?plotHor for more details.') }

  }
}


#' Plot visible path of celestial object on horizon
#'
#' This function adds the visible path of a celestial object to a horizon plot.
#' @param orbit Object of \emph{skyscapeR.orbit} format.
#' @param col String with colour to plot path in. Defaults to red.
#' @seealso \code{\link{plot.skyscapeR.horizon}}, \code{\link{orbit}}
#' @import utils stats graphics
#' @noRd
plotOrb<- function(orbit, col) {

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
#' ss <- star.phases('Aldebaran',-4000, c(35,-8))
#' plotPhases(ss)
#' }
plotPhases = function(starphase, ...) {
  if (class(starphase)!='skyscapeR.starphase') { stop('No skyscapeR.starphase object found.') }

  col <- RColorBrewer::brewer.pal(4,'Accent')
  seasons <- c("RS","R","S","")
  par(mar=c(2,1,1,1))

  plot(-100,-100, xlim=c(1,365), ylim=c(0,1), main=paste0(starphase$star$name,' @ ', starphase$year), xlab="", ylab="", axes=F, ...)
  axis(1, at=c(0,91,182,273,365), labels=c('dS','Eq','jS','Eq','dS'), lwd.ticks=2)
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
        text(mean(ind.i),0.5,seasons[i])
      }
    }

  }

  events <- c()
  for (i in 1:NROW(starphase$phase)) {
    if (starphase$phase[i] == 'Curtailed Passage') { events <- c(events,'AR','HS') }
    if (starphase$phase[i] == 'Arising and Lying Hidden') { events <- c(events,'AS','HR') }
  }

  for (i in 1:(NROW(starphase$phase)*2)) {
    ind <- starphase$raw$events[[i]]
    ind.bb <- split(ind,cumsum(c(1,abs(diff(ind))>3)))

    for (j in 1:NROW(ind.bb)) {
      ind.i <- ind.bb[[j]]
      x.poly <- c(ind.i[1]-.5,tail(ind.i,1)+.5,tail(ind.i,1)+.5,ind.i[1]-.5)
      y.poly <- c(0,0,1,1)
      polygon(x.poly, y.poly, col=MESS::col.alpha(col[4], alpha=.3), border=NA)
      text(mean(ind.i),0.5,events[i], cex=0.7, font=2)
    }
  }
}
