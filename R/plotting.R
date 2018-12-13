# copyright for plots
cp <- list(
  x = -0.05,
  y = 1.05,
  text = paste0('skyscapeR ', packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'),
  font = list ( size = 7 ),
  xref = "paper",
  xanchor = "left",
  yref = "paper",
  yanchor = "top",
  showarrow = FALSE
)


#' Polar plot of orientations (azimuths)
#'
#' This function creates a polar plot of azimuthal data using \emph{plotly}.
#' @param az Array of azimuths or data frame with column named \emph{True.Azimuth}. Values
#' outside the [0, 360] range will be ignored.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the azimuths of celestial objects. Beware that this assumes a single
#' location (given by parameter loc) and a flat horizon of zero degrees.
#' @param loc (Optional) This can be either a vector with the latitude and longitude of the
#' location, or a \emph{skyscapeR.horizon} object. Only necessary for plotting potential
#' celestial targets.
#' @param col (Optional) Single colour or colour pallete to use for plotting measurements.
#' @param lwd (Optional) Line width to plot measurements. Defaults to 1.
#' @param lty (Optional) Line type to plot measurements. Defaults to 1.
#' @param len (Optional) Length of line to plot measurements.
#' @param legend (Optional) Show legend. Defaults to FALSE.
#' @export
#' @import utils stats plotly
#' @examples
#' # Plot some azimuth data:
#' az <- c(120, 100, 93, 97, 88, 115, 112, 67)
#' plotAz(az)
#'
#' # To visualise this data against the common solar and lunar targets:
#' tt <- sky.objects(c('sun','moon'), epoch=-2000, lty=1, lwd=1.5, col=c('black','red'))
#' plotAz(az, obj=tt, loc=c(35,-8), legend=T)
#'
#' # To visualise the same data with shorter lines
#' plotAz(az, obj=tt, loc=c(35,-8), len=0.2)
#'
#' # To display only celestial objects
#' plotAz(az=NULL, obj=tt, loc=c(35,-8), legend=T)
plotAz = function(az, obj, loc, col='blue', lwd=1, lty=1, len=1, legend=F) {
  if (class(az)=='data.frame') {
    az <- az$True.Azimuth
    names <- az$Name
  } else {
    names <- as.character(1:NROW(az))
  }
  if (length(col)==1) { col <- rep(col, length(az)) }
  len <- 1-len

  # Setup
  p <- plotly::plot_ly(
    type = 'scatterpolar',
    mode = 'lines'
  )

  if (!missing(obj)) {
    # Add celestial objects
    for (i in 1:obj$n) {
      rise <- c(); set <- c();
      orb <- orbit(obj$decs[i], loc, refraction=F)
      forb <- splinefun(orb$az, orb$alt)
      rise <- uniroot(forb, interval=c(0, 180))$root
      set <- uniroot(forb, interval=c(180, 360))$root
      tt <- round(c(rise, set),1)

      sl <- T
      for (j in 1:length(tt)) {
        p <- plotly::add_trace(p,
                               r = seq(len,1,by=.2),
                               theta = rep(tt[j],length(seq(len,1,by=.2))),
                               line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
                               text = paste0(colnames(obj$decs)[i],'\nAzimuth: ',round(tt[j],2),'º\nEpoch: ', yr2epoch(obj$epoch)),
                               hoverinfo = 'text',
                               name = colnames(obj$decs)[i],
                               showlegend = sl
        )
        sl <- F
      }
    }
  }

  # Add measurements
  sl <- T
  for (i in 1:length(az)) {
    p <- plotly::add_trace(p,
                           r = seq(len,1,by=.2),
                           theta = rep(az[i],length(seq(len,1,by=.2))),
                           line = list( color = col[i], width = lwd, dash = lty2dash(lty) ),
                           text = paste0('Site: ',names[i],'\nAzimuth: ',round(az[i],2),'º'),
                           hoverinfo = 'text',
                           name='Data',
                           showlegend = sl
    )
    sl <- F
  }

  # Layout
  p <- plotly::layout(p,
                      polar = list(
                        radialaxis = list( visible = F ),
                        angularaxis = list(
                          direction = 'clockwise',
                          tickmode = 'array',
                          tickvals = seq(0,359,45),
                          ticktext = c('N','NE','E','SE','S','SW','W','NW'),
                          tickfont = list( size = 14 )
                        )
                      ),
                      font = list(
                        family = 'Arial',
                        size = 11,
                        color = '#000'
                      ),
                      annotations = cp,
                      showlegend = legend
  )


  p
}


#' Bar plot of orientations (declinations)
#'
#' This function creates a bar plot of orientation data.
#' @param val Array of declinations.
#' @param unc Single value or array of declination uncertainty
#' @param names (Optional) Array of names of measurements in \emph{val}
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declinations of celestial objects.
#' @param obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the polar plot. Defaults to \emph{TRUE}.
#' @param col (Optional) Colour to plot measurments in. Defaults to \emph{blue}.
#' @param shade (Optional) Boolean to control whether to shade a polygon of measurements. Defaults
#' to \emph{TRUE}
#' @param mark (Optional) Boolean to control whether to mark the declination value. Defaults to
#' \emph{TRUE}
#' @param sort (Optional) Boolean to control whether to sort the measurements by their declination
#' value. Defaults to \emph{FALSE}
#' @param xlim (Optional) Array of limits for x-axis.
#' @param ylim (Optional) Array of limits for y-axis.
#' @export
#' @import utils stats graphics
#' @seealso \code{\link{sky.objects}}
#' @examples
#' # Plot some declination data:
#' decs <- c(10, 12, -5, 4)
#' plotBars(decs, unc=5)
#'
#' # To visualize this data against the common solar and lunar targets:
#' tt <- sky.objects(c('sun','moon'), epoch=-2000, lty=c(2,3))
#' plotBars(decs, unc=5, obj=tt)
plotBars <- function(val, unc, names, obj, obj.label=T, col='blue', shade=T, mark=T, sort=F, xlim, ylim) {
  if (min(val)< -90 | max(val)>90) { stop('It appears that val includes azimuth values This function can only be used for declination values.')}
  if (NROW(unc)==1) { unc <- rep(unc, NROW(val)) }
  if (sort) {
    ind <- sort(val, decreasing=T, index.return=T)$ix
    val <- val[ind]
    unc <- unc[ind]
    if (!missing(names)) { names <- names[ind] }
  }

  if (!missing(names)) { par(mar=c(4, 9, 2, 2) + 0.1) } else { par(mar=c(4, 2, 2, 2) + 0.1) }
  if (missing(xlim)) { xlim <- c(min(val-unc)-5, max(val+unc)+5) }
  if (missing(ylim)) { ylim <- c(0.5, NROW(val)+0.5) }
  plot.default(-100,-100, xlab='Declination (º)', ylab='', xlim=xlim, ylim=ylim, axes=F, yaxs='i')
  axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
  axis(1, at=0, labels = 0)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 2) { axis(1, at=seq(-90,360,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 1) { axis(1, at=seq(-90,360,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,360,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,360,10), lwd=0.5, labels=F) }
  if (!missing(names)) { axis(2, at=1:NROW(names), labels=names, las=2, cex=0.7) } else { axis(2, at=1:NROW(val), las=2) }

  box()
  mtext(paste0('skyscapeR v',packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  if (shade) { border <- NA } else { border <- col; col <- NA  }

  for (i in 1:NROW(val)) {
    xp <- c(val[i]-unc[i], val[i]+unc[i], val[i]+unc[i], val[i]-unc[i])
    yp <- c(i-0.4,i-0.4,i+0.4,i+0.4)
    polygon(xp, yp, col=MESS::col.alpha(col,0.5), border=border)
    if (mark) { lines(c(val[i],val[i]), c(i-0.4,i+0.4), col=col, lwd=2) }
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



#' Plot a histogram
#'
#' This function creates a plot of a histogram
#' @param hh Object of \emph{skyscapeR.hist} format, created using \code{\link{histogram}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param col (Optional) Colour to plot the histogram in. Defaults to blue.
#' @param xrange Array of two values restricting the horizontal range of the plot.
#' @param legend (Optional) Show legend. Defaults to FALSE.
#' @import utils stats graphics
#' @export
#' @seealso \code{\link{histogram}}, \code{\link{curvigram}}, \code{\link{sky.objects}}
#' @examples
#' # Plot the histogram of Recumbent Stone Circles:
#' data(RugglesRSC)
#' hist <- histogram(RugglesRSC$Dec, unc=2)
#' plot(hist, xrange=c(-40,0))
#'
#' # Redo the plot to include lunar extreme declinations:
#' LEx <- sky.objects('moon', -2000, col='red', lty=2)
#' plot(hist, obj=LEx, xrange=c(-40,0))
plot.skyscapeR.hist <- function(hh, obj, col='blue', xrange, legend=F){

  xx <- c(); yy <- c()
  for (i in 1:NROW(hh$data$density)) {
    xx <- c(xx, hh$data$dec[i], hh$data$dec[i+1])
    yy <- c(yy, rep(hh$data$density[i],2))
  }

  p <- plotly::plot_ly(type = 'scatter', mode = 'lines')
  p <- plotly::add_trace(p, x = ~xx, y = ~yy, fill = 'tozeroy',
                         hoverinfo='skip',
                         name='Data',
                         line = list( color = MESS::col.alpha(col, 0.7)),
                         fillcolor = MESS::col.alpha(col, 0.5))
  p <- plotly::add_markers(p, x = ~(hh$data$dec+diff(hh$data$dec)[1]/2)[1:NROW(hh$data$density)], y = ~hh$data$density,
                           text = paste0('Dec Range: [', hh$data$dec[1:NROW(hh$data$density)], 'º;', hh$data$dec[2:NROW(hh$data$dec)],'º]\nDensity: ',hh$data$density),
                           hoverinfo = 'text',
                           marker = list(size=2, color=MESS::col.alpha(col, 0.7)),
                           showlegend=F)


  if (missing(xrange)) { xrange <- range(hh$data$dec) }

  # objects
  if (!missing(obj)) {
    for (i in 1:obj$n) {
      if (sum(obj$decs[,i] >= xrange[1] & obj$decs[,i] <= xrange[2])) {
        p <- plotly::add_trace(p, x = rep(obj$decs[1,i],NROW(seq(0,max(hh$data$density), length.out = 20))), y = seq(0,max(hh$data$density), length.out = 20),
                               line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
                               text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[1,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[1])),
                               hoverinfo = 'text',
                               name = colnames(obj$decs)[i] )
        if (length(obj$epoch)==2) {
          p <- plotly::add_trace(p, x = rep(obj$decs[2,i],NROW(seq(0,max(hh$data$density), length.out = 20))), y = seq(0,max(hh$data$density), length.out = 20),
                                 line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
                                 text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[2,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[2])),
                                 hoverinfo = 'text',
                                 showlegend = F,
                                 fill = 'tonexty')
        }
      }
    }
  }

  p <- plotly::layout(p,
                      xaxis = list(title = 'Declination (º)', range = xrange),
                      yaxis = list(title = 'Density'),
                      annotations = cp,
                      showlegend = legend)
  p
}





#' Plot a curvigram
#'
#' This function creates a plot of a curvigram
#' @param cc Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param col (Optional) Colour to plot the curvigram in. Defaults to blue..
#' @param xlim Array of two values restricting the horizontal range of the plot.
#' @param legend (Optional) Show legend. Defaults to FALSE.
#' @import utils stats graphics
#' @export
#' @seealso \code{\link{curvigram}}, \code{\link{histogram}}, \code{\link{sky.objects}}
#' @examples
#' # Plot the curvigram of Recumbent Stone Circles:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, unc=2)
#' plot(curv)
#'
#' # Redo the plot to include lunar extreme declinations:
#' LEx <- sky.objects('moon', -2000, col='red', lty=2)
#' plot(curv, obj=LEx, legend=T)
plot.skyscapeR.curv <- function(cc, obj, col='blue', xlim, legend=F){

  xx <- cc$data$dec
  yy <- cc$data$density

  p <- plotly::plot_ly(type = 'scatter', mode = 'lines')
  p <- plotly::add_trace(p, x = ~xx, y = ~yy, fill = 'tozeroy',
                         line = list( color = MESS::col.alpha(col, 0.7)),
                         fillcolor = MESS::col.alpha(col, 0.5),
                         text = paste0('Dec: ', round(cc$data$dec,1), 'º\nDensity: ',round(cc$data$density,3)),
                         hoverinfo = 'text',
                         name = 'Data')

  if (missing(xlim)) { xlim <- range(cc$data$dec) }

  # objects
  if (!missing(obj)) {
    for (i in 1:obj$n) {
      if (sum(obj$decs[,i] >= xlim[1] & obj$decs[,i] <= xlim[2])) {
        p <- plotly::add_trace(p, x = rep(obj$decs[1,i],NROW(seq(0,max(cc$data$density), length.out = 20))), y = seq(0,max(cc$data$density), length.out = 20),
                               line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
                               text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[1,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[1])),
                               hoverinfo = 'text',
                               name = colnames(obj$decs)[i] )
        if (length(obj$epoch)==2) {
          p <- plotly::add_trace(p, x = rep(obj$decs[2,i],NROW(seq(0,max(cc$data$density), length.out = 20))), y = seq(0,max(cc$data$density), length.out = 20),
                                 line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
                                 text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[2,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[2])),
                                 hoverinfo = 'text',
                                 showlegend = F,
                                 fill = 'tonexty')
        }
      }
    }
  }

  p <- plotly::layout(p,
                      xaxis = list(title = 'Declination (º)', range = xlim),
                      yaxis = list(title = 'Density'),
                      annotations = cp,
                      showlegend = legend)
  p
}



#' Plot significance test results
#'
#' This function creates a plot of the results of significance test
#' @param sig Object of \emph{skyscapeR.sigTest} format, created using \code{\link{sigTest}}.
#' @param xlim (Optional) Array of two values restricting the horizontal range of the plot.
#' @param show.pval (Optional) Boolean to control whether the global p-value is written.
#' Default is True.
#' @param show.local (Optional) Boolean to control whether the local regions of significance
#' are plotted. Default is False.
#' @import utils stats graphics
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, sd=2)
#' sig <- sigTest(curv, null.hyp=nh.Uniform(c(57,2)))
#'
#' plot(sig, show.local=T)
#' }
plot.skyscapeR.sigTest <- function(sig, xlim, show.pval=T, show.local=F) {
  # empirical
  emp <- sig$data$empirical
  if (missing(xlim)) { xlim <- sort(emp$data$dec[c(min(which(emp$data$density >= 1e-12)), max(which(emp$data$density >= 1e-12)))]) }
  if (show.local) { ylim <- c(-max(emp$data$density)*.05,max(emp$data$density)) } else { ylim <- c(0,max(emp$data$density))}

  if (sig$metadata$type=='hist') {
    xx <- c(); yy <- c()
    for (i in 1:NROW(emp$data$density)) {
      xx <- c(xx, emp$data$dec[i], emp$data$dec[i+1])
      yy <- c(yy, rep(emp$data$density[i],2))
    }
  } else {
    xx <- emp$data$dec
    yy <- emp$data$density
  }
  plot.default(xx, yy, type='l', main='', xlab='Declination', ylab='Density', lwd=2, col='blue', xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', axes=F); box()
  xp <- c(xx, rev(xx))
  yp <- c(yy, rep(0, length(xx)))
  polygon(xp, yp, col=MESS::col.alpha('blue',0.5), border=NA)
  axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
  axis(1, at=0, labels = 0)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 2) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 1) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0.5, labels=F) }
  axis(2)
  mtext(paste0('skyscapeR v',packageVersion('skyscapeR'),' (', substr(packageDescription('skyscapeR')$Date,1,4),')'),3, adj=0, cex=0.5)

  # sigtest
  if (sig$metadata$type=='hist') {
    yy1 <- c(); yy2 <- c(); yy3 <- c()
    for (i in 1:NROW(sig$CE.mean)) {
      yy1 <- c(yy1, rep(sig$data$CE.mean[i],2))
      yy2 <- c(yy2, rep(sig$data$CE.upper[i],2))
      yy3 <- c(yy3, rep(sig$data$CE.lower[i],2))
    }
  } else {
    yy1 <- sig$data$CE.mean
    yy2 <- sig$data$CE.upper
    yy3 <- sig$data$CE.lower
  }
  lines(xx, yy1, col='grey')
  xp <- c(xx, rev(xx))
  yp <- c(yy2, rev(yy3))
  polygon(xp, yp, col=MESS::col.alpha('grey', 0.5), border=NA)

  # global p-value
  if (show.pval) {
    if (sig$metadata$global.p.value == 0 ) {
      pval <- paste0("global p-value < ", round(1/(sig$metadata$nsims+1),4))
    } else {
      pval <- paste0('global p-value = ', sig$metadata$global.p.value)
    }
    text(par('usr')[2], abs(diff(par('usr')[3:4]))*.90, pos=2, pval, cex=0.7)
  }

  # regions of significance
  if (show.local) {
    abline(0,0, lwd=1)
    aux <- as.matrix(sig$metadata$local[,1:3])
    for (i in 1:NROW(aux)) {
      if (sig$metadata$local[i,4] == '+') { col <- 'darkgreen' } else { col <- 'red' }
      xp <- c(aux[i,1], aux[i,1], aux[i,2], aux[i,2])
      yp <- c(ylim[1], 0, 0, ylim[1])
      polygon(xp, yp, col=MESS::col.alpha(col, .4), border=NA)
      # arrows(aux[i,1], ylim[1]/2, aux[i,2], ylim[1]/2, col='red', pch=3, angle=90, code=3, lwd=0.7, length=0.05)
      text(mean(aux[i,1:2]), ylim[1]/2, labels=stars.pval(aux[i,3]), col='black', cex=0.9)
    }
  }
}


#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param hor Object of \emph{skyscapeR.horizon} format.
#' @param show.unc (Optional) Boolean that controls whether to display uncertainty in altitude.
#'  Default is \emph{FALSE}.
#' @param show.axes (Optional) Boolean that controls whether to display azimuth values on horizontal
#' axis and altitude on vertical axis. Default is \emph{FALSE}.
#' @param max.alt Maximum altitude to display. Defaults to 45 degrees.
#' @param az0 Leftmost azimuth of plot. Defaults to 0 degrees, i.e. North at the left.
#' @param zoom Boolean that controls whether to provide a zoomed-in view of 100 degrees in
#' azimuth and 5 degrees of altitude above the horizon line. Default is \emph{FALSE}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the paths of celestial objects.
#' @param measure (Optional) A data.frame object with columns \emph{True.Azimuth} and
#' \emph{Altitude} such as the ones produced by \code{\link{reduct.compass}} or
#' \code{\link{reduct.theodolite}}. If no column named \emph{Altitude} is found then it will
#' plot all azimuths are zero degrees altitude.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @import utils stats graphics grDevices
#' @seealso \code{\link{downloadHWT}}, \code{\link{sky.objects}}
#' @examples
#' # Plot a horizon retrieved from HeyWhatsThat:
#' hor <- downloadHWT('HIFVTBGK')
#' plot(hor)
#'
#' # Add the paths of the solstices and equinoxes sun in the year 1999 BC:
#' tt <- sky.objects('sun', -2000, 'blue')
#' plot(hor, objects=tt)
plot.skyscapeR.horizon <- function(hor, show.unc=F, show.axes=F, max.alt, az0 = 0, zoom=F, obj, measure, ...) {

  # rejiggle so plot starts at given start point
  if (az0 < -360) { az0 <- az0 + 360 }
  ind <- which(hor$data$az < az0)
  if (NROW(ind)>0) {
    hor$data$az <- c(hor$data$az, hor$data$az[ind] + 360)
    hor$data$alt <- c(hor$data$alt, hor$data$alt[ind])
  }

  yl <- floor(min(hor$data$alt, na.rm=T))-5
  if (missing(max.alt)) { yl[2] <- 45 } else { yl[2] <- max.alt }

  if (zoom) {
    xl <- c(az0, az0+100)
    yl <- c(-5,5)
    show.axes <- T
  } else {
    xl <- c(az0, az0 + 360)
  }

  if (show.unc) { show.axes <- T }

  if (show.axes) {
    xx <- "azimuth"
    ll <- seq(-360,720,by=10)
    at <- ll
    yy <- 'altitude'
    par(mar=c(4,4,1,1))
  } else {
    xx <- ""
    ll <- c("N","E","S","W","N","E","S","W","N","E","S","W","N")
    at <- seq(-360,720,by=90)
    yy <- ''
    par(mar=c(2,1,1,1))
  }


  plot(-99999,-99999, xlab = xx, ylab = yy, yaxs='i', xaxs='i', ylim = yl, xlim = xl, axes=F, lwd=5, ...)
  if (show.axes == T) {
    axis(1, at = at, labels = ll)
    axis(2)
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
        plot(orb, obj$col[i])
      } else {
        orb1 <- orbit(obj$decs[1,i], hor, res=0.5)
        orb2 <- orbit(obj$decs[2,i], hor, res=0.5)
        plot(orb1, obj$col[i])
        plot(orb2, obj$col[i])
      }
    }
  }

  # Horizon line
  lines(hor$data$az, hor$data$alt)

  # Altitude uncertainty
  if (show.unc) {
    x <- c(hor$data$az, rev(hor$data$az))
    y <- c(hor$data$alt-hor$data$alt.unc, rev(hor$data$alt+hor$data$alt.unc))
    polygon(x, y, col='grey', border='grey')
    lines(hor$data$az, hor$data$alt, col='black', lwd=1.2)
  } else {
    x <- c(hor$data$az, rev(hor$data$az))
    y <- c(hor$data$alt, rep(-20,NROW(hor$data$az)))
    polygon(x ,y, col=rgb(217/255,95/255,14/255,1))
  }

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
plot.skyscapeR.orbit<- function(orbit, col) {

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
#' plot(ss)
#' }
plot.skyscapeR.starphase = function(starphase, ...) {

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
