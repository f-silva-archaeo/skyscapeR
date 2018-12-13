## TODO
# (1) do other plotting functions (sigTest, starphase)
# (2) deal with epoch ranges on objects
# (3) standaerd inpout parameters like xrange rather than xlim
# (4) finish plotBars_ly and horizon
# (5) implement optional x-axis dimensions (declinaton, solar dates in gregorian or julian)

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




########################################
















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

# decs <- c(10, 12, -5, 4)
# plotBars(decs, unc=5)
# plotBars_ly(decs, unc=5)
# val <- decs
# unc <- 5

plotBars_ly <- function(val, unc, names, obj, obj.label=T, col='blue', shade=T, mark=T, sort=F, xlim, ylim) {
  if (min(val)< -90 | max(val)>90) { stop('It appears that val includes azimuth values This function can only be used for declination values.')}
  if (NROW(unc)==1) { unc <- rep(unc, NROW(val)) }
  if (sort) {
    ind <- sort(val, decreasing=T, index.return=T)$ix
    val <- val[ind]
    unc <- unc[ind]
    if (!missing(names)) { names <- names[ind] }
  }



  p <- plotly::plot_ly()

  for(i in 1:(length(val))){
    p <- plotly::add_trace(p,
                   x = c(val[i] - unc[i], val[i] + unc[i]),  # x0, x1
                   y = c(i, i),  # y0, y1
                   mode = "lines",
                   line = list(color = MESS::col.alpha('blue',.4), width = 40),
                   showlegend = F,
                   # hoverinfo = "text",

                   # Create custom hover text

                   # text = paste("Task: ", df$Task[i], "<br>",
                                # "Duration: ", df$Duration[i], "days<br>",
                                # "Resource: ", df$Resource[i]),

                   evaluate = T  # needed to avoid lazy loading
    )
  }

  p





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





#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param hor Object of \emph{skyscapeR.horizon} format.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the paths of celestial objects.
#' @param data (Optional) A data.frame object with columns \emph{True.Azimuth} and
#' @param show.unc (Optional) Boolean that controls whether to display uncertainty in altitude.
#'  Default is \emph{FALSE}.
#' @param show.ground (Optional) Boolean that controls whether to dispaly the ground. Default is \emph{TRUE}.
#' @param show.axes (Optional) Boolean that controls whether to display azimuth values on horizontal
#' axis and altitude on vertical axis. Default is \emph{FALSE}.
#' @param xrange Range of azimuth axis. Defaults to c(0, 360).
#' @param yrange Range of altitude axis. Defaults to c(-10, 45).
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
#' plot(hor, obj=tt)
plot_ly.skyscapeR.horizon <- function(hor, obj, data, show.unc=F, show.ground=T, show.axes=F, xrange, yrange) {
  if(missing(xrange)) { xrange <- c(0,360) }
  if(missing(yrange)) { yrange <- c(-10,45) }

  xx <- c(hor$data$az-360, hor$data$az, hor$data$az+360)
  yy <- rep(hor$data$alt, 3)

  p <- plotly::plot_ly(x = ~xx, y = ~yy, type = 'scatter', mode = 'lines',
                       line = list(color='black'),
                       hoverinfo = 'skip')

  ## TODO
  # fix colours
  # add orbits of objects (why not all being displayed????)
  # add epoch on orbits of
  # add things for input parameter show.axes

  ## objects
  if (!missing(obj)) {
    for (i in 1:obj$n) {
      assign(paste0('orb',i), orbit(obj$decs[i], hor, res=0.1))
      # ind <- which(orb$alt >= yrange[1]-10); orb$alt <- orb$alt[ind]; orb$az <- orb$az[ind]
      # p <- plotly::add_trace(p, x = ~get(paste0('orb',i))$az-360, y = ~get(paste0('orb',i))$alt,
      #                        line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
      #                        text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[1,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[1]),'\nAzimuth: ', round(get(paste0('orb',i))$az,1),'º\nAltitude: ', round(get(paste0('orb',i))$alt,2)-360,'º'),
      #                        hoverinfo = 'text',
      #                        showlegend = F)
      p <- plotly::add_lines(p, x = ~get(paste0('orb',i))$az, y = ~get(paste0('orb',i))$alt,
                             line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
                             text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[1,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[1]),'\nAzimuth: ', round(get(paste0('orb',i))$az,1),'º\nAltitude: ', round(get(paste0('orb',i))$alt,2),'º'),
                             hoverinfo = 'text',
                             name = colnames(obj$decs)[i] )
      # p <- plotly::add_trace(p, x = ~get(paste0('orb',i))$az+360, y = ~get(paste0('orb',i))$alt,
      #                        line = list( color = obj$col[i], width = obj$lwd[i], dash = lty2dash(obj$lty[i]) ),
      #                        text = paste0(colnames(obj$decs)[i],'\nDeclination: ',round(obj$decs[1,i],1),'º\nEpoch: ', yr2epoch(obj$epoch[1]),'\nAzimuth: ', round(get(paste0('orb',i))$az,1),'º\nAltitude: ', round(get(paste0('orb',i))$alt,2)+360,'º'),
      #                        hoverinfo = 'text',
      #                        showlegend = F)
    }
  }

  ## ground
  if (show.ground) {
    p <- plotly::add_trace(p, x = ~xx, y = ~yy, showlegend = F, hoverinfo = 'none')
    p <- plotly::add_trace(p, x = ~xx, y = ~rep(-20,length(yy)), fill = 'tonexty', fillcolor = 'rgba(209,156,31,1)', showlegend = F, hoverinfo = 'none')
  }


  ## altitude uncertainty
  if (show.unc) {
    p <- plotly::add_ribbons(p, x = ~xx,
                             ymin = ~rep(hor$data$alt - hor$data$alt.unc, 3),
                             ymax = ~rep(hor$data$alt + hor$data$alt.unc, 3),
                             line = list(color = MESS::col.alpha('grey80',0.5)),
                             fillcolor = MESS::col.alpha('grey80',0.5),
                             name = "Horizon Altitude Uncertainty",
                             hoverinfo = 'none')
  }

  ## Horizon line
  p <- plotly::add_trace(p, x = ~xx, y = ~yy, type = 'scatter', mode = 'lines',
                       line = list(color='black'),
                       text = paste0('Horizon\nAzimuth: ',round(xx,1),'º\nAltitude: ', round(yy,2),'º'),
                       hoverinfo = 'text',
                       name = 'Horizon')

  p <- plotly::layout(p,
                      xaxis = list(title = 'Azimuth',
                                   range = xrange,
                                   tickvals = seq(-360,360+360,45),
                                   ticktext = c(rep(c('N','NE','E','SE','S','SW','W','NW'),3),'N'),
                                   showgrid = F,
                                   zeroline = F),
                      yaxis = list(title = 'Altitude',
                                   range = yrange,
                                   showgrid = F,
                                   zeroline = F),
                      annotations = cp,
                      showlegend = F)
  p
}
