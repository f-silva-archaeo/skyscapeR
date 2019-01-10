#' Create a simplistic sketchview of the sky at a given moment in time
#'
#' @param time Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. Defaults to GMT.
#' @param cal (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. Defaults to Gregorian.
#' @param xrange Range of azimuths to display, preferably no larger than 120 degrees.
#' @param yrange Range of altitudes to display, preferably no larger than 120 degrees.
#' @param sun (Optional) Boolean on whether the sun should be displayed. Defaults to \emph{TRUE}.
#' @param moon (Optional) Boolean on whether the moon should be displayed. Defaults to \emph{TRUE}.
#' @param planets (Optional) Boolean on whether the visible planets should be displayed. Defaults to \emph{FALSE}.
#' @param exagerate (Optional) Boolean on whether the size of the sun, moon and planets should be exagerated,
#' which can be usedful when attempting wider viewing angles. Defaults to \emph{TRUE}.
#' @param nstars (Optional) Number of stars to display. Defaults to \emph{1000}.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param atm (Optional) Atmospheric pressure (in mbar). Default is 1013.25 mbar.
#' @param temp (Optional) Atmospheric temperature (in Celsius). Default is 15 degrees.
#' @export
#' @examples
#' sky.sketch(time='2019-01-10 16:51', loc=c(35,-8))
sky.sketch <- function(time, timezone='GMT', cal='G', xrange=c(30,150), yrange=c(-45,45), sun=T, moon=T, planets=F, exagerate=T, nstars=1000, loc, atm=1013.25, temp=15) {

  if (abs(diff(xrange)) > 120) { message('xrange is larger than 120 degrees. This will cause a distorted view of the sky.')}
  if (abs(diff(yrange)) > 120) { message('yrange is larger than 120 degrees. This will cause a distorted view of the sky.')}

  ## prep stars info
  fpath <- system.file("ephemeris", "sefstars.txt", package="swephR")
  cnames <- c('traditional name','nomenclature name','equinox','RA hr','RA min', 'RA sec', 'Dec deg', 'Dec min', 'Dec sec', 'pm RA', 'pm Dec', 'rad vel', 'ann plx', 'mag V', 'DM zone', 'DM number')
  sefstars <- read.csv(fpath, as.is=T, header=F, comment.char='#', col.names=cnames)
  sefstars[,1] <- c()
  sefstars <- unique(sefstars)
  ind <- sort(sefstars$mag.V, index.return=T)$ix
  sefstars <- sefstars[ind,]

  ## proc
  if (class(time)=='character') { jd <- time2jd(time, timezone, cal) } else { jd <- time }

  plot(-999,-999, xlim=xrange, ylim=yrange, xlab='', ylab='')

  if (nstars > 0) {
    for (i in 1:nstars) {
      star <- star(paste0(',',sefstars$nomenclature.name[i]), year=as.numeric(substr(time,1,4)))
      mag <- star$app.mag
      star <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(star$RA, star$Dec))$xaz[c(1,3)]
      star[1] <- star[1]+180; if (star[1]>360) { star[1] <- star[1]-360 }
      points(star[1], star[2], pch=19, cex=mag2size(mag), col='black')
      ## TODO colour?
    }
  }

  if (sun) {
    sun <- body.position('sun', jd)
    sun <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(sun$RA, sun$Dec))$xaz[c(1,3)]
    sun[1] <- sun[1]+180; if (sun[1]>360) { sun[1] <- sun[1]-360 }
    plotrix::draw.circle(sun[1], sun[2], radius=0.5+3*exagerate, border='black', col='yellow', lwd=.5)
  }

  if (moon) {
    moon <- body.position('moon', jd)
    moon <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(moon$RA, moon$Dec))$xaz[c(1,3)]
    moon[1] <- moon[1]+180; if (moon[1]>360) { moon[1] <- moon[1]-360 }
    moonphase <- moonphase(jd)
    plotrix::draw.circle(moon[1], moon[2], radius=0.5+3*exagerate, border='black', col='white', lwd=.5)
    ## TODO phase
  }

  if (planets) {
    planets <- c('mercury','venus','mars', 'jupiter', 'saturn')
    for (i in 1:5) {
      planet <- body.position(planets[i], jd)
      planet <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(planet$RA, planet$Dec))$xaz[c(1,3)]
      planet[1] <- planet[1]+180; if (planet[1]>360) { planet[1] <- planet[1]-360 }
      points(planet[1], planet[2], pch=19, cex=.5, col='darkred')
    }

  }

  polygon(c(0,360,360,0),c(0,0,-90,-90), border=NA, col=MESS::col.alpha('grey',0.7))
  abline(h=0, lwd=2)
}
