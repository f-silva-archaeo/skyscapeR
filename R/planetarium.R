#' Create a simplistic sketch of the sky at a given moment in time
#'
#' @param time Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. #' If not given the value set
#' by \code{\link{skyscapeR.vars}} will be used instead.
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
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
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temprature for erfraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @export
#' @examples
#' sky.sketch(time='2019/01/10 18:51', loc=c(35,-8,100))
sky.sketch <- function(time, timezone, calendar, xrange=c(30,150), yrange=c(-45,45), sun=T, moon=T, planets=F, exagerate=T, nstars=1000, loc, atm, temp) {

  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (abs(diff(xrange)) > 120) { message('xrange is larger than 120 degrees. This will cause a distorted view of the sky.')}
  if (abs(diff(yrange)) > 120) { message('yrange is larger than 120 degrees. This will cause a distorted view of the sky.')}

  ## prep stars info
  fpath <- system.file("ephemeris", "sefstars.txt", package="swephR")
  cnames <- c('traditional name','nomenclature name','equinox','RA hr','RA min', 'RA sec', 'Dec deg', 'Dec min', 'Dec sec', 'pm RA', 'pm Dec', 'rad vel', 'ann plx', 'mag V', 'DM zone', 'DM number')
  sefstars <- read.csv(fpath, as.is=T, header=F, comment.char='#', col.names=cnames)
  sefstars <- sefstars[-which(sefstars$mag.V==0),]
  ind <- sort(sefstars$mag.V, index.return=T)$ix
  sefstars <- sefstars[ind,]
  df <- sefstars[,-1]
  ind <- which(duplicated(df))
  sefstars <- sefstars[-ind,]; rm(df)

  ## proc
  if (class(time)=='character') { jd <- time2jd(time, timezone, calendar) } else { jd <- time }

  plot(-999,-999, xlim=xrange, ylim=yrange, xlab='', ylab='')

  if (nstars > 0) {
    for (i in 1:nstars) {
      star <- star(paste0(',',sefstars$nomenclature.name[i]), year=as.numeric(substr(time,1,4)))
      mag <- star$app.mag
      star <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(star$coord$RA, star$coord$Dec))$xaz[c(1,3)]
      star[1] <- star[1]+180; if (star[1]>360) { star[1] <- star[1]-360 }
      points(star[1], star[2], pch=19, cex=mag2size(mag), col='black')
      ## TODO add colour
    }
  }

  if (sun) {
    sun <- as.numeric(body.position('sun', jd, loc = loc, refraction = T, atm = atm, temp = temp, verbose = F)$horizontal)
    plotrix::draw.circle(sun[1], sun[2], radius=0.5+1*exagerate, border='black', col='yellow', lwd=.5)
  }

  if (moon) {
    moon <- as.numeric(body.position('moon', jd, loc = loc, refraction = T, atm = atm, temp = temp, verbose = F)$horizontal)
    moonphase <- moonphase(jd)
    plotrix::draw.circle(moon[1], moon[2], radius=0.5+1*exagerate, border='black', col='white', lwd=.5)
    ## TODO add lunar phase
  }

  if (planets) {
    planets <- c('mercury','venus','mars', 'jupiter', 'saturn')
    for (i in 1:5) {
      planet <- body.position(planets[i], jd, loc = loc, refraction = T, atm = atm, temp = temp, verbose = F)$horizontal
      # planet <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(planet$RA, planet$Dec))$xaz[c(1,3)]
      # planet[1] <- planet[1]+180; if (planet[1]>360) { planet[1] <- planet[1]-360 }
      points(planet[1], planet[2], pch=19, cex=.5, col='darkred')
    }

  }

  polygon(c(0,360,360,0),c(0,0,-90,-90), border=NA, col=MESS::col.alpha('grey',0.7))
  abline(h=0, lwd=2)
}
