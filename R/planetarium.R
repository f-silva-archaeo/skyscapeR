#' Create a simplistic sketch of the sky at a given moment in time
#'
#' @param time Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (e.g. "GMT", "CET") or
#' a string with continent followed by country capital (e.g. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. If not given the value set
#' by \code{\link{skyscapeR.vars}} will be used instead.
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param xrange Range of azimuths to display, preferably no larger than 120 degrees.
#' @param yrange Range of altitudes to display, preferably no larger than 120 degrees.
#' @param sun (Optional) Boolean on whether the sun should be displayed. Defaults to \emph{TRUE}.
#' @param moon (Optional) Boolean on whether the moon should be displayed. Defaults to \emph{TRUE}.
#' @param planets (Optional) Boolean on whether the visible planets should be displayed. Defaults to \emph{FALSE}.
#' @param exagerate (Optional) Boolean on whether the size of the sun, moon and planets should be exaggerated,
#' which can be useful when attempting wider viewing angles. Defaults to \emph{TRUE}.
#' @param max.mag (Optional) Maximum magnitude of stars to consider. Default is 6.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temperature for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @export
#' @examples
#' sky.sketch(time='2019/01/10 18:51', loc=c(35,-8,100))
sky.sketch <- function(time, timezone, calendar, xrange=c(30,150), yrange=c(-45,45), sun=T, moon=T, planets=F, exagerate=T, max.mag=6, loc, atm, temp) {

  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (abs(diff(xrange)) > 120) { message('xrange is larger than 120 degrees. This will cause a distorted view of the sky.')}
  if (abs(diff(yrange)) > 120) { message('yrange is larger than 120 degrees. This will cause a distorted view of the sky.')}

  ## prep stars info
  ss <- skyscapeR.env$sefstars
  ind <- which(ss$mag.V>max.mag); ss <- ss[-ind,] # magnitude filter
  nstars <- NROW(ss)

  ## proc
  if (class(time)=='character') { jd <- time2jd(time, timezone, calendar) } else { jd <- time }

  par(mar=c(2,2,1,1))
  plot(-999,-999, xlim=xrange, ylim=yrange, xlab='', ylab='')

  if (nstars > 0) {
    for (i in 1:nstars) {
      aux <- try(star(ss$nomenclature.name[i], year=as.numeric(substr(time,1,4))), silent=T)
      if (class(aux) == 'try-error') next
      mag <- aux$app.mag
      star <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(aux$coord$RA, aux$coord$Dec))$xaz[c(1,3)]
      star[1] <- star[1]+180; if (star[1]>360) { star[1] <- star[1]-360 }
      points(star[1], star[2], pch=19, cex=mag2size(mag), col='black')
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
      points(planet[1], planet[2], pch=19, cex=.5, col='darkred')
    }

  }

  polygon(c(0,360,360,0),c(0,0,-90,-90), border=NA, col=MESS::col.alpha('grey',0.7))
  abline(h=0, lwd=2)
}
