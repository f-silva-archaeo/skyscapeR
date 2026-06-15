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
#' @param const (Optional) A string or vector of strings specifying constellations to display (e.g., "Ori", "modern").
#' @param show.az Logical. Whether to display azimuth labels (default: FALSE).
#' @param col.ground Character. The colour of the ground in the sketch (default: "darkgreen").
#'
#' @return A plot displaying the sky at a given time and location.
#' @export
#' @examples
#' sky.sketch(time='2019/01/10 18:51', loc=c(35,-8,100))
#'
#' # showing one constellation
#' sky.sketch(time='2019/01/10 18:51', loc=c(35,-8,100), const='Ori')
#'
#' # using a skyscapeR.horizon and three constellations
#' az <- c(0,90,180,270,360); alt <- c(0,5,5,0,0)
#' hor <- createHor(az, alt, 0.1, c(40.1,-8), 'Test')
#' sky.sketch(time='2019/01/10 18:51', loc=hor, const=c('Ori','Tau','Gem'))
#'
#' # showing all modern constellations
#' sky.sketch(time='2019/01/10 18:51', loc=hor, const='modern')
sky.sketch <- function(time, loc, sun=TRUE, moon=TRUE, planets=FALSE, const, xrange=c(30,150), yrange=c(-15,75), timezone, calendar, exagerate=TRUE, max.mag=6, show.az=FALSE, col.ground='darkgreen', atm, temp) {

  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }
  if (inherits(loc, "skyscapeR.horizon")) {
    hor <- loc
    loc <- hor$metadata$georef
  } else { hor <- NULL }

  if (abs(diff(xrange)) > 120) { message('xrange is larger than 120 degrees. This will cause a distorted view of the sky.')}
  if (abs(diff(yrange)) > 120) { message('yrange is larger than 120 degrees. This will cause a distorted view of the sky.')}

  ## prep stars info
  ss <- read.csv(paste0(system.file('ephemeris',package='skyscapeR'),'/skyscapeR-sefstars.txt'), header=FALSE)
  ind <- which(ss$V14>max.mag); if (length(ind)>0) { ss <- ss[-ind,] } # magnitude filter
  ss <- ss[-which(is.na(ss$V14)),]
  nstars <- NROW(ss)

  ## proc
  if (is(time,'character')) { jd <- time2jd(time, timezone, calendar) } else { jd <- time }

  par(mar=c(2,2,1,1))
  plot(-999,-999, xlim=xrange, ylim=yrange, xlab='', ylab='', axes=FALSE)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-40,360+40,0.1), lwd=0.2, labels=FALSE) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-40,360+40,0.5), lwd=0.2, labels=FALSE) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-40,360+40,1), lwd=0.5, labels=FALSE) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-40,360+40,5), lwd=0.5, labels=FALSE) }
  if (scale < 90 & scale >= 10) { axis(1, at=seq(-40,360+40,10), lwd=0.5, labels=FALSE) }

  if (show.az == TRUE) {
    if (scale >= 10 & scale < 45) { scale <- 10 }
    if (scale >= 45 & scale < 90) { scale <- 45 }
    if (scale >= 90) { scale <- 90 }
    axis(1, at = seq(-90,360+90,scale), labels = seq(-90,360+90,scale), lwd=0)

  } else {

    ll <- c("N","NE","E","SE","S","SW","W","NW","N","NE","E","SE","S","SW","W","NW","N","NE","E","SE","S","SW","W","NW","N")
    axis(1, at = seq(-360,720,by=45), labels = ll, lwd=0.5)
  }

  if (nstars > 0) {
    for (i in 1:nstars) {
      aux <- try(star(ss$V2[i], year=as.numeric(substr(time,1,4))), silent=TRUE)
      if (is(aux,'try-error')) next
      mag <- aux$app.mag
      star <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(aux$coord$RA, aux$coord$Dec))$xaz[c(1,3)]
      star[1] <- star[1]+180; if (star[1]>360) { star[1] <- star[1]-360 }
      points(star[1], star[2], pch=19, cex=mag2size(mag), col='black')
      points(star[1]-360, star[2], pch=19, cex=mag2size(mag), col='black')
      points(star[1]+360, star[2], pch=19, cex=mag2size(mag), col='black')
    }

    if (!missing(const)) {
      # constellation lines
      id <- suppressWarnings(constLines(const))

      for (j in 1:NROW(id)){
        for (i in 1:NROW(id[[j]])) {
          s1 <- star(id[[j]][i,1], year=as.numeric(substr(time,1,4)))
          s2 <- star(id[[j]][i,2], year=as.numeric(substr(time,1,4)))
          s1 <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(s1$coord$RA, s1$coord$Dec))$xaz[c(1,3)]
          s1[1] <- s1[1]+180; if (s1[1]>360) { s1[1] <- s1[1]-360 }
          s2 <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(s2$coord$RA, s2$coord$Dec))$xaz[c(1,3)]
          s2[1] <- s2[1]+180; if (s2[1]>360) { s2[1] <- s2[1]-360 }
          if (s2[1]-s1[1] > 350) {
            s2[1] <- s2[1]-360 }
          if (s2[1]-s1[1] < -350) {
            s2[1] <- s2[1]+360 }
          lines(c(s1[1],s2[1]), c(s1[2], s2[2]), lwd=0.5)
          lines(c(s1[1]-360,s2[1]-360), c(s1[2], s2[2]), lwd=0.5)
          lines(c(s1[1]+360,s2[1]+360), c(s1[2], s2[2]), lwd=0.5)
        }
      }
    }

  }

  if (sun) {
    sun <- as.numeric(body.position('sun', jd, loc = loc, refraction = TRUE, atm = atm, temp = temp, verbose = FALSE)$horizontal)
    plotrix::draw.circle(sun[1], sun[2], radius=0.5+1*exagerate, border='black', col='yellow', lwd=.5)
  }

  if (moon) {
    moon <- as.numeric(body.position('moon', jd, loc = loc, refraction = TRUE, atm = atm, temp = temp, verbose = FALSE)$horizontal)
    moonphase <- moonphase(jd)
    plotrix::draw.circle(moon[1], moon[2], radius=0.5+1*exagerate, border='black', col='white', lwd=.5)
    ## TODO add lunar phase
  }

  if (planets) {
    planets <- c('mercury','venus','mars', 'jupiter', 'saturn')
    cols <- c('darkblue','darkgrey','darkred', 'darkgreen', 'darkorange')
    for (i in 1:5) {
      planet <- body.position(planets[i], jd, loc = loc, refraction = TRUE, atm = atm, temp = temp, verbose = FALSE)$horizontal
      points(planet[1], planet[2], pch=19, cex=.5, col=cols[i])
    }

  }

  if (inherits(hor, "skyscapeR.horizon")) {
    line(hor$data$az, hor$data$alt)
    xx <- c(hor$data$az, rev(hor$data$az))
    yy <- c(hor$data$alt, rep(-90,NROW(hor$data$az)))
    polygon(xx, yy, col=col.ground, lwd=1.5)
  } else {
    polygon(c(-360,360+360,360+360,-360),c(0,0,-90,-90), border=NA, col=MESS::col.alpha('grey',0.7))
    abline(h=0, lwd=2)
  }
  box()
}


#' Create a planisphere view of the sky for a given epoch and location
#'
#' This function generates a planisphere, which is a circular star chart
#' that represents the visible sky for a given epoch and observer's latitude.
#'
#' @param epoch Numerical.The epoch (year) for which the planisphere is calculated.
#' @param lat Numerical.Latitude of the observer in degrees.
#' @param max.mag  Numerical.Maximum magnitude of stars to be displayed (default: 6).
#' @param show.ecliptic  Optional.Whether to display the ecliptic (TRUE or FALSE).
#' @param col.ecliptic Character.Colour of the ecliptic line (default: "orange").
#' @param show.equator Optional.Whether to display the celestial equator (True or False).
#' @param col.equator  Character.Colour of the celestial equator (default: "red").
#' @param star.col  Character.Colour of the stars (default: "black").
#'
#' @return A polar plot of the visible sky, showing stars, the celestial equator,
#'         and the ecliptic for the specified epoch and location.
#' @export
#' @examples
#' planisphere(2023, lat=50, max.mag=6, show.ecliptic=TRUE, show.equator=TRUE)
planisphere <- function(epoch, lat, max.mag=6, show.ecliptic=TRUE, col.ecliptic='orange', show.equator=TRUE, col.equator='red', star.col='black') {

  ss <- read.csv(paste0(system.file('ephemeris',package = 'swephR'),'/sefstars.txt'), header=FALSE)
  colnames(ss) <- c('name', 'identifier', 'ICRS', 'RA.1', 'RA.2', 'RA.3', 'Dec.1', 'Dec.2', 'Dec.3', 'pm.1', 'pm.2', 'radvel', 'plx', 'magV')
  ind <- which(ss$magV > max.mag | is.na(ss$magV)); ss <- ss[-ind,] # magnitude filter

  ra <- c(); dec <- c(); mag <- c(); type <- NA
  for (i in 1:NROW(ss)) {
    aux <- star(ss$identifier[i], year=epoch)
    ra[i] <- aux$coord$RA
    dec[i] <- aux$coord$Dec
    mag[i] <- aux$app.mag
  }

  ind <- which(dec < -90+lat)
  dec <- dec[-ind]
  ra <- ra[-ind]
  mag <- mag[-ind]

  labels <- c(seq(12,24,1),seq(1,11,1))
  label.pos <- (180 + 360/24*seq(0,length(labels)))/180*pi

  DD <- FALSE
  if (show.equator) {
    # draws equator (red path)
    plotrix::polar.plot(rep(90,360), seq(0,359,1), radial.lim=seq(0,90-(-90+lat),10),
                        start=-270, clockwise = TRUE, rp.type='p', line.col=col.equator,
                        labels = labels, label.pos = label.pos,
                        show.grid.labels=FALSE, main=ifelse(epoch>0,paste(epoch,'AD'),paste(abs(epoch)+1,'BC')))
    DD <- TRUE
  }

  if (show.ecliptic) {
    # draw ecliptic (yellow path)
    dat <- paste0(epoch,'/01/01 00:00:00')
    jd0 <- time2jd(dat)
    jd <- seq(jd0, jd0+366, 1)

    sun.pos <- body.position('sun', jd, loc=c(lat, 0, 0))$equatorial

    plotrix::polar.plot(90-sun.pos$Dec, sun.pos$RA, rp.type='p', line.col=col.ecliptic,
                        start=-270, clockwise = TRUE, radial.lim=seq(0,90-(-90+lat),10),
                        labels = labels, label.pos = label.pos,
                        show.grid.labels=FALSE, main=ifelse(epoch>0,paste(epoch,'AD'),paste(abs(epoch)+1,'BC')),
                        add=DD)
    DD <- TRUE
  }

  # draw stars
  plotrix::polar.plot(90-dec, ra, rp.type='s', point.symbols=19, cex=sapply(mag, mag2size),
                      start=-270, clockwise = TRUE, radial.lim=seq(0,90-(-90+lat),10),
                      labels = labels, label.pos = label.pos,
                      show.grid.labels=FALSE, main=ifelse(epoch>0,paste(epoch,'AD'),paste(abs(epoch)+1,'BC')),
                      add=DD, point.col=star.col)}

