#' @noRd
events <- function(name) {
  return(switch(name,
                'December Solstice' = 'dS',
                'June Solstice' = 'jS',
                'Equinox' = 'eq',
                'Major Lunar Extreme (southern)' = 'sMjLX',
                'Minor Lunar Extreme (southern)' = 'smnLX',
                'Major Lunar Extreme (northern)' = 'nMjLX',
                'Minor Lunar Extreme (northern)' = 'nmnLX',
                name))
}


#' Creates a \emph{skyscapeR.object} for plotting of celestial objects at given epoch
#'
#' This function creates an object containing all the necessary information to
#' plot celestial objects/events unto the many plotting functions of \emph{skyscapeR}
#' package.
#' @param names The name(s) of the celestial object(s) or event(s) of interest.
#' These can be one of the following soli-lunar events: \emph{jS}, \emph{dS}, \emph{eq},
#'  \emph{nmnLX}, \emph{nMjLX},
#' \emph{smnLX}, \emph{sMjLX}, or the name of any star in the database. As shorthand, the names
#' \emph{sun} and \emph{moon} can be used to represent all the above solar and lunar events,
#' respectively. Alternatively, custom declination values can also be used.
#' @param epoch The year or year range (as an array) one is interested in.
#' @param loc (Optional) This can be either a vector with the latitude and longitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param col (Optional) The color for plotting, and differentiating these objects.
#' Defaults to red for all objects.
#' @param lty (Optional) Line type (see \code{\link{par}}) used for differentiation.
#' Only activated for single year epochs.
#' @param lwd (Optional) Line width (see \code{\link{par}}) used for differentiation.
#' Only activated for single year epochs.
#' @export
#' @examples
#' # Create a object with solar targets for epoch range 4000-2000 BC:
#' tt <- sky.objects('solar extremes', c(-4000,-2000))
#'
#' # Create an object with a few stars for same epoch:
#' tt <- sky.objects(c('Sirius', 'Betelgeuse', 'Antares'), c(-4000,-2000),
#' col=c('white', 'red', 'orange'))
#'
#' # Create an object with solstices and a custom declination value:
#' tt <- sky.objects(c('December Solstice','June Solstice', -13), c(-4000,-2000))
sky.objects = function(names, epoch, loc=FALSE, col = 'red', lty = 1, lwd = 1) {
  if (length(names) > 1 & length(col)==1) { col <- rep(col, length(names)) }
  if (length(names) > 1 & length(lty)==1) { lty <- rep(lty, length(names)) }
  if (length(names) > 1 & length(lwd)==1) { lwd <- rep(lwd, length(names)) }

  i <- which(names=='solar extremes')
  if (length(i)>0) {
    names <- names[-i]
    names <- c('December Solstice','June Solstice',names)
    aux <- col[i]; col <- col[-i]; col <- c(rep(aux,2),col)
    aux <- lty[i]; lty <- lty[-i]; lty <- c(rep(aux,2),lty)
    aux <- lwd[i]; lwd <- lwd[-i]; lwd <- c(rep(aux,2),lwd)
  }

  i <- which(names=='lunar extremes')
  if (length(i)>0) {
    names <- names[-i]
    names <- c('Major Lunar Extreme (southern)', 'Minor Lunar Extreme (southern)', 'Minor Lunar Extreme (northern)', 'Major Lunar Extreme (northern)',names)
    aux <- col[i]; col <- col[-i]; col <- c(rep(aux,4),col)
    aux <- lty[i]; lty <- lty[-i]; lty <- c(rep(aux,4),lty)
    aux <- lwd[i]; lwd <- lwd[-i]; lwd <- c(rep(aux,4),lwd)
  }

  N <- NROW(names)

  if (NROW(col)==1) { col <- rep(col,N) }
  if (NROW(lty)==1) { lty <- rep(lty,N) }
  if (NROW(lwd)==1) { lwd <- rep(lwd,N) }

  tt <- c()
  tt.col <- c()
  tt.lty <- c()
  tt.lwd <- c()

  for (i in 1:N) {
    # stars
    if (swephR::swe_fixstar2_ut(names[i], swephR::swe_julday(2000, 1, 1, 12.,1), 2048+16384+16)$serr=="") {
      aux <- array(NA, c(NROW(epoch),1))
      for (j in 1:NROW(epoch)) {
        aux[j,] <- star(names[i], epoch[j])$coord$Dec
      }
      colnames(aux) <- names[i]
      tt <- cbind(tt, aux)
      tt.col <- c(tt.col, col[i])
      tt.lty <- c(tt.lty, lty[i])
      tt.lwd <- c(tt.lwd, lwd[i])
      next
    } else
      # custom dec
      if (!is.na(suppressWarnings(as.numeric(names[i])))) {
        aux <- array(NA, c(NROW(epoch),1))
        aux[,1] <- rep(as.numeric(names[i]),NROW(epoch))
        colnames(aux) <- 'Custom Dec'
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, col[i])
        tt.lty <- c(tt.lty, lty[i])
        tt.lwd <- c(tt.lwd, lwd[i])
        next
      } else {

        # try calling the functions
        aux <- array(NA, c(NROW(epoch),1))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- do.call(events(names[i]), list(year=epoch[j], loc=loc))
        }
        colnames(aux) <- events(names[i])
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, col[i])
        tt.lty <- c(tt.lty, lty[i])
        tt.lwd <- c(tt.lwd, lwd[i])
      }

  }
  rownames(tt) <- epoch

  # check min and max decs
  if (NROW(epoch)==2) {
    ind <- which(colnames(tt) != 'Custom Dec')
    for (i in ind) { tt[,i] <- minmaxdec(events(colnames(tt)[i]), from=min(epoch), to=max(epoch)) }
    rownames(tt) <- c('min','max')
  }

  # return result
  object <- c()
  object$n <- NCOL(tt)
  object$decs <- tt
  object$epoch <- epoch
  object$loc <- loc
  object$col <- tt.col
  if (NROW(epoch)==1) {
    object$lty <- tt.lty
    object$lwd <- tt.lwd
  }
  class(object) <- "skyscapeR.object"
  return(object)
}


#' Declination of December Solstice for a given year
#'
#' This function calculates the declination of the sun
#' at December Solstice for a given year, based upon
#' obliquity estimation and corrected average parallax.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date}.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the sun.
#'  Defaults to 0.00224.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{jS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{spatial.equinox}},
#' \code{\link{parallax.corr}}
#' @examples
#' # December Solstice geocentric declination for year 4001 BC:
#' dS(-4000)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' dS(-4000, loc=50)
dS = function(year = skyscapeR.env$cur.year, loc = FALSE, parallax = 0.00224, altitude = 0, verbose=TRUE) {
  aux <- -obliquity(year)

  if (class(loc)=='logical') {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Declination of June Solstice for a given year
#'
#' This function calculates the declination of the sun
#' at June Solstice for a given year, based upon
#' obliquity estimation and corrected average parallax.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date}.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the sun.
#'  Defaults to 0.00224.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{spatial.equinox}},
#' \code{\link{parallax.corr}}
#' @examples
#' # June Solstice geocentric declination for year 4001 BC:
#' jS(-4000)
#'
#' # Topocentric declination for same year and latitude of 50 degres N:
#' jS(-4000, loc=50)
jS = function(year = skyscapeR.env$cur.year, loc = FALSE, parallax = 0.00224, altitude = 0, verbose=TRUE) {
  aux <- obliquity(year)

  if (class(loc)=='logical') {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Declination of the sun at the Equinox
#'
#' This function calculates the declination of the sun
#' at the astronomical equinox with corrected average parallax.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the sun.
#'  Defaults to 0.00224.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{jS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{spatial.equinox}}, \code{\link{parallax.corr}}
#' @examples
#' # Equinoctial geocentric declination:
#' eq()
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' eq(loc=50)
eq = function(loc=FALSE, parallax = 0.00224, altitude = 0, verbose=TRUE) {
  aux <- 0

  if (class(loc)=='logical' & verbose) {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Declination of the zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the zenith for a given location with
#' corrected average parallax. If this phenomena does
#' not occur at given location (i.e. if location is
#' outside the tropical band) the function returns a
#' \emph{NULL} value.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param parallax (Optional) Average parallax value for the sun.
#'  Defaults to 0.00224.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{antizenith}}, \code{\link{spatial.equinox}}, \code{\link{parallax.corr}}
#' @examples
#' # Zenith sun declination for Mexico City:
#' zenith(19.419)
#'
#' # There is no zenith sun phenomenon in London:
#' zenith(51.507)
zenith = function(loc, parallax = 0.00224, altitude = 0) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$metadata$georef[1]
  } else { lat <- loc[1] }

  if (lat > jS(loc=loc) | lat < dS(loc=loc)) {
    return(NULL)
  } else {
    aux <- lat
    aux <- aux - parallax.corr(parallax, loc, altitude)
    return(aux)
  }
}

#' Declination of the anti-zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the anti-zenith for a given location with
#' corrected average parallax. If this phenomena does
#' not occur at given location (i.e. if location is
#' outside the tropical band) the function returns a
#' \emph{NULL} value.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param parallax (Optional) Average parallax value for the sun.
#'  Defaults to 0.00224.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{spatial.equinox}}, \code{\link{parallax.corr}}
#' @examples
#' # Anti-zenith sun declination for Mexico City:
#' antizenith(19.419)
#'
#' # There is no anti-zenith sun phenomena in London:
#' antizenith(51.507)
antizenith = function(loc, parallax = 0.00224, altitude = 0) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$metaata$georef[1]
  } else { lat <- loc[1] }

  if (lat > jS(loc=loc) | lat < dS(loc=loc)) {
    return(NULL)
  } else {
    aux <- -lat
    aux <- aux - parallax.corr(parallax, loc, altitude)
    return(aux)
  }
}

#' Declination of the spatial equinox for a given location
#'
#' @param hor This should be a \emph{skyscapeR.horizon} object.
#' @param parallax (Optional) Average parallax value for the sun.
#'  Defaults to 0.00224.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}},
#' \code{\link{zenith}}, \code{\link{antizenith}}, \code{\link{parallax.corr}}
#' @examples
#' \dontrun{
#' hor <- createHWT(34.174051531543405, 110.81299818694872, name='Xipo, Lingbao')
#' spatial.equinox(hor)
#' }
spatial.equinox = function(hor, parallax = 0.00224) {
  if (class(hor)!='skyscapeR.horizon') stop('Full horizon profile required for spatial equinox calculation.')

  options(warn=-1)
  obj <- sky.objects('solar extremes', 2000, loc=hor)
  rise <- c(); set <- c();
  for (i in 1:2) {
    orb <- orbit(obj$decs[i], obj$loc, refraction=FALSE)
    forb <- approxfun(orb$az, orb$alt, yleft=-9999, yright=-9999)
    riset <- approxfun(hor$data$az,forb(hor$data$az)-hor$data$alt, yleft=-9999, yright=-9999)
    rise[i] <- uniroot(riset, interval=c(0, 180))$root
    set[i] <- uniroot(riset, interval=c(180, 360))$root
  }
  options(warn=0)
  out <- c()
  out$rise <- az2dec(mean(rise), loc=hor) - parallax.corr(parallax, hor, hor2alt(hor, mean(rise)))
  out$set <- az2dec(mean(set), loc=hor) - parallax.corr(parallax, hor, hor2alt(hor, mean(set)))

  return(out)
}

#' Declination of northern minor Lunar Extreme for a given year
#'
#' This function calculates the declination of the northern
#' minor Lunar Extreme for a given year, based upon
#' obliquity estimation and corrected average parallax.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date}.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the moon
#'  Defaults to 0.952.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{smnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Northern minor Lunar Extreme geocentric declination for year 2501 BC:
#' nmnLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' nmnLX(-2500, loc=50)
nmnLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  aux <- obliquity(year) - (5.145+0.145)

  if (class(loc)=='logical') {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}


#' Declination of southern minor Lunar Extreme for a given year
#'
#' This function calculates the declination of the southern
#' minor Lunar Extreme for a given year, based upon
#' obliquity estimation and corrected average parallax.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date}.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the moon
#'  Defaults to 0.952.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Southern minor Lunar Extreme geocentric declination for year 2501 BC:
#' smnLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' smnLX(-2500, loc=50)
smnLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  aux <- -(obliquity(year) - (5.145+0.145))

  if (class(loc)=='logical') {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}


#' Declination of northern major Lunar Extreme for a given year
#'
#' This function calculates the declination of the northern
#' major Lunar Extreme for a given year, based upon
#' obliquity estimation and corrected average parallax.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date}.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the moon
#'  Defaults to 0.952.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{sMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Northern major Lunar Extreme geocentric declination for year 2501 BC:
#' nMjLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' nMjLX(-2500, loc=50)
nMjLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  aux <- obliquity(year) + (5.145+0.145)

  if (class(loc)=='logical') {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}


#' Declination of southern major Lunar Extreme for a given year
#'
#' This function calculates the declination of the southern
#' major Lunar Extreme for a given year, based upon
#' obliquity estimation and corrected average parallax.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date}.
#' @param loc (Optional) This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object. If missing or \emph{FALSE},
#' function will output geocentric declination.
#' @param parallax (Optional) Average parallax value for the moon
#'  Defaults to 0.952.
#' @param altitude (Optional) Altitude of the sun. Defaults to 0 degrees.
#' @param verbose (Optional) Boolean to control output of warnings and messages.
#' Defaults to TRUE.
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{nMjLX}},
#' \code{\link{parallax.corr}}
#' @examples
#' # Southern major Lunar Extreme geocentric declination for year 2501 BC:
#' sMjLX(-2500)
#'
#' # Topocentric declination for same year and latitude of 50 degrees N:
#' sMjLX(-2500, loc=50)
sMjLX = function(year = skyscapeR.env$cur.year, loc=FALSE, parallax = 0.952, altitude = 0, verbose = TRUE) {
  aux <- -(obliquity(year) + (5.145+0.145))

  if (class(loc)=='logical') {
    if (verbose) { cat('No location given, therefore outputting geocentric declination.\n') }
  } else {
    aux <- aux - parallax.corr(parallax, loc, altitude)
  }
  return(aux)
}

#' Equinoctial Full Moons
#'
#' This function calculates the date, rise/set times, azimuth and declination
#' for sun and moon on the days of the Spring Full Moon (SFM) and Autumn Full
#' Moon (AFM), for a given year and location.
#' @param season (Optional) Either 'spring' or 'autumn'. Default is 'spring.
#' @param rise (Optional) Boolean to choose whether to calculate Equinoctial Full
#' Moon rises or sets. Defaults to \emph{TRUE}.
#' @param year Epoch(s) for which to do calculations. Can be either a single value (the year),
#' two values (range of years), or a vector of years.
#' @param loc This can be either a \emph{skyscapeR.horizon} object, or a vector with the
#' latitude, longitude and elevation of the site, in this order.
#' @param min.phase (Optional) Minimum lunar phase (between 0 and 1) for which a
#' moon is considered to be full. Defaults to 0.99.
#' @param refraction (Optional) Boolean for whether or not atmospheric refraction should be taken into account.
#' Defaults to \emph{TRUE}.
#' @param atm (Optional) Atmospheric pressure (in mbar). Only needed if \emph{refraction} is set to \emph{TRUE}.
#' Default is 1013.25 mbar.
#' @param temp (Optional) Atmospheric temperature (in Celsius). Only needed if \emph{refraction} is set to \emph{TRUE}.
#' Default is 15 degrees.
#' @param calendar (Optional) Calendar used to output dates. G for gregorian and J for julian. Defaults to \emph{Gregorian}.
#' @param timezone (Optional) Timezone for output of rising and setting time either as a known acronym
#' (e.g. "GMT", "CET") or a string with continent followed by country capital (e.g. "Europe/London"). See
#' \link{timezones} for details. Defaults to system timezone.
#' @import swephR parallel foreach doParallel
#' @export
#' @examples
#' # Spring Full Moon from a location in Portugal in the year 2018
#' EFM(year=2018, loc=c(35,-8,100))
#'
#' # Autumn Full Moons in the last three years
#' \dontrun{
#' EFM(season='autumn', year=c(2019,2021), loc=c(35,-8,100))
#' }
EFM <- function(season='spring', rise=T, year, loc, min.phase=.99, refraction, atm, temp, timezone, calendar) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }
  checkYear(year)

  if (length(year)==2) { year <- seq(year[1], year[2], 1) }

  suns <- data.frame(year=NA, date=NA, time=NA, azimuth=NA, declination=NA, stringsAsFactors = F); moons <- suns
  if (length(year) > 1) { pb <- txtProgressBar(max = length(year), style=3) }

  for (k in 1:length(year)) {
    jd0 <- time2jd(paste(year[k],'01','01',sep='/'), timezone, calendar)

    # find celestial crossovers (i.e. declination-based)
    jd <- seq(jd0, jd0+365, 1)
    sun <- body.position('sun', jd, loc=loc, refraction = refraction, atm=atm, temp=temp)$equatorial[,2]
    moon <- body.position('moon', jd, loc=loc, refraction = refraction, atm=atm, temp=temp)$equatorial[,2]
    phase <- moonphase(jd, timezone, calendar)

    ind <- which(c(0,diff(sign(sun[phase >= min.phase] - moon[phase >= min.phase]))) != 0)
    crossovers <- jd[phase >= min.phase][ind]

    if (season=='spring') { i <- which(month(jd2time(crossovers)) < 6) } else { i <- which(month(jd2time(crossovers)) > 6)  }

    # find visible crossovers (i.e. rising/setting azimuth-based)
    jd <- seq(crossovers[i]-33, crossovers[i]+33, 1/4)
    phase <- moonphase(jd, timezone, calendar); excl <- which(phase < min.phase); jd <- jd[-excl]
    jd <- time2jd(unique(substr(jd2time(jd, timezone, calendar),1,which(strsplit(jd2time(jd, timezone, calendar), "")[[1]]==" ")-1)), timezone, calendar)

    ss <- c()
    for (j in 1:length(jd)) {  ### TODO parallelise this
      if (rise) {
        sun <- riseset('sun', jd=jd[j], loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$rise
        moon <- riseset('moon', jd=jd[j], loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$rise
      } else {
        sun <- riseset('sun', jd=jd[j], loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$set
        moon <- riseset('moon', jd=jd[j], loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$set
      }
      ss[j] <- sign(sun$azimuth-moon$azimuth)
      jd[j] <- time2jd(paste(substr(jd2time(jd[j], timezone, calendar),1,which(strsplit(jd2time(jd[j], timezone, calendar), "")[[1]]==" ")-1), moon$time), timezone, calendar)
    }
    phase <- moonphase(jd, timezone, calendar); excl <- which(phase < min.phase); if (length(excl)>0) { jd <- jd[-excl]; ss <- ss[-excl] }

    if (sum(c(0,diff(ss))!=0)) { ind <- which(c(0,diff(ss))!=0) } else { stop('Error: No crossover found. Check with package maintainer.') }
    jd <- jd[ind]
    if (rise) {
      sun <- riseset('sun', jd=jd, loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$rise
      moon <- riseset('moon', jd=jd, loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$rise
    } else {
      sun <- riseset('sun', jd=jd, loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$set
      moon <- riseset('moon', jd=jd, loc=loc, calendar=calendar, timezone=timezone, refraction=refraction, atm=atm, temp=temp)$set
    }
    date <- jd2time(jd)

    dates <- substr(date,1,which(strsplit(date, "")[[1]]==" ")-1)
    suns[k,] <- c(year[k], dates, sun$time, sun$azimuth, sun$declination)
    moons[k,] <- c(year[k], dates, moon$time, moon$azimuth, moon$declination)

    if (length(year) > 1) { setTxtProgressBar(pb, k) }
  }

  suns$year <- as.numeric(suns$year); suns$azimuth <- as.numeric(suns$azimuth); suns$declination <- as.numeric(suns$declination)
  moons$year <- as.numeric(moons$year); moons$azimuth <- as.numeric(moons$azimuth); moons$declination <- as.numeric(moons$declination)
  out <- list(Event = paste0(season, ' full moon'), Moon= moons, Sun= suns)
  return(out)
}
