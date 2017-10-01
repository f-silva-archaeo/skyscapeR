devtools::use_package("astrolibR")
devtools::use_package("palinsol")
devtools::use_package("pracma")
jd <- astrolibR::jdcnv(2000, 1, 1, 0.)  # J2000.0
cur.year <- as.numeric(format(Sys.Date(), "%Y")) # current year



#' Calculates declination from azimuth and altitude measurements
#'
#' This function calculates the declination corresponding to an
#' orientation , i.e. an azimuth. The altitude can either be given
#'  or, alternatively, if a skyscapeR.horizon object is provided,
#'  the corresponding horizon altitude will be automatically retrieved.
#' This function is a wrapper for function \code{\link[astrolibR]{hor2eq}}
#' of package 'astrolibR'.
#' @param az Azimuth(s) for which to calculate declination(s). See examples below.
#' @param loc Location, can be either a skyscapeR.horizon object or, alternatively,
#' a latitude.
#' @param alt Altitude of orientation. Optional, if left empty and a skyscapeR.object
#' is provided then this is will automaticallty retrieved from the horizon data.
#' @param ... Any other parameters to be passed unto  \code{\link[astrolibR]{hor2eq}}.
#' @export
#' @seealso \code{\link[astrolibR]{hor2eq}}
#' @examples
#' hor <- download.HWT('HIFVTBGK')
#'
#' dec <- az2dec(92, hor)
#' dec <- az2dec(92, hor, alt=4)
#'
#' # Can also be used for an array of azimuths:
#' decs <- az2dec( c(87,92,110), hor )
az2dec = function(az, loc, alt, ...){
  if (class(loc) != 'skyscapeR.horizon') { hor <- c(); hor$georef <- c(loc, 0) } else { hor <- loc }
  if (missing(alt) & class(loc) == 'skyscapeR.horizon') {
    hh <- splinefun(hor$az, hor$alt)
    alt <- hh(az)
  }
  dec <- round( astrolibR::hor2eq(alt, az, jd, hor$georef[1], hor$georef[2], precess_ = F, ...)$dec, 2); names <- ''
  return(dec)
}



#' Computes obliquity based on Laskar (2004) tables
#'
#' This function calculates the obliquity for a given year,
#' by interpolating the tables provided by Laskar (2004).
#' It is a wrapper for function \code{\link[palinsol]{la04}}
#' of package 'palinsol'.
#' @param year Year for which to calculate the obliquity.
#' Defaults to present year as given by Sys.Date()
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, Astron.
#'  Astroph., 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @export
#' @seealso \code{\link[palinsol]{astro}}, \code{\link{dS}}, \code{\link{jS}}
#' @examples
#' #' # Obliquity for year 3999 BC:
#' obliquity(-4000)
obliquity = function(year = cur.year) {
  aux <- palinsol::la04(year-1950, degree=T)
  names(aux) <- NULL
  return(aux[1])
}


#' Declination of December Solstice for a given year
#'
#' This function calculates the declination of the sun
#' at December Solstice for a given year, based upon
#' obliquity estimation.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{jS}}, \code{\link{eq}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' # December Solstice declination for year 3999 BC:
#' dS(-4000)
dS = function(year = cur.year) {
  aux <- obliquity(year)
  return(-aux)
}


#' Declination of June Solstice for a given year
#'
#' This function calculates the declination of the sun
#' at June Solstice for a given year, based upon
#' obliquity estimation.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' # June Solstice declination for year 3999 BC:
#' jS(-4000)
jS = function(year = cur.year) {
  aux <- obliquity(year)
  return(aux)
}


#' Declination of northerm minor Lunar Extreme for a given year
#'
#' This function calculates the declination of the northern
#' minor Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{smnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}}
#' @examples
#' # Northern minor Lunar Extreme declination for year 2499 BC:
#' nmnLX(-2500)
nmnLX = function(year = cur.year) {
  aux <- obliquity(year) - 5.145
  return(aux)
}


#' Declination of southern minor Lunar Extreme for a given year
#'
#' This function calculates the declination of the southern
#' minor Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}}
#' @examples
#' # Southern minor Lunar Extreme declination for year 2499 BC:
#' smnLX(-2500)
smnLX = function(year = cur.year) {
  aux <- obliquity(year) - 5.145
  return(-aux)
}


#' Declination of northern major Lunar Extreme for a given year
#'
#' This function calculates the declination of the northern
#' major Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{sMjLX}}
#' @examples
#' # Northern major Lunar Extreme declination for year 2499 BC:
#' nMjLX(-2500)
nMjLX = function(year = cur.year) {
  aux <- obliquity(year) + 5.145
  return(aux)
}


#' Declination of southern major Lunar Extreme for a given year
#'
#' This function calculates the declination of the southern
#' major Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{nMjLX}}, \code{\link{smnLX}}
#' @examples
#' # Southern major Lunar Extreme declination for year 2499 BC:
#' sMjLX(-2500)
sMjLX = function(year = cur.year) {
  aux <- obliquity(year) + 5.145
  return(-aux)
}

#' Declination of sun at the equinoxes
#'
#' This function always returns a value of zero, which is the
#' declination of the sun on the day of the (astronomical)
#' equinoxes.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' eq()
eq = function() {
  return(0)
}

#' Declination of the zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the zenith for a given location. If
#'  this phenomena does not occur at given location
#'  (i.e. if location is outside the tropical band)
#'  the function returns a NULL value.
#' @param loc This can be either the latitude of the l
#' ocation, or a skyscapeR.horizon object
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{antizenith}}
#' @examples
#' # Zenith sun declination for Mexico City:
#' zenith(19.419)
#'
#' # There is no zenith sun phenomena in London:
#' zenith(51.507)
zenith = function(loc) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
  } else { lat <- loc }

  if (lat > jS() | lat < dS()) {
    return(NULL)
  } else { return(lat) }
}

#' Declination of the anti-zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the anti-zenith, or nadir, for a given
#' location. If this phenomena does not occur at given
#' location (i.e. if location is outside the tropical
#' band) the function returns a NULL value.
#' @param loc This can be either the latitude of the l
#' ocation, or a skyscapeR.horizon object
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{zenith}}
#' @examples
#' # Anti-zenith sun declination for Mexico City:
#' antizenith(19.419)
#'
#' # There is no anti-zenith sun phenomena in London:
#' antizenith(51.507)
antizenith = function(loc) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
  } else { lat <- loc }

  if (lat > jS() | lat < dS()) {
    return(NULL)
  } else { return(-lat) }
}




#' Time shift object of class 'skyscapeR.star'
#'
#' This function calculates the coordinates (RA and DEC)
#' of a star for any time in the past. It uses functions
#' \code{\link[astrolibR]{precess}}, \code{\link[astrolibR]{co_nutate}}
#' and \code{\link[astrolibR]{co_aberration}}.
#' @param star Object created using \code{\link{star}}
#' @param year Year for which to calculate the coordinates.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link[astrolibR]{precess}}, \code{\link[astrolibR]{co_nutate}},
#' \code{\link[astrolibR]{co_aberration}}
#' @examples
#' # Time shift data on Sirius to 2500 BC:
#' Sirius <- star('Sirius')
#' Sirius.2500BC <- palaeo.star(Sirius, -2501)
palaeo.star = function(star, year = cur.year) {
  if (class(star)=="star") {
    # precession
    trash <- capture.output(coor <- astrolibR::precess(star$ra,star$dec, 2000, year), file=NULL)

    # proper motion
    adj.pm <- star$proper.motion*(year-2000)/1000
    coor$ra <- coor$ra + ten(0,0,adj.pm[1])*sign(adj.pm[1])
    coor$dec <- coor$dec + ten(0,0,adj.pm[2])*sign(adj.pm[2])

    # nutation
    adj.nu <- astrolibR::co_nutate(jdcnv(year,1,1,0), star$ra, star$dec)
    coor$ra <- coor$ra + ten(0,0,adj.nu$d_ra)*sign(adj.nu$d_ra)
    coor$dec <- coor$dec + ten(0,0,adj.nu$d_dec)*sign(adj.nu$d_dec)

    # aberration
    adj.ab <- astrolibR::co_aberration(jdcnv(year,1,1,0), star$ra, star$dec, deg2rad(obliquity(year)))
    coor$ra <- coor$ra + ten(0,0,adj.ab$d_ra)*sign(adj.ab$d_ra)
    coor$dec <- coor$dec + ten(0,0,adj.ab$d_dec)*sign(adj.ab$d_dec)

    # output
    star$ra <- coor$ra
    star$dec <- coor$dec
    star$epoch <- as.character(year)

    return(star)
  } else { stop('No object of class skyscaper.star detected. Please check man pages: ?palaeosky')}

}


#' Calculate visible path of celestial object at given location
#'
#' This function calculates the visibile path of a celestial
#' object from any location on earth. It outputs a 'skyscapeR.orbit'
#' object, which includes AZ and ALT information.
#' @param dec Declination of object.
#' @param loc Locations, either a skyscapeR.object or a vector
#' containing the latitude and longitude of location, in this order.
#' @export
#' @examples
#' # Visible path of sun on June Solstice on year 3999 BC from London:
#' sun.dec <- jS(-4000)
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon )
#' path <- orbit(sun.dec, loc)
#' plot(path$az, path$alt, ylim=c(0,90), type='l', xlab='AZ', ylab='ALT', col='red', lwd=2)
orbit = function(dec, loc) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
    lon <- loc$georef[2]
  } else {
    lat <- loc[1]
    lon <- loc[2]
  }

  ra <- seq(0,360, by=0.1)
  aux <- array(NA,c(NROW(ra),2))

  for (i in 1:NROW(ra)) {
    trash <- capture.output(tmp <- astrolibR::eq2hor(ra[i],dec,jd,lat,lon),file=NULL)
    aux[i,] <- c(tmp$az,tmp$alt)
  }

  # return result
  orbit <- c()
  orbit$az <- aux[,1]
  orbit$alt <- aux[,2]
  orbit$dec <- dec
  orbit$georef <- c(lat,lon)
  class(orbit) <- "skyscapeR.orbit"
  return(orbit)
}
