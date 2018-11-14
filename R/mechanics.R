jd <- swephR::swe_julday(2000,1,1,12,1) # J2000.0
cur.year <- as.numeric(format(Sys.Date(), "%Y")) # current year
swephR::swe_set_ephe_path(system.file("ephemeris", "", package = "swephRdata"))


#' Compute astronomical parameters in the past or future based
#' on \emph{Laskar et al (2004)} tables
#'
#' This function calculates the obliquity, eccentricity and
#' longitude of the perihelion for a given year,
#' by interpolating the tables provided by \emph{Laskar et al (2004)}.
#' It is a modified version of function \code{\link[palinsol]{la04}}
#' of package \emph{palinsol}.
#' @param t Year for which to calculate the obliquity,
#' in years before of after 1950.
#' @param degree (Optional) Boolean that controls whether results
#' are output in degrees or radians. Defaults to TRUE.
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @references Michel Crucifix (2016). palinsol: Insolation for
#'  Palaeoclimate Studies. R package version 0.93.
#'  [CRAN page](https://CRAN.R-project.org/package=palinsol)
#' @noRd
#' @seealso \code{\link[palinsol]{astro}}, \code{\link{obliquity}}
#' @examples
#' Laskar04(-4000)
Laskar04 = function(t, degree=TRUE) {
  # Modified version of palinsol::la04
  data(Laskar04, envir=environment())

  tka = t/1000.
  if (tka>0)
  {
      F <-  floor(tka)
      ORB <- Laskar04$la04future[F+1, ]
      if (! (tka == F)) {
        D  <- tka - floor(tka)
        diff <- Laskar04$la04future[F+2, ] - ORB
        # note : if the diff in varpi is greater than pi,
        # this probably means that we have skipped 2*pi,
        # so we need to correct accordingly
        if (diff$varpi > pi) diff$varpi = diff$varpi - 2*pi
        if (diff$varpi < -pi) diff$varpi = diff$varpi + 2*pi
        #
        ORB <- ORB + D*diff
      }
  } else {
      F <-  floor(tka)
      ORB <- Laskar04$la04past[-F+1, ]
      if (! (tka == F)) {
        D  <- tka - F
        diff <- Laskar04$la04past[-F, ] - ORB
        # note : if the diff in varpi is greater than pi,
        # this probably means that we have skipped 2*pi,
        # so we need to correct accordingly
        if (diff$varpi > pi) diff$varpi = diff$varpi - 2*pi
        if (diff$varpi < -pi) diff$varpi = diff$varpi + 2*pi
        #
        ORB <- ORB + D*diff
      }
  }
  if (degree) {rad2deg <- 180/pi
  ORB['eps'] <- ORB['eps']*rad2deg
  ORB['varpi'] <- ORB['varpi']*rad2deg}

  # must return an array (0.92 -> 0.93)
  names <- c('eps','ecc','varpi')
  OUT = as.numeric(ORB[names]) ; names(OUT) <- names
  OUT
}




#' Computes obliquity of the ecliptic
#'
#' This function calculates the obliquity for a given epoch. It is a
#' wrapper for function \code{\link[swephR]{swe_calc}} of package \emph{swephR}.
#' @param year Year for which to calculate the obliquity.
#' Defaults to present year as given by Sys.Date()
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc}}
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @examples
#' #' # Obliquity for year 3999 BC:
#' obliquity(-4000)
obliquity = function(year = cur.year) {
  jd <- swephR::swe_julday(year,1,1,12,1) # J2000.0
  aux <- swephR::swe_calc(jd, -1, 0)$xx[1]
  return(aux)
}



#' Calculate visible path of celestial object at given location
#'
#' This function calculates the visible path of a celestial
#' object from any location on earth. It outputs a \emph{skyscapeR.orbit}
#' object, which includes AZ and ALT information.
#' @param dec Declination of object.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param res The resolution (in degrees of RA) with which to calculate the path.
#' @param refraction (Optional) If set to TRUE it will calculate apparent rather
#' than true altitudes. Default is TRUE.
#' @param ...  Any other parameters to be passed unto \code{\link[swephR]{swe_azalt}}.
#' @import swephR
#' @export
#' @examples
#' # Visible path of sun on June Solstice on year 3999 BC from London:
#' sun.dec <- jS(-4000)
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon )
#' path <- orbit(sun.dec, loc)
#' plot(path$az, path$alt, ylim=c(0,90), type='l', xlab='AZ', ylab='ALT', col='red', lwd=2)
orbit = function(dec, loc, res=0.5, refraction=T, ...) {
  if (class(loc)=='skyscapeR.horizon') { loc <- loc$metadata$georef }

  ra <- seq(0, 360, by=res)
  aux <- array(NA, c(NROW(ra),2))

  for (i in 1:NROW(ra)) {
    tmp <- eq2horFS(ra[i], dec, jd, loc, refraction, ...)

    aux[i,] <- c(tmp$az,tmp$alt)
  }

  # return result
  orbit <- c()
  orbit$az <- aux[,1]
  orbit$alt <- aux[,2]
  orbit$dec <- dec
  orbit$georef <- loc
  class(orbit) <- "skyscapeR.orbit"
  return(orbit)
}



#' Returns the azimuth of the sun at a given time from a specific location
#'
#' This function returns the azimuth of the sun at a given time and location,
#' useful for data reduction of theodolite measurements using the sunsight
#' technique (\code{\link{reduct.theodolite}}).
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param time String containing the date and time in the following format:
#' "YYYY-MM-DD HH:MM:SS"
#' @param timezone Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London").
#' @param limb (Optional) Measured limb of the sun. Options are \emph{left}, \emph{right}.
#' If missing the centre of the sun will be output.
#' @import swephR
#' @export
#' @seealso \code{\link{reduct.theodolite}}
#' @examples
#' sunAz(c(52,-3), '2017-10-04 12:32:14', 'Europe/London')
sunAz = function(loc, time, timezone, limb) {
  if (class(loc)=='skyscapeR.horizon') { loc <- loc$metadata$georef }
  if (is.null(dim(loc))) { dim(loc) <- c(1, NROW(loc)) }

  az <- c()
  for (i in 1:NROW(loc)) {
    pb.date <- as.POSIXct(time[i], timezone[i])
    UT <- format(pb.date, tz="UTC",usetz=TRUE)
    UT <- as.POSIXlt(UT, 'UTC')

    swephR::swe_set_topo(loc[i,1], loc[i,2], 0)
    jd <- swephR::swe_julday(UT$year+1900, UT$mon+1, UT$mday, UT$hour+UT$min/60+UT$sec/3600,1)
    ss <- swephR::swe_calc_ut(jd, 0, 32*1024+2048)
    xin <- ss$xx[1:2]
    az[i] <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), 1013.25, 15, xin)$xaz[1]+180
    if (az[i] > 360) { az[i] <- az[i]-360 }

    if (!missing(limb)) {
      if (limb=="left") { az[i] <- az[i] - 32/60/2 }
      if (limb=="right") { az[i] <- az[i] + 32/60/2 }
    }
  }

  return(az)
}
