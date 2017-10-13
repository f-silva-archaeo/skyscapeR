jd <- astrolibR::jdcnv(2000, 1, 1, 0.)  # J2000.0
cur.year <- as.numeric(format(Sys.Date(), "%Y")) # current year


#' Computes obliquity based on \emph{Laskar (2004)} tables
#'
#' This function calculates the obliquity for a given year,
#' by interpolating the tables provided by \emph{Laskar (2004)}.
#' It is a wrapper for function \code{\link[palinsol]{la04}}
#' of package 'palinsol'.
#' @param year Year for which to calculate the obliquity.
#' Defaults to present year as given by Sys.Date()
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @export
#' @import palinsol
#' @seealso \code{\link[palinsol]{astro}}, \code{\link{dS}}, \code{\link{jS}}
#' @examples
#' #' # Obliquity for year 3999 BC:
#' obliquity(-4000)
obliquity = function(year = cur.year) {
  # loadNamespace("palinsol")
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
#' Defaults to present year as given by \emph{Sys.Date()}.
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
#' Defaults to present year as given by \emph{Sys.Date()}.
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
#' Defaults to present year as given by \emph{Sys.Date()}.
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
#' Defaults to present year as given by \emph{Sys.Date()}.
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
#' Defaults to present year as given by \emph{Sys.Date()}.
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
#' Defaults to present year as given by \emph{Sys.Date()}.
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
#' @param bh \emph{NULL} parameter. Can be left empty.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' eq()
eq = function(bh=NULL) {
  return(0)
}

#' Declination of the zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the zenith for a given location. If
#'  this phenomena does not occur at given location
#'  (i.e. if location is outside the tropical band)
#'  the function returns a \emph{NULL} value.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
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
#' band) the function returns a \emph{NULL} value.
#' @param loc This can be either the latitude of the l
#' ocation, or a \emph{skyscapeR.horizon} object.
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




#' Time shift object of class \emph{skyscapeR.star}
#'
#' This function calculates the coordinates (RA and DEC)
#' of a star for any time in the past. It uses functions
#' \code{\link[astrolibR]{precess}}, \code{\link[astrolibR]{co_nutate}}
#' and \code{\link[astrolibR]{co_aberration}}.
#' @param star Object created using \code{\link{star}}.
#' @param year Year for which to calculate the coordinates.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link[astrolibR]{precess}}, \code{\link[astrolibR]{co_nutate}},
#' \code{\link[astrolibR]{co_aberration}}
#' @examples
#' # Time shift data on Sirius to 2500 BC:
#' Sirius <- star('Sirius')
#' Sirius.2500BC <- palaeo.star(Sirius, -2501)
palaeo.star = function(star, year = cur.year) {
  if (class(star)=="skyscapeR.star") {
    # precession
    trash <- capture.output(coor <- astrolibR::precess(star$ra,star$dec, 2000, year), file=NULL)

    # proper motion
    adj.pm <- star$proper.motion*(year-2000)/1000
    coor$ra <- coor$ra + astrolibR::ten(0,0,adj.pm[1])*sign(adj.pm[1])
    coor$dec <- coor$dec + astrolibR::ten(0,0,adj.pm[2])*sign(adj.pm[2])

    # nutation
    adj.nu <- astrolibR::co_nutate(astrolibR::jdcnv(year,1,1,0), star$ra, star$dec)
    coor$ra <- coor$ra + astrolibR::ten(0,0,adj.nu$d_ra)*sign(adj.nu$d_ra)
    coor$dec <- coor$dec + astrolibR::ten(0,0,adj.nu$d_dec)*sign(adj.nu$d_dec)

    # aberration
    adj.ab <- astrolibR::co_aberration(astrolibR::jdcnv(year,1,1,0), star$ra, star$dec, pracma::deg2rad(obliquity(year)))
    coor$ra <- coor$ra + astrolibR::ten(0,0,adj.ab$d_ra)*sign(adj.ab$d_ra)
    coor$dec <- coor$dec + astrolibR::ten(0,0,adj.ab$d_dec)*sign(adj.ab$d_dec)

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
#' object from any location on earth. It outputs a \emph{skyscapeR.orbit}
#' object, which includes AZ and ALT information.
#' @param dec Declination of object.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param res The resolution (in degrees of RA) with which to calculate the path.
#' @param ...  Any other parameters to be passed unto \code{\link[astrolibR]{eq2hor}}.
#' @export
#' @examples
#' # Visible path of sun on June Solstice on year 3999 BC from London:
#' sun.dec <- jS(-4000)
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon )
#' path <- orbit(sun.dec, loc)
#' plot(path$az, path$alt, ylim=c(0,90), type='l', xlab='AZ', ylab='ALT', col='red', lwd=2)
orbit = function(dec, loc, res=0.5, ...) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
    lon <- loc$georef[2]
  } else {
    lat <- loc[1]
    lon <- loc[2]
  }

  ra <- seq(0, 360, by=res)
  aux <- array(NA, c(NROW(ra),2))

  for (i in 1:NROW(ra)) {
    tmp <- eq2horFS(ra[i], dec, jd, lat, lon, precess_=F, ...)
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


#' Creates a \emph{skyscapeR.object} for plotting of celestial objects at given epoch
#'
#' This function creates an object containing all the necessary information to
#' plot celestial objects/events unto the many plotting functions of \emph{skyscapeR}
#' package.
#' @param names The name(s) of the celestial object(s) or event(s) of interest.
#' These can be one of the following soli-lunar events: \emph{jS}, \emph{dS}, \emph{eq}, \emph{nmnLX}, \emph{nMjLX},
#' \emph{smnLX}, \emph{sMjLX}, or the name of any star in the database. As shorthand, the names
#' \emph{sun} and \emph{moon} can be used to represent all the above solar and lunar events,
#' respectively. Alternatively, custom declination values can also be used.
#' @param epoch The year or year range (as an array) one is interested in.
#' @param col (Optional) The colour for plotting, and differentiating these objects.
#' Defaults to red for all objects.
#' @param lty (Optional) Line type (see \code{\link{par}}) used for differentiation.
#' Only activated for single year epochs.
#' @param lwd (Optional) Line width (see \code{\link{par}}) used for differentiation.
#' Only activated for single year epochs.
#' @export
#' @examples
#' # Create a object with solar targets for epoch range 4000-2000 BC:
#' tt <- object('sun', c(-4000,-2000))
#'
#' # Create an object with a few stars for same epoch:
#' tt <- object(c('Sirius', 'Betelgeuse', 'Antares'), c(-4000,-2000), col=c('white', 'red', 'orange'))
#'
#' # Create an object with solstices and a custom declination value:
#' tt <- object(c('dS','jS', -13), c(-4000,-2000))
object = function(names, epoch, col = 'red', lty = 1, lwd = 1) {
  N <- NROW(names)

  if (NROW(col)==1) { col <- rep(col,N) }
  if (NROW(lty)==1) { lty <- rep(lty,N) }
  if (NROW(lwd)==1) { lwd <- rep(lwd,N) }

  tt <- c()
  tt.col <- c()
  tt.lty <- c()
  tt.lwd <- c()

  options(warn=-2)
  for (i in 1:N) {
    if (is.na(as.numeric(names[i]))) {
      # sun and moon shorthand
      if (pracma::strcmp(names[i], 'sun')) {
        aux <- array(NA, c(NROW(epoch),3))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- c(dS(epoch[j]), jS(epoch[j]), eq(epoch[j]))
        }
        colnames(aux) <- c('dS','jS','eq')
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, rep(col[i], 3))
        tt.lty <- c(tt.lty, rep(lty[i], 3))
        tt.lwd <- c(tt.lwd, rep(lwd[i], 3))
        next
      }

      if (pracma::strcmp(names[i], 'moon')) {
        aux <- array(NA, c(NROW(epoch),4))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- c(sMjLX(epoch[j]), smnLX(epoch[j]), nmnLX(epoch[j]), nMjLX(epoch[j]))
        }
        colnames(aux) <- c('sMjLX', 'smnLX', 'nmnLX', 'nMjLX')
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, rep(col[i], 4))
        tt.lty <- c(tt.lty, rep(lty[i], 4))
        tt.lwd <- c(tt.lwd, rep(lwd[i], 4))
        next
      }

      if (pracma::strcmp(names[i], 'sunandmoon')) {
        aux <- array(NA, c(NROW(epoch),7))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- c(dS(epoch[j]), jS(epoch[j]), eq(epoch[j]), sMjLX(epoch[j]), smnLX(epoch[j]), nmnLX(epoch[j]), nMjLX(epoch[j]))
        }
        colnames(aux) <- c('dS','jS','eq','sMjLX', 'smnLX', 'nmnLX', 'nMjLX')
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, rep(col[i], 7))
        tt.lty <- c(tt.lty, rep(lty[i], 7))
        tt.lwd <- c(tt.lwd, rep(lwd[i], 7))
        next
      }

      ####### TO DO ::: include individual solar and lunar targets

      # stars  :::: NEEDS CHANGING TO OUTPUT MIN AND MAX DEC
      data(stars, envir=environment())
      if (sum(sapply(as.character(stars$NAME), pracma::strcmp, s2=names[i]))) {
        aux <- array(NA, c(NROW(epoch),1))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- palaeo.star(do.call(star, list(names[i])), epoch[j])$dec
        }
        colnames(aux) <- names[i]
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, col[i])
        tt.lty <- c(tt.lty, lty[i])
        tt.lwd <- c(tt.lwd, lwd[i])
        next
      }
    }

    # custom dec
    if (!is.na(as.numeric(names[i]))) {
      aux <- array(NA, c(NROW(epoch),1))
      aux[,1] <- rep(as.numeric(names[i]),NROW(epoch))
      colnames(aux) <- paste0('Custom Dec:',names[i])
      tt <- cbind(tt, aux)
      tt.col <- c(tt.col, col[i])
      tt.lty <- c(tt.lty, lty[i])
      tt.lwd <- c(tt.lwd, lwd[i])
      next
    }

    # if none of the above fit the bill try calling the functions
    aux <- array(NA, c(NROW(epoch),1))
    for (j in 1:NROW(epoch)) {
      aux[j,] <- do.call(names[i], list(epoch[j]))
    }
    colnames(aux) <- names[i]
    tt <- cbind(tt, aux)
    tt.col <- c(tt.col, col[i])
    tt.lty <- c(tt.lty, lty[i])
    tt.lwd <- c(tt.lwd, lwd[i])
  }
  options(warn=0)
  rownames(tt) <- epoch


  # return result
  object <- c()
  object$n <- NCOL(tt)
  object$decs <- tt
  object$epoch <- epoch
  object$col <- tt.col
  if (NROW(epoch)==1) {
    object$lty <- tt.lty
    object$lwd <- tt.lwd
  }
  class(object) <- "skyscapeR.object"
  return(object)
}

#' Returns the azimuth of the sun at a given time from a specific location
#'
#' This function returns the azimuth of the sun at a given time and location,
#' useful for data reduction of theodolite mesaurements using the sunsight
#' technique (\code{\link{reduction.theod}}).
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param time String containing the date and time in the following format:
#' "YYYY-MM-DD HH:MM:SS"
#' @param timezone Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London").
#' @export
#' @seealso \code{\link{reduction.theod}}
#' @examples
#' sunAz(c(52,-3), '2017-10-04 12:32:14', 'Europe/London')
sunAz = function(loc, time, timezone) {
  if (class(loc)=='skyscapeR.horizon') { loc <- loc$georef }

  pb.date <- as.POSIXct(time, timezone)
  UT <- format(pb.date, tz="UTC",usetz=TRUE)
  UT <- as.POSIXlt(UT, 'UTC')
  jd <- astrolibR::jdcnv(UT$year+1900, UT$mon+1, UT$mday, UT$hour+UT$min/60+UT$sec/3600)
  ss <- astrolibR::sunpos(jd)
  az <- eq2horFS(ss$ra, ss$dec, jd, loc[1], loc[2], precess_ = F)$az

  return(az)
}
