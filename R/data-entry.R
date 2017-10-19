jd <- astrolibR::jdcnv(2000, 1, 1, 0.)  # J2000.0

#' Calculates declination from azimuth and altitude measurements
#'
#' This function calculates the declination corresponding to an
#' orientation , i.e. an azimuth. The altitude can either be given
#'  or, alternatively, if a \emph{skyscapeR.horizon} object is provided,
#'  the corresponding horizon altitude will be automatically retrieved.
#' This function is a wrapper for function \code{\link[astrolibR]{hor2eq}}
#' of package \emph{astrolibR}.
#' @param az Azimuth(s) for which to calculate declination(s). See examples below.
#' @param loc Location, can be either a \emph{skyscapeR.horizon} object or, alternatively,
#' a latitude.
#' @param alt Altitude of orientation. Optional, if left empty and a skyscapeR.object
#' is provided then this is will automatically retrieved from the horizon data via \code{\link{hor2alt}}
#' @param ... Any other parameters to be passed unto  \code{\link[astrolibR]{hor2eq}}.
#' @import astrolibR
#' @export
#' @seealso \code{\link[astrolibR]{hor2eq}}, \code{\link{hor2alt}}
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
  if (missing(alt) & class(loc) == 'skyscapeR.horizon') { alt <- hor2alt(hor, az) }

  prec <- max(nchar(sub('.*\\.', '', as.character(az))))
  dec <- round( astrolibR::hor2eq(alt, az, jd, hor$georef[1], hor$georef[2], precess_ = F, ...)$dec, prec)
  return(dec)
}

#' Retrieves horizon altitude for a given azimuth from a given horizon profile
#'
#' This function retrieves the horizon altitude for a given azimuth from
#' a previously created \emph{skyscapeR.horizon} object via spline interpolation.
#' @param hor A \emph{skyscapeR.horizon} object from which to retrieve horizon altitude.
#' @param az Array of azimuth(s) for which to retrieve horizon altitude(s).
#' @export
#' @import stats
#' @seealso \code{\link{createHor}}, \code{\link{download.HWT}}
#' @examples
#' hor <- download.HWT('HIFVTBGK')
#' hor2alt(hor, 90)
hor2alt = function(hor, az) {
  hh <- splinefun(hor$az, hor$alt)
  alt <- round(hh(az), 2)
  return(alt)
}

#' Estimates magnetics declination (diff. between true and magnetic
#' north) based on IGRF 12th gen model
#'
#' This function estimates the magnetic declination at a given location
#' and moment in time, using the \emph{12th generation International
#' Geomagnetic Reference Field (IGRF)} model. This function is a wrapper
#' for function \code{\link[oce]{magneticField}} of package \emph{oce}.
#' @param loc Location, can be either a \emph{skyscapeR.horizon} object or, alternatively,
#' a latitude.
#' @param date Date for which to calculate magnetic declination in the format: 'YYYY/MM/DD'
#' @export
#' @seealso \code{\link[oce]{magneticField}}
#' @examples
#' # Magnetic Declination for London on April 1st 2016:
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon )
#' mag.dec( loc, "2016/04/01" )
mag.dec = function(loc, date) {
  if (class(loc) == 'skyscapeR.horizon') { loc <- loc$georef }
  if (is.null(dim(loc))) { dim(loc) <- c(1,NROW(loc)) }
  aux <- oce::magneticField(loc[,2], loc[,1], as.POSIXlt(date, format="%Y/%m/%d"))$declination
  return(aux)
}
