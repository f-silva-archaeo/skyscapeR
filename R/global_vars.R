################################################
skyscapeR.env <- new.env(parent = emptyenv())

# Timezone and Calendar
skyscapeR.env$timezone <- 'Europe/London'
skyscapeR.env$calendar <- 'Gregorian'

# Atmospheric Refraction
skyscapeR.env$refraction <- TRUE
skyscapeR.env$atm <- 1013.25
skyscapeR.env$temp <- 15

# Frame of Reference
skyscapeR.env$dec <- 'topo'

# star ephemeris source
skyscapeR.env$stars <- 'skyscapeR'

# Current Year (cannot be changed)
skyscapeR.env$cur.year <- as.numeric(format(Sys.Date(), "%Y"))

################################################

#' See and change the global variables used by skyscapeR
#'
#' @param timezone Timezone of input either as a known acronym (e.g. "GMT", "CET") or
#' a string with continent followed by country capital (e.g. "Europe/London"). See
#' \code{\link{timezones}} for details. Default is the system timezone
#' @param calendar Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Defaults to \emph{Gregorian}.
#' @param refraction  Whether atmospheric refraction is to be taken into account. Default is TRUE.
#' @param atm Atmospheric pressure for refraction calculation. Default is 1013.25 mbar.
#' @param temp Atmospheric temperature for refraction calculation. Default is 15 degrees.
#' @param dec Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. Defaults to topocentric.
#' @export
skyscapeR.vars = function(timezone, calendar, refraction, atm, temp, dec) {
  if (!missing(timezone)) { skyscapeR.env$timezone <- timezone }
  if (!missing(calendar)) { skyscapeR.env$calendar <- calendar }
  if (!missing(refraction)) { skyscapeR.env$refraction <- refraction }
  if (!missing(atm)) { skyscapeR.env$atm <- atm }
  if (!missing(temp)) { skyscapeR.env$temp <- temp }
  if (!missing(dec)) { skyscapeR.env$dec <- dec }

  aux <- ls(skyscapeR.env)
  return(mget(aux, skyscapeR.env))
}
#' Swap between swephR and skyscapeR versions of stellar ephemeris
#'
#' @param source Package name for the source of the stellar ephemeris data. Can be
#' either \emph{swephR} or \emph{skyscapeR}. Defaults to the latter.
#' @export
swapStars <- function(source='skyscapeR') {

  swephR_file <- system.file('ephemeris', 'sefstars.txt', package='swephR')
  swephRdata_file <- system.file('ephemeris', 'sefstars.txt', package='swephRdata')

  if (source == 'skyscapeR') {
    source_file <- system.file('ephemeris', 'skyscapeR-sefstars.txt', package='skyscapeR')
  } else if (source == 'swephR') {
    source_file <- system.file('ephemeris', 'swephR-sefstars.txt', package='skyscapeR')
  } else {
    stop('Source not recognised.')
  }

  # check if change is needed
  if (!identical(readLines(swephR_file, n=1), readLines(source_file, n=1))) {
      file.copy(source_file, swephR_file, overwrite = TRUE)
      file.copy(source_file, swephRdata_file, overwrite = TRUE)
      packageStartupMessage('Replaced swephR stellar ephemeris file.')
  }
}

globalVariables('star.names')



