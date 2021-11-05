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
  folder <- .libPaths()
  if (source=='skyscapeR') {
    # copy original swephR file with backup in skyscapeR folder
    file.copy(paste0(system.file('ephemeris',package='swephR'),'/sefstars.txt'), paste0(system.file('ephemeris',package='swephR'),'/sefstars.bkp'), overwrite=T)
    file.copy(paste0(system.file('ephemeris',package='swephRdata'),'/sefstars.txt'), paste0(system.file('ephemeris',package='swephRdata'),'/sefstars.bkp'), overwrite=T)
    file.copy(paste0(system.file('ephemeris',package='swephR'),'/sefstars.txt'), paste0(system.file('ephemeris',package='skyscapeR'),'/sefstars.bkp'), overwrite=T)

    # replace with new one
    file.copy(paste0(system.file('ephemeris',package='skyscapeR'),'/skyscapeR-sefstars.txt'), paste0(system.file('ephemeris',package='swephR'),'/sefstars.txt'), overwrite=T)
    file.copy(paste0(system.file('ephemeris',package='skyscapeR'),'/skyscapeR-sefstars.txt'), paste0(system.file('ephemeris',package='swephRdata'),'/sefstars.txt'), overwrite=T)
    cat('Replaced swephR stellar ephemeris file with skyscapeR version.')

  } else if (source=='swephR') {
    if (file.exists(paste0(system.file('ephemeris',package='swephR'),'/sefstars.bkp'))) {
      file.copy(paste0(system.file('ephemeris',package='swephR'),'/sefstars.bkp'), paste0(system.file('ephemeris',package='swephR'),'/sefstars.txt'), overwrite=T)
      file.copy(paste0(system.file('ephemeris',package='swephRdata'),'/sefstars.bkp'), paste0(system.file('ephemeris',package='swephRdata'),'/sefstars.txt'), overwrite=T)
    } else if (file.exists(paste0(system.file('ephemeris',package='skyscapeR'),'/sefstars.bkp'))) {
      file.copy(paste0(system.file('ephemeris',package='skyscapeR'),'/sefstars.bkp'), paste0(system.file('ephemeris',package='swephR'),'/sefstars.txt'), overwrite=T)
      file.copy(paste0(system.file('ephemeris',package='skyscapeR'),'/sefstars.bkp'), paste0(system.file('ephemeris',package='swephRdata'),'/sefstars.txt'), overwrite=T)
    } else { stop('Original swephR version of stellar ephemeris already in place.') }
    cat('Replaced stellar ephemeris file with original swephR version.')
  } else { stop('Source not recognised.')}
}

globalVariables('star.names')
