skyscapeR.env <- new.env(parent = emptyenv())

# Timezone and Calendar
skyscapeR.env$timezone <- ''
skyscapeR.env$calendar <- 'Gregorian'
skyscapeR.env$year <- 'BC/AD'

# Atmospheric Refraction
skyscapeR.env$refraction <- TRUE
skyscapeR.env$atm <- 1013.25
skyscapeR.env$temp <- 15

# Frame of Reference
skyscapeR.env$dec <- 'topo'


# if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
# if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
# if (missing(dec)) { dec <- skyscapeR.env$dec }
# if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
# if (missing(atm)) { atm <- skyscapeR.env$atm }
# if (missing(temp)) { temp <- skyscapeR.env$temp }


#' See and change the global variables used by skyscapeR
#'
#' @param timezone Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \code{\link{timezones}} for details. Default is the system timezone
#' @param calendar Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Defaults to \emph{Gregorian}.
#' @param year How to interpret years in time input. Either 'BC/AD' or '+/-', the difference
#' being that in 'BC/AD' there is no year zero, whereas in '+/-' year 0 corresponds to 1 BC.
#' Default is 'BC/AD'
#' @param refraction  Whether atmospheric refraction is to be taken into account. Default is TRUE.
#' @param atm Atmospheric pressure for refraction calculation. Default is 1013.25 mbar.
#' @param temp Atmospheric temprature for erfraction calculation. Default is 15 degrees.
#' @param dec Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. Defaults to topocentric.
#' @export
#' @examples
#' # Julian date at noon GMT on Christmas day 2018
#' time2jd('2018-12-25 12:00:00', 'GMT')
skyscapeR.vars = function(timezone, calendar, year, refraction, atm, temp, dec) {
  if (!missing(timezone)) { skyscapeR.env$timezone <- timezone }
  if (!missing(calendar)) { skyscapeR.env$calendar <- calendar }
  if (!missing(year)) { skyscapeR.env$year <- year }
  if (!missing(refraction)) { skyscapeR.env$refraction <- refraction }
  if (!missing(atm)) { skyscapeR.env$atm <- atm }
  if (!missing(temp)) { skyscapeR.env$temp <- temp }
  if (!missing(dec)) { skyscapeR.env$dec <- dec }

  return(mget(ls(skyscapeR.env), skyscapeR.env))
}
