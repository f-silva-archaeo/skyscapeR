################################################
skyscapeR.env <- new.env(parent = emptyenv())

# Timezone and Calendar
skyscapeR.env$timezone <- ''
skyscapeR.env$calendar <- 'Gregorian'

# Atmospheric Refraction
skyscapeR.env$refraction <- TRUE
skyscapeR.env$atm <- 1013.25
skyscapeR.env$temp <- 15

# Frame of Reference
skyscapeR.env$dec <- 'topo'

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
#' @examples
#' # Julian date at noon GMT on Christmas day 2018
#' time2jd('2018-12-25 12:00:00', 'GMT')
skyscapeR.vars = function(timezone, calendar, refraction, atm, temp, dec) {
  if (!missing(timezone)) { skyscapeR.env$timezone <- timezone }
  if (!missing(calendar)) { skyscapeR.env$calendar <- calendar }
  if (!missing(refraction)) { skyscapeR.env$refraction <- refraction }
  if (!missing(atm)) { skyscapeR.env$atm <- atm }
  if (!missing(temp)) { skyscapeR.env$temp <- temp }
  if (!missing(dec)) { skyscapeR.env$dec <- dec }

  aux <- ls(skyscapeR.env)
  return(mget(aux[-which(aux=='sefstars')], skyscapeR.env))
}


################################################
### Clean-up and prep sefstars file from swephR
fpath <- system.file("ephemeris", "sefstars.txt", package="swephR")
cnames <- c('traditional name','nomenclature name','equinox','RA hr','RA min', 'RA sec', 'Dec deg', 'Dec min', 'Dec sec', 'pm RA', 'pm Dec', 'rad vel', 'ann plx', 'mag V', 'DM zone', 'DM number')
sefstars <- read.csv(fpath, as.is=T, header=F, comment.char='#', col.names=cnames, strip.white=T)
sefstars <- sefstars[-which(sefstars$mag.V==0),]
ind <- sort(sefstars$mag.V, index.return=T)$ix
sefstars <- sefstars[ind,]

# removes duplicated entries (double stars, etc)
ind <- which(duplicated(sefstars$nomenclature.name)); sefstars <- sefstars[-ind,]
ind <- which(sefstars$nomenclature.name=='80Uma'); sefstars <- sefstars[-ind,]

# removes stars without name as swephR::swe_fixstar2_ut cannot pick them up anyway
# ind <- which(ss$traditional.name==""); ss <- ss[-ind,]

skyscapeR.env$sefstars <- sefstars
rm(fpath, cnames, ind, sefstars)
################################################

