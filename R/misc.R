#swephR::swe_set_ephe_path(system.file("ephemeris", "", package = "swephRdata"))
jd <- swephR::swe_julday(2000,1,1,12,1) # J2000.0

#' @noRd
extractLatitude <- function(loc){
  if (class(loc)=='skyscapeR.horizon') { latitude <- loc$metadata$georef[1] }
  if (class(loc)=='numeric') { latitude <- loc[1]}
  # check
  if (latitude > 90 | latitude < -90) { stop('Latitude must be between -90 and 90 degrees.')}
  return(latitude)
}


#' @noRd
stars.pval <- function(p.value) {
if (class(p.value)=='character') { p.value <- as.numeric(substr(p.value, 3,nchar(p.value))) }
out <- symnum(p.value, corr = FALSE, na = FALSE,
              cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
              symbols = c("***", "**", "*", "+", "ns"), legend=F)
return(as.character(out))
}

#' @noRd
dnorm <- function(x, mean, sd) {
  return(exp(-(x-mean)^2/(2*sd^2)) / sqrt(2*pi*sd^2)) }  ## replaces R function (twice as fast)

#' @noRd
tWS <- function(days, year) {
  WS <- findWS(year)$ind
  aux <- days-WS
  aux[aux<1] <- aux[aux<1]+365
  return(aux)
}


#' @noRd
eq2hor <- function(ra, dec, loc, refraction, atm, temp) {
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  xx <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],loc[3]), atm, temp, c(ra,dec))$xaz

  if (refraction) { alt <- xx[3] } else { alt <- xx[2] }
  az <- xx[1]-180
  if (az < 0) { az <- az + 360 }
  if (az > 360) { az <- az - 360 }
  return(list(alt = alt, az = az))
}




#' @noRd
minmaxdec = function(name, from, to, loc=FALSE) {
  xx <- seq(from, to, 100)
  dd <- c()

  # Stars
  fpath <- system.file("ephemeris", "sefstars.txt", package="swephR")
  cnames <- c('traditional name','nomenclature name','equinox','RA hr','RA min', 'RA sec', 'Dec deg', 'Dec min', 'Dec sec', 'pm RA', 'pm Dec', 'rad vel', 'ann plx', 'mag V', 'DM zone', 'DM number')
  sefstars <- read.csv(fpath, as.is=T, header=F, comment.char='#', col.names=cnames, strip.white=T)
  if (sum(as.character(sefstars$traditional.name) == name)) {
    for (i in 1: NROW(xx)) {
      dd[i] <- star(name, xx[i])$coord$Dec
    }
  } else {
    # Solar-Lunar targets
    for (i in 1: NROW(xx)) {
      dd[i] <- do.call(name, list(year=xx[i], loc=loc, verbose=F))
    }
  }
  ff <- splinefun(xx,dd)
  xxz <- seq(from,to,0.01)
  dd <- ff(xxz)
  mm <- c(min(dd),max(dd))

  return(mm)
}

#' Retrieves horizon altitude for a given azimuth from a given horizon profile
#'
#' This function retrieves the horizon altitude for a given azimuth from
#' a previously created \emph{skyscapeR.horizon} object via spline interpolation.
#' @param hor A \emph{skyscapeR.horizon} object from which to retrieve horizon altitude.
#' @param az Array of azimuth(s) for which to retrieve horizon altitude(s).
#' @export
#' @import stats
#' @seealso \code{\link{createHor}}, \code{\link{downloadHWT}}
#' @examples
#' hor <- downloadHWT('HIFVTBGK')
#' hor2alt(hor, 90)
hor2alt = function(hor, az) {
  alt <- approx(c(hor$data$az-360,hor$data$az,hor$data$az+360), rep(hor$data$alt,3), xout=az)$y
  alt <- round(alt, 2)

  if (!is.null(hor$data$alt.unc)) {
    hh <- approx(c(hor$data$az-360,hor$data$az,hor$data$az+360), rep(hor$data$alt.unc,3), xout=az)$y
    alt.unc <- round(hh, 2)

    alt <- data.frame(alt=alt, alt.unc=alt.unc)
  }
  return(alt)
}

#' Converts degree measurements in deg-min-sec (º ' ") format into decimal-point degree format.
#'
#' @param dd Degree
#' @param mm (Optional) Arcminutes
#' @param ss (Optional) Arcseconds
#' @export
#' @examples
#' deg <- ten(24, 52, 16)
ten <- function (dd, mm = 0, ss = 0) {
  np = nargs()
  testneg = NULL
  sign = 1
  if ((np == 1 && is.character(dd))) {
    temp = gsub(":", " ", dd)
    testneg = grep("-", temp)
    if (length(testneg) == 1)
      temp = sub("-", "", temp)
    if (length(testneg) == 1)
      sign = -1
    vector = as.double(strsplit(temp, " ")[[1]])
  }
  else {
    vector = as.double(c(dd, mm, ss))
  }
  fac = c(1, 60, 3600)
  vector = abs(vector)
  return(sign * sum(vector/fac))
}

#' @noRd
lty2dash <- function(lty){
  if (lty==1) { return('solid') }
  if (lty==2) { return('dash') }
  if (lty==3) { return('dot') }
  if (lty==4) { return('dashdot') }
  if (lty==5) { return('longdash') }
  if (lty==6) { return('longdashdot') }
}

#' @noRd
epoch2yr <- function(epoch){
  aux <- as.numeric(sub(" .*", "", epoch))
  sign <- sub(".* ", "", epoch)
  if (sign=='BCE') { aux <- -aux}
  return(aux)
}


#' @noRd
checkbody <- function(body){
  body <- toupper(body)
  SE <- get(data('SE', package='swephR', envir = environment())); nn <- names(SE[3:13]); nn[1] <- 'ECLIPTIC'
  ind <- which(nn==body)+2
  if (length(ind)==0) { stop('Solar System body not recognised. Please check name.')}
  body <- as.numeric(SE[ind])
  return(body)
}



#' @noRd
mag2size <- function(mag) {
  if (mag == 0 | is.na(mag) | mag > 6) { size <- 0 }
  if (mag <= 6) { size <- .1 }
  if (mag <= 5) { size <- .2 }
  if (mag <= 4) { size <- .3 }
  if (mag <= 3) { size <- .4 }
  if (mag <= 2) { size <- .5 }
  if (mag <= 1) { size <- .6 }
  if (mag < 0) { size <- .7 }
  return(size)
}


#' Calculates declination from azimuth and altitude measurements
#'
#' This function calculates the declination corresponding to an
#' orientation , i.e. an azimuth. The altitude can either be given
#'  or, alternatively, if a \emph{skyscapeR.horizon} object is provided,
#'  the corresponding horizon altitude will be automatically retrieved.
#' This function is a wrapper for function \code{\link[swephR]{swe_azalt_rev}}
#' of package \emph{swephR}.
#' @param az Azimuth(s) for which to calculate declination(s). See examples below.
#' @param loc Location, can be either a \emph{skyscapeR.horizon} object or, alternatively,
#' an array of latitude values.
#' @param alt Altitude of orientation. Optional, if left empty and a skyscapeR.object
#' is provided then this is will automatically retrieved from the horizon data via \code{\link{hor2alt}}
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_azalt_rev}}, \code{\link{hor2alt}}
#' @examples
#' hor <- downloadHWT('HIFVTBGK')
#'
#' dec <- az2dec(92, hor)
#' dec <- az2dec(92, hor, alt=4)
#'
#' # Can also be used for an array of azimuths:
#' decs <- az2dec( c(87,92,110), hor )
az2dec = function(az, loc, alt){
  if (missing(alt) & class(loc) == 'skyscapeR.horizon') { alt <- hor2alt(hor, az)[,1] }
  if (class(loc) != 'skyscapeR.horizon') {
    hor <- c()
    if (length(loc) < length(az) & length(loc)==3) { hor$metadata$georef <- matrix(rep(loc,NROW(az)), ncol=3, byrow=T)}
    if (length(loc) < length(az) & length(loc)==1) { hor$metadata$georef <- matrix(rep(c(loc,0,0),NROW(az)), ncol=3, byrow=T)}
    if (length(loc) == length(az)) { hor$metadata$georef <- cbind(loc, 0) }
    if (length(loc) == 3*NROW(az)) { hor$metadata$georef <- loc; dim(hor$metadata$georef) <- c(NROW(az),3) }
  } else { hor <- loc }

  if (length(alt) == 1) { alt <- rep(alt, NROW(az)) }

  prec <- max(nchar(sub('.*\\.', '', as.character(az))))

  dec <- c()
  for (i in 1:NROW(az)) {
    dec[i] <- round( swephR::swe_azalt_rev(jd, 1, c(hor$metadata$georef[2],hor$metadata$georef[1],hor$metadata$georef[3]), c(az[i]-180, alt[i]))$xout[2], prec)
  }

  return(dec)
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
  if (class(loc) == 'skyscapeR.horizon') { loc <- loc$metadata$georef }
  if (is.null(dim(loc))) { dim(loc) <- c(1,NROW(loc)) }
  aux <- oce::magneticField(loc[,2], loc[,1], as.POSIXlt(date, format="%Y/%m/%d"))$declination
  return(aux)
}
