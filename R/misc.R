#swephR::swe_set_ephe_path(system.file("ephemeris", "", package = "swephRdata"))
jd <- swephR::swe_julday(2000,1,1,12,1) # J2000.0

#' @noRd
checkYear <- function(year){
  if (sum(year < -13201 | year > 17191)>0) { stop('Year(s) outside of Swiss Ephemeris range (13201 BCE to 17191 CE).')}
}

#' @noRd
extractLatitude <- function(loc){
  if (class(loc)[1]=='skyscapeR.horizon') { latitude <- loc$metadata$georef[1] }
  if (class(loc)[1]=='numeric') { latitude <- loc[1]}
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
dnorm <- function(x, mean, sd) { return(exp(-(x-mean)^2/(2*sd^2)) / sqrt(2*pi*sd^2)) }  ## replaces R function (twice as fast)

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

  jd <- swephR::swe_julday(2000,1,1,12,1)
  xx <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],loc[3]), atm, temp, c(ra,dec))$xaz

  if (refraction) { alt <- xx[3] } else { alt <- xx[2] }
  az <- xx[1]-180
  if (az < 0) { az <- az + 360 }
  if (az > 360) { az <- az - 360 }

  return(list(alt = alt, az = az))
}

#' @noRd
veq2hor <- Vectorize(eq2hor, 'ra', SIMPLIFY=T)

#' @noRd
minmaxdec = function(name, from, to, loc=FALSE) {
  xx <- seq(from, to, 100)
  dd <- c()

  # Stars
  sefstars <- skyscapeR.env$sefstars
  if (sum(as.character(sefstars$traditional.name) == name) | sum(as.character(sefstars$nomenclature.name) == name)) {
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



#' Estimates magnetic declination (difference between true and magnetic
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
#' loc <- c( 51.5074, -0.1278 )
#' mag.dec( loc, "2016/04/01" )
mag.dec = function(loc, date) {
  if (class(loc)[1] == 'skyscapeR.horizon') { loc <- loc$metadata$georef }
  if (length(loc) < length(date)) { loc <- matrix(loc,NROW(date),3, byrow=T) }
  if (NROW(date) < NROW(loc)) { date <- matrix(date,NROW(loc),1, byrow=T) }

  if (is.null(dim(loc))) { dim(loc) <- c(1,NROW(loc)) }
  aux <- rep(NA, NROW(loc))
  for (i in 1:NROW(loc)) {
    aux[i] <- oce::magneticField(loc[i,2], loc[i,1], as.POSIXlt(date[i], format="%Y/%m/%d"))$declination
  }
  return(aux)
}

#' @noRd
subset <- function(x, subset, ...){
  if (class(x) == 'skyscapeR.pdf') {
    dd <- x
    dd$metadata$name <- dd$metadata$name[subset]
    dd$metadata$az <- dd$metadata$az[subset]
    dd$metadata$unc <- dd$metadata$unc[subset]
    dd$metadata$horizon <- dd$metadata$horizon[subset]
    dd$data <- dd$data[subset]
    return(dd)
  } else { return(subset(x, subset, ...)) }
}


#' @noRd
sampleList <- function(ll) {
  n <- NROW(ll)
  out <- c()
  for (i in 1:n){
    out[i] <- sample(ll[[i]][,1], 1, prob=ll[[i]][,2])
  }
  return(out)
}

#' Returns the high-density region of a probability distribution
#'
#' @param x A \emph{skyscapeR.spd} or \emph{skyscapeR.pdf} object.
#' @param mass (Optional) Probability mass of the region. Default is 0.954.
#' @export
hpdi <- function(x, mass=0.954){
  if (class(x) == 'skyscapeR.spd') { grd <- x$data }
  if (class(x) == 'skyscapeR.pdf' | class(x) == 'list' | class(x) == 'data.frame') { grd <- x }
  sorted <- sort(grd$y, decreasing=TRUE)
  heightIdx = min( which( cumsum( sorted) >= sum(grd$y, na.rm=T) * mass ) )
  height = sorted[heightIdx]
  indices = which( grd$y >= height )
  gaps <- which(diff(indices) > 1)
  starts <- indices[c(1, gaps + 1)]
  ends <- indices[c(gaps, length(indices))]
  result <- cbind(start = grd$x[ends], end = grd$x[starts])
  if (result[2] < result[1]) { aux <- result[1]; result[1] <- result[2]; result[2] <- aux }
  return(result)
}


#' Converts p-value into symbol
#'
#' @param p.value p-value
#' @export
pval2stars <- function(p.value) {
  if (class(p.value)=='character') { p.value <- as.numeric(substr(p.value, 3,nchar(p.value))) }
  out <- symnum(p.value, corr = FALSE, na = FALSE,
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols = c("***", "**", "*", "+", "ns"), legend=F)
  return(as.character(out))
}
