swephR::swe_set_ephe_path(system.file("ephemeris", "", package = "swephRdata"))

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
eq.to.Hor <- function(xx,lat,lon) {
  rr <- eq2horFS(xx[1],xx[2],xx[3],cbind(lat,lon))$alt
  return(rr)
}

#' @noRd
findWS <- function(year) {
  jd0 <- swephR::swe_julday(year, 1, 1, 12, 12)
  yy <- seq(jd0, length.out = 365)
  ss <- swephR::swe_calc_ut(yy, 0, 2048)$xx[,2]
  ind <- which.min(ss)

  return(list(ind=ind, month=dd.to.DD(ind)[1], day=dd.to.DD(ind)[2]))
}

#' @noRd
calWS <- function(WS) {
  days <- seq(1,365,1)
  days <- 365-WS$ind+days
  days[days>365] <- days[days>365]-365
  return(days)
}

#' @noRd
dd.to.DD <- function(day, char=F) {
  if (day>=1 & day<=365) {
    month.name
    months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    summon <- cumsum(months)

    mm <- tail(which(summon < day),1) + 1
    dd <- day - summon[tail(which(summon <= day),1)]
    if (NROW(mm) < 1) {
      mm <- 1
      dd <- day
    }
    if (dd == 0) { dd <- months[12] }

    if (char) {
      return( paste0(month.name[mm], ' ',dd) )
    } else { return(c(mm,dd)) }
  }
}

#' Fixed eq2hor function
#' @noRd
eq2horFS <- function(ra, dec, jd=jd, loc, refraction=F, atm=1013.25, temp=15) {
  xx <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, c(ra,dec))$xaz

  if (refraction) { alt <- xx[3] } else { alt <- xx[2] }
  az <- xx[1]-180
  if (az < 0) { az <- az + 360 }
  if (az > 360) { az <- az - 360 }
  return(list(alt = alt, az = az))
}


#' @noRd
minmaxdec = function(name, from, to) {
  xx <- seq(from, to, 100)
  dd <- c()

  # Stars
  data(stars, envir=environment())
  if (sum(as.character(stars$NAME) == name)) {
    for (i in 1: NROW(xx)) {
      dd[i] <- star(name, xx[i])$dec
    }
  } else {
    # Solar-Lunar targets
    for (i in 1: NROW(xx)) {
      dd[i] <- do.call(name, list(xx[i]))
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
yr2epoch <- function(year){
  aux <- as.character(abs(year))
  if (year > 0) { aux <- paste(aux, 'CE') }
  if (year < 0)  { aux <- paste(aux, 'BCE') }
  if (year == 0) { aux <- '1 CE' }

  return(aux)
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
  data(SE, package='swephR', envir = environment()); nn <- names(SE[3:13]); nn[1] <- 'ECLIPTIC'
  ind <- which(nn==body)+2
  if (length(ind)==0) { stop('Solar System body name not recognised. Please recheck.')}
  body <- as.numeric(SE[ind])
  return(body)
}


#' Converts date and time (in any timezone) to Julian date
#'
#' @param time String containing the date and time in the format "YYYY-MM-DD HH:MM:SS".
#' Use \code{\link{timestring}} if needed.
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \code{\link{timezones}} for details. Default is GMT.
#' @param cal (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Defaults to Gregorian.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_julday}}, \code{\link{as.POSIXlt}}, \code{\link{timezones}},
#' \code{\link{timestring}}
#' @examples
#' # Julian date at noon GMT on Christmas day 2018
#' time2jd('2018-12-25 12:00:00', 'GMT')
time2jd <- function(time, timezone = 'GMT', cal='G') {
  pb.date <- as.POSIXct(time, timezone)
  UT <- format(pb.date, tz="UTC",usetz=TRUE)
  UT <- as.POSIXlt(UT, 'UTC')

  if (cal=='G') { cal <- 1 }
  if (cal=='J') { cal <- 0 }

  jd <- swephR::swe_julday(UT$year+1900, UT$mon+1, UT$mday, UT$hour+UT$min/60+UT$sec/3600, cal)
  return(jd)
}

#' Converts Julian date and time (in any timezone) to julian date
#'
#' @param jd Julian date in numeric format
#' @param timezone (Optional) Desired timezone for output either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \code{\link{timezones}} for details. Default is GMT.
#' @param cal (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. Defaults to Gregorian.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_julday}}, \code{\link{as.POSIXlt}}, \code{\link{timezones}}
#' @examples
#' jd <- time2jd('2018-12-25 12:00:00', 'GMT') # Julian date at noon GMT on Christmas day 2018
#' jd2time(jd, 'CET') # converts julian date to Central European timezone
#'
jd2time <- function(jd, timezone='GMT', cal='G') {
  if (cal=='G') { cal <- 1 }
  if (cal=='J') { cal <- 0 }

  time <- swe_revjul(jd, cal)

  ts <- timestring(time$year, time$month, time$day, floor(time$hour), floor((time$hour %% 1)*60), floor(((((time$hour %% 1)*60)) %% 1) *60))
  ts <- as.POSIXct(ts, 'UTC')
  ts <- format(ts, tz=timezone, usetz=TRUE)

  return(substr(ts,1, 19))
}

#' Converts date and time numeric values to a single string
#'
#' @param year Year
#' @param month Month
#' @param day Day
#' @param hour Hour
#' @param minute Minute
#' @param second Second
#' @export
#' @examples
#' timestring(2018, 12, 25, 2, 34)
timestring <- function(year, month, day, hour=12, minute=0, second=0) {
  month <- as.character(month); if (nchar(month) < 2) { month <- paste0('0', month) }
  day <- as.character(day); if (nchar(day) < 2) { day <- paste0('0', day) }
  hour <- as.character(hour); if (nchar(hour) < 2) { hour <- paste0('0', hour) }
  minute <- as.character(minute); if (nchar(minute) < 2) { minute <- paste0('0', minute) }
  second <- as.character(second); if (nchar(second) < 2) { second <- paste0('0', second) }

  ts <- paste(paste(year, month, day, sep='-'), paste(hour, minute, second, sep=':'))
  return(ts)
}


#' Converts day-month in 'MM-DD' format to a more readable format
#'
#' @param date date in 'MM-DD' format
#' @export
#' @examples
#' long.date('01-01')
#' long.date('08-23')
long.date <- function(date){
  day <- as.numeric(substr(date,4,5))
  month <- month.abb[as.numeric(substr(date,1,2))]

  return( paste(day, month) )
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

