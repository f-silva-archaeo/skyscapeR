

#' Converts date and time (in any timezone) to Julian date
#'
#' @param time String containing the date and time in the format "YYYY/MM/DD HH:MM:SS". BCE dates
#' should use negative sign. Use \code{\link{timestring}} if needed.
#' @param timezone (Optional) Timezone of input either as a known acronym (e.g. "GMT", "CET") or
#' a string with continent followed by country capital (e.g. "Europe/London"). See
#' \code{\link{timezones}} for details. Default is the system timezone
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Defaults to \emph{Gregorian}.
#' @param verbose (Optional) Controls whether messages should be displayed. Default is \emph{FALSE}.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_julday}}, \code{\link{as.POSIXlt}}, \code{\link{timezones}},
#' \code{\link{timestring}}
#' @examples
#' # Julian date at noon GMT on Christmas day 2018
#' time2jd('2018/12/25 12:00:00', 'GMT')
time2jd <- function(time, timezone, calendar, verbose=F) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  cal <- switch(calendar,  G = 1, Gregorian = 1, J = 0, Julian = 0, -1)

  out <- rep(NA, length(time))
  for (i in 1:length(time)) {
    UT <- c()
    ind <- which(strsplit(time[i], "")[[1]]=="/")
    UT$year <- as.numeric(substr(time[i], 1, ind[1]-1))
    UT$mon <- as.numeric(substr(time[i], ind[1]+1, ind[2]-1))
    UT$day <- as.numeric(substr(time[i], ind[2]+1, ind[2]+2))
    UT$hour <- 0; UT$min <- 0; UT$sec <- 0
    ind <- which(strsplit(time[i], "")[[1]]==":")
    if (length(ind)>0) {
      UT$hour <- as.numeric(substr(time[i], ind[1]-2, ind[1]-1))
      UT$min <- as.numeric(substr(time[i], ind[1]+1, ind[1]+2))
      if (length(ind)>1) { UT$sec <- as.numeric(substr(time[i], ind[2]+1, ind[2]+2)) } else { UT$sec <- 0 }
    }

    if (timezone != '' & timezone != 'UTC' & timezone != 'GMT') {
      if (verbose) { message('Converting to UTC.') }
      pb.date <- as.POSIXct(timestring(2000, UT$mon, UT$day, UT$hour, UT$min, UT$sec), timezone)
      aux <- format(pb.date, tz="UTC", usetz=TRUE)
      aux <- as.POSIXlt(aux, 'UTC')
      UT$hour <- aux$hour
      UT$min <- aux$min
      UT$sec <- aux$sec
    }

    out[i] <- swephR::swe_julday(UT$year, UT$mon, UT$day, UT$hour+UT$min/60+UT$sec/3600, cal)
  }

  return(out)
}



#' Converts Julian date and time (in any timezone) to julian date
#'
#' @param jd Julian date in numeric format
#' @param timezone (Optional) Desired timezone for output either as a known acronym (e.g. "GMT", "CET") or
#' a string with continent followed by country capital (e.g. "Europe/London"). See
#' \code{\link{timezones}} for details. Default is system timezone.
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. Defaults to \emph{Gregorian}.
#' @param verbose (Optional) Controls whether messages should be displayed. Default is \emph{FALSE}.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_julday}}, \code{\link{as.POSIXlt}}, \code{\link{timezones}}
#' @examples
#' jd <- time2jd('2018/12/25 12:00:00', 'GMT') # Julian date at noon GMT on Christmas day 2018
#' jd2time(jd, 'CET') # converts julian date to Central European timezone
jd2time <- function(jd, timezone, calendar, verbose=F) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  cal <- switch(calendar,  G = 1, Gregorian = 1, J = 0, Julian = 0, -1)

  out <- rep(NA, length(jd))
  for (i in 1:length(jd)) {
    time <- swephR::swe_revjul(jd[i], cal)
    time$minute <- floor((time$hour %% 1)*60)
    time$second <- floor(((((time$hour %% 1)*60)) %% 1) *60)
    time$hour <- floor(time$hour)

    if (timezone != '' & timezone != 'UTC' & timezone != 'GMT') {
      if (verbose) { message('Converting from UTC.') }
      ts <- as.POSIXct(timestring(2000, time$month, time$day, time$hour, time$minute, time$second), 'UTC')
      ts <- format(ts, tz=timezone, usetz=TRUE)
      aux <- as.POSIXlt(ts, 'UTC')
      time$hour <- aux$hour
      time$minute <- aux$min
      time$second <- aux$sec
    }
    out[i] <- timestring(time$year, time$month, time$day, time$hour, time$minute, time$second)
  }
  return(out)
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

  ts <- paste(paste(year, month, day, sep='/'), paste(hour, minute, second, sep=':'))
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
findWS <- function(year, calendar) {
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  jd0 <- swephR::swe_julday(year, 1, 1, 12, 1)
  jd <- seq(jd0, jd0+365, length.out = 365*24)
  dec <- swephR::swe_calc_ut(jd, 0, 2048)$xx[,2]
  ind <- which.min(dec)

  aux <- jd2time(jd[ind], calendar=calendar)
  ind2 <- which(strsplit(aux, "")[[1]]==" ")

  return(list(ind=round(ind/24,0)+1, date=substr(aux,1,ind2-1), jd=jd[ind]))
}



#' @noRd
calWS <- function(WS) {
  days <- seq(1,365,1)
  days <- 365-WS$ind+days
  days[days>365] <- days[days>365]-365
  return(days)
}



#' @noRd
dd.to.DD_unvec <- function(day, char=F, WS=F) {
  if (WS==T) { day <- day-11 }
  if (day<=0) { day <- day+ 365 }
  if (char==F) {
    out <- matrix(NA, nrow=length(day), ncol=2)
  } else { out <- c() }

  for (i in 1:length(day)) {
    if (day[i] >= 1 & day[i] <= 365) {
      month.name
      months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
      summon <- cumsum(months)

      mm <- tail(which(summon < day[i]),1) + 1
      dd <- day[i] - summon[tail(which(summon <= day[i]),1)]
      if (NROW(mm) < 1) {
        mm <- 1
        dd <- day[i]
      }
      if (dd == 0) { dd <- months[12] }

      if (char) {
        out[i] <- paste0(month.name[mm], ' ',dd)
      } else { out[i,] <- c(mm,dd) }
    }
  }
  return(out)
}
#' @noRd
dd.to.DD <- Vectorize(dd.to.DD_unvec, 'day')


#' Converts year number (epoch) to calendar year
#'
#' @param year year number
#' @export
#' @examples
#' BC.AD(100)
#' BC.AD(0)
#' BC.AD(-1)
#' BC.AD(-99)
#' BC.AD(-100)
BC.AD <- function(year) {
  out <- year
  ind <- which(year<=0); out[ind] <- paste(abs(year[ind]-1), 'BC')
  ind <- which(year>0); out[ind] <- paste(abs(year[ind]), 'AD')
  return(out)
}


#' @noRd
day <- function(date, as.char=F) {
  ind1 <- which(strsplit(date, "")[[1]]=="/")
  ind2 <- which(strsplit(date, "")[[1]]==" "); if (length(ind2)==0) { ind2 <- length(date)+1 }
  aux <- substr(date,ind1[2]+1,ind2-1)
  if (!as.char) { aux <- as.numeric(aux) }
  return(aux)
}


#' @noRd
month <- function(date, as.char=F) {
  ind <- which(strsplit(date, "")[[1]]=="/")
  aux <- substr(date,ind[1]+1,ind[2]-1)
  if (!as.char) { aux <- as.numeric(aux) }
  return(aux)
}


#' @noRd
year <- function(date, as.char=F) {
  ind <- which(strsplit(date, "")[[1]]=="/")
  aux <- substr(date,1,ind[1]-1)
  if (!as.char) { aux <- as.numeric(aux) }
  return(aux)
}

