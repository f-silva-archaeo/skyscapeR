#' Data reduction for theodolite measurements using the sun-sight method
#'
#' This function calculates the true azimuth of a structure measured with
#' a theodolite using the sunsight technique.
#' @param loc Location, can be either a \emph{skyscapeR.horizon} object or,
#' alternatively, a latitude.
#' @param az Array of azimuths. Use \code{\link[astrolibR]{ten}} to convert to
#' decimal point format if necessary.
#' @param date Date of measurements as a string in the format: 'YYYY/MM/DD'
#' @param time Time of sun-sight measurement in the format: 'HH:MM:SS'
#' @param tz Timezone of input wither as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London").
#' @param az.sun (Optional) Measured azimuth of the sun. Defaults to zero.
#' @param alt (Optional) Altitude, necessary for automatic declination calculation.
#' If missing and \emph{loc} is a \emph{skuscapeR.horizon} object then the altitude
#' will be automatically read from the horizon profile.
#' @export
#' @seealso \code{\link{sunAz}}, \code{\link[astrolibR]{ten}}, \code{\link{sixty}}
#' @references Ruggles, C.L.N. (1999). \emph{Astronomy in Prehistoric Britain and Ireland}.
#' Yale University Press.
#' @examples
#' lat <- ten(35,50,37.8)
#' lon <- ten(14,34,6.4)
#' az <- c( ten(298,24,10), ten(302,20,40))
#' az.sun <- ten(327,29,50)
#' date <- "2016/02/20"
#' time <- "11:07:17"
#'
#' data <- reduct.theodolite(c(lat,lon), az, date , time, tz= "Europe/Malta", az.sun)
#'
#' # Declination will be automatically calculated if the altitude is also given:
#' data <- reduct.theodolite(c(lat,lon), az, date , time, tz= "Europe/Malta", az.sun, alt=c(2,5))
#'
#' # Alternatively, the altitude can be automatically retrieved from a horizon profile:
#' hor <- download.HWT('UFXERSLQ')
#' data <- reduct.theodolite(hor, az, date, time, tz= "Europe/Malta", az.sun)
reduct.theodolite = function(loc, az, date, time, tz, az.sun = 0, alt) {
  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- loc$georef } else { hor <- NULL }

  date <- as.Date(date, "%Y/%m/%d")
  time <- paste(date, time)

  diff <- az - rep(az.sun, NROW(az))
  ind <- which(abs(diff)>180); if (length(ind)>0) { diff[ind] <- az[ind] - rep(az.sun-360, NROW(ind)) }

  az.sun.corr <- sunAz(loc, time, tz)
  az.corr <- az.sun.corr + diff

  df <- data.frame(Latitude=loc[1], Longitude=loc[2], Uncorrected.Az=az, Date.Time=time, Sun.Az=az.sun.corr, Corrected.Az=az.corr)

  if (!missing(alt)) {
    dec <- az2dec(az.corr, loc, alt)
    df$Altitude = alt
    df$Declination <- dec
  } else if (class(hor)=='skyscapeR.horizon') {
    dec <- az2dec(az.corr, hor)
    df$Altitude <- hor2alt(hor, az.corr)
    df$Declination <- dec
  }

  return(df)
}


#' Data reduction for compass measurements
#'
#' This function calculates the true azimuth of a structure measured with
#' a compass.
#' @param loc Location, can be either a \emph{skyscapeR.horizon} object or,
#' alternatively, a latitude.
#' @param mag.az Array of magnetic azimuth measurements.
#' @param date (Optional) Date of measurements as a string in the format: 'YYYY/MM/DD'.
#' Only necessary is \emph{magdec} is not given.
#' @param magdec (Optional) Magnetic declination, if known.
#' @param alt (Optional) Altitude, necessary for automatic declination calculation.
#' If missing and \emph{loc} is a \emph{skuscapeR.horizon} object then the altitude
#' will be automatically read from the horizon profile.
#' @export
#' @seealso \code{\link{mag.dec}}, \code{\link{az2dec}}, \code{\link{hor2alt}}
#' @examples
#' loc <- c(35,-7)
#' mag.az <- c(89.5, 105, 109.5)
#' data <- reduct.compass(loc, mag.az, "2016/04/02")
#'
#' # Declination will be automatically calculated if the altitude is also given:
#' data <- reduct.compass(loc, mag.az, "2016/04/02", alt=c(1,2,0))
#'
#' # Alternatively, the altitude can be automatically retrieved from a horizon profile:
#' hor <- download.HWT('NML6GMSX')
#' data <- reduct.compass(hor, mag.az, "2016/04/02")
reduct.compass = function(loc, mag.az, date, magdec, alt) {
  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- loc$georef } else { hor <- NULL }

  if (missing(magdec) & !missing(date)) {
    magdec <- mag.dec(loc, date)
  }
  true.az <- mag.az + magdec

  df <- data.frame(Latitude=loc[1], Longitude=loc[2], Magnetic.Az=mag.az, Date=date, Mag.Dec=magdec, True.Az=true.az)

  if (!missing(alt)) {
    dec <- az2dec(true.az, loc, alt)
    df$Altitude = alt
    df$Declination <- dec
  } else if (class(hor)=='skyscapeR.horizon') {
    dec <- az2dec(true.az, hor)
    df$Altitude <- hor2alt(hor, true.az)
    df$Declination <- dec
  }

  return(df)
}
