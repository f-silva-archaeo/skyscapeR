#' Data reduction for theodolite measurements using the sun-sight method
#'
#' This function calculates the true azimuth of a structure measured with
#' a theodolite using the sunsight technique.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param az Array of azimuths. Use \code{\link[astrolibR]{ten}} to convert to
#' decimal point format if necessary.
#' @param date Date of measurements as a string in the format: 'YYYY/MM/DD'
#' @param time Time of sun-sight measurement in the format: 'HH:MM:SS'
#' @param tz Timezone of input wither as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London").
#' @param az.sun (Optional) Measured azimuth of the sun. Defaults to zero.
#' @param alt (Optional) Altitude, necessary for automatic declination calculation.
#' If missing and \emph{loc} is a \emph{skyscapeR.horizon} object then the altitude
#' will be automatically read from the horizon profile.
#' @param name (Optional) Names or labels to identify each measurement.
#' @param ID (Optional) IDs or codes to identify each measurement.
#' @param HWT.ID (Optional) HeyWhatsThat IDs relating to a previously generated
#' horizon profile for measurement.
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
reduct.theodolite = function(loc, az, date, time, tz, az.sun = 0, alt, name, ID, HWT.ID) {
  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- loc$georef } else { hor <- NULL }

  if (NROW(loc) < NROW(az)*2) { loc <- matrix(loc,NROW(az),2, byrow=T) }
  if (NROW(date) < NROW(az)) { date <- matrix(date,NROW(az),1, byrow=T) }
  if (NROW(time) < NROW(az)) { time <- matrix(time,NROW(az),1, byrow=T) }
  if (NROW(tz) < NROW(az)) { tz <- matrix(tz,NROW(az),1, byrow=T) }
  if (NROW(az.sun) < NROW(az)) { az.sun <- matrix(az.sun,NROW(az),1, byrow=T) }
  if (!missing(alt)) { if (NROW(alt) < NROW(az)) { alt <- matrix(alt,NROW(az),1, byrow=T) } }

  date <- as.Date(date, "%Y/%m/%d")
  time <- paste(date, time)

  diff <- az - az.sun
  ind <- which(abs(diff)>180); if (length(ind)>0) { diff[ind] <- az[ind] - az.sun[ind]-360 }

  prec <- max(nchar(sub('.*\\.', '', as.character(az))))
  az.sun.corr <- round(sunAz(loc, time, tz), prec)
  az.corr <- az.sun.corr + diff

  df <- c()
  if (!missing(ID)) { df$ID <- ID }
  if (!missing(name)) { df$Name <- name }
  df$Latitude=loc[,1]
  df$Longitude=loc[,2]
  if (!missing(HWT.ID)) { df$HWT.ID <- HWT.ID }
  df$Uncorrected.Azimuth=az
  df$Date.Time=time
  df$Sun.Az=az.sun.corr
  df$True.Azimuth=az.corr

  if (!missing(alt)) {
    message('Altitude values found. Calculating declination...')
    dec <- az2dec(az.corr, loc, alt)
    df$Altitude = alt
    df$Declination <- dec
  } else if (class(hor)=='skyscapeR.horizon') {
    message('Horizon profile found. Obtaining altitude values and calculating declination...')
    dec <- az2dec(az.corr, hor)
    df$Altitude <- hor2alt(hor, az.corr)
    df$Declination <- dec
  } else { message('No altitude values or horizon profile found. Declination values were not calculated.') }

  df <- as.data.frame(df)
  return(df)
}


#' Data reduction for compass measurements
#'
#' This function calculates the true azimuth of a structure measured with
#' a compass.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param mag.az Array of magnetic azimuth measurements.
#' @param date (Optional) Date of measurements as a string in the format: 'YYYY/MM/DD'.
#' Only necessary is \emph{magdec} is not given.
#' @param magdec (Optional) Magnetic declination, if known.
#' @param alt (Optional) Altitude, necessary for automatic declination calculation.
#' If missing and \emph{loc} is a \emph{skyscapeR.horizon} object then the altitude
#' will be automatically read from the horizon profile.
#' @param name (Optional) Names or labels to identify each measurement.
#' @param ID (Optional) IDs or codes to identify each measurement.
#' @param HWT.ID (Optional) HeyWhatsThat IDs relating to a previously generated
#' horizon profile for measurement.
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
reduct.compass = function(loc, mag.az, date, magdec, alt, name, ID, HWT.ID) {
  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- loc$georef } else { hor <- NULL }

  if (NROW(loc) < NROW(mag.az)*2) { loc <- matrix(loc,NROW(mag.az),2, byrow=T) }
  if (!missing(date) & (NROW(date) < NROW(mag.az))) { date <- matrix(date,NROW(mag.az),1, byrow=T) }
  if (!missing(alt)) { if (NROW(alt) < NROW(mag.az)) { alt <- matrix(alt,NROW(mag.az),1, byrow=T) } }

  if (missing(magdec) & !missing(date)) {
    magdec <- mag.dec(loc, date)
  }
  prec <- max(nchar(sub('.*\\.', '', as.character(mag.az))))
  true.az <- round(mag.az + magdec, prec)
  df <- c()
  if (!missing(ID)) { df$ID <- ID }
  if (!missing(name)) { df$Name <- name }
  df$Latitude=loc[,1]
  df$Longitude=loc[,2]

  df$Magnetic.Azimuth=mag.az
  df$Date=date
  df$Mag.Dec=magdec
  df$True.Azimuth=true.az
  if (!missing(HWT.ID)) { df$HWT.ID <- HWT.ID }
  if (!missing(alt)) {
    message('Altitude values found. Calculating declination...')
    dec <- az2dec(true.az, loc, alt)
    df$Altitude = alt
    df$Declination <- dec
  } else if (class(hor)=='skyscapeR.horizon') {
    message('Horizon profile found. Obtaining altitude values and calculating declination...')
    dec <- az2dec(true.az, hor)
    df$Altitude <- hor2alt(hor, true.az)
    df$Declination <- dec
  } else { message('No altitude values or horizon profile found. Declination values were not calculated.') }

  df <- as.data.frame(df)
  return(df)
}
