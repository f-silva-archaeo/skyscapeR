#' Data reduction for theodolite measurements using the sun-sight method
#'
#' This function calculates the true azimuth of a structure measured with
#' a theodolite using the sunsight technique.
#' @param loc Location, can be either a \emph{skyscapeR.horizon} object or, alternatively,
#' a latitude.
#' @param az Array of azimuths. Use \code{\link[astrolibR]{ten}} to convert to decimal point format if necessary.
#' @param date Date of measurements as a string in the format: 'YYYY-MM-DD'
#' @param time Time of sun-sight measurement in the format: 'HH:MM:SS'
#' @param tz Timezone of input wither as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London").
#' @param az.sun (Optional) Measured azimuth of the sun. Defaults to zero.
#' @export
#' @seealso \code{\link{sunAz}}, \code{\link[astrolibR]{ten}}, \code{\link{sixty}}
#' @references Ruggles, C.L.N. (1999). \emph{Astronomy in Prehistoric Britain and Ireland}. Yale University Press.
#' @examples
#' lat <- ten(35,50,37.8)
#' lon <- ten(14,34,6.4)
#' az <- c( ten(298,24,10), ten(302,20,40))
#' az.sun <- ten(327,29,50)
#' date <- "2016/02/20"
#' time <- "11:07:17"
#'
#' true.az <- reduction.theod(c(lat,lon), az, date , time, tz= "Europe/Malta", az.sun)
reduction.theod = function(loc, az, date, time, tz, az.sun = 0) {
  if (class(loc)=='skyscapeR.horizon') { loc <- loc$georef }

  date <- as.Date(date, "%Y/%m/%d")
  time <- paste(date, time)

  diff <- az - rep(az.sun, NROW(az))
  ind <- which(abs(diff)>180); if (length(ind)>0) { diff[ind] <- az[ind] - rep(az.sun-360, NROW(ind)) }

  az.sun.corr <- sunAz(loc, time, tz)

  return(az.sun.corr + diff)
}


