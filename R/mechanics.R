cur.year <- as.numeric(format(Sys.Date(), "%Y")) # current year

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
    if (length(loc) < length(az) & length(loc)==2) { hor$metadata$georef <- matrix(rep(loc,NROW(az)), ncol=3, byrow=T)}
    if (length(loc) < length(az) & length(loc)==1) { hor$metadata$georef <- matrix(rep(c(loc,0),NROW(az)), ncol=3, byrow=T)}
    if (length(loc) == length(az)) { hor$metadata$georef <- cbind(loc, 0) }
    if (length(loc) == 3*NROW(az)) { hor$metadata$georef <- loc; dim(hor$metadata$georef) <- c(NROW(az),3) }
  } else { hor <- loc }

  if (length(alt) == 1 & length(az) > 1) { alt <- rep(alt, NROW(az)) }

  prec <- max(nchar(sub('.*\\.', '', as.character(az))))

  dec <- c()
  for (i in 1:NROW(az)) {
    dec[i] <- round( swephR::swe_azalt_rev(jd, 1, c(hor$metadata$georef[2],hor$metadata$georef[1],0), c(az[i]-180, alt[i]))$xout[2], prec)
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



#' Computes obliquity of the ecliptic
#'
#' This function calculates the obliquity for a given epoch. It is a
#' wrapper for function \code{\link[swephR]{swe_calc_ut}} of package \emph{swephR}.
#' @param year Year for which to calculate the obliquity.
#' Defaults to present year as given by Sys.Date()
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc_ut}}
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @examples
#' #' # Obliquity for year 3999 BC:
#' obliquity(-4000)
obliquity = function(year = cur.year) {
  aux <- c()
  for (i in 1:length(year)) {
    jd <- swephR::swe_julday(year[i],1,1,12,1)
    aux[i] <- swephR::swe_calc_ut(jd, -1, 0)$xx[1]
  }
  return(aux)
}

#' @noRd
body.position.unvec = function(obj='sun', time, timezone, calendar, dec, loc = NULL, refraction, atm, temp, verbose=T) {

  if (is.null(loc)) {
    if (verbose) { cat('No location given. Forcing output to equatorial coordinates in the geocentric frame of reference.\n') }
    dec <- 'geo'
    nohoriz <- T
  } else { nohoriz <- F }

  body <- checkbody(obj)

  coords.eq <- data.frame(RA=NA, Dec=NA)
  if (!nohoriz) { coords.hor <- data.frame(az=NA, alt=NA) }

  if (class(time)=='character') { jd <- time2jd(time, timezone, calendar) } else { jd <- time }
  if (dec == 'geo') {
    coords.eq[1,] <- swephR::swe_calc_ut(jd, body, 2048)$xx[1:2]

  } else if (dec == 'topo') {
    if (class(loc)=='skyscapeR.horizon') { loc <- loc$georef }
    swephR::swe_set_topo(loc[2],loc[1],loc[3])
    coords.eq[1,] <- swephR::swe_calc_ut(jd, body, 2048+32*1024+16)$xx[1:2]
  }

  if (!nohoriz) {
    aux <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],loc[3]), atm, temp, as.numeric(coords.eq))$xaz

    if (refraction) {
      coords.hor[1,] <- aux[c(1,3)] ## apparent altitude
      coords.eq[1,] <- swephR::swe_azalt_rev(jd, 1, c(loc[2],loc[1],loc[3]), aux)$xout[1:2] ## apparent declination
    } else { coords.hor <- aux[1:2] }
    coords.hor[1,1] <- coords.hor[1,1] - 180
    if (coords.hor[1,1] > 360) { coords.hor[1,1] <- coords.hor[1,1] - 360 }
    if (coords.hor[1,1] < 0) { coords.hor[1,1] <- coords.hor[1,1] + 360 }
  }

  out <- c()
  out$equatorial <- coords.eq
  if (!nohoriz) { out$horizontal <- coords.hor }
  return(out)
}

#' @noRd
body.position.vec <- Vectorize(body.position.unvec, 'time', SIMPLIFY = FALSE)


#' Computes position of Solar System bodies in equatorial coordinates
#'
#' This function calculates the geocentric or topocentric declination and right ascension of solar
#'  system bodies at a given time. It is a wrapper for function
#'  \code{\link[swephR]{swe_calc_ut}} of package \emph{swephR}.
#' @param obj (Optional) String containing name of the solar system body of interest. Can be
#' any of the planets (inc. Pluto), the Moon, the Sun or the Ecliptic. Defaults to 'sun'.
#' @param time Either a string containing the date and time in the format "YYYY/MM/DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. #' If not given the value set
#' by \code{\link{skyscapeR.vars}} will be used instead.
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param dec (Optional) Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param loc (Optional) Location, only needed if output is in topocentric declination.
#' @param refraction (Optional) Whether atmospheric refraction is to be taken into account.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temprature for erfraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link{timestring}}, \code{\link{time2jd}}
#' @examples
#' # Position of the sun at noon GMT on Christmas day 2018:
#' body.position('sun', '2018/12/25 12:00:00', timezone='GMT')
#'
#' # Declination of the moon at same time
#' body.position('moon', '2018/12/25 12:00:00', timezone='GMT')$equatorial$Dec
body.position <- function(obj='sun', time, timezone, calendar, dec, loc=NULL, refraction, atm, temp, verbose=T) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(dec)) { dec <- skyscapeR.env$dec }
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (length(time)==1) {
    return(body.position.unvec(obj, time, timezone, calendar, dec, loc, refraction, atm, temp, verbose))
  } else {
    aux <- body.position.vec(obj, time, timezone, calendar, dec, loc, refraction, atm, temp, verbose)
    out <- c()
    out$equatorial <- as.data.frame(matrix(unlist(lapply(aux,'[[', 'equatorial')), ncol=2, byrow=T)); names(out$equatorial) <- names(aux[[1]]$equatorial)
    if (!is.null(loc)) { out$horizontal <- as.data.frame(matrix(unlist(lapply(aux,'[[', 'horizontal')), ncol=2, byrow=T)); names(out$horizontal) <- names(aux[[1]]$horizontal) }
    return(out)
  }
}



#' Computes the phase of the moon
#'
#' This function calculates the moon phase, in percentage of full. It is a wrapper
#' for function \code{\link[swephR]{swe_pheno_ut}} of package \emph{swephR}.
#' @param time Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. Defaults to system timezone.
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. Defaults to Gregorian.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_pheno_ut}}
#' @examples
#' # Moonphase at noon GMT on Christmas day 2018:
#' moonphase('2018/12/25 12:00:00', 'GMT')
moonphase <- function(time, timezone, calendar) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  aux <- c()
  for (i in 1:length(time)) {
    if (class(time)=='character') { jd <- time2jd(time[i], timezone, calendar) } else { jd <- time[i] }
    aux[i] <- swephR::swe_pheno_ut(jd, 1, 0)$attr[2]
  }

  return(aux)
}


#' Computes the rising and setting azimuth, declination and time of a Solar System object
#' for a given location and day
#'
#' @param obj (Optional) String containing name of the solar system body of interest. Can be
#' any of the planets (inc. Pluto), the Moon, the Sun or the Ecliptic. Defaults to 'sun'.
#' @param date Either a string containing the date in the format "YYYY/MM/DD"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#'  Can also be a single year ("YYYY") or a month ("YYYY/MM") to calculate risings and settings for every day in the
#'  year or month, respectively. Not necessary if \emph{jd} is given.
#' @param jd (Optional) A numeric containing the julian date (see \code{\link{time2jd}}) for which to
#' calculate rising and settings. Only needed if \emph{date} is not given.
#' @param calendar (Optional) Calendar used in parameter \emph{date}. G for gregorian and J for julian.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param alt (Optional) The altitude of the horizon to consider. Defaults to zero degrees.
#' @param timezone (Optional) Timezone for output of rising and setting time either as a known acronym
#' (eg. "GMT", "CET") or a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param dec (Optional) Output declination: \emph{geo} for the geocentric, or \emph{topo} for the topocentric
#' frame of reference. If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param refraction (Optional) Whether atmospheric refraction is to be taken into account.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temprature for erfraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link[swephR]{swe_azalt}}
#' @examples
#' # Rising and setting of the sun on june solstice 2018, from the location of London
#' riseset('sun', '2018/06/21', loc=c(51.5, 0.11, 100))
#'
#' # Rising ans setting of the moon on june solstice 2018, using a horizon profile
#' hor <- downloadHWT('HIFVTBGK') # Liverpool cathedral
#' riseset('moon', '2018/06/21', loc=hor)
#'
#' # Rising and setting of the sun throughout February 1999, from the location of London
#' riseset('sun', '1999/02', loc=c(51.5, 0.11, 100))
#'
#' # Rising and setting of the sun throughout 3000 BC, from the location of London
#' riseset('sun', -3000, loc=c(51.5, 0.11, 100))
riseset <- function(obj = 'sun', date, jd, alt=0, loc, calendar, timezone, dec, refraction, atm, temp, verbose=T) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }
  if (missing(dec)) { dec <- skyscapeR.env$dec }
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  out <- c()
  out$object <- obj
  body <- checkbody(obj) ## TODO add stars

  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- c(hor$metadata$georef, hor$metadata$elevation) }

  if (!missing(jd) & !missing(date)) { stop('Both date and JD were found, please provide only one.') }

  if (!missing(jd)) {
    date <- jd2time(jd, timezone, calendar)
    date <- substr(date,1,which(strsplit(date, "")[[1]]==" ")-1)
    date <- paste(date, '00:00:00')
    jd0 <- time2jd(date, timezone, calendar)
  }

  if (nchar(date)==10) {
    date <- paste(date, '00:00:00')
    jd0 <- time2jd(date, timezone, calendar)
  }

  if (class(date)=='numeric') {
    dat <- paste(date,'/01/01 00:00:00')
    jd0 <- time2jd(dat, timezone, calendar)
    jd <- seq(jd0, jd0+366, 1)
    ind <- which(substr(jd2time(jd, timezone, calendar), 1, which(strsplit(jd2time(jd, timezone, calendar), "")[[1]]=="/")[1]-1)==as.character(date))
    jd0 <- jd[ind]
  }
  if (nchar(date)==7) {
    dat <- paste0(date,'/01 00:00:00')
    jd0 <- time2jd(dat, timezone, calendar)
    jd <- seq(jd0, jd0+31, 1)
    ind <- which(substr(jd2time(jd, timezone, calendar), 1, which(strsplit(jd2time(jd, timezone, calendar), "")[[1]]=="/")[2]-1)==date)
    jd0 <- jd[ind]
  }

  rises <- data.frame(date=NA, time=NA, azimuth=NA, declination=NA, stringsAsFactors = F); sets <- rises
  if (length(jd0) > 1 & verbose) { pb <- txtProgressBar(max = length(jd0), style=3) }

  for (k in 1:length(jd0)) {
    rise <- swe_rise_trans_true_hor(jd0[k], body, '', 0, 1+256+ifelse(refraction,512,0), c(loc[2],loc[1],loc[3]), atm, temp, alt)$tret
    aux <- body.position(obj, rise, '', '', dec, loc, refraction, atm, temp, verbose=F)
    aux1 <- data.frame(azimuth = aux$horizontal$az, declination = aux$equatorial$Dec, time = jd2time(rise, timezone, calendar), stringsAsFactors = F)

    set <- swe_rise_trans_true_hor(jd0[k], body, '', 0, 2+256+ifelse(refraction,512,0), c(loc[2],loc[1],loc[3]), atm, temp, alt)$tret
    aux <- body.position(obj, set, '', '', dec, loc, refraction, atm, temp, verbose=F)
    aux2 <- data.frame(azimuth = aux$horizontal$az, declination = aux$equatorial$Dec, time = jd2time(set, timezone, calendar), stringsAsFactors = F)

    date <- substr(jd2time(jd0[k], timezone, calendar),1,which(strsplit(jd2time(jd0[k], timezone, calendar), "")[[1]]==" ")-1)
    rises[k,] <- c(date, substr(aux1$time,which(strsplit(aux1$time, "")[[1]]==" ")+1,nchar(aux1$time)), aux1$azimuth, aux1$declination)
    sets[k,] <- c(date, substr(aux2$time,which(strsplit(aux2$time, "")[[1]]==" ")+1,nchar(aux2$time)), aux2$azimuth, aux2$declination)

    if (length(jd0) > 1 & verbose) { setTxtProgressBar(pb, k) }
  }

  rises$azimuth <- as.numeric(rises$azimuth); rises$declination <- as.numeric(rises$declination)
  sets$azimuth <- as.numeric(sets$azimuth); sets$declination <- as.numeric(sets$declination)

  out$rise = rises
  out$set = sets
  return(out)
}




#' Calculate visible path of celestial object at given location
#'
#' This function calculates the visible path of a celestial
#' object from any location on earth. It outputs a \emph{skyscapeR.orbit}
#' object, which includes AZ and ALT information.
#' @param dec Declination of object.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param res The resolution (in degrees of RA) with which to calculate the path.
#' @param refraction (Optional) Whether atmospheric refraction is to be taken into account.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temprature for erfraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @import swephR
#' @export
#' @examples
#' # Visible path of sun on June Solstice on year 3999 BC from London:
#' sun.dec <- jS(-4000)
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon, 0 )
#' path <- orbit(sun.dec, loc)
#' plot(path$az, path$alt, ylim=c(0,90), type='l', xlab='AZ', ylab='ALT', col='red', lwd=2)
orbit = function(dec, loc, res=0.25, refraction, atm, temp) {
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (class(loc)=='skyscapeR.horizon') { loc <- loc$metadata$georef }

  ra <- seq(0, 360, by=res)
  aux <- array(NA, c(NROW(ra),2))

  for (i in 1:NROW(ra)) {
    tmp <- eq2hor(ra[i], dec, loc, refraction, atm, temp)

    aux[i,] <- c(tmp$az,tmp$alt)
  }

  # trim off and sort by aazimuth
  ind <- which (aux[,2] > -20)
  aux <- aux[ind,]

  ind <- sort(aux[,1], index.return=T)$ix
  aux <- aux[ind,]

  # return result
  orbit <- c()
  orbit$az <- aux[,1]
  orbit$alt <- aux[,2]
  orbit$dec <- dec
  orbit$georef <- loc
  class(orbit) <- "skyscapeR.orbit"
  return(orbit)
}



#' Returns the azimuth of the sun at a given time from a specific location
#'
#' This function returns the azimuth of the sun at a given time and location,
#' useful for data reduction of theodolite measurements using the sunsight
#' technique (\code{\link{reduct.theodolite}}).
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude, longitude and elevation of location, in this order.
#' @param time String containing the date and time in the following format:
#' "YYYY-MM-DD HH:MM:SS"
#' @param timezone Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London").
#' @param limb (Optional) Measured limb of the sun. Options are \emph{left}, \emph{right}.
#' If missing the centre of the sun will be output.
#' @param alt (Optional) Boolean that triggers output of altitude of the sun at exact time.
#' Default is FALSE.
#' @import swephR
#' @export
#' @seealso \code{\link{reduct.theodolite}}
#' @examples
#' sunAz(c(52,-3,100), '2017-10-04 12:32:14', 'Europe/London')
sunAz = function(loc, time, timezone, limb, alt=F) {
  if (missing(timezone)) { timezone <- skyscapeR.env$timezone }

  if (class(loc)=='skyscapeR.horizon') { loc <- c(loc$metadata$georef, loc$metadata$elevation) }
  if (is.null(dim(loc))) { dim(loc) <- c(1, NROW(loc)) }

  az <- c(); at <- c()
  for (i in 1:NROW(loc)) {
    pb.date <- as.POSIXct(time[i], timezone[i])
    UT <- format(pb.date, tz="UTC",usetz=TRUE)
    UT <- as.POSIXlt(UT, 'UTC')

    swephR::swe_set_topo(loc[i,2], loc[i,1], loc[i,3])
    jd <- swephR::swe_julday(UT$year+1900, UT$mon+1, UT$mday, UT$hour+UT$min/60+UT$sec/3600,1)
    ss <- swephR::swe_calc_ut(jd, 0, 32*1024+2048)
    xin <- ss$xx[1:2]
    aux <- swephR::swe_azalt(jd, 1, c(loc[i,2],loc[i,1],loc[i,3]), 1013.25, 15, xin)
    az[i] <- aux$xaz[1]+180
    if (az[i] > 360) { az[i] <- az[i]-360 }

    if (!missing(limb)) {
      if (limb=="left") { az[i] <- az[i] - 32/60/2 }
      if (limb=="right") { az[i] <- az[i] + 32/60/2 }
    }
    if (alt) { at[i] <- aux$xaz[2] }
  }
  if (alt) {
    df <- data.frame(az = az, alt = at)
    return(df)
  } else {
    return(az)
  }
}


#' Solar Date
#'
#' Returns the calendar date when the sun has the same declination as the input declinations.
#' @param dec Single value or array of declinations.
#' @param year Year for which to do calculations.
#' @param calendar (Optional) Calendar used for output. G for gregorian and J for julian. Defaults to Gregorian.
#' @import swephR
#' @export
#' @examples
#' solar.date(-23, 2018)
#' solar.date(-23, 1200, calendar='G')
#' solar.date(-23, 1200, calendar='J')
solar.date <- function(dec, year, calendar){
  if (missing(calendar)) { calendar <- skyscapeR.env$calendar }

  jd0 <- time2jd(timestring(year,1,1,12), timezone='UTC', calendar)

  out <- c(); cdate <- c()
  for (i in 1:366) {
    cdate[i] <- substr(jd2time(jd0+i-1, calendar=calendar),6,10)
    out[i] <- body.position(obj='sun', jd0+i-1, timezone='UTC', calendar=calendar, verbose = FALSE)$equatorial$Dec
  }
  rr <- range(out)

  ind <- which(dec >= rr[1] & dec <= rr[2])
  if (length(ind)< length(dec)) { message(paste(length(dec)-length(ind),' value(s) were outside of solar range and have been excluded.')) }
  dec <- dec[ind]
  aux <- matrix(NA, ncol=length(dec), nrow=3)
  for (i in 1:length(dec)) {
    ff <- approxfun(1:366, out-dec[i])
    dd <- round(rootSolve::uniroot.all(ff, interval=c(1,366)))
    tt <- rbind(long.date(cdate[dd])); if (length(tt)==1) { tt <- c(tt,NA) }
    aux[,i] <- c(dec[i], tt)
  }

  rownames(aux) <- c('Dec','Date1', 'Date2')
  return(aux)
}



#' Corrected parallax for a given location and object altitude
#'
#' Given the average parallax, this fucntion corrects this value for a given latitude of the observer
#' and for the altitude of the celestial object.
#' @param parallax Average parallax to correct (e.g. 0.00224 for the Sun, or 0.952 for the Moon)
#' @param loc This can be either the latitude of the location, or a \emph{skyscapeR.horizon} object.
#' @param altitude (Optional) Altitude of the celestial object.. Defaults to 0 degrees.
#' @import swephR
#' @export
#' @examples
#' # Parallax correction for the moon, as seen from latitude 50ºN and at 0º altitude
#' parallax.corr(0.952, 50, 0)
parallax.corr <- function(parallax, loc, altitude=0) {
  latitude <- extractLatitude(loc)
  return(parallax*cos(altitude/180*pi)*sin(latitude/180*pi)) # V Reijs/SE formula
  # return(parallax * cos(altitude/180*pi) * (1-sin(latitude/180*pi)^2/298.3)) # Nautical formula
}
