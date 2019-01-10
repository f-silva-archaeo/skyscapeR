#jd <- swephR::swe_julday(2000,1,1,12,1) # J2000.0
cur.year <- as.numeric(format(Sys.Date(), "%Y")) # current year
swephR::swe_set_ephe_path(system.file("ephemeris", "", package = "swephRdata"))


#' Computes obliquity of the ecliptic
#'
#' This function calculates the obliquity for a given epoch. It is a
#' wrapper for function \code{\link[swephR]{swe_calc}} of package \emph{swephR}.
#' @param year Year for which to calculate the obliquity.
#' Defaults to present year as given by Sys.Date()
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc}}
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @examples
#' #' # Obliquity for year 3999 BC:
#' obliquity(-4000)
obliquity = function(year = cur.year) {
  jd <- swephR::swe_julday(year,1,1,12,1) # J2000.0
  aux <- swephR::swe_calc(jd, -1, 0)$xx[1]
  return(aux)
}



#' Computes position of Solar System bodies in geocentric equatorial coordinates
#'
#' This function calculates the geocentric declination and right ascension of solar
#'  system bodies at a given time. It is a wrapper for function
#'  \code{\link[swephR]{swe_calc}} of package \emph{swephR}.
#' @param body (Optional) String containing name of the solar system body of interest. Can be
#' any of the planets (inc. Pluto), the Moon, the Sun or the Ecliptic. Defaults to 'sun'.
#' @param time Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. Defaults to GMT.
#' @param cal (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. Defaults to Gregorian.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc}}, \code{\link{timestring}}, \code{\link{time2jd}}
#' @examples
#' # Position of the sun at noon GMT on Christmas day 2018:
#' body.position('sun', '2018-12-25 12:00:00', timezone='GMT')
#'
#' # Declination of the moon at same time
#' body.position('moon', '2018-12-25 12:00:00', timezone='GMT')$Dec
body.position = function(body='sun', time, timezone='Europe/London', cal='G') {
  body <- checkbody(body)

  if (class(time)=='character') { jd <- time2jd(time, timezone, cal) } else { jd <- time }
  aux <- swephR::swe_calc(jd, body, 2048)

  out <- c()
  out$RA <- aux$xx[1]
  out$Dec <- aux$xx[2]
  return(out)
}


#' Computes the phase of the moon
#'
#' This function calculates the moon phase, in percentage of full. It is a wrapper
#' for function \code{\link[swephR]{swe_pheno_ut}} of package \emph{swephR}.
#' @param time Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param timezone (Optional) Timezone of input either as a known acronym (eg. "GMT", "CET") or
#' a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Only needed if \emph{time} is a string. Defaults to GMT.
#' @param cal (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' Only needed if \emph{time} is a string. Defaults to Gregorian.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_pheno_ut}}
#' @examples
#' # Moonphase at noon GMT on Christmas day 2018:
#' moonphase('2018-12-25 12:00:00', 'GMT')
moonphase <- function(time, timezone='Europe/London', cal='G') {
  if (class(time)=='character') { jd <- time2jd(time, timezone, cal) } else { jd <- time }
  return(swephR::swe_pheno_ut(jd, 1, 0)$attr[2])
}


#' @import swephR
#' @noRd
vecAzAlt <- function(jd, body, loc, refraction=T, atm=1013.25, temp=15) {
  swephR::swe_set_topo(loc[2], loc[1],0)
  xin <- swephR::swe_calc_ut(jd, body, 32*1024+2048)$xx[1:2]
  aux <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), atm, temp, xin)$xaz
  if (refraction) { aux <- aux[c(1,3)] } else { aux <- aux[1:2] }
  return(cbind(xin,aux))
}



#' Computes the rising and setting azimuth, declination and time of a Solar System object
#' for a given location and day
#'
#' @param body (Optional) String containing name of the solar system body of interest. Can be
#' any of the planets (inc. Pluto), the Moon, the Sun or the Ecliptic. Defaults to 'sun'.
#' @param date Either a string containing the date in the format "YYYY-MM-DD"
#'  (see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).
#' @param cal (Optional) Calendar used in parameter \emph{date}. G for gregorian and J for julian. Defaults to Gregorian.
#' @param timezone (Optional) Timezone for output of rising and setting time either as a known acronym
#' (eg. "GMT", "CET") or a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Defaults to GMT.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param refraction (Optional) Boolean for whether or not atmospheric refraction should be taken into account.
#' Defaults to _TRUE_.
#' @param atm (Optional) Atmospheric pressure (in mbar). Only needed if \emph{refraction} is set to _TRUE_. Default is 1013.25 mbar.
#' @param temp (Optional) Atmospheric temperature (in Celsius). Only needed if \emph{refraction} is set to _TRUE_. Default is 15 degrees.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link[swephR]{swe_azalt}}
#' @examples
#' # Rising and setting of the sun on june solstice 2018, from the location of London
#' riseset('sun', '2018-06-21', loc=c(51.5, 0.11))
#'
#' # Rising ans setting of the moon on june solstice 2018, using a horizon profile
#' hor <- downloadHWT('HIFVTBGK') # Liverpool cathedral
#' riseset('moon', '2018-06-21', loc=hor)
riseset <- function(body = 'sun', date, cal='G', timezone='GMT', loc, refraction=T, atm=1013.25, temp=15) {
  body <- checkbody(body)

  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- hor$metadata$georef }

  if (class(date)=='character') {
    if (nchar(date)==10) { date <- paste(date, '00:00:00') }
    jd0 <- time2jd(date, timezone='GMT', cal)
    } else { jd0 <- date }

  jd <- seq(jd0, jd0+1, 0.001)
  ss <- sapply(jd, vecAzAlt, body, loc=loc, refraction=refraction, atm=atm, temp=temp)

  # azimuth
  if(exists('hor')) {
    f <- approxfun(ss[3,], ss[4,] - hor2alt(hor, ss[3,])$alt) ## shifted by visible horizon, when available
    message('Using given horizon profile.')
  } else {
    f <- approxfun(ss[3,], ss[4,])
  }
  az <- rootSolve::uniroot.all(f, interval=c(0,360))

  # declination
  dd <- approxfun(ss[3,], ss[2,])
  dec <- dd(az)

  # time
  tt <- approxfun(ss[3,], jd)
  time <- tt(az)
  aux <- c()
  aux[1] <- substr(jd2time(time[1], timezone),12,21)
  aux[2] <- substr(jd2time(time[2], timezone),12,21)


  # swephR calculates azimuths clockwise from south
  az <- az + 180
  az[az>=360] <- az[az>=360]-360

  out <- c()
  ind <- which(az<180)
  out$rise <- data.frame(azimuth = az[ind], declination = dec[ind], time = aux[ind])

  ind <- which(az>180)
  out$set <- data.frame(azimuth = az[ind], declination = dec[ind], time = aux[ind])
  return(out)
}

#' Computes the daily rising and setting azimuth, declination and time of a Solar System object
#' for a given location (and horizon) for a full year
#'
#' @param body (Optional) String containing name of the solar system body of interest. Can be
#' any of the planets (inc. Pluto), the Moon, the Sun or the Ecliptic. Defaults to 'sun'.
#' @param year Year for which to calculate rising and settings.
#' @param cal (Optional) Calendar used for output. G for gregorian and J for julian. Defaults to Gregorian.
#' @param timezone (Optional) Timezone for output of rising and setting time either as a known acronym
#' (eg. "GMT", "CET") or a string with continent followed by country capital (eg. "Europe/London"). See
#' \link{timezones} for details. Defaults to GMT.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param refraction (Optional) Boolean for whether or not atmospheric refraction should be taken into account.
#' Defaults to _TRUE_.
#' @param atm (Optional) Atmospheric pressure (in mbar). Only needed if \emph{refraction} is set to _TRUE_. Default is 1013.25 mbar.
#' @param temp (Optional) Atmospheric temperature (in Celsius). Only needed if \emph{refraction} is set to _TRUE_. Default is 15 degrees.
#' @param ncores (Optional) Number of processing cores to use for parallelisation. Defaults to the number of
#' available cores minus 1.
#' @param verbose (Optional) Boolean to decide whether or not
#' @import parallel foreach doParallel swephR
#' @export
#' @seealso \code{\link[swephR]{swe_calc_ut}}, \code{\link[swephR]{swe_azalt}}
#' @examples
#' # daily rising and setting of the moon in 2006, from the location of London
#' riseset.year('moon', 2006, loc=c(51.5, 0.11))
riseset.year <- function(body = 'sun', year, cal='G', timezone='GMT', loc, refraction=T, atm=1013.25, temp=15, ncores=parallel::detectCores()-1, verbose=T) {
  body <- checkbody(body)

  if (class(loc)=='skyscapeR.horizon') { hor <- loc; loc <- hor$metadata$georef }

  jd0 <- time2jd(paste(year,'01','01',sep='-'), timezone=timezone, cal)

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, library(skyscapeR))
  parallel::clusterEvalQ(cl, library(swephR))
  if (verbose) { message(paste0('Running on ', ncores, ' processing cores.')) }

  out <- matrix(NA, nrow=365, ncol=7)
  out <- foreach::foreach(i = 1:365, .combine=rbind, .inorder = F) %dopar% {
    jd <- seq(jd0+i-1, jd0+i, 0.001)
    ss <- sapply(jd, vecAzAlt, body, loc=loc, refraction=refraction, atm=atm, temp=temp)

    # azimuth
    if(exists('hor')) {
      f <- approxfun(ss[3,], ss[4,] - hor2alt(hor, ss[3,])$alt) ## shifted by visible horizon, when available
      message('Using given horizon profile.')
    } else {
      f <- approxfun(ss[3,], ss[4,])
    }
    az <- rootSolve::uniroot.all(f, interval=c(0,360))

    # declination
    dd <- approxfun(ss[3,], ss[2,])
    dec <- dd(az)

    # time
    tt <- approxfun(ss[3,], jd)
    time <- tt(az)
    aux <- c()
    aux[1] <- substr(jd2time(time[1], timezone),12,21)
    aux[2] <- substr(jd2time(time[2], timezone),12,21)

    # swephR calculates azimuths clockwise from south
    az <- az + 180
    az[az>=360] <- az[az>=360]-360

    c(i, az, dec, aux)
  }
  aux <- out
  ind <- sort(as.numeric(aux[,1]), index.return=T)$ix
  aux <- aux[ind,]
  parallel::stopCluster(cl)
  if (verbose) { message('Done.') }

  out <- c()
  ind <- which(as.numeric(aux[1,2:3])<180)
  out$rise <- data.frame(day = as.numeric(aux[,1]), azimuth = as.numeric(aux[,1+ind]), declination = as.numeric(aux[,3+ind]), time = aux[,5+ind])
  row.names(out$rise) <- NULL

  ind <- which(as.numeric(aux[1,2:3])>180)
  out$set <- data.frame(day = as.numeric(aux[,1]), azimuth = as.numeric(aux[,1+ind]), declination = as.numeric(aux[,3+ind]), time = aux[,5+ind])
  row.names(out$set) <- NULL

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
#' @param refraction (Optional) If set to TRUE it will calculate apparent rather
#' than true altitudes. Default is TRUE.
#' @param ...  Any other parameters to be passed unto \code{\link[swephR]{swe_azalt}}.
#' @import swephR
#' @export
#' @examples
#' # Visible path of sun on June Solstice on year 3999 BC from London:
#' sun.dec <- jS(-4000)
#' london.lat <- 51.5074 #N
#' london.lon <- -0.1278 #W
#' loc <- c( london.lat, london.lon )
#' path <- orbit(sun.dec, loc)
#' plot(path$az, path$alt, ylim=c(0,90), type='l', xlab='AZ', ylab='ALT', col='red', lwd=2)
orbit = function(dec, loc, res=0.25, refraction=T, ...) {
  if (class(loc)=='skyscapeR.horizon') { loc <- loc$metadata$georef }

  ra <- seq(0, 360, by=res)
  aux <- array(NA, c(NROW(ra),2))

  for (i in 1:NROW(ra)) {
    tmp <- eq2horFS(ra[i], dec, time2jd('2000-01-01 12:00'), loc, refraction, ...)

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
#' containing the latitude and longitude of location, in this order.
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
#' sunAz(c(52,-3), '2017-10-04 12:32:14', 'Europe/London')
sunAz = function(loc, time, timezone = 'Europe/London', limb, alt=F) {
  if (class(loc)=='skyscapeR.horizon') { loc <- loc$metadata$georef }
  if (is.null(dim(loc))) { dim(loc) <- c(1, NROW(loc)) }

  az <- c(); at <- c()
  for (i in 1:NROW(loc)) {
    pb.date <- as.POSIXct(time[i], timezone[i])
    UT <- format(pb.date, tz="UTC",usetz=TRUE)
    UT <- as.POSIXlt(UT, 'UTC')

    swephR::swe_set_topo(loc[i,1], loc[i,2], 0)
    jd <- swephR::swe_julday(UT$year+1900, UT$mon+1, UT$mday, UT$hour+UT$min/60+UT$sec/3600,1)
    ss <- swephR::swe_calc_ut(jd, 0, 32*1024+2048)
    xin <- ss$xx[1:2]
    aux <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],0), 1013.25, 15, xin)
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
#' @param cal (Optional) Calendar used for output. G for gregorian and J for julian. Defaults to Gregorian.
#' @import swephR
#' @export
#' @examples
#' solar.date(-23, 2018)
#' solar.date(-23, 1200, cal='G')
#' solar.date(-23, 1200, cal='J')
solar.date <- function(dec, year, cal='G'){

  jd0 <- time2jd(timestring(year,1,1,12), timezone='UTC', cal)

  out <- c(); cdate <- c()
  for (i in 1:366) {
    cdate[i] <- substr(jd2time(jd0+i-1, cal=cal),6,10)
    out[i] <- body.position(body='sun', jd0+i-1, timezone='UTC', cal=cal)$Dec
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
