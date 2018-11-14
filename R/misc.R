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
  if (sum(sapply(as.character(stars$NAME), pracma::strcmp, s2=name))) {
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
