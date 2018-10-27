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
  rr <- eq2horFS(xx[1],xx[2],xx[3],lat,lon, precess_=F)$alt
  return(rr)
}

#' @noRd
findWS <- function(year) {
  jd0 <- astrolibR::juldate(c(year,1,1,12)) + 2400000
  yy <- seq(jd0, length.out = 365)
  ss <- astrolibR::sunpos(yy)$dec
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

#' Fixed eq2hor function from astrolibR package
#' @noRd
eq2horFS = function (ra, dec, jd, lat = 43.0783, lon = -89.865, ws = F,
                     obsname, b1950, precess_ = TRUE, nutate_ = TRUE, refract_ = TRUE,
                     aberration_ = TRUE, altitude = 0, ...)
{
  d2r = pi/180
  h2r = pi/12
  h2e = 15
  if (!missing(b1950))
    s_now = "   (j1950)"
  else s_now = "   (j2000)"
  j_now = (jd - 2451545)/365.25 + 2000
  if (precess_) {
    if (!missing(b1950)) {
      for (i in 1:length(jd)) {
        tmp = astrolibR::precess(ra[i], dec[i], 1950, j_now[i],
                                 fk4 = TRUE)
        ra[i] = tmp$ra
        dec[i] = tmp$dec
      }
    }
    else {
      for (i in 1:length(jd)) {
        tmp = astrolibR::precess(ra[i], dec[i], 2000, j_now[i])
        ra[i] = tmp$ra
        dec[i] = tmp$dec
      }
    }
  }
  tmp = astrolibR::co_nutate(jd, ra, dec)
  dra1 = tmp$d_ra
  ddec1 = tmp$d_dec
  eps = tmp$eps
  d_psi = tmp$d_psi
  tmp = astrolibR::co_aberration(jd, ra, dec, eps)
  dra2 = tmp$d_ra
  ddec2 = tmp$d_dec
  eps = tmp$eps
  ra = ra + (dra1 * nutate_ + dra2 * aberration_)/3600
  dec = dec + (ddec1 * nutate_ + ddec2 * aberration_)/3600
  lmst = astrolibR::ct2lst(lon, 0, jd)
  lmst = lmst * h2e
  last = lmst + d_psi * cos(eps)/3600
  ha = last - ra
  w = (ha < 0)
  ha[w] = ha[w] + 360
  ha = ha%%360
  tmp = astrolibR::hadec2altaz(ha, dec, lat, ws = ws)
  alt = tmp$alt
  az = tmp$az
  if (refract_)
    alt = astrolibR::co_refract(alt, altitude = altitude, ..., to_observed = TRUE)

  return(list(alt = alt, az = az, ha = ha))
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
