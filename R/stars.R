#' Create \emph{skyscapeR.star} object
#'
#' This function retrieves information for a given star
#' and saves it in the \emph{skyscapeR.star} format ready to be
#' used by other skyscapeR package function.
#' @param string This can be either the traditional name for the star or its Bayer designation.
#' @param year Year for which to calculate the coordinates.
#' Defaults to current year.
#' @import swephR
#' @export
#' @seealso \code{\link[swephR]{swe_fixstar2_ut}}, \code{\link[swephR]{swe_fixstar2_mag}}
#' @examples
#' # Retrieve data for Aldebaran:
#' Aldeb <- star('Aldebaran')
#'
#' # Retrieve data for Aldebaran on 2999 BC:
#' ss <- star('Aldebaran', -3000)
star <- function(string, year=skyscapeR.env$cur.year) {
  aux <- star.names$identifier[which(star.names$Western==string)]
  if (length(aux)==1) {
    name <- string
    string <-  aux
  } else {
    name <- star.names$Western[which(star.names$identifier==string)]
  }

  string <- paste0(',',string)
  aux <- swephR::swe_fixstar2_ut(string, swephR::swe_julday(2000, 1, 1, 12,1), 2048+16384+16)
  if (aux$serr!="") stop('Star not recognised.')

  aux <- data.frame(year=NA, RA=NA, Dec=NA)
  for (i in 1:length(year)) {
    info <- swephR::swe_fixstar2_ut(string, swephR::swe_julday(year[i], 1, 1, 12,1), 2048+16384+16)
    aux[i,] <- c(year[i], info$xx[1], info$xx[2])
  }

  if (length(name)==0) { name <- ''}

  star <- c()
  star$name <- name
  star$Bayer <- strsplit(info$star,',')[[1]][2]
  star$app.mag <- swephR::swe_fixstar2_mag(paste(string))$mag
  star$coord <- aux
  class(star) <- "skyscapeR.star"

  return(star)
}


#' @noRd
ff <- function(x, loc, atm, temp, pos) { return(swephR::swe_azalt(x, 1, c(loc[2],loc[1],loc[3]), atm, temp, c(pos[1], pos[2]))$xaz) }

#' @noRd
calc.phase <- function(year, day, sun.coord, star, loc, arcus_visionis, alt.hor, k, limit, alt.rs, res, refraction, atm, temp, pb, pbi) {
  phase <- rep("", length(day))

  for (i in 1:length(day)) {
    jd0 <- time2jd('2000/01/01 00:00:00')
    jd <- seq(jd0, jd0 + 1, res)
    # Sun
    Sun.alt <- sapply(jd, ff, loc, atm, temp, sun.coord[day[i],])[3,]
    Sun.alt <- approx(jd, Sun.alt, xout=seq(jd0, jd0 + 1, 1/24/60))$y
    # Star
    Star.alt <- sapply(jd, ff, loc, atm, temp, as.numeric(star$coord[2:3]))[3,]
    Star.alt <- approx(jd, Star.alt, xout=seq(jd0, jd0 + 1, 1/24/60))$y

    # atmospheric extinction
    airmass <- 1/(cos((90-Star.alt)/180*pi) + 0.025*exp(-11*cos((90-Star.alt)/180*pi)))
    mag <- star$app.mag + airmass*k; mag[mag<0] <- NA

    # Q1: can the star be seen?
    ind <- which(Sun.alt < alt.hor & Star.alt > alt.hor & Star.alt-Sun.alt >= arcus_visionis & mag <= limit)

    if (length(ind) > 1) {
      alt <- Star.alt[ind]
      ind2 <- which(alt >= alt.hor & alt <= alt.rs)
      if (length(ind2) > 1) {
        test <- split(ind2, cumsum(c(1,abs(diff(ind[ind2])) > 1)))
        if (length(which(lengths(test)<2))>0) { test <- test[-which(lengths(test)<2)] }
        # Q2: can the star be seen rising or setting?
        for (j in 1:NROW(test)) {
          if (mean(diff(alt[test[[j]]])) > 0) {
            phase[i] <- paste0(phase[i], 'R')
          } else {
            phase[i] <- paste0(phase[i], 'S')
          }
        }
      } else { phase[i] <- 'V' }
    } else {
      # invisibility
      phase[i] <- 'I'
    }

    setTxtProgressBar(pb, pbi+i)
  }
  return(phase)
}

#' Calculate the seasons and phase type of a star
#'
#' This function calculates the seasons (Rising, Setting, etc.)
#' and phase types (Arising and Lying Hidden, Curtailed Passage) of a
#' star for a given location and epoch. This functions uses the
#' \emph{arcus visionis} approximation of Purrington (1988) and
#' the atmospheric extinction approximation of Schaefer (1989). For
#' the nomenclature used, and description of star phase types, see Brady (2015).
#' @param star Either the star name or a \emph{skyscapeR.star} object.
#' @param year The year of interest.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param alt.hor (Optional) The altitude of the horizon to consider.
#' Defaults to zero degrees.
#' @param k (Optional) Extinction coefficient (see Schaefer 1989).
#' Defaults to 0.2, corresponding to a poor night on mountain top
#' or best night at a dry sea level site.
#' @param limit (Optional) The maximum magnitude of a star that can be
#' visible with the naked eye. Defaults to 6.
#' @param alt.rs (Optional) The maximum altitude of a star's first or last
#' visibility for it to still be considered to be as rising or setting.
#' Defaults to ten degrees.
#' @param res (Optional) Resolution of calculation. The smaller this
#' figure the slower the computation. Defaults to 1/24/6, i.e. every 10 minutes.
#' @param refraction (Optional) Whether atmospheric refraction is to be taken into account.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temperature for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @export
#' @seealso \code{\link{plot.skyscapeR.starphases}}
#' @references Purrington, Robert D. (1988) Heliacal Rising and Setting:
#' Quantitative Aspects, \emph{Journal for the History of Astronomy
#' (Archaeoastronomy Supplement 12)} 19, S72-S84. Available online at
#' [SAO/NASA ADS Astronomy Abstract Service](http://adsabs.harvard.edu/abs/1988JHAS...19...72P)
#' @references Brady, Bernadette (2015) Star Phases: the Naked-eye Astronomy of the Old Kingdom
#' Pyramid Texts. In F Silva and N Campion (eds) \emph{Skyscapes: The Role and Importance of
#' the Sky in Archaeology}. Oxford: Oxbow Books, pp. 76-86.
#' @examples
#' ss1 <- star.phases('Aldebaran',-4000, c(35,-8,200))
#'
#' # One can then look at the star's phase type:
#' ss1$metadata$type
#'
#' # Date range of seasons:
#' ss1$metadata$seasons
#'
#' # Date range of phase-type events:
#' ss1$metadata$events
#'
#' # And plot them:
#' plot(ss1)
#'
#' # You can play with the parameters and see how predictions change:
#' ss1 <- star.phases('Aldebaran',-4000, c(35,-8,200), alt.hor=2, alt.rs=5)
#' plot(ss1)
star.phases <- function(star, year, loc, alt.hor = 0, k = 0.2, limit = 6, alt.rs = 10, res = 1/24/6, refraction, atm, temp) {
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$metadata$georef[1]
    lon <- loc$metadata$georef[2]
    elev <- loc$metadata$georef[3]
  } else {
    lat <- loc[1]
    lon <- loc[2]
    elev <- loc[3]
  }
  swephR::swe_set_topo(lon, lat, elev)

  # load star data
  if (class(star)!='skyscapeR.star') { star <- star(star, year) }
  arcus_visionis <- 2.1*star$app.mag + 10

  # check if circumpolar or always invisible
  if (lat>=0) {
    if (star$coord$Dec < -(90-lat+alt.hor)) stop('Star is always below the horizon at this location and epoch.')
    if (star$coord$Dec >= (90-lat+alt.hor)) stop('Star is always above the horizon at this location and epoch.')
  }
  if (lat<0) {
    if (star$coord$Dec > (90-lat)) stop('Star is always below the horizon at this location and epoch.')
    if (star$coord$Dec <= -(90-lat)) stop('Star is always above the horizon at this location and epoch.')
  }

  # pre-calculations of solar position
  ob <- obliquity(year)
  ll <- seq(0,359.9,0.1)
  Dec <- asin(sin(ob/180*pi)*sin(ll/180*pi))/pi*180
  RA <- atan2(sin(ll/180*pi)*cos(ob/180*pi), cos(ll/180*pi))/pi*180
  RA[RA<0] <- RA[RA<0]+360
  lon <- seq(0,359.99,by=360/365)
  dd <- approx(seq(0,359.9,0.1), Dec, xout=lon)$y
  rr <- approx(seq(0,359.9,0.1), RA, xout=lon)$y
  sun.pos <- cbind(rr,dd)

  # calculations
  phase <- rep('', 365)
  star.pos <- as.numeric(star$coord[2:3])

  pb <- txtProgressBar(max=12+2*4+4*7, style=3)
  # step 1: once a month
  day <- seq(1,365,30)
  aux <- calc.phase(year, day, sun.pos, star, loc, arcus_visionis, alt.hor, k, limit, alt.rs, res=1/24/2, refraction, atm, temp, pb, 0)
  phase[day] <- aux

  # step 2: fill in blanks and identify months where star season changes
  ind <- which(phase!=""); pp <- phase[ind]
  change <- c(0,which(pp!=c(pp[-1], pp[1])))
  for (i in 1:(length(change)-1)) {
    phase[ind[change[i]+1]:ind[change[i+1]]] <- phase[ind[change][i]]
  }

  # step 3: re-run on months identified before
  ind <- which(phase=="")
  ind <- split(ind, cumsum(seq_along(ind) %in% (which(diff(ind)>1)+1)))
  ind <- ind[which(sapply(ind, length)>=14)]

  pb0 <- 12
  for (i in 1:NROW(ind)) {
    days <- ind[[i]][c(7,14)]
    aux <- calc.phase(year, days, sun.pos, star, loc, arcus_visionis, alt.hor, k, limit, alt.rs, res=1/24/2, refraction, atm, temp, pb, pb0)
    phase[days] <- aux
    pb0 <- pb0+2
  }

  # step 4: fill in blanks and identify weeks where star season changes
  ind <- which(phase!=""); pp <- phase[ind]
  change <- c(0,which(pp!=c(pp[-1], pp[1])))
  for (i in 1:(length(change)-1)) {
    phase[ind[change[i]+1]:ind[change[i+1]]] <- phase[ind[change][i]]
  }

  # step 5: re-run on weeks identified before with higher resolution
  ind <- which(phase=="")
  ind <- split(ind, cumsum(seq_along(ind) %in% (which(diff(ind)>1)+1)))

  pb0 <- 20
  for (i in 1:NROW(ind)) {
    days <- ind[[i]]
    aux <- calc.phase(year, days, sun.pos, star, loc, arcus_visionis, alt.hor, k, limit, alt.rs, res, refraction, atm, temp, pb, pb0)
    phase[days] <- aux
    pb0 <- pb0 + length(ind[[i]])
  }

  # find december solstice date and pin calendar to it
  days <- seq(1,365)-which.min(dd)
  days[days<0] <- days[days<0]+365

  ind <- sort(days, index.return=T)$ix
  days <- days[ind]; phase <- phase[ind]

  uu <- unique(phase)
  if (length(uu)==1) {
    if (uu=='V') { type <- 'circumpolar'}
    if (uu=='I') { type <- 'invisible' }
  } else {
    if ('V' %in% uu) { type <- 'curtailed passage' }
    if ('I' %in% uu) { type <- 'arising and lying hidden' }
    if ('I' %in% uu & 'V' %in% uu) { type <- 'dual phase'}
  }

  # identify events
  events <- data.frame(event=NA, day=NA)
  ttt <- which(phase=='I')
  if (length(ttt)>0) {
    events[nrow(events)+1,] <- data.frame(event='acronycal setting', day=days[min(ttt)-1])
    events[nrow(events)+1,] <- data.frame(event='heliacal rising', day=days[max(ttt)+1])
  }
  ttt <- which(phase=='V')
  if (length(ttt)>0) {
    events[nrow(events)+1,] <- data.frame(event='acronycal rising', day=days[min(ttt)-1])
    events[nrow(events)+1,] <- data.frame(event='heliacal setting', day=days[max(ttt)+1])
  }
  events <- events[-1,]
  events <- events[sort(events$day, index.return=T)$ix,]
  rownames(events) <- NULL
  events$day <- dd.to.DD(events$day, char=T, WS=T)

  # cleanup
  phase <- sapply(phase, function(x) { paste(unique(strsplit(x, "")[[1]]), collapse='')}, USE.NAMES =F)

  # identify seasons
  seasons <- data.frame(season=NA, begin=NA, end=NA, length=NA)
  ttt <- which(phase=='R')
  if (length(ttt)>0) {
    if (sum(diff(ttt)>1)>0) {
      db <- days[ttt[which(diff(ttt)>1)+1]]; de <- days[ttt[which(diff(ttt)>1)]]
    } else {
      db <- days[min(ttt)]; de <- days[max(ttt)]
    }
    seasons[nrow(seasons)+1,] <- data.frame(season='rise only', begin=db, end=de)
  }

  ttt <- which(phase=='SR'); phase[ttt] <- 'RS'
  ttt <- which(phase=='RS')
  if (length(ttt)>0) {
    if (sum(diff(ttt)>1)>0) {
      db <- days[ttt[which(diff(ttt)>1)+1]]; de <- days[ttt[which(diff(ttt)>1)]]
    } else {
      db <- days[min(ttt)]; de <- days[max(ttt)]
    }
    seasons[nrow(seasons)+1,] <- data.frame(season='rise and set', begin=db, end=de)
  }

  ttt <- which(phase=='S')
  if (length(ttt)>0) {
    if (sum(diff(ttt)>1)>0) {
      db <- days[ttt[which(diff(ttt)>1)+1]]; de <- days[ttt[which(diff(ttt)>1)]]
    } else {
      db <- days[min(ttt)]; de <- days[max(ttt)]
    }
    seasons[nrow(seasons)+1,] <- data.frame(season='set only', begin=db, end=de)
  }

  ttt <- which(phase=='I')
  if (length(ttt)>0) {
    if (sum(diff(ttt)>1)>0) {
      db <- days[ttt[which(diff(ttt)>1)+1]]; de <- days[ttt[which(diff(ttt)>1)]]
    } else {
      db <- days[min(ttt)]; de <- days[max(ttt)]
    }
    seasons[nrow(seasons)+1,] <- data.frame(season='arising and lying hidden', begin=db, end=de)
  }

  ttt <- which(phase=='V')
  if (length(ttt)>0) {
    if (sum(diff(ttt)>1)>0) {
      db <- days[ttt[which(diff(ttt)>1)+1]]; de <- days[ttt[which(diff(ttt)>1)]]
    } else {
      db <- days[min(ttt)]; de <- days[max(ttt)]
    }
    seasons[nrow(seasons)+1,] <- data.frame(season='curtailed passage', begin=db, end=de)
  }

  seasons <- seasons[-1,]
  seasons <- seasons[sort(seasons$begin, index.return=T)$ix,]
  rownames(seasons) <- NULL
  seasons$length <- abs(seasons$end-seasons$begin)
  seasons$begin <- dd.to.DD(seasons$begin, char=T, WS=T)
  seasons$end <- dd.to.DD(seasons$end, char=T, WS=T)
  seasons

  # results
  out <- c()
  out$metadata$star <- star
  out$metadata$year <- year
  out$metadata$loc <- loc
  out$metadata$type <- type
  out$metadata$events <- events
  out$metadata$seasons <- seasons
  out$data <- data.frame(day = days, phase = phase)
  class(out) <- 'skyscapeR.starphases'
  return(out)
}




