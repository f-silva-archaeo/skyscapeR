cur.year <- as.numeric(format(Sys.Date(), "%Y")) # current year

#' Create \emph{skyscapeR.star} object
#'
#' This function retrieves information for a given star
#' and saves it in the \emph{skyscapeR.star} format ready to be
#' used by other skyscapeR package function.
#' @param string This can be either: 1) the common name for the star,
#'  2) the Hipparchos catalogue number (if string begins with HIP),
#'  3) the Bright Star Catalogue number (if string begins with HR),
#'  4) the Messier catalogue number (if string begins with M), or
#'  5) the Bayer designation (if string has exactly 5 characters,
#'   the second of which is a space).
#' @param year Year for which to calculate the coordinates.
#' Defaults to J2000.0 epoch.
#' @export
#' @seealso \code{\link[astrolibR]{precess}}, \code{\link[astrolibR]{co_nutate}},
#' \code{\link[astrolibR]{co_aberration}}
#' @examples
#' # Retrieve data for Aldebaran:
#' Aldeb <- star('Aldebaran')
#'
#' # Retrieve data for the Pleiades (M45):
#' P1 <- star('Pleiades')
#' P2 <- star('M45')
#'
#' # Retrieve data for Sirius on 2999 BC:
#' ss <- star('Sirius', -3000)
star <- function(string, year) {
  data(stars, envir = environment())
  ind <- which(stars$NAME == string)
  if (substr(string,1,3)=="HIP") { ind <- which(stars$HIP.ID == string)}
  if (substr(string,1,2)=="HR") { ind <- which(stars$HR.ID == string)}
  if (substr(string,1,1)=="M" & nchar(string)==3) { ind <- which(stars$MESSIER.ID == string)}
  if (nchar(string)==5 & substr(string,2,2)==" ") { ind <- which(stars$BAYER.DESIGNATION == string)}

  if (length(ind) == 0) { stop('No star with that name or designation has been found.') }

  star <- c()
  star$name <- as.character(stars$NAME[ind])
  star$constellation <- as.character(stars$CONSTELLATION[ind])
  star$colour <- as.character(stars$COLOUR[ind])
  star$app.mag <- stars$VMAG[ind]
  star$ra <- stars$RA[ind]
  star$dec <- stars$DEC[ind]
  star$proper.motion <- c(stars$PM_RA[ind],stars$PM_DEC[ind])
  star$epoch <- "J2000.0"
  class(star) <- "skyscapeR.star"

  # Time shift
  if (!missing(year)) {
    # precession
    trash <- capture.output(coor <- astrolibR::precess(star$ra,star$dec, 2000, year), file=NULL)

    # proper motion
    adj.pm <- star$proper.motion*(year-2000)/1000
    coor$ra <- coor$ra + astrolibR::ten(0,0,adj.pm[1])*sign(adj.pm[1])
    coor$dec <- coor$dec + astrolibR::ten(0,0,adj.pm[2])*sign(adj.pm[2])

    # nutation
    adj.nu <- astrolibR::co_nutate(astrolibR::jdcnv(year,1,1,0), star$ra, star$dec)
    coor$ra <- coor$ra + astrolibR::ten(0,0,adj.nu$d_ra)*sign(adj.nu$d_ra)
    coor$dec <- coor$dec + astrolibR::ten(0,0,adj.nu$d_dec)*sign(adj.nu$d_dec)

    # aberration
    adj.ab <- astrolibR::co_aberration(astrolibR::jdcnv(year,1,1,0), star$ra, star$dec, pracma::deg2rad(obliquity(year)))
    coor$ra <- coor$ra + astrolibR::ten(0,0,adj.ab$d_ra)*sign(adj.ab$d_ra)
    coor$dec <- coor$dec + astrolibR::ten(0,0,adj.ab$d_dec)*sign(adj.ab$d_dec)

    # output
    star$ra <- coor$ra
    star$dec <- coor$dec
    star$epoch <- as.character(year)
  }

  return(star)
}

#' Calculate the seasons and phases of a star
#'
#' This function calculates the seasons (Rising, Setting, etc.)
#' and phases (Arising and Lying Hidden, Curtailed Passage) of a
#' star for a given location and epoch. This functions uses the
#' \emph{arcus visionis} approximation of Purrington (1988). For
#' the nomenclature used, and description of star phases, see Brady (2015).
#' @param star Either the star name or a \emph{skyscapeR.star} object.
#' @param year The year of interest.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude and longitude of location, in this order.
#' @param alt.hor (Optional) The altitude of the horizon to consider.
#' Defaults to zero degrees.
#' @param alt.rs (Optional) The maximum altitude of a star's first or last
#' visibility for it to still be considered to be as rising or setting.
#' Defaults to ten degrees.
#' @param res (Optional) Resolution of calculation. The smaller this
#' figure the slower the computation. Defaults to 24/3600 = 1 sec.
#' @param ncores (Optional) Number of processing cores to use for parallelisation.
#' Defaults to the number of available cores minus 1.
#' @export
#' @import parallel numDeriv
#' @seealso \code{\link{plotPhases}}
#' @references Purrington, Robert D. (1988) Heliacal Rising and Setting:
#' Quantitative Aspects, \emph{Journal for the History of Astronomy
#' (Archaeoastronomy Supplement 12)} 19, S72-S84. Available online at
#' [SAO/NASA ADS Astronomy Abstract Service](http://adsabs.harvard.edu/abs/1988JHAS...19...72P)
#' @references Brady, Bernadette (2015) Star Phases: the Naked-eye Astronomy of the Old Kingdom
#' Pyramid Texts. In F Silva and N Campion (eds) \emph{Skyscapes: The Role and Importance of
#' the Sky in Archaeology}. Oxford: Oxbow Books, pp. 76-86.
#' @examples
#' \dontrun{
#' ss1 <- star.phases('Aldebaran',-4000, c(35,-8))
#'
#' # One can then look at the star's phase:
#' ss1$phase
#'
#' # Date range of seasons:
#' ss1$seasons
#'
#' # Date range of events/phases:
#' ss1$events
#'
#' # And plot them:
#' plotPhases(ss1)
#'
#' # You can play with the parameters and see how predictions change:
#' ss1 <- star.phases('Aldebaran',-4000, c(35,-8), alt.hor=2, alt.rs=5)
#' plotPhases(ss1)
#' }
star.phases <- function(star, year, loc, alt.hor = 0, alt.rs = 10, res = 24/3600, ncores) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
    lon <- loc$georef[2]
  } else {
    lat <- loc[1]
    lon <- loc[2]
  }

  # load star data
  if (class(star)!='skyscapeR.star') { star <- star(star, year) }
  arcus_visionis <- 2.1*star$app.mag + 10

  # init parallel cluster
  if (missing(ncores)) { ncores <-  parallel::detectCores()-1 }
  cl <- parallel::makeCluster(ncores, type = "PSOCK")
  parallel::clusterEvalQ(cl, library(skyscapeR))
  parallel::clusterExport(cl, list("lat", "lon"), envir=environment())
  message(paste0('Running calculations on ', ncores, ' processing cores. This may take a while...'))

  # calculations
  jd0 <- astrolibR::juldate(c(year,1,1,0)) + 2400000
  jd.t <- seq(jd0, jd0+365, res)
  tx <- jd.t - jd0+1

  # Sun
  Sun <- astrolibR::sunpos(jd.t)
  X <- cbind(Sun$ra,Sun$dec); X <- cbind(X,jd.t)
  Sun.alt <- parallel::parApply(cl, X, 1, eq.to.Hor, lat=lat, lon=lon)

  # Star
  X <- cbind(rep(star$ra, NROW(jd.t)), rep(star$dec, NROW(jd.t))); X <- cbind(X,jd.t)
  Star.alt <- parallel::parApply(cl, X, 1, eq.to.Hor, lat=lat, lon=lon)

  parallel::stopCluster(cl)

  day <- which(Sun.alt >= alt.hor)
  civ <- which(Sun.alt >= -6 & Sun.alt < alt.hor)
  naut <- which(Sun.alt >= -12 & Sun.alt < -6)
  astr <- which(Sun.alt >= -18 & Sun.alt < -12)
  night <- which(Sun.alt < -18)

  # Seasons
  nightime <- c(civ,naut,astr,night)
  ind <- which( Star.alt[nightime] >= alt.hor & Star.alt[nightime] <= alt.rs & Star.alt[nightime]-Sun.alt[nightime] >= arcus_visionis)
  ss.ind <- nightime[ind]

  ff <- splinefun(tx, Star.alt)
  dd <- numDeriv::grad(ff, tx)

  RISE <- tx[ss.ind[which(dd[ss.ind]>0)]]; RISE <- unique(trunc(RISE,0)); RISE <- sort(RISE)
  SET <- tx[ss.ind[which(dd[ss.ind]<0)]]; SET <- unique(trunc(SET,0)); SET <- sort(SET)
  RISESET <- intersect(sort(RISE), sort(SET)) # days of rise/set season
  RISE <- setdiff(RISE, RISESET) # days of rise season
  SET <- setdiff(SET, RISESET) # days of set season
  days <- seq(1,365,1); NOHOR <- setdiff(days, sort(c(RISE,SET,RISESET)))

  # Events
  ff.sun <- splinefun(tx, Sun.alt)
  dd.sun <- numDeriv::grad(ff.sun, tx)

  twilight <- c(civ,naut,astr)
  ind <- which( Star.alt[twilight] >= alt.hor & Star.alt[twilight] <= alt.rs & (Star.alt[twilight]-Sun.alt[twilight]) >= arcus_visionis)
  ss.ind <- twilight[ind]
  MR <- tx[ss.ind[which(dd[ss.ind]>0 & dd.sun[ss.ind]>0) ]]; MR <- unique(trunc(MR,0)); MR <- sort(MR) # days of morning rising
  MS <- tx[ss.ind[which(dd[ss.ind]<0 & dd.sun[ss.ind]>0) ]]; MS <- unique(trunc(MS,0)); MS <- sort(MS) # days of morning setting
  ER <- tx[ss.ind[which(dd[ss.ind]>0 & dd.sun[ss.ind]<0) ]]; ER <- unique(trunc(ER,0)); ER <- sort(ER) # days of evening rising
  ES <- tx[ss.ind[which(dd[ss.ind]<0 & dd.sun[ss.ind]<0) ]]; ES <- unique(trunc(ES,0)); ES <- sort(ES) # days of evening setting

  # convert dates to proleptic gregorian ##### TO DO CLEANUP
  RISESET <- sort(tWS(RISESET, year))
  RISE <- sort(tWS(RISE, year))
  SET <- sort(tWS(SET, year))
  NOHOR <- sort(tWS(NOHOR, year))
  MR <- sort(tWS(MR, year))
  MS <- sort(tWS(MS, year))
  ER <- sort(tWS(ER, year))
  ES <- sort(tWS(ES, year))

  # stellar event detection
  res.events <- list(); phase <- c()
  bb <- split(NOHOR,cumsum(c(1,diff(NOHOR)>3)))
  for (j in 1:NROW(bb)) {
    if (length(which(RISE==min(bb[[j]])-4)) > 0) {
      ## CP-phase
      phase <- c(phase,"Curtailed Passage")
      res.events$AR <- ER
      res.events$HS <- MS
    }
    if (length(which(SET==min(bb[[j]])-4)) > 0) {
      ## ALH-phase
      phase <- c(phase,"Arising and Lying Hidden")
      res.events$AS <- ES
      res.events$HR <- MR
    }
  }

  # prettify output of seasons
  ind.bb <- split(RISE,cumsum(c(1,diff(RISE)>3))); if (NROW(ind.bb)==2) { RISE <- c(ind.bb$`2`,ind.bb$`1`) }
  sRISE <- paste0(dd.to.DD(RISE[1], char=T), ' - ', dd.to.DD(tail(RISE,1), char=T) )

  ind.bb <- split(SET,cumsum(c(1,diff(SET)>3))); if (NROW(ind.bb)==2) { SET <- c(ind.bb$`2`,ind.bb$`1`) }
  sSET <- paste0(dd.to.DD(SET[1], char=T), ' - ', dd.to.DD(tail(SET,1), char=T) )

  if (NROW(phase)==1) {
    ind.bb <- split(RISESET,cumsum(c(1,diff(RISESET)>3))); if (NROW(ind.bb)==2) { RISESET <- c(ind.bb$`2`,ind.bb$`1`) }
    sRS <- paste0(dd.to.DD(RISESET[1], char=T), ' - ', dd.to.DD(tail(RISESET,1), char=T) )

    ind.bb <- split(NOHOR,cumsum(c(1,diff(NOHOR)>3))); if (NROW(ind.bb)==2) { NOHOR <- c(ind.bb$`2`,ind.bb$`1`) }
    sNO <- paste0(dd.to.DD(NOHOR[1], char=T), ' - ', dd.to.DD(tail(NOHOR,1), char=T) )
  } else {
    # dual-phase star
    ind.bb <- split(NOHOR,cumsum(c(1,diff(NOHOR)>3))); NOHOR1 <- ind.bb$`1`; NOHOR2 <- ind.bb$`2`
    sNO1 <- paste0(dd.to.DD(NOHOR1[1], char=T), ' - ', dd.to.DD(tail(NOHOR1,1), char=T) )
    sNO2 <- paste0(dd.to.DD(NOHOR2[1], char=T), ' - ', dd.to.DD(tail(NOHOR2,1), char=T) )
  }

  res.seasons <- list()
  if(NROW(phase)==1 & phase=="Arising and Lying Hidden") {
    seasons <- rbind(sRISE, sRS, sSET, sNO); rownames(seasons) <- c('Rising','Rising and Setting','Setting', phase)
    res.seasons$R <- RISE; res.seasons$RS <- RISESET; res.seasons$S <- SET; res.seasons$ALH <- NOHOR
  } else if(NROW(phase)==1 & phase=="Curtailed Passage") {
    seasons <- rbind(sSET, sRS, sRISE, sNO); rownames(seasons) <- c('Setting','Rising and Setting','Rising', phase)
    res.seasons$S <- SET; res.seasons$RS <- RISESET; res.seasons$R <- RISE; res.seasons$CP <- NOHOR
  } else if(NROW(phase)==2) {
    if(phase[1]=="Arising and Lying Hidden") {
      seasons <- rbind(sSET, sNO1, sRISE, sNO2); rownames(seasons) <- c('Setting',phase[1],'Rising',phase[2])
      res.seasons$S <- SET; res.seasons$ALH <- NOHOR1; res.seasons$R <- RISE; res.seasons$CP <- NOHOR2
    } else {
      seasons <- rbind(sRISE, sNO1, sSET, sNO2); rownames(seasons) <- c('Rising',phase[1],'Setting',phase[2])
      res.seasons$R <- RISE; res.seasons$CP <- NOHOR1; res.seasons$S <- SET; res.seasons$ALH <- NOHOR2
    }
  }

  # prettify output of events
  stellar.events <- data.frame(short=c('HR','HS','AR','AS'), long=c('Heliacal Rising', 'Heliacal Setting','Achronycal Rising', 'Achronycal Setting'))
  events <- c(); ss <- c()
  for (i in 1:NROW(res.events)) {
    ind.bb <- split(res.events[[i]],cumsum(c(1,diff(res.events[[i]])>3))); if (NROW(ind.bb)==2) { aux <- c(ind.bb$`2`,ind.bb$`1`) } else { aux <- res.events[[i]] }
    events <- rbind(events, paste0(dd.to.DD(aux[1], char=T), ' - ', dd.to.DD(tail(aux,1), char=T) ))
    ss[i] <- as.character(stellar.events[which(stellar.events$short == names(res.events[i])),2])
  }
  rownames(events) <- ss

  # results
  out <- c()
  out$star <- star
  out$year <- year
  out$loc <- loc
  out$phase <- phase
  out$seasons <- seasons
  out$events <- events
  out$alt$hor <- alt.hor
  out$alt$rs <- alt.rs
  out$raw$seasons <- list(RS <- RISESET, R <- RISE, S <- SET, NH <- NOHOR )
  out$raw$events <- res.events
  class(out) <- 'skyscapeR.starphase'
  return(out)
}
