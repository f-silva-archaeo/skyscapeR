

#' Find celestial targets within declination and time ranges
#'
#' @param decrange Range of declination to consider.
#' @param timerange Temporal range to consider
#' @param max.mag (Optional) Maximum magnitude of stars to consider. Defaults to 2.5
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude, longitude and elevation of location, in this order. Defaults
#' to FALSE, thus checking only geocentric declination.
#' @param dates (Optional) Whether to display time period when identified celestial targets
#' occur within declination range. Default is FALSE.
#' @param calendar (Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @export
#' @examples
#' \dontrun{
#' findTargets(c(-25,-17.5), c(-2500,-1750))
#'
#' # if a location is given then the zenith and anti-zenith sun will also be looked at:
#' findTargets(c(3,12), c(-2500,-1750), loc=c(8.6, 7.3, 200))
#'
#' # if a horizon profile is given then the spatial equinox will also be looked at:
#' hor <- downloadHWT('J657KVEV')
#' findTargets(c(-7,2), c(-2500,-1750), loc=hor)
#' }
findTargets <- function(decrange, timerange, max.mag=2.5, loc=FALSE, dates=FALSE, calendar=skyscapeR.env$calendar) {
  targets <- c()

  ## solar
  targets$solar <- data.frame(name='Test', min.dec=12, max.dec=12, date1='1 Jan', date2='1 Feb', stringsAsFactors=F)
  targets$solar[1,] <- c('june solstice', range(jS(seq(min(timerange), max(timerange),1), loc, verbose=FALSE)), '21 Jun', NA)
  targets$solar[2,] <- c('december solstice', range(dS(seq(min(timerange), max(timerange),1), loc, verbose=FALSE)), '21 Dec', NA)
  targets$solar[3,] <- c('astronomical equinox', 0, NA, '21 Mar', '21 Sep')
  if (class(loc)[1] == 'skyscapeR.horizon') {
    aux <- solar.date(sort(as.numeric(spatial.equinox(loc)$declination)), mean(timerange), calendar, verbose=F)
    targets$solar[4,] <- c('spatial equinox', as.numeric(aux[1,]), paste(aux[2,1],'/',aux[2,2]), paste(aux[3,2],'/',aux[3,1]))
  }
  if(class(loc)[1] != 'logical') {
    if (!is.null(zenith(loc))) {
      aux <- solar.date(sort(as.numeric(zenith(loc))), mean(timerange), calendar, verbose=F)
      targets$solar[5,] <- c('zenith sun', as.numeric(aux[1,]), NA, aux[2,], aux[3,])
      aux <- solar.date(sort(as.numeric(antizenith(loc))), mean(timerange), calendar, verbose=F)
      targets$solar[6,] <- c('antizenith sun', as.numeric(aux[1,]), NA, aux[2,], aux[3,])
    }
  }

  aux1 <- try(solar.date(seq(min(decrange),max(decrange),0.01), min(timerange), calendar, verbose=F), silent=T)
  aux2 <- try(solar.date(seq(min(decrange),max(decrange),0.01), max(timerange), calendar, verbose=F), silent=T)

  ttt <- c()
  if (class(aux1)[1] != 'try-error') { ttt <- c(ttt, aux1[2,],aux1[3,]) }
  if (class(aux2)[1] != 'try-error') { ttt <- c(ttt, aux2[2,],aux2[3,]) }
  if (!is.null(ttt)) {
    dd <- as.numeric(substr(ttt,1,2))
    mm <- substr(ttt,3,6); mm <- gsub(" ","",mm)
    mm <- match(mm, month.abb)

    months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    summon <- c(0,cumsum(months))
    xx <- sort(unique(summon[mm]+dd))

    if (sum(xx==365) & sum(xx==1)) {
      xx <- c(xx,xx+365)
      xx <- xx[-which(xx<365/2)]
      xx <- xx[-which(xx>365+365/2)]
    }
    ind <- which(diff(xx)>2)
    dts <- c(xx[1], xx[ind], xx[ind+1], xx[length(xx)])
    dts[dts>365] <- dts[dts>365]-365
    dts <- dd.to.DD(dts)
    dts <- rbind(sprintf("%02d",dts[1,]), sprintf("%02d",dts[2,]))

    targets$solar[7,] <- c('sunrise/set', decrange[1], decrange[2], paste(long.date(paste(dts[,1], collapse = '-')),'-',long.date(paste(dts[,2], collapse = '-'))), NA)
    if (NCOL(dts)>2) { targets$solar[7,5] <- paste(long.date(paste(dts[,3], collapse = '-')),'-',long.date(paste(dts[,4], collapse = '-'))) }
  }

  targets$solar[,2] <- as.numeric(targets$solar[,2])
  targets$solar[,3] <- as.numeric(targets$solar[,3])

  targets$solar$min.dec <- round(targets$solar$min.dec,2)
  targets$solar$max.dec <- round(targets$solar$max.dec,2)

  targets$solar <- targets$solar[-which(is.na(targets$solar[,1])),]
  for (i in 1:NROW(targets$solar)) {
    if ((min(decrange) > max(targets$solar[i,c(2,3)], na.rm=T))  | (max(decrange) < min(targets$solar[i,c(2,3)], na.rm=T))) {
      targets$solar[i,] <- c(NA,NA,NA,NA,NA)
    }
  }
  targets$solar <- targets$solar[-which(is.na(targets$solar[,1])),]
  rownames(targets$solar) <- NULL
  if (NROW(targets$solar) ==0) { targets$solar <- c() }



  ## lunar
  targets$lunar <- data.frame(name='Test', min.dec=12, max.dec=12, stringsAsFactors=F)

  time <- seq(min(timerange), max(timerange),1)
  aux <- sMjLX(time, loc)


  targets$lunar[1,] <- c('southern major lunar extreme', range(sMjLX(seq(min(timerange), max(timerange),1), loc, verbose=FALSE)))
  targets$lunar[2,] <- c('southern minor lunar extreme', range(smnLX(seq(min(timerange), max(timerange),1), loc, verbose=FALSE)))
  targets$lunar[3,] <- c('northern minor lunar extreme', range(nmnLX(seq(min(timerange), max(timerange),1), loc, verbose=FALSE)))
  targets$lunar[4,] <- c('northern major lunar extreme', range(nMjLX(seq(min(timerange), max(timerange),1), loc, verbose=FALSE)))
  ## TODO add EFM distribution

  targets$lunar[,2] <- as.numeric(targets$lunar[,2])
  targets$lunar[,3] <- as.numeric(targets$lunar[,3])

  targets$lunar$min.dec <- round(targets$lunar$min.dec,2)
  targets$lunar$max.dec <- round(targets$lunar$max.dec,2)

  for (i in 1:NROW(targets$lunar)) {
    if ((targets$lunar$min.dec[i] < min(decrange) & targets$lunar$max.dec[i] < min(decrange)) | (targets$lunar$min.dec[i] > max(decrange) & targets$lunar$max.dec[i] > max(decrange))) {
      targets$lunar[i,] <- c(NA,NA,NA)
    }
  }
  targets$lunar <- targets$lunar[-which(is.na(targets$lunar[,1])),]
  rownames(targets$lunar) <- NULL
  if (NROW(targets$lunar) ==0) { targets$lunar <- c() }



  ## stellar
  ss <- read.csv(paste0(system.file('ephemeris',package = 'swephR'),'/sefstars.txt'), header=F)
  colnames(ss) <- c('name', 'identifier', 'ICRS', 'RA.1', 'RA.2', 'RA.3', 'Dec.1', 'Dec.2', 'Dec.3', 'pm.1', 'pm.2', 'radvel', 'plx', 'magV')
  ind <- which(ss$magV>max.mag); ss <- ss[-ind,] # magnitude filter

  targets$stellar <- data.frame(name='Test', identifier='Test', magV=1.2, min.dec=12, max.dec=12, daterange=NA, stringsAsFactors=F)
  for (i in 1:NROW(ss)) {
    targets$stellar[i,] <- c(ss$name[i], ss$identifier[i], as.character(ss$magV[i]), minmaxdec(ss$identifier[i], timerange[1], timerange[2]), dating(ss$identifier[i], decrange, decrange[2], timerange[1], timerange[2]))
  }

  targets$stellar[,4] <- as.numeric(targets$stellar[,4])
  targets$stellar[,5] <- as.numeric(targets$stellar[,5])

  targets$stellar$min.dec <- round(targets$stellar$min.dec,2)
  targets$stellar$max.dec <- round(targets$stellar$max.dec,2)

  for (i in 1:NROW(targets$stellar)) {
    if ((targets$stellar$min.dec[i] < min(decrange) & targets$stellar$max.dec[i] < min(decrange)) | (targets$stellar$min.dec[i] > max(decrange) & targets$stellar$max.dec[i] > max(decrange))) {
      targets$stellar[i,] <- c(NA,NA,NA,NA,NA,NA)
    }
  }

  targets$stellar <- targets$stellar[-which(is.na(targets$stellar[,2])),]
  for (i in 1:NROW(targets$stellar)) {
    aux <- star.names$Western[which(star.names$identifier == targets$stellar[i,2])]
    if (length(aux)>0) { targets$stellar[i,1] <- aux } else { targets$stellar[i,1] <- '' }
  }

  rownames(targets$stellar) <- NULL
  if (NROW(targets$stellar) ==0) { targets$stellar <- c() }

  ## for "security" reasons
  if (dates==FALSE) targets$stellar$daterange <- NULL

  return(targets)
}
