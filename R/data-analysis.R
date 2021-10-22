

#' Find celestial targets within declination and time ranges
#'
#' @param decrange Range of declination to consider.
#' @param timerange Temporal range to consider
#' @param max.mag (Optional) Maximum magnitude of stars to consider. Defaults to 2.5
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude, longitude and elevation of location, in this order. Defaults
#' to FALSE, thus checking only geocentric declination.
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
findTargets <- function(decrange, timerange, max.mag=2.5, loc=FALSE, calendar=skyscapeR.env$calendar) {
  targets <- c()

  ## solar
  targets$solar <- data.frame(name='Test', min.dec=12, max.dec=12, date1='1 Jan', date2='1 Feb', stringsAsFactors=F)
  targets$solar[1,] <- c('june solstice', sort(jS(timerange, loc, verbose=FALSE)), '21 Jun', NA)
  targets$solar[2,] <- c('december solstice', sort(dS(timerange, loc, verbose=FALSE)), '21 Dec', NA)
  targets$solar[3,] <- c('astronomical equinox', 0, NA, '21 Mar', '21 Sep')
  if (class(loc)[1] == 'skyscapeR.horizon') {
    aux <- solar.date(sort(as.numeric(spatial.equinox(loc))), mean(timerange), calendar, verbose=F)
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

  aux1 <- try(solar.date(seq(decrange[1],decrange[2],0.1), min(timerange), calendar, verbose=F), silent=T)
  aux2 <- try(solar.date(seq(decrange[1],decrange[2],0.1), max(timerange), calendar, verbose=F), silent=T)
  if (class(aux1)[1] != 'try-error' & class(aux2)[1] == 'try-error') { targets$solar[7,] <- c('sunrise/set', decrange[1], decrange[2], paste(aux1[2,1],'-',aux1[2,NCOL(aux1)]), paste(aux1[3,NCOL(aux1)],'-',aux1[3,1])) }
  if (class(aux1)[1] == 'try-error' & class(aux2)[1] != 'try-error') { targets$solar[7,] <- c('sunrise/set', decrange[1], decrange[2], paste(aux2[2,1],'-',aux2[2,NCOL(aux2)]), paste(aux2[3,NCOL(aux2)],'-',aux2[3,1])) }
  if (class(aux1)[1] != 'try-error' & class(aux2)[1] != 'try-error') { targets$solar[7,] <- c('sunrise/set', decrange[1], decrange[2], paste(aux1[2,1],'-',aux2[2,NCOL(aux2)]), paste(aux1[3,NCOL(aux1)],'-',aux2[3,1])) }

  targets$solar[,2] <- as.numeric(targets$solar[,2])
  targets$solar[,3] <- as.numeric(targets$solar[,3])

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
  targets$lunar[1,] <- c('southern major lunar extreme', sort(sMjLX(timerange, loc, verbose=FALSE)))
  targets$lunar[2,] <- c('southern minor lunar extreme', sort(smnLX(timerange, loc, verbose=FALSE)))
  targets$lunar[3,] <- c('northern minor lunar extreme', sort(nmnLX(timerange, loc, verbose=FALSE)))
  targets$lunar[4,] <- c('northern major lunar extreme', sort(nMjLX(timerange, loc, verbose=FALSE)))
  ## TODO add EFM distribution

  targets$lunar[,2] <- as.numeric(targets$lunar[,2])
  targets$lunar[,3] <- as.numeric(targets$lunar[,3])

  for (i in 1:NROW(targets$lunar)) {
    if ((targets$lunar$min.dec[i] < min(decrange) & targets$lunar$max.dec[i] < min(decrange)) | (targets$lunar$min.dec[i] > max(decrange) & targets$lunar$max.dec[i] > max(decrange))) {
      targets$lunar[i,] <- c(NA,NA,NA)
    }
  }
  targets$lunar <- targets$lunar[-which(is.na(targets$lunar[,1])),]
  rownames(targets$lunar) <- NULL
  if (NROW(targets$lunar) ==0) { targets$lunar <- c() }

  ## stellar
  ss <- skyscapeR.env$sefstars
  ind <- which(ss$mag.V>max.mag); ss <- ss[-ind,] # magnitude filter

  targets$stellar <- data.frame(name='Test', Bayer='Test', vmag=1.2, min.dec=12, max.dec=12, stringsAsFactors=F)
  for (i in 1:NROW(ss)) {
    if (nchar(ss$traditional.name[i])>0) { aux <- ss$traditional.name[i] } else { aux <- ss$nomenclature.name[i] }
    targets$stellar[i,] <- c(ss$traditional.name[i], ss$nomenclature.name[i], as.character(ss$mag.V[i]), minmaxdec(aux, timerange[1], timerange[2]))
  }

  targets$stellar[,4] <- as.numeric(targets$stellar[,4])
  targets$stellar[,5] <- as.numeric(targets$stellar[,5])

  for (i in 1:NROW(targets$stellar)) {
    if ((targets$stellar$min.dec[i] < min(decrange) & targets$stellar$max.dec[i] < min(decrange)) | (targets$stellar$min.dec[i] > max(decrange) & targets$stellar$max.dec[i] > max(decrange))) {
      targets$stellar[i,] <- c(NA,NA,NA,NA,NA)
    }
  }
  targets$stellar <- targets$stellar[-which(is.na(targets$stellar[,1])),]
  rownames(targets$stellar) <- NULL
  if (NROW(targets$stellar) ==0) { targets$stellar <- c() }

  return(targets)
}
