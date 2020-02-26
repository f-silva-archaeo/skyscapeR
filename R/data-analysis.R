#' Computes declination curvigram
#'
#' This function computes the curvigram of declinations,
#' using provided measurement uncertainty and a Gaussian
#' kernel, i.e. using the method of Silva (2017). When all
#' measurements have the same associated
#' uncertainty this function wraps \code{\link{density}} and
#' accepts the same input for \emph{bw} in \emph{unc}.
#' @param dec Array of declination values.
#' @param unc (Optional) Uncertainty in declination values.
#' Either a single value or string to be applied
#' to all measurements (see \code{\link{bw.nrd}}), or an
#' array of values of the same length as \emph{dec}. Defaults
#' to 2 degrees.
#' @param norm (Optional) Boolean specifying whether the resulting curvigram
#' should be normalized to unity. Defaults to \emph{FALSE}.
#' @param cut (Optional) Number of uncertainties beyond the extremes
#' of the data at which to trim the curvigram. Defaults to 4. See \code{\link{density}}.
#' @param range (Optional) As an alternative to \emph{cut} you can
#' stipulate the range of declination values to output as an array of two values.
#' See \emph{from, to} in \code{\link{density}}.
#' @seealso \code{\link{density}}
#' @export
#' @import stats
#' @references Silva, Fabio (2017) Inferring Alignments I: Exploring the Accuracy
#' and Precision of Two Statistical Approaches, \emph{Journal of Skyscape Archaeology}
#' 3(1), 93-111. DOI: 10.1558/jsa.31958
#' @examples
#' # Curvigram of Ruggles' Recumbent Stone Circle data:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, 2)
#' plot(curv)
curvigram <- function(dec, unc = 2, norm = F, cut = 4, range) {
  if (missing(range)) { range <- c(min(dec) - cut*max(unc), max(dec) + cut*max(unc)) }

  if (length(unc)==1) {
    dens <- density(dec, bw=unc, from=range[1], to=range[2], n=1024)

    xx <- dens$x
    spd <- dens$y
    sd <- dens$bw

  } else {
    xx <- seq(-90,90, by = 0.01)
    spd <- array(0, NROW(xx))
    for (i in 1:NROW(dec)) {
      spd <- spd + dnorm(xx, dec[i], unc[i])
    }
    ind <- which(xx >= range[1] & xx <= range[2])
    xx <- xx[ind]
    spd <- spd[ind]
  }

  if (norm) { spd <- spd/max(spd) }

  result <- c()
  result$metadata$unc <- unc
  result$metadata$norm <- norm
  result$metadata$range <- range
  result$data$dec <- xx
  result$data$density <- spd
  class(result) <- "skyscapeR.curv"
  return(result)
}



#' Function to find celestial targets within declination and time ranges
#'
#' @param decrange Range of declination to consider.
#' @param timerange Temporal range to consider
#' @param nstars (Optional) NUmber of stars to consider i ndescending order of apparent
#' magnitude. Defaults to 100.
#' @param loc Location, either a \emph{skyscapeR.object} or a vector
#' containing the latitude, longitude and elevation of location, in this order. Defaults
#' to FALSE, thus checking only geocentric declinations.
#' @export
#' @examples
#' findTargets(c(-35,-32), c(-2500,-1750))
findTargets <- function(decrange, timerange, nstars=100, loc=FALSE) {
  targets <- data.frame(name='Test', constellation='Test', vmag=1.2, dec1=12, dec2=12, stringsAsFactors=F)

  ## solar extremes
  targets[1,] <- c('june solstice', NA, NA, jS(timerange, loc, verbose=FALSE))
  targets[2,] <- c('december solstice', NA, NA, dS(timerange, loc, verbose=FALSE))

  ## lunar extremes
  targets[3,] <- c('southern major lunar extreme', NA, NA, sMjLX(timerange, loc, verbose=FALSE))
  targets[4,] <- c('southern minor lunar extreme', NA, NA, smnLX(timerange, loc, verbose=FALSE))
  targets[5,] <- c('northern minor lunar extreme', NA, NA, nmnLX(timerange, loc, verbose=FALSE))
  targets[6,] <- c('northern major lunar extreme', NA, NA, nMjLX(timerange, loc, verbose=FALSE))

  ## stars
  #identify top nstars in order of magnitude
  fpath <- system.file("ephemeris", "sefstars.txt", package="swephR")
  cnames <- c('traditional name','nomenclature name','equinox','RA hr','RA min', 'RA sec', 'Dec deg', 'Dec min', 'Dec sec', 'pm RA', 'pm Dec', 'rad vel', 'ann plx', 'mag V', 'DM zone', 'DM number')
  sefstars <- read.csv(fpath, as.is=T, header=F, comment.char='#', col.names=cnames, strip.white=T)
  sefstars <- sefstars[-which(sefstars$mag.V==0),]
  ind <- sort(sefstars$mag.V, index.return=T)$ix
  sefstars <- sefstars[ind,]
  df <- sefstars[,-1]
  ind <- which(duplicated(df))
  ss <- sefstars[-ind,]; rm(sefstars, df)

  for (i in 1:nstars) {
    targets[6+i,] <- c(as.character(ss$traditional.name[i]), as.character(ss$nomenclature.name[i]), as.character(ss$mag.V[i]), minmaxdec(as.character(ss$traditional.name[i]), timerange[1], timerange[2]))
  }

  targets[,4] <- as.numeric(targets[,4])
  targets[,5] <- as.numeric(targets[,5])

  ## cleanup
  for (i in 1:NROW(targets)) {
    if ((targets$dec1[i] < min(decrange) & targets$dec2[i] < min(decrange)) | (targets$dec1[i] > max(decrange) & targets$dec2[i] > max(decrange))) {
      targets[i,] <- c(NA,NA,NA,NA,NA)
    }
  }
  targets <- targets[-which(is.na(targets[,1])),]
  return(targets)
}
