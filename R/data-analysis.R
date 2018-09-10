#' Computes declination histograms
#'
#' This function computes the histogram of declinations,
#' using provided measurement uncertainty and a bin size.
#' @param dec Array of declination values.
#' @param unc (Optional) Uncertainty in declination values,
#' which will be used as bin size. Must be a single value.
#' @param norm (Optional) Boolean specifying whether the resulting curvigram
#' should be normalized to unity. Defaults to \emph{FALSE}.
#' @param cut (Optional) Number of uncertainties beyond the extremes
#' of the data at which to trim the curvigram. Defaults to 4. See \code{\link{density}}.
#' @param range (Optional) As an alternative to \emph{cut} you can
#' stipulate the range of declination values to output as an array of two values.
#' See \emph{from, to} in \code{\link{density}}.
#' @seealso \code{\link{curvigram}}
#' @export
#' @import stats
#' @examples
#' # Histogram of Ruggles' Recumbent Stone Circle data:
#' data(RugglesRSC)
#' hist <- histogram(RugglesRSC$Dec, 2)
#' plot(hist)
histogram <- function(dec, unc = 2, norm = F, cut = 4, range) {
  if (missing(range)) {
    range <- c(min(dec) - cut * max(unc), max(dec) + cut *
                 max(unc))
  }

  if (length(unc) > 1) {
    stop('Cannot use function histogram() with multiple measurement uncertainties. Suggest curvigram() instead.')
  }

  breaks <- seq(-90, 90, unc)

  dens <- hist(dec, breaks = breaks, plot=F)
  xx <- breaks
  if (norm) { hh <- dens$density } else { hh <-  dens$counts}

  ind <- which(breaks[1:(NROW(breaks)-1)] >= range[1] & breaks[1:(NROW(breaks)-1)] <= range[2])
  xx <- c(xx[ind], xx[tail(ind,1)+1])
  hh <- hh[ind]

  result <- c()
  result$metadata$unc <- unc
  result$metadata$norm <- norm
  result$metadata$range <- range
  result$data$dec <- xx
  result$data$density <- hh
  class(result) <- "skyscapeR.hist"
  return(result)
}


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
    xx <- seq(-90,90, length.out = 0.01)
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
