#' Computes declination curvigram
#'
#' This function computes the curvigram of declinations,
#' using provided measurement uncertainty and a Gaussian
#' kernel.
#' @param dec Array of declination values
#' @param sd Either an array of declination values (same length as dec),
#' or a single value to be applied to all declinations. Defaults to 2.
#' @param norm Boolean specifying whether the resulting curvigram
#' should be normalized to unity. Defaults to TRUE.
#' @export
#' @examples
#' # Curvigram of Ruggles' Recumbent Stone Circle data:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, 2)
#' plot(curv$dec, curv$density, type='l', lwd=2, col='blue', xlab='DEC', ylab='')
curvigram <- function(dec, sd = 2, norm = T) {
  if (length(sd)==1) { sd <- rep(sd, NROW(dec)) }

  xx <- seq(-90,90,by=0.01)
  spd <- array(0, NROW(xx))
  for (i in 1:NROW(dec)) {
    spd <- spd + dnorm(xx, dec[i], sd[i])
  }

  if (norm) { spd <- spd/max(spd) }

  result <- c()
  result$dec <- xx
  result$density <- spd
  class(result) <- "skyscapeR.curv"
  return(result)
}
