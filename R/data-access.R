#' Create \emph{skyscapeR.star} object
#'
#' This function retrieves information for a given star
#' and saves it in the skyscapeR.star format ready to be
#' used by other skyscapeR package function.
#' @param string This can be either: 1) the common name for the star,
#'  2) the Hipparchos catalogue number (if string begins with HIP),
#'  3) the Bright Star Catalogue number (if string begins with HR),
#'  4) the Messier catalogue number (if string begins with M), or
#'  5) the Bayer designation (if string has exactly 5 characters,
#'   the second of which is a space).
#' @export
#' @examples
#' # Retrieve data for Aldebaran:
#' Aldeb <- star('Aldebaran')
#'
#' # Retrieve data for the Pleiades (M45):
#' P1 <- star('Pleiades')
#' P2 <- star('M45')
star <- function(string) {
  data(stars, envir = environment())
  ind <- which(stars$NAME == string)
  if (substr(string,1,3)=="HIP") { ind <- which(stars$HIP.ID == string)}
  if (substr(string,1,2)=="HR") { ind <- which(stars$HR.ID == string)}
  if (substr(string,1,1)=="M" & nchar(string)==3) { ind <- which(stars$MESSIER.ID == string)}
  if (nchar(string)==5 & substr(string,2,2)==" ") { ind <- which(stars$BAYER.DESIGNATION == string)}

  if (length(ind) == 0) {
    stop('No star with that name or designation has been found.')
  } else {
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
    return(star)
  }
}
