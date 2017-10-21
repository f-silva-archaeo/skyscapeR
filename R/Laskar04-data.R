#' Astronomical elements from Laskar et al. (2004)
#'
#' Astronomical elements (longitude of perihelion,
#' obliquity and eccentricity) by step of 1ka, from
#' -51M BP to present (la04past) and from present
#' to + 21M BP (la04future).
#'
#' @docType data
#' @usage data(Laskar04)
#' @format A list of two data frames with 4 variables:
#'  \describe{
#'    \item{time}{Time in years before or after J1950.0}
#'    \item{ecc}{Eccentricity}
#'    \item{eps}{Obliquity}
#'    \item{varpi}{Longitude of Perihelion}
#'    }
#' @keywords datasets
#' @references Laskar, J. et al. (2004), A long-term numerical
#' solution for the insolation quantities of the Earth, \emph{Astron.
#'  Astroph.}, 428, 261-285, doi:10.1051/0004-6361:20041335.
#' @references Michel Crucifix (2016). palinsol: Insolation for
#'  Palaeoclimate Studies. R package version 0.93.
#'  [CRAN page](https://CRAN.R-project.org/package=palinsol)
#' @examples
#' data(Laskar04)
"Laskar04"
