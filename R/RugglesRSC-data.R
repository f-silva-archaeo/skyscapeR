#' Recumbent Stone Circle Data
#'
#' Declination data from C.L.N. Ruggles' fieldwork on the Scottish Recumbent Stone Circles.
#'
#' @docType data
#' @usage data(RugglesRSC)
#' @format A data frame with 37 rows and 2 variables:
#'  \describe{
#'    \item{Dec}{Declination}
#'    \item{ID}{Site ID}
#'    }
#' @keywords datasets
#' @references Ruggles, C.L.N. (1999). \emph{Astronomy in Prehistoric Britain and Ireland}. Yale University Press.
#' @examples
#' data(RugglesRSC)
#' kde <- density(RugglesRSC$Dec, 2)
#' plot(kde)
"RugglesRSC"
