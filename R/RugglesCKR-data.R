#' Cork and Kerry Stone Row Data
#'
#' Data from C.L.N. Ruggles' fieldwork on the Stone Rows of Cork and Kerry.
#'
#' @docType data
#' @usage data(RugglesCKR)
#' @format A data frame with 41 rows and 5 variables:
#'  \describe{
#'    \item{Ref}{Site Ref}
#'    \item{NE.SW}{String indicating whether they are towards the NE or SW}
#'    \item{Az.Hill}{Azimuth of hill towards which the Stone Rows are pointing}
#'    \item{Alt.Hill}{Altitude of hill towards which the Stone Rows are pointing}
#'    \item{Dec.Hill}{Declination of hill towards which the Stone Rows are pointing}
#'    }
#' @keywords datasets
#' @references Ruggles, C.L.N. (1999). \emph{Astronomy in Prehistoric Britain and Ireland}. Yale University Press.
#' @examples
#' data(RugglesCKR)
#' kde <- density(RugglesCKR$Dec.Hill, 2)
#' plot(kde)
"RugglesCKR"
