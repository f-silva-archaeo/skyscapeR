utils::globalVariables(c("constellations"))
#' Constellations
#'
#' Constellations database. Currently only IAU / Modern western constellations
#' included, but can be expanded to include those of other societies.
#'
#' @docType data
#' @usage data(constellations)
#' @format A data frame with 300 rows and 2 variables:
#'  \describe{
#'    \item{set}{Object identifier used for ephemeris}
#'    \item{name1}{Full name of constellation}
#'    \item{name2}{Shorter name for constellation}
#'    \item{name3}{Alternative name for constellation}
#'    \item{name4}{Alternative name for constellation}
#'    \item{HIP}{String with HIP catalogue number of stars to join together.
#'    This should be in the same format that Stellarium uses for its sky cultures.}
#'    }
#' @keywords datasets
"constellations"
