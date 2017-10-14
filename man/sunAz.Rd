% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mechanics.R
\name{sunAz}
\alias{sunAz}
\title{Returns the azimuth of the sun at a given time from a specific location}
\usage{
sunAz(loc, time, timezone)
}
\arguments{
\item{loc}{Location, either a \emph{skyscapeR.object} or a vector
containing the latitude and longitude of location, in this order.}

\item{time}{String containing the date and time in the following format:
"YYYY-MM-DD HH:MM:SS"}

\item{timezone}{Timezone of input either as a known acronym (eg. "GMT", "CET") or
a string with continent followed by country capital (eg. "Europe/London").}
}
\description{
This function returns the azimuth of the sun at a given time and location,
useful for data reduction of theodolite mesaurements using the sunsight
technique (\code{\link{reduct.theodolite}}).
}
\examples{
sunAz(c(52,-3), '2017-10-04 12:32:14', 'Europe/London')
}
\seealso{
\code{\link{reduct.theodolite}}
}