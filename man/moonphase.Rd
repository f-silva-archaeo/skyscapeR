% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mechanics.R
\name{moonphase}
\alias{moonphase}
\title{Computes the phase of the moon}
\usage{
moonphase(time, timezone, calendar)
}
\arguments{
\item{time}{Either a string containing the date and time in the format "YYYY-MM-DD HH:MM:SS"
(see \code{\link{timestring}}), or a numeric containing the julian date (see \code{\link{time2jd}}).}

\item{timezone}{(Optional) Timezone of input either as a known acronym (e.g. "GMT", "CET") or
a string with continent followed by country capital (e.g. "Europe/London"). See
\link{timezones} for details. Only needed if \emph{time} is a string. Defaults to system timezone.}

\item{calendar}{(Optional) Calendar used in parameter \emph{time}. G for gregorian and J for julian.
Only needed if \emph{time} is a string. Defaults to Gregorian.}
}
\description{
This function calculates the moon phase, in percentage of full. It is a wrapper
for function \code{\link[swephR]{swe_pheno_ut}} of package \emph{swephR}.
}
\examples{
# Moonphase at noon GMT on Christmas day 2018:
moonphase('2018/12/25 12:00:00', 'GMT')
}
\seealso{
\code{\link[swephR]{swe_pheno_ut}}
}
