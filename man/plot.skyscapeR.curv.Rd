% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot.skyscapeR.curv}
\alias{plot.skyscapeR.curv}
\title{Plot a curvigram}
\usage{
\method{plot}{skyscapeR.curv}(cc, col = "blue", shading = T, xrange, yrange, obj, show.obj.label = T)
}
\arguments{
\item{cc}{Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.}

\item{col}{(Optional) Colour to plot the curvigram in. Defaults to blue.}

\item{shading}{(Optional) Whether to shade the curvigram. Defaults to TRUE.}

\item{xrange}{(Optional) Array of two values restricting the horizontal range of the plot.}

\item{yrange}{(Optional) Array of two values restricting the vertical range of the plot.}

\item{obj}{(Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
for displaying the declination of celestial objects.}

\item{show.obj.label}{(Optional) Boolean to control whether to label the celestial objects in
the polar plot. Defaults to \emph{TRUE}.}
}
\description{
This function creates a plot of a curvigram
}
\examples{
# Plot the curvigram of Recumbent Stone Circles:
data(RugglesRSC)
curv <- curvigram(RugglesRSC$Dec, unc=2)
plot(curv)

# Redo the plot to include lunar extreme declinations:
LEx <- sky.objects('lunar extremes', -2000, col='red', lty=2)
plot(curv, obj=LEx)
}
\seealso{
\code{\link{curvigram}}, \code{\link{histogram}}, \code{\link{sky.objects}}
}
