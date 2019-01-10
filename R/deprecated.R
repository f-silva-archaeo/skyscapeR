#' Declination distribution corresponding to a uniform azimuthal
#' distribution for null hypothesis significance testing
#'
#' This function returns the declination distribution that
#' corresponds to a uniform distribution in azimuths (i.e. random
#'  orientation) for use in significance testing.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @param alt (Optional) The horizon altitude to use in
#' \code{\link{az2dec}} conversion. Defaults to 0 degrees.
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' aux <- nh.Uniform(loc=c(52,-2), alt=2)
#' plot(aux$dec, aux$density, type='l')
#' }
#' @details This function is deprecated. Please see \code{\link{sigTest}} instead.
nh.Uniform = function(loc, alt=0) {
  .Defunct('sigTest')
}


#' Plot visible path of celestial object on horizon
#'
#' This function adds the visible path of a celestial object to a horizon plot.
#' @param orbit Object of \emph{skyscapeR.orbit} format.
#' @param col String with colour to plot path in. Defaults to red.
#' @seealso \code{\link{plot.skyscapeR.horizon}}, \code{\link{orbit}}
#' @import utils stats graphics
#' @noRd
plot.skyscapeR.orbit<- function(orbit, col) {
  .Defunct('plot.skyscapeR.horizon')
}


#' Declination distribution of the Sun throughout the year
#'  for null hypothesis significance testing
#'
#' This function returns the declination distribution of
#' the the Sun throughout the year for use in significance testing.
#' @param year Year for which to calculate the distribution.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' aux <- nh.SolarRange(-4000)
#' plot(aux$dec, aux$density, type='l')
#' }
#' @details This function is deprecated. Please see \code{\link{sigTest}} instead.
nh.SolarRange = function(year = cur.year) {
  .Defunct('sigTest')
}


#' Declination distribution of the Moon throughout the year
#'  for null hypothesis significance testing
#'
#' This function returns the declination distribution of
#' the the Sun throughout the year for use in significance testing.
#' @param year Year for which to calculate the distribution.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' aux <- nh.SolarRange(-4000)
#' plot(aux$dec, aux$density, type='l')
#' }
#' @details This function is deprecated. Please see \code{\link{sigTest}} instead.
nh.LunarRange = function(year = cur.year) {
  .Defunct('sigTest')
}


#' Declination distribution of the Summer Full Moon
#'  for null hypothesis significance testing
#'
#' This function returns the declination distribution of
#' the Summer Full Moon for use in significance testing.
#' @param min.phase (Optional) This should be the minimum
#'  lunar phase (i.e. percentage illumination) for the moon
#'   to be considered full. The value should range between
#'   0 (dark moon) and 1 (full moon). Defaults to 0.99.
#' @param min.sundec (Optional) This should be the minimum
#' solar declination for the moon to be considered a
#' \emph{summmer} full moon. Defaults to 20 degrees,
#' corresponding to a month before or after june solstice.
#' @param year Year for which to calculate the obliquity.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' aux <- nh.SummerFM(.99, 20, -4000)
#' plot(aux$dec, aux$density, type='l')
#' }
#' @details This function is deprecated. Please see \code{\link{sigTest}} instead.
nh.SummerFM = function(min.phase = .99, min.sundec = 20, year = cur.year) {
  .Defunct('sigTest')
}


#' Plot a curvigram
#'
#' This function creates a plot of a curvigram.
#' @param curv Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the curvigram. Defaults to \emph{TRUE}.
#' @param xlim Array of two values restricting the horizontal range of the plot.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @import utils stats graphics
#' @export
#' @seealso \code{\link{curvigram}}, \code{\link{sky.objects}}, \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' # Plot the curvigram of Recumbent Stone Circles:
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, unc=2)
#' plotCurv(curv, xlim=c(-40,0))
#'
#' # Redo the plot to include lunar extreme declinations:
#' LEx <- sky.objects('moon', -2000, col='red', lty=2)
#' plotCurv(curv, obj=LEx, xlim=c(-40,0))
#' }
#' @details This function is deprecated. Please see \code{\link{plot.skyscapeR.curv}} instead.

plotCurv = function(curv, obj, obj.label=T, xlim=NULL, ...) {
  .Defunct('plot.skyscapeR.curv')
  if (class(curv)!='skyscapeR.curv') { stop('No skyscapeR.curv object found.') }
  plot(curv, obj, obj.label, xlim=xlim, ...)
}


#' Plot horizon data
#'
#' This function creates a plot of horizon data.
#' @param hor Object of \emph{skyscapeR.horizon} format.
#' @param show.az Boolean that controls whether to display azimuth values on horizontal axis.
#'  Defaults to \emph{FALSE}.
#' @param max.alt Maximum altitude to display. Defaults to 45 degrees.
#' @param az0 Leftmost azimuth of plot. Defaults to 0 degrees, i.e. North at the left.
#' @param zoom Boolean that controls whether to provide a zoomed-in view of 100 degrees in
#' azimuth and 5 degrees of altitude above the horizon line. Defaults to \emph{FALSE}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the paths of celestial objects.
#' @param measure (Optional) A data.frame object with columns \emph{True.Azimuth} and
#' \emph{Altitude} such as the ones produced by \code{\link{reduct.compass}} or
#' \code{\link{reduct.theodolite}}. If no column named \emph{Altitude} is found then it will
#' plot all azimuths are zero degrees altitude.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @import utils stats graphics grDevices
#' @seealso \code{\link{downloadHWT}}, \code{\link{sky.objects}}
#' @examples
#' \dontrun{
#' # Plot a horizon retrieved from HeyWhatsThat:
#' hor <- downloadHWT('HIFVTBGK')
#' plotHor(hor)
#'
#' # Add the paths of the solstices and equinoxes sun in the year 1999 BC:
#' tt <- sky.objects('sun', -2000, 'blue')
#' plotHor(hor, objects=tt)
#' }
#' @details This function is deprecated. Please see \code{\link{plot.skyscapeR.horizon}} instead.
plotHor <- function(hor, show.az=F, max.alt, az0 = 0, zoom=F, obj, measure, ...) {
  .Defunct('plot.skyscapeR.horizon')
  if (class(hor)!='skyscapeR.horizon') { stop('No skyscapeR.horizon object found.') }
  plot(hor, show.unc=F, show.az, max.alt, az0, zoom, obj, measure, ...)
}


#' Plot stellar phase and seasonality
#'
#' This function creates a plot of stellar seasonality and phases/events.
#' @param starphase Object of \emph{skyscapeR.starphase} format.
#' @param ... Any other parameters to be passed unto \code{\link{plot.default}}.
#' @export
#' @import MESS RColorBrewer
#' @seealso \code{\link{star.phases}}
#' @examples
#' # Plot the seasonality of Aldebaran for 3999 BCE:
#' \dontrun{
#' ss <- star.phases('Aldebaran',-4000, c(35,-8))
#' plotPhases(ss)
#' }
#' @details This function is deprecated. Please see \code{\link{plot.skyscapeR.starphase}} instead.
plotPhases = function(starphase, ...) {
  .Defunct('plot.skyscapeR.starphase')
  if (class(starphase)!='skyscapeR.starphase') { stop('No skyscapeR.starphase object found.') }
  plot(starphase, ...)
}


#' Plot visible path of celestial object on horizon
#'
#' This function adds the visible path of a celestial object to a horizon plot.
#' @param orbit Object of \emph{skyscapeR.orbit} format.
#' @param col String with colour to plot path in. Defaults to red.
#' @seealso \code{\link{plot.skyscapeR.horizon}}, \code{\link{orbit}}
#' @import utils stats graphics
#' @noRd
plotOrb<- function(orbit, col) {
  .Defunct('plot.skyscapeR.orbit')
  if (class(orbit)!='skyscapeR.orbit') { stop('No skyscapeR.orbit object found.') }
  plot(orbit, col)
}


#' Plot a z-score transformed curvigram
#'
#' This function creates a plot of a z-score transformed curvigram, which is to say the
#' curvigram transformed into sigma units, based on a previously generated significance test.
#' @param signif A \emph{skyscapeR.sig} object created with \code{\link{sigTest}}.
#' @param obj (Optional) A \emph{skyscapeR.object} object created with \code{\link{sky.objects}}
#' for displaying the declination of celestial objects.
#' @param obj.label (Optional) Boolean to control whether to label the celestial objects in
#' the curvigram. Defaults to \emph{TRUE}.
#' @param xlim Array of two values restricting the horizontal range of the plot.
#' @export
#' @seealso \code{\link{sigTest}}
#' @examples
#' \dontrun{
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, unc=2)
#' sig <- sigTest(curv, nh.Uniform(c(57,2)))
#'
#' plotZscore(sig)
#' }
#' @details This function is deprecated. Please see \code{\link{plot.skyscapeR.sigTest}} instead.
plotZscore = function(signif, obj, obj.label=T, xlim=NULL) {
  .Defunct('plot.skyscapeR.sigTest')
}


#' Download horizon data from \emph{HeyWhatsThat}
#'
#' This function downloads previously created horizon data
#' from \emph{HeyWhatsThat}, given its ID, and saves it as
#' a \emph{skyscapeR.horizon} object.
#' @param HWTID This is the 8 character ID attributed by
#' \emph{HeyWhatsThat.com}
#' @export
#' @import utils
#' @references \href{http://heywhatsthat.com/}{HeyWhatsThat.com}
#' @seealso \code{\link{createHWT}}
#' @examples
#' \dontrun{
#' # Retrieve horizon data for \href{https://www.heywhatsthat.com/?view=HIFVTBGK}{Liverpool Cathedral}:
#' hor <- download.HWT('HIFVTBGK')
#' }
download.HWT <- function(HWTID) {
  .Defunct('downloadHWT')
}
