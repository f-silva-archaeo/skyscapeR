#' Creates a \emph{skyscapeR.object} for plotting of celestial objects at given epoch
#'
#' This function creates an object containing all the necessary information to
#' plot celestial objects/events unto the many plotting functions of \emph{skyscapeR}
#' package.
#' @param names The name(s) of the celestial object(s) or event(s) of interest.
#' These can be one of the following soli-lunar events: \emph{jS}, \emph{dS}, \emph{eq},
#'  \emph{nmnLX}, \emph{nMjLX},
#' \emph{smnLX}, \emph{sMjLX}, or the name of any star in the database. As shorthand, the names
#' \emph{sun} and \emph{moon} can be used to represent all the above solar and lunar events,
#' respectively. Alternatively, custom declination values can also be used.
#' @param epoch The year or year range (as an array) one is interested in.
#' @param col (Optional) The colour for plotting, and differentiating these objects.
#' Defaults to red for all objects.
#' @param lty (Optional) Line type (see \code{\link{par}}) used for differentiation.
#' Only activated for single year epochs.
#' @param lwd (Optional) Line width (see \code{\link{par}}) used for differentiation.
#' Only activated for single year epochs.
#' @export
#' @examples
#' \dontrun{
#' # Create a object with solar targets for epoch range 4000-2000 BC:
#' tt <- sky.objects('sun', c(-4000,-2000))
#'
#' # Create an object with a few stars for same epoch:
#' tt <- sky.objects(c('Sirius', 'Betelgeuse', 'Antares'), c(-4000,-2000),
#' col=c('white', 'red', 'orange'))
#'
#' # Create an object with solstices and a custom declination value:
#' tt <- sky.objects(c('dS','jS', -13), c(-4000,-2000))
#' }
sky.objects = function(names, epoch, col = 'red', lty = 1, lwd = 1) {
  N <- NROW(names)

  if (NROW(col)==1) { col <- rep(col,N) }
  if (NROW(lty)==1) { lty <- rep(lty,N) }
  if (NROW(lwd)==1) { lwd <- rep(lwd,N) }

  tt <- c()
  tt.col <- c()
  tt.lty <- c()
  tt.lwd <- c()

  options(warn=-2)
  for (i in 1:N) {
    if (is.na(as.numeric(names[i]))) {
      # sun and moon shorthand
      if (pracma::strcmp(names[i], 'sun')) {
        aux <- array(NA, c(NROW(epoch),3))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- c(dS(epoch[j]), jS(epoch[j]), eq(epoch[j]))
        }
        colnames(aux) <- c('dS','jS','eq')
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, rep(col[i], 3))
        tt.lty <- c(tt.lty, rep(lty[i], 3))
        tt.lwd <- c(tt.lwd, rep(lwd[i], 3))
        next
      }

      if (pracma::strcmp(names[i], 'moon')) {
        aux <- array(NA, c(NROW(epoch),4))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- c(sMjLX(epoch[j]), smnLX(epoch[j]), nmnLX(epoch[j]), nMjLX(epoch[j]))
        }
        colnames(aux) <- c('sMjLX', 'smnLX', 'nmnLX', 'nMjLX')
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, rep(col[i], 4))
        tt.lty <- c(tt.lty, rep(lty[i], 4))
        tt.lwd <- c(tt.lwd, rep(lwd[i], 4))
        next
      }

      if (pracma::strcmp(names[i], 'sunandmoon')) {
        aux <- array(NA, c(NROW(epoch),7))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- c(dS(epoch[j]), jS(epoch[j]), eq(epoch[j]), sMjLX(epoch[j]), smnLX(epoch[j]), nmnLX(epoch[j]), nMjLX(epoch[j]))
        }
        colnames(aux) <- c('dS','jS','eq','sMjLX', 'smnLX', 'nmnLX', 'nMjLX')
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, rep(col[i], 7))
        tt.lty <- c(tt.lty, rep(lty[i], 7))
        tt.lwd <- c(tt.lwd, rep(lwd[i], 7))
        next
      }

      # stars
      data(stars, envir=environment())
      if (sum(sapply(as.character(stars$NAME), pracma::strcmp, s2=names[i]))) {
        aux <- array(NA, c(NROW(epoch),1))
        for (j in 1:NROW(epoch)) {
          aux[j,] <- star(names[i], epoch[j])$dec
        }
        colnames(aux) <- names[i]
        tt <- cbind(tt, aux)
        tt.col <- c(tt.col, col[i])
        tt.lty <- c(tt.lty, lty[i])
        tt.lwd <- c(tt.lwd, lwd[i])
        next
      }

      # if none of the above fit the bill try calling the functions
      aux <- array(NA, c(NROW(epoch),1))
      for (j in 1:NROW(epoch)) {
        aux[j,] <- do.call(names[i], list(epoch[j]))
      }
      colnames(aux) <- names[i]
      tt <- cbind(tt, aux)
      tt.col <- c(tt.col, col[i])
      tt.lty <- c(tt.lty, lty[i])
      tt.lwd <- c(tt.lwd, lwd[i])
    }

    # custom dec
    if (!is.na(as.numeric(names[i]))) {
      aux <- array(NA, c(NROW(epoch),1))
      aux[,1] <- rep(as.numeric(names[i]),NROW(epoch))
      colnames(aux) <- paste0('Custom Dec:',names[i])
      tt <- cbind(tt, aux)
      tt.col <- c(tt.col, col[i])
      tt.lty <- c(tt.lty, lty[i])
      tt.lwd <- c(tt.lwd, lwd[i])
      next
    }

  }
  options(warn=0)
  rownames(tt) <- epoch

  # check min and max decs
  if (NROW(epoch)==2) {
    ind <- which(substr(colnames(tt),1,6) != 'Custom')
    dec.range <- sapply(colnames(tt[,ind]), minmaxdec, from=min(epoch), to=max(epoch))
    aux <- matrix(NA, 2, NCOL(tt)); aux[,ind] <- dec.range; aux[,-ind] <- tt[,-ind]; rownames(aux) <- c('min','max')
    tt <- rbind(tt,aux)
  }

  # return result
  object <- c()
  object$n <- NCOL(tt)
  object$decs <- tt
  object$epoch <- epoch
  object$col <- tt.col
  if (NROW(epoch)==1) {
    object$lty <- tt.lty
    object$lwd <- tt.lwd
  }
  class(object) <- "skyscapeR.object"
  return(object)
}


#' Declination of December Solstice for a given year
#'
#' This function calculates the declination of the sun
#' at December Solstice for a given year, based upon
#' obliquity estimation.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{jS}}, \code{\link{eq}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' # December Solstice declination for year 3999 BC:
#' dS(-4000)
dS = function(year = cur.year) {
  aux <- obliquity(year)
  return(-aux)
}


#' Declination of June Solstice for a given year
#'
#' This function calculates the declination of the sun
#' at June Solstice for a given year, based upon
#' obliquity estimation.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link{obliquity}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' # June Solstice declination for year 3999 BC:
#' jS(-4000)
jS = function(year = cur.year) {
  aux <- obliquity(year)
  return(aux)
}


#' Declination of northern minor Lunar Extreme for a given year
#'
#' This function calculates the declination of the northern
#' minor Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link{smnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}}
#' @examples
#' # Northern minor Lunar Extreme declination for year 2499 BC:
#' nmnLX(-2500)
nmnLX = function(year = cur.year) {
  aux <- obliquity(year) - 5.145
  return(aux)
}


#' Declination of southern minor Lunar Extreme for a given year
#'
#' This function calculates the declination of the southern
#' minor Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{nMjLX}}, \code{\link{sMjLX}}
#' @examples
#' # Southern minor Lunar Extreme declination for year 2499 BC:
#' smnLX(-2500)
smnLX = function(year = cur.year) {
  aux <- obliquity(year) - 5.145
  return(-aux)
}


#' Declination of northern major Lunar Extreme for a given year
#'
#' This function calculates the declination of the northern
#' major Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{smnLX}}, \code{\link{sMjLX}}
#' @examples
#' # Northern major Lunar Extreme declination for year 2499 BC:
#' nMjLX(-2500)
nMjLX = function(year = cur.year) {
  aux <- obliquity(year) + 5.145
  return(aux)
}


#' Declination of southern major Lunar Extreme for a given year
#'
#' This function calculates the declination of the southern
#' major Lunar Extreme for a given year, by simple addition
#' of obliquity with maximum lunar inclination.
#' @param year Year for which to calculate the declination.
#' Defaults to present year as given by \emph{Sys.Date()}.
#' @export
#' @seealso \code{\link{nmnLX}}, \code{\link{nMjLX}}, \code{\link{smnLX}}
#' @examples
#' # Southern major Lunar Extreme declination for year 2499 BC:
#' sMjLX(-2500)
sMjLX = function(year = cur.year) {
  aux <- obliquity(year) + 5.145
  return(-aux)
}

#' Declination of sun at the equinoxes
#'
#' This function always returns a value of zero, which is the
#' declination of the sun on the day of the (astronomical)
#' equinoxes.
#' @param bh \emph{NULL} parameter. Can be left empty.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{zenith}}, \code{\link{antizenith}}
#' @examples
#' eq()
eq = function(bh=NULL) {
  return(0)
}

#' Declination of the zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the zenith for a given location. If
#'  this phenomena does not occur at given location
#'  (i.e. if location is outside the tropical band)
#'  the function returns a \emph{NULL} value.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{antizenith}}
#' @examples
#' # Zenith sun declination for Mexico City:
#' zenith(19.419)
#'
#' # There is no zenith sun phenomena in London:
#' zenith(51.507)
zenith = function(loc) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
  } else { lat <- loc }

  if (lat > jS() | lat < dS()) {
    return(NULL)
  } else { return(lat) }
}

#' Declination of the anti-zenith sun for a given location
#'
#' This function returns the declination of the sun
#' when it is at the anti-zenith, or nadir, for a given
#' location. If this phenomena does not occur at given
#' location (i.e. if location is outside the tropical
#' band) the function returns a \emph{NULL} value.
#' @param loc This can be either the latitude of the
#' location, or a \emph{skyscapeR.horizon} object.
#' @export
#' @seealso \code{\link{jS}}, \code{\link{dS}}, \code{\link{eq}}, \code{\link{zenith}}
#' @examples
#' # Anti-zenith sun declination for Mexico City:
#' antizenith(19.419)
#'
#' # There is no anti-zenith sun phenomena in London:
#' antizenith(51.507)
antizenith = function(loc) {
  if (class(loc)=='skyscapeR.horizon') {
    lat <- loc$georef[1]
  } else { lat <- loc }

  if (lat > jS() | lat < dS()) {
    return(NULL)
  } else { return(-lat) }
}
