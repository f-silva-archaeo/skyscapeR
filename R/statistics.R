
#' @noRd
bernoulli <- function(n,p,r) {
  return( factorial(n) / (factorial(r) * factorial(n-r)) * p^r * (1-p)^(n-r) )
}


#' Execute a (series of) Bernoulli trial(s)
#'
#' This function allows one to calculate the probability of having r structures
#' out of n, orientated towards a target with probability p.
#' @param n Total number of structures
#' @param p Probability of target (e.g. ratio of azimuths)
#' @param r Number of structures orientated towards target (hits)
#' @param type (Optional) Type of probability to output. Possibilities are: (a) \emph{single}
#' in which case the result of a single Bernoulli trial is reported; or (b) \emph{tail} in
#' which case it calculates Bernoulli trials for all hit values between 1 and (r-1) and then
#' outputs 1 minus the calculated probability. The latter is effectively a p-value. Default is \emph{tail}
#' @export
#' @references Ruggles, C (1999) \emph{Astronomy in Prehistoric Britain and Ireland}. Yale University Press.
#' @examples
#' # probability of having at least 10 out of 30 structures
#' # aligned to targets covering 20% of the horizon
#' bernoulli.trial(30, 0.2, 10)
bernoulli.trial <- function(n,p,r, type='tail') {
  if (type=='single') { return(bernoulli(n,p,r)) }

  if (type=='tail') {
    out <- 0
    for (i in 1:r) {
      out <- out + bernoulli(n,p,i-1)
    }
    return(1-out)
  }
}


#' Convert discrete azimuth measurements into probability distributions
#'
#' @param pdf (Optional) String describing the probability distribution to be used. At the
#' moment only \emph{normal} and \emph{uniform} are supported. Default is \emph{normal}
#' @param az An array of azimuths
#' @param unc Azimuth uncertainties as either an array of the same length as \emph{az} or a single
#' value to be applied to all measurements
#' @param name (Optional) An array of names to identify each measurement
#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @param .cutoff (Optional) Value of probability distribution(s) at which point it will be cutoff to save on memory. Default is 1e-4
#' @param .res (Optional) Azimuth resolution with which to output probability distribution(s). Default is 0.01 degrees.
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. <doi:10.1016/j.jas.2020.105138>
#' @examples
#' test <- az.pdf(az=c(87,93,90,110), unc=3)
#' plot(test)
az.pdf <- function(pdf='normal', az, unc, name, verbose=T, .cutoff=1e-4, .res=0.01) {

  ## probability distribution
  if (is(pdf,'character')) {

    code <- switch(pdf, N = 0, Normal = 0, n = 0, normal = 0, U = 1, Uniform = 1, u = 1, uniform = 1, -1)

    if (code == 0) {
      if (verbose) { cat('Normal probability distribution(s) chosen\n') }
      azpdf <- function(x,i) dnorm(x, az[i], unc[i])
    } else if (code == 1) {
      if (verbose) { cat('Uniform probability distribution(s) chosen\n') }
      azpdf <- function(x,i) dunif(x, az[i]-unc[i], az[i]+unc[i])
    } else { stop('Probability density function not recognized. Please use only normal or uniform.')}

    if (missing(az) | missing (unc)) {
      stop('Missing azimuths and/or uncertainties.')
    }

  } else if (is(pdf, 'data.frame') & sum(colnames(pdf) == c('az','dens'))) {
    if (verbose) { stop('Custom probability distribution found. This feature is not yet functional.') }

  } else if (is(pdf, 'list') & is(pdf[[1]], 'data.frame') & sum(colnames(pdf[[1]]) == c('az','dens'))) {
    if (verbose) { stop('Custom probability distribution(s) found. This feature is not yet functional.') }

  } else {
    stop('Probability density function not recognized. Please use only normal or uniform.')# or a list of data.frames with columns "az" and "dens".')
  }

  ## init parameters
  n <- length(az)
  if (n > 1 & length(unc) == 1) {
    if (verbose) { cat('Single uncertainty value found. Using it for all measurements.\n') }
    unc <- rep(unc, length(az))
  }

  xx <- seq(0-90, 360+90, .res)

  out <- list()
  if (verbose & length(az)>1) { pb <- txtProgressBar(max = length(az), style = 3) }

  for (k in 1:n) {
    yy <- azpdf(xx,k)

    ## clean-up due to circularity of azimuths
    ind <- which(xx < 0)
    yy[ind + 360/.res] <- yy[ind + 360/.res] + yy[ind]
    x <- xx[-ind]; yy <- yy[-ind]
    ind <- which(x > 359.99)
    yy[ind - 360/.res] <- yy[ind - 360/.res] + yy[ind]
    x <- x[-ind]; yy <- yy[-ind]

    ind <- which(yy >= .cutoff)
    out[[k]] <- data.frame(x = x[ind], y = yy[ind])

    if (verbose & length(az)>1) { setTxtProgressBar(pb, k) }
  }

  if (missing(name)) { name <- paste('Measurement',seq(1,length(az))) }

  mtdta <- list(
    coord = 'az',
    pdf = pdf,
    name = name,
    az = az,
    unc = unc)
  mtdta$param <- list(.cutoff=.cutoff, .res=.res)
  out <- list(metadata = mtdta, data = out)
  class(out) <- 'skyscapeR.pdf'

  if (verbose) { cat('\nDone.') }
  return(out)
}


#' @noRd
coordtrans_unc <- function(pdf, hor, xrange, refraction, atm, temp, verbose=T, .res=0.1) {

  if ((is(pdf, 'skyscapeR.pdf') & pdf$metadata$coord != 'az') & (!is(pdf, 'skyscapeR.spd'))) stop('Azimuthal probability density function(s) not recognized or not in correct format.')
  if (!is(hor, 'skyscapeR.horizon') & !is(hor, 'list')) stop('Horizon(s) not recognized or not in correct format.')

  ## init parameters
  if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
  if (missing(atm)) { atm <- skyscapeR.env$atm }
  if (missing(temp)) { temp <- skyscapeR.env$temp }

  if (is(pdf, 'skyscapeR.spd')) {
    n <- 1
    pdf$data <- list(pdf$data)
  } else {
    n <- length(pdf$metadata$az)
  }
  .cutoff <- pdf$metadata$param$.cutoff


  if (is(hor, "skyscapeR.horizon")) {
    if (verbose) { cat('Single horizon profile found. Using it for all measurements.\n') }
    j <- 1; hh <- list(hor)
  } else {
    pointr::ptr('j','k')
    hh <- hor
  }

  xx <- seq(-90, 360+90, .1)

  out <- list()
  if (verbose & n>1) { pb <- txtProgressBar(max = n, style = 3) }

  for (k in 1:n) {
    azs <- pdf$data[[k]]$x

    fhor <- approx(hh[[j]]$data$az, hh[[j]]$data$alt, xout=azs)$y
    fhor.unc <- approx(hh[[j]]$data$az, hh[[j]]$data$alt.unc, xout=azs)$y

    aux1 <- az2dec(azs, hh[[j]], fhor-4*fhor.unc, refraction=refraction, atm=atm, temp=temp)
    aux2 <- az2dec(azs, hh[[j]], fhor+4*fhor.unc, refraction=refraction, atm=atm, temp=temp)

    if (missing(xrange)) {
      min.dec <- min(aux1, aux2, na.rm=T); max.dec <- max(aux1, aux2, na.rm=T)
    } else {
      min.dec <- xrange[1]; max.dec <- xrange[2]
    }

    dec <- seq(min.dec, max.dec, by = .res)
    mean.dec <- (aux2+aux1)/2
    sd.dec <- (aux2-aux1)/8

    dens <- rep(NA, length(dec))

    fn <- function(x, mean.dec, sd.dec, azs, pdf, xx) {
      d <- dnorm(x, mean.dec, sd.dec)
      f <- approxfun(azs, d, yleft=0, yright=0)

      # do convolution manually
      conv <- f(xx) * approx(pdf$data[[k]]$x, pdf$data[[k]]$y, xout=xx, yleft=0, yright=0)$y
      dens <- MESS::auc(xx, conv)
      return(dens)
    }
    dens <- sapply(dec, fn, mean.dec=mean.dec, sd.dec=sd.dec, azs=azs, pdf=pdf, xx=xx)

    dens <- spline(dec, dens, xout=seq(min.dec, max.dec, .01))
    dec <- dens$x; dens <- dens$y

    if (!is(pdf, 'skyscapeR.spd')) { dens <- dens/MESS::auc(dec,dens) } # normalise

    out[[k]] <- data.frame(x=dec, y=dens)
    if (verbose & n>1) { setTxtProgressBar(pb, k) }
  }

  mtdta <- list(
    coord = 'dec',
    name = pdf$metadata$name
  )
  mtdta$param <- list(refraction = refraction, atm=atm, temp=temp, .cutoff=.cutoff, .res=.res)
  mtdta$horizon <- deparse(substitute(hor))
  mtdta$az.pdf <- deparse(substitute(pdf))
  out <- list(metadata = mtdta, data = out)
  if (is(pdf, 'skyscapeR.spd')) {
    class(out) <- 'skyscapeR.spd'
    out$data <- out$data[[1]]
  } else {
    class(out) <- 'skyscapeR.pdf'
  }
  if (verbose) { cat('\nDone.') }
  return(out)
}

#' Coordinate-transform azimuth prob distributions into declination prob distributions
#'
#' @param pdf A \emph{skyscapeR.pdf} object created with \code{\link{az.pdf}} or a \emph{skyscapeR.spd} object created with \code{\link{spd}}
#' @param hor A \emph{skyscapeR.horizon} object created with \code{\link{createHor}} or \code{\link{downloadHWT}}
#' @param xrange (Optional) Array of values (min and max) for SPD when transforming a \emph{skyscapeR.spd} object.
#' @param refraction (Optional) Whether atmospheric refraction is to be taken into account.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temperature for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @param .res (Optional) Declination resolution with which to output probability distribution(s). Default is 0.1 degrees.
#' @import pointr
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. <doi:10.1016/j.jas.2020.105138>
#' @examples
#' Az <- az.pdf(az=c(87,93,90,110), unc=3)
#' hor <- createHor(az=c(0,360), alt=c(0,0), loc=c(35,-8,25)) # flat horizon with 0 degrees of altitude
#' Dec <- coordtrans(Az, hor)
#' plot(Dec)
coordtrans <- compiler::cmpfun(coordtrans_unc, options=list(enableJIT = 3))


#' Summed probability density (SPD)
#'
#' @param pdf A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param normalise (Optional) Boolean to control whether to normalize the SPD. Default is FALSE
#' @param xrange (Optional) Array of values (min and max) for SPD if different from range of \emph{pdf}. If given,
#' it overrides \emph{.cutoff}
#' @param .cutoff (Optional) Value of SPD at which point it will be cutoff to save on memory. Default is 1e-5
#' @param .res (Optional) Resolution with which to output SPD. Default is 0.01 degrees.
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. <doi:10.1016/j.jas.2020.105138>
#' @examples
#' # SPD of azimuths
#' Az <- az.pdf(az=c(87,93,90,110), unc=3)
#' s1 <- spd(Az)
#' plot(s1)
#'
#' # SPD of declinations
#' hor <- createHor(az=c(0,360), alt=c(0,0), loc=c(35,-8,25)) # flat horizon with 0 degrees of altitude
#' Dec <- coordtrans(Az, hor)
#' s2 <- spd(Dec)
#' plot(s2)
spd <- function(pdf, normalise = F, xrange, .cutoff = 1e-5, .res = 0.01) {
  if (!is(pdf, 'skyscapeR.pdf')) { stop('Please provide a valid skyscapeR.pdf object.') }

  aux <- pdf$data
  x <- lapply(aux, "[[", 'x')
  y <- lapply(aux, "[[", 'y')
  if (missing(xrange)) {
    xrange <- range(x)
    cutoff <- T
  } else { cutoff <- F }

  xx <- seq(xrange[1], xrange[2], .res)
  spd <- rep(0, length(xx))
  for (i in 1:NROW(x)) {
    spd <- spd + approx(x[[i]], y[[i]], xout=xx, yleft=0, yright=0)$y
  }

  if (cutoff) { spd[spd < .cutoff] <- 0 }

  if (normalise) { spd <- spd/sum(spd, na.rm = T) }
  out <- c()
  out$metadata <- c()
  out$metadata$coord <- pdf$metadata$coord
  out$metadata$xrange <- xrange
  if (cutoff) { out$metadata$param$.cutoff <- .cutoff }
  out$metadata$param$.res <- .res
  out$data <- data.frame(x = xx, y = spd)

  class(out) <- 'skyscapeR.spd'
  return(out)
}



#' @noRd
randomTest_unc <- function(pdf, .res=0.1, nsims=1000, conf=.95, tails=2, normalise=F, ncores=parallel::detectCores()-1, verbose=T, hor, refraction, atm, temp) {
  quick <- F
  if (pdf$metadata$coord == 'dec' | !missing(hor) ) {
    cat('Running statistical significance test on declination.\n')
    coord <- 'dec'
    xrange <- c(-90,90)
    if (pdf$metadata$coord == 'dec') {
      az <- get(pdf$metadata$az.pdf, pos = 1)
      hor <- get(pdf$metadata$horizon, pos = 1)
      refraction <- pdf$metadata$param$refraction
      atm <- pdf$metadata$param$atm
      temp <- pdf$metadata$param$temp
      empirical <- spd(pdf, xrange=xrange, normalise = normalise, .res=.res)
    } else {
      az <- pdf
      if (!is(hor,"skyscapeR.horizon")) { stop('Please include a valid skyscapeR.horizon object.') }
      if (verbose) { cat('Single horizon profile found. Using it for all measurements.\n') }
      if (missing(refraction)) { refraction <- skyscapeR.env$refraction }
      if (missing(atm)) { atm <- skyscapeR.env$atm }
      if (missing(temp)) { temp <- skyscapeR.env$temp }
      quick <- T
      empirical <- spd(pdf, .res=.res)
      empirical <- coordtrans(empirical, hor, xrange=xrange, refraction=refraction, atm=atm, temp=temp, verbose=F, .res=.res)
    }
  } else {
    cat('Running statistical significance test on azimuth.\n')
    coord <- 'az'
    xrange <- c(0,360)
    az <- pdf
    empirical <- spd(pdf, normalise = normalise, xrange=xrange, .res=.res)
    ncores <- 1 ## Forces this due to strange error
  }

  if (ncores > 1) {
    cl <- snow::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    if (verbose) {
      cat(paste0('Running ', nsims,' simulations on ', ncores, ' processing cores.\n'))
      pb <- txtProgressBar(max=nsims, style=3)
      progress <- function(i) setTxtProgressBar(pb, i)
      opts <- list(progress = progress)
    }
    res <- matrix(NA, nsims, length(empirical$data$y))
    res <- foreach (i = 1:nsims, .combine=rbind, .inorder = F, .options.snow = opts) %dopar% {
      simAz <- sample(seq(0, 360, az$metadata$param$.res), length(az$metadata$name), replace=T)
      simPDF <- az.pdf(pdf=az$metadata$pdf, az=simAz, unc=az$metadata$unc, verbose=F)

      if (quick==T) {
        simPDF <- spd(simPDF)
        coordtrans(simPDF, hor, xrange=xrange, refraction=refraction, atm=atm, temp=temp, verbose=F, .res=.res)$data$y
      } else {
        if (coord == 'dec') {
          simPDF <- coordtrans(simPDF, hor, refraction=refraction, atm=atm, temp=temp, verbose=F, .res=.res)
        }
        spd(simPDF, xrange=xrange, normalise=normalise, .res=.res)$data$y
      }
    }
    snow::stopCluster(cl)
    if (verbose) { cat('\n') }

  } else {
    if (verbose) {
      cat(paste0('Running ', nsims,' simulations on a single processing core.\n'))
      pb <- txtProgressBar(max=nsims, style=3)
    }
    res <- matrix(NA, nsims, length(empirical$data$y))
    for (i in 1:nsims) {
      simAz <- sample(seq(0, 360, az$metadata$param$.res), length(az$metadata$name), replace=T)
      simPDF <- az.pdf(pdf=az$metadata$pdf, az=simAz, unc=az$metadata$unc, verbose=F)

      if (quick==T) {
        simPDF <- spd(simPDF)
        res[i,] <- coordtrans(simPDF, hor, xrange=xrange, refraction=refraction, atm=atm, temp=temp, verbose=F, .res=.res)$data$y
      } else {
        if (coord == 'dec') {
          simPDF <- coordtrans(simPDF, hor, refraction=refraction, atm=atm, temp=temp, verbose=F, .res=.res)
        }
        res[i,] <- spd(simPDF, xrange=xrange, normalise=normalise, .res=.res)$data$y
      }

      if (verbose) { setTxtProgressBar(pb, i) }
    }
    if (verbose) { cat('\n') }
  }

  ## confidence envelope
  zScore.sim <- matrix(0, nrow=nsims, ncol=length(empirical$data$y))
  zMean <- colMeans(res)
  zStd <- apply(res, 2, sd)

  zScore.emp <- (empirical$data$y - zMean)/zStd
  zScore.sim <- apply(res, 2, scale)

  if (verbose) { cat(paste0('Performing a ',tails,'-tailed test at the ', conf*100, '% significance level.\n')) }
  if (tails==2) {
    lvl.up <- 1-(1-conf)/2; lvl.dn <- (1-conf)/2
  } else if (tails==1) {
    lvl.up <- conf; lvl.dn <- 0
  } else { stop() }

  upper <- apply(zScore.sim, 2, quantile, probs = lvl.up, na.rm = T)
  upCI <- zMean + upper*zStd
  upCI[is.na(upCI)] <- 0
  if (tails==2) {
    lower <- apply(zScore.sim, 2, quantile, probs = lvl.dn, na.rm = T)
    loCI <- zMean + lower*zStd
    loCI[is.na(loCI)] <- 0
  }

  ## global p-value
  above <- which(zScore.emp > upper); emp.stat <- sum(zScore.emp[above] - upper[above])
  if (tails==2) { below <- which(zScore.emp < lower); emp.stat <- emp.stat + sum(lower[below] - zScore.emp[below]) }

  sim.stat <- abs(apply(zScore.sim, 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper))
  if (tails==2) { sim.stat <- sim.stat + abs(apply(zScore.sim, 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower)) }
  global.p <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)

  ## local p-value
  ind <- split(above, cumsum(c(1,diff(above) > 1))); ind <- ind[which(lengths(ind) > 1)]
  local <- data.frame(type=NA, start=0, end=0, p.value=0); j <- 0
  if (length(ind)>0) {
    for (j in 1:NROW(ind)) {
      emp.stat <- sum(zScore.emp[ind[[j]]] - upper[ind[[j]]])
      sim.stat <- abs(apply(zScore.sim[,ind[[j]]], 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper[ind[[j]]]))

      local[j,]$type <- '+'
      local[j,]$start <- min(empirical$data$x[ind[[j]]])
      local[j,]$end <- max(empirical$data$x[ind[[j]]])
      local[j,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
    }
  }

  if (tails==2) {
    ind <- split(below, cumsum(c(1,diff(below) > 1))); ind <- ind[which(lengths(ind) > 1)]
    if (length(ind)>0) {
      for (k in 1:NROW(ind)) {
        emp.stat <- sum(lower[ind[[k]]] - zScore.emp[ind[[k]]])
        sim.stat <- abs(apply(zScore.sim[,ind[[k]]], 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower[ind[[k]]]))

        local[j+k,]$type <- '-'
        local[j+k,]$start <- min(empirical$data$x[ind[[k]]])
        local[j+k,]$end <- max(empirical$data$x[ind[[k]]])
        local[j+k,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
      }
    }
  }

  ## cleanup
  rownames(local) <- c()
  aux <- apply(local[,c(2,3)], 1, diff)
  ind <- which(aux <= 10*pdf$metadata$param$.res + 1000*.Machine$double.eps)
  if (length(ind)>0) { local <- local[-ind,] }
  rownames(local) <- c()

  ## output
  out <- c()
  out$metadata$coord <- coord
  out$metadata$nsims <- nsims
  out$metadata$conf <- conf
  out$metadata$tails <- tails
  out$metadata$normalise <- normalise
  out$metadata$global.pval <- global.p
  out$metadata$local.pval <- local

  out$data$empirical <- empirical$data
  out$data$null.hyp <- list(x = empirical$data$x, CE.mean = zMean, CE.upper = upCI)
  if (tails==2) { out$data$null.hyp$CE.lower <- loCI }
  class(out) <- 'skyscapeR.sigTest'

  return(out)
}


#' Significance test against the null hypothesis of random orientation
#'
#' @param pdf A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param nsims (Optional) Number of simulations to run. Default is 1000.
#' @param conf (Optional) Confidence envelope (in percentage) of the null model to calculate. Default is .95
#' @param tails (Optional) Whether to calculate 1-tailed p-value (greater than) or 2-tailed p-value (smaller than or greater than).
#' Default is 2.
#' @param normalise (Optional) Boolean to control whether to normalize SPDs. Default is FALSE
#' @param ncores (Optional) Number of CPU cores to use. Default is the number of available cores minus 1.
#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @param .res (Optional) Resolution with which to output probability distribution(s). Default is 0.1 degrees.
#' @param hor (Optional) A single \emph{skyscapeR.horizon} object created with \code{\link{createHor}} or \code{\link{downloadHWT}}.
#' This and the following arguments are only needed if \emph{pdf} contains only azimuths. If provided, the
#' code will assume that all azimuthal measurements are for the same site with horizon given in \emph{hor} and
#' that the user want to run the significance test in declination. Coordinate transformation will therefore be done automatically.
#' @param refraction (Optional) Whether atmospheric refraction is to be taken into account.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param atm (Optional) Atmospheric pressure for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @param temp (Optional) Atmospheric temperature for refraction calculation.
#' If not given the value set by \code{\link{skyscapeR.vars}} will be used instead.
#' @import snow doSNOW foreach
#' @export
#' @references Silva, F (2020) A probabilistic framework and significance test for the analysis of structural orientations
#'  in skyscape archaeology \emph{Journal of Archaeological Science} 118, 105138. <doi:10.1016/j.jas.2020.105138>
#' @examples
#' \dontrun{
#' # significance test for azimuth
#' Az <- az.pdf(az=c(87,93,90,110), unc=3)
#' st1 <- randomTest(Az, nsims=1000)
#' plot(st1)
#'
#' # significance test for declination
#' hor <- createHor(az=c(0,360), alt=c(0,0), loc=c(35,-8,25)) # flat horizon with 0 degrees of altitude
#' Dec <- coordtrans(Az, hor)
#' st2 <- randomTest(Dec, nsims=1000)
#' plot(st2)
#' }
randomTest <- compiler::cmpfun(randomTest_unc, options=list(enableJIT = 3))



#' @noRd
modelTest_unc <- function(pdf, model, nsims=1000, conf=.95, tails=2, normalise=F, ncores=parallel::detectCores()-1, save.sim=F, verbose=T) {
  if (missing(model)) { stop('model is necessary to run this test.') }
  if (class(model) != 'data.frame' & class(model) != 'list') { stop('model must be given as a data.frame or list of data.frames.') }
  if (class(model) == 'list' & NROW(model) != NROW(pdf$data)) { stop('model must be a list of data.frames with equal number to the empirical observations.') }
  if (class(model) == 'data.frame') {
    model <- list(model)
    model <- rep(model, NROW(pdf$data))
  }
  if (save.sim & ncores>1) { stop('Cannot save simulated values when running on more than one core.') }
  if (pdf$metadata$coord == 'dec') {
    stop('Must implement dec version.')
  }
  if (pdf$metadata$coord == 'az') {
    xrange <- c(0,360)
    az <- pdf
  }

  ## empirical SPD
  if (verbose) { cat('Creating Empirical SPD...') }
  empirical <- spd(pdf, xrange=xrange, normalise = normalise)
  if (verbose) { cat('Done.\n') }


  if (ncores > 1) {
    ## bootstrapping
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, source("src v1.3.R"))
    if (verbose) { cat(paste0('Running ', nsims,' simulations on ', ncores, ' processing cores. This may take a while...')) }

    res <- matrix(NA, nsims, length(empirical$data$y))
    res <- foreach (i = 1:nsims, .combine=rbind, .inorder = F) %dopar% {
      simAz <- sampleList(model)
      simPDF <- az.pdf(pdf = pdf$metadata$pdf, az=simAz, unc=pdf$metadata$unc, verbose=F)

      if (pdf$metadata$coord == 'dec') {
        simPDF <- coordtrans(simPDF, hor, refraction=pdf$metadata$param$refraction, atm=pdf$metadata$param$atm, temp=pdf$metadata$param$temp, verbose=F, .prec=pdf$metadata$param$.prec)
      }
      spd(simPDF, xrange=xrange, normalise = normalise)$data$y
    }
    parallel::stopCluster(cl)
    if (verbose) { cat('Done.\n') }

  } else {
    if (verbose) { cat(paste0('Running ', nsims,' simulations on a single processing core. This may take a while...')) }
    if (save.sim) { sim.az <- matrix(NA, nsims, length(az)) }
    res <- matrix(NA, nsims, length(empirical$data$y))
    if (verbose) { pb <- txtProgressBar(max=nsims, style=3) }
    for (i in 1:nsims) {
      simAz <- sampleList(model)
      simPDF <- az.pdf(pdf = pdf$metadata$pdf, az=simAz, unc=pdf$metadata$unc, verbose=F)

      if (pdf$metadata$coord == 'dec') {
        simPDF <- coordtrans(simPDF, hor, refraction=pdf$metadata$param$refraction, atm=pdf$metadata$param$atm, temp=pdf$metadata$param$temp, verbose=F, .prec=pdf$metadata$param$.prec)
      }
      res[i,] <- spd(simPDF, xrange=xrange, normalise = normalise)$data$y

      if (save.sim) { sim.az[i,] <- simAz}
      if (verbose) { setTxtProgressBar(pb, i) }
    }
    if (verbose) { cat('Done.\n') }
  }

  ## confidence envelope
  zScore.sim <- matrix(0, nrow=nsims, ncol=length(empirical$data$y))
  zMean <- colMeans(res)
  zStd <- apply(res, 2, sd)

  zScore.emp <- (empirical$data$y - zMean)/zStd
  zScore.sim <- apply(res, 2, scale)

  ## corrections for null models with zeros (to avoid infinity values)
  zScore.emp[!is.finite(zScore.emp)] <- 9999
  zScore.emp[which(zStd == 0 & empirical$data$y < 0.001)] <- NA
  zScore.sim[zScore.sim=='NaN'] <- 0

  if (verbose) { cat(paste0('Performing a ',tails,'-tailed test at the ', conf*100, '% significance level.\n')) }
  if (tails==2) {
    lvl.up <- 1-(1-conf)/2; lvl.dn <- (1-conf)/2
  } else if (tails==1) {
    lvl.up <- conf; lvl.dn <- 0
  } else { stop() }

  upper <- apply(zScore.sim, 2, quantile, probs = lvl.up, na.rm = T)
  upCI <- zMean + upper*zStd
  # upCI[is.na(upCI)] <- 0
  if (tails==2) {
    lower <- apply(zScore.sim, 2, quantile, probs = lvl.dn, na.rm = T)
    loCI <- zMean + lower*zStd
    # loCI[is.na(loCI)] <- 0
  }

  ## global p-value
  above <- which(zScore.emp > upper); emp.stat <- sum(zScore.emp[above] - upper[above])
  if (tails==2) { below <- which(zScore.emp < lower); emp.stat <- emp.stat + sum(lower[below] - zScore.emp[below]) }

  sim.stat <- abs(apply(zScore.sim, 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper))
  if (tails==2) { sim.stat <- sim.stat + abs(apply(zScore.sim, 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower)) }
  global.p <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)

  ## local p-value
  ind <- split(above, cumsum(c(1,diff(above) > 1))); ind <- ind[which(lengths(ind) > 1)]
  local <- data.frame(type=NA, start=0, end=0, p.value=0); j <- 0
  if (length(ind)>0) {
    for (j in 1:NROW(ind)) {
      emp.stat <- sum(zScore.emp[ind[[j]]] - upper[ind[[j]]])
      sim.stat <- abs(apply(zScore.sim[,ind[[j]]], 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper[ind[[j]]]))

      local[j,]$type <- '+'
      local[j,]$start <- min(empirical$data$x[ind[[j]]])
      local[j,]$end <- max(empirical$data$x[ind[[j]]])
      local[j,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
    }
  }

  if (tails==2) {
    ind <- split(below, cumsum(c(1,diff(below) > 1))); ind <- ind[which(lengths(ind) > 1)]
    if (length(ind)>0) {
      for (k in 1:NROW(ind)) {
        emp.stat <- sum(lower[ind[[k]]] - zScore.emp[ind[[k]]])
        sim.stat <- abs(apply(zScore.sim[,ind[[k]]], 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower[ind[[k]]]))

        local[j+k,]$type <- '-'
        local[j+k,]$start <- min(empirical$data$x[ind[[k]]])
        local[j+k,]$end <- max(empirical$data$x[ind[[k]]])
        local[j+k,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
      }
    }
  }

  ## cleanup
  rownames(local) <- c()
  aux <- apply(local[,c(2,3)], 1, diff)
  ind <- which(aux <= 5*pdf$metadata$param$.prec + 1000*.Machine$double.eps)
  if (length(ind)>0) { local <- local[-ind,] }
  rownames(local) <- c()

  ## output
  out <- c()
  out$metadata$coord <- pdf$metadata$coord
  out$metadata$nsims <- nsims
  out$metadata$conf <- conf
  out$metadata$tails <- tails
  out$metadata$normalise <- normalise
  out$metadata$global.pval <- global.p
  out$metadata$local.pval <- local

  out$data$empirical <- empirical$data
  out$data$null.hyp <- list(x = empirical$data$x, CE.mean = zMean, CE.upper = upCI)
  if (tails==2) { out$data$null.hyp$CE.lower <- loCI }
  if (save.sim) { out$data$simulated <- sim.az}
  class(out) <- 'skyscapeR.sigTest'

  return(out)
}

#' Significance test against a null model
#'
#' @param pdf A \emph{skyscapeR.pdf} object created with either \code{\link{az.pdf}} or \code{\link{coordtrans}}
#' @param model A \emph{data.frame} or list of \emph{data.frames} with equal number to the empirical observations.
#' See example for structure
#' @param nsims (Optional) Number of simulations to run. Default is 1000.
#' @param conf (Optional) Confidence envelope (in percentage) of the null model to calculate. Default is .95
#' @param tails (Optional) Whether to calculate 1-tailed p-value (greater than) or 2-tailed p-value (smaller than or greater than).
#' Default is 2.
#' @param normalise (Optional) Boolean to control whether to normalize SPDs. Default is FALSE
#' @param ncores (Optional) Number of CPU cores to use. Default is the number of available cores minus 1.
#' @param verbose (Optional) Boolean to control whether or not to display text. Default is TRUE.
#' @param .res (Optional) Resolution with which to output probability distribution(s). Default is 0.1 degrees.
#' @import snow doSNOW foreach
#' @export
modelTest <- compiler::cmpfun(modelTest_unc, options=list(enableJIT = 3))
