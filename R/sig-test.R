#' Perform a null hypothesis significance test of a given curvigram
#'
#' This function performs a null hypothesis significance test, for a given curvigram
#' and null hypothesis and outputs a p-value as well as all the information needed for
#' ancillary plotting.
#' @param curv Object of \emph{skyscapeR.curv} format, created using \code{\link{curvigram}}
#' @param null.hyp Object of \emph{skyscapeR.nh} format, created with one of the Null Hypothesis models
#' of \emph{skyscapeR} (see See Also section below).
#' @param level (Optional) Level of confidence for p-value calculation and output. Defaults to 0.95,
#' i.e. a 95\% confidence envelope.
#' @param type (Optional) Whether the test is to be '1-tailed' or '2-tailed'. Defaults to '2-tailed'.
#' @param nsims (Optional) Number of simulations to run. The higher this number the slower this process will
#' be, but the lower it is the less power the method has. Defaults to 2000 as a base minimum to test for
#' significance at the p=0.0005 level, but the recommended value is 10,000.
#' @param ncores (Optional) Number of processing cores to use for parallelisation. Defaults to the number of
#' available cores minus 1.
#' @export
#' @import parallel foreach numDeriv doParallel
#' @importFrom rootSolve uniroot.all
#' @seealso \code{\link{nh.Uniform}}, \code{\link{nh.SummerFM}}, \code{\link{plotCurv}}, \code{\link{plotZscore}}
#' @examples
#' \dontrun{
#' data(RugglesRSC)
#' curv <- curvigram(RugglesRSC$Dec, sd=2)
#' sig <- sigTest(curv, null.hyp=nh.Uniform(c(57,2)))
#'
#' plotCurv(curv, signif=sig)
#' }
sigTest = function(curv, null.hyp, level=.95, type='2-tailed', nsims=2000, ncores) {
  requireNamespace('foreach')
  # redo curvigram with fixed range
  N <- length(curv$mes)
  raw <- curv$mes; mes.unc <- curv$mes.unc
  range <- c(min(null.hyp$dec)-10, max(null.hyp$dec)+10) ######  change number of points as well?
  curv <- curvigram(raw, unc=mes.unc, range=range)

  # monte carlo resampling
  if (missing(ncores)) { ncores <-  parallel::detectCores()-1 }
  cl <- parallel::makeCluster(ncores, type = "PSOCK")
  parallel::clusterEvalQ(cl, library(skyscapeR))
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers()
  message(paste0('Running calculations on ', ncores, ' processing cores. This may take a while...'))

  res <- foreach (i = 1:(1.2*nsims), .combine=rbind, .inorder = F, .errorhandling = 'remove') %dopar% {
    simData <- sample(null.hyp$dec, N, prob=null.hyp$density, replace=T)
    simCurv <- curvigram(simData, unc=mes.unc, range=range)

    simCurv$density
  }
  parallel::stopCluster(cl)

  if (NROW(res) > nsims) { res <- res[1:nsims,] }
  if (NROW(res) < nsims) { stop('oops')  }

  # z-score transformation
  zScore.sim <- matrix(0, nrow=nsims, ncol=length(curv$density))
  zMean <- colMeans(res)
  zStd <- apply(res, 2, sd)

  zScore.emp <- (curv$density - zMean)/zStd
  for (i in 1:nsims) {
    zScore.sim[i,] <- (res[i,] - zMean)/zStd
  }

  # confidence interval
  if (type=='2-tailed') {
    lvl.up <- 1-(1-level)/2
    lvl.dn <- (1-level)/2
    message(paste0('Performing a 2-tailed test at the ',level*100, '% significance level.'))
  } else if (type=='1-tailed') {
    lvl.up <- level
    lvl.dn <- 0
    message(paste0('Performing a 1-tailed test at the ',level*100, '% significance level.'))
  }
  upper <- apply(zScore.sim, 2, quantile, lvl.up, na.rm=T)
  lower <- apply(zScore.sim, 2, quantile, lvl.dn, na.rm=T)
  upCI <- zMean + upper*zStd
  loCI <- zMean + lower*zStd

  # overall p-value
  area.up <- zScore.emp - upper; area.up[area.up<0] <- NA; stat.emp <- sum(area.up, na.rm=T)
  if (type=='2-tailed') { area.dn <- lower - zScore.emp ; area.dn[area.dn<0] <- NA; stat.emp <- stat.emp + sum(area.dn, na.rm=T)}

  np.sims <- 0
  for (i in 1:nsims) {
    area.up <- zScore.sim[i,] - upper; area.up[area.up<0] <- NA; stat.sim <- sum(area.up, na.rm=T)
    if (type=='2-tailed') { area.dn <- lower - zScore.sim[i,] ; area.dn[area.dn<0] <- NA; stat.sim <- stat.sim + sum(area.dn, na.rm=T)}
    if (stat.sim >= stat.emp) {np.sims <- np.sims + 1}
  }
  pval <- np.sims / nsims

  # zScores of peaks
  func0 <- splinefun(curv$dec, zScore.emp)
  fd <- numDeriv::grad(func0,curv$dec); func1 <- splinefun(curv$dec,fd)
  sd <- numDeriv::grad(func1,curv$dec); func2 <- splinefun(curv$dec,sd)
  roots <- rootSolve::uniroot.all(func1, interval= range)
  ind <- which(func2(roots) < 0); maxima <- roots[ind]
  ind <- which(is.finite(zScore.emp))
  funcZ <- splinefun(curv$dec[ind], zScore.emp[ind])
  sigma <- funcZ(maxima)

  datarange <- c(min(curv$mes) - 4*max(curv$mes.unc), max(curv$mes) + 4*max(curv$mes.unc))
  ind <- which(maxima <= datarange[1] | maxima >= datarange[2]); if (length(ind) > 0) { maxima <- maxima[-ind]; sigma <- sigma[-ind] }

  # output
  out <- c()
  out$p.value <- pval
  out$level <- level
  out$type <- type
  out$maxima <- rbind(maxima,sigma); rownames(out$maxima) <- c('dec','zScore')
  out$nsims <- nsims
  out$null.hyp <- rbind(curv$dec,loCI,zMean,upCI); rownames(out$null.hyp) <- c('dec',paste0(as.character(lvl.dn*100),'%'),'mean',paste0(as.character(lvl.up*100),'%'))
  out$null.hyp.z <- rbind(curv$dec,lower,zScore.emp,upper); rownames(out$null.hyp.z) <- c('dec',paste0(as.character(lvl.dn*100),'%'),'zScore',paste0(as.character(lvl.up*100),'%'))
  out$data.range <- datarange
  class(out) <- 'skyscapeR.sig'
  return(out)
}




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
#' @seealso \code{\link{sigTest}}, \code{\link{nh.SummerFM}}, \code{\link{nh.SolarRange}}
#' @examples
#' \dontrun{
#' aux <- nh.Uniform(loc=c(52,-2), alt=2)
#' plot(aux$dec, aux$density, type='l')
#' }
nh.Uniform = function(loc, alt=0) {
  if (class(loc) == 'skyscapeR.horizon') { loc <- loc$georef }

  xx <- seq(0, 360, .01)
  ddd <- sapply(xx, az2dec, loc, alt)
  aux <- density(ddd, bw=0.1)

  out <- c()
  out$type <- 'nh.Uniform'
  out$param$loc <- loc
  out$param$alt <- alt
  out$dec <- aux$x
  out$density <- aux$y
  class(out) <- 'skyscapeR.nh'

  return(out)
}


#' Declination distribution of the Sun throughout the year
#'  for null hypothesis significance testing
#'
#' This function returns the declination distribution of
#' the the Sun throughout the year for use in significance testing.
#' @param year Year for which to calculate the distribution.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{sigTest}}, \code{\link{nh.Uniform}}, \code{\link{nh.SummerFM}}
#' @examples
#' \dontrun{
#' aux <- nh.SolarRange(-4000)
#' plot(aux$dec, aux$density, type='l')
#' }
nh.SolarRange = function(year = cur.year) {
  options(warn=-1)
  jd0 <- astrolibR::juldate(c(year-100,1,1,12)) + 2400000
  xx <- seq(jd0, jd0+365*50,1)
  sundec <- astrolibR::sunpos(xx)$dec
  aux <- density(sundec, 0.1)
  options(warn=0)

  out <- c()
  out$type <- 'nh.SolarRange'
  out$param$year <- year
  out$dec <- aux$x
  out$density <- aux$y
  class(out) <- 'skyscapeR.nh'

  return(out)
}


#' Declination distribution of the Moon throughout the year
#'  for null hypothesis significance testing
#'
#' This function returns the declination distribution of
#' the the Sun throughout the year for use in significance testing.
#' @param year Year for which to calculate the distribution.
#' Defaults to present year as given by Sys.Date()
#' @export
#' @seealso \code{\link{sigTest}}, \code{\link{nh.Uniform}}, \code{\link{nh.SummerFM}}
#' @examples
#' \dontrun{
#' aux <- nh.SolarRange(-4000)
#' plot(aux$dec, aux$density, type='l')
#' }
nh.LunarRange = function(year = cur.year) {
  options(warn=-1)
  jd0 <- astrolibR::juldate(c(year-100,1,1,12)) + 2400000
  xx <- seq(jd0, jd0+365*200,1)
  moondec <- astrolibR::moonpos(xx)$dec
  aux <- density(moondec)
  options(warn=0)

  out <- c()
  out$type <- 'nh.MoonRange'
  out$param$year <- year
  out$dec <- aux$x
  out$density <- aux$y
  class(out) <- 'skyscapeR.nh'

  return(out)
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
#' @seealso \code{\link{sigTest}}, \code{\link{nh.Uniform}}, \code{\link{nh.SolarRange}}
#' @examples
#' \dontrun{
#' aux <- nh.SummerFM(.99, 20, -4000)
#' plot(aux$dec, aux$density, type='l')
#' }
nh.SummerFM = function(min.phase = .99, min.sundec = 20, year = cur.year) {
  options(warn=-1)
  jd0 <- astrolibR::juldate(c(year-100,1,1,12)) + 2400000
  xx <- seq(jd0, jd0+365*200,1)
  mp <- astrolibR::mphase(xx)
  sundec <- astrolibR::sunpos(xx)$dec
  ind <- which(mp >= min.phase & sundec >= min.sundec)
  dec <- astrolibR::moonpos(xx[ind])$dec
  # diff <- obliquity(1900)- obliquity(year); dec <- dec + diff        # another way to time-shift
  aux <- density(dec)
  options(warn=0)


  out <- c()
  out$type <- 'nh.SummerFM'
  out$param$min.phase <- min.phase
  out$param$min.sundec <- min.sundec
  out$param$year <- year
  out$dec <- aux$x
  out$density <- aux$y
  class(out) <- 'skyscapeR.nh'

  return(out)
}


