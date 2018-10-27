#' Perform a null hypothesis significance test of a given curvigram
#'
#' This function performs a null hypothesis significance test, for a given dataset
#' and outputs a p-value as well as all the information needed for  plotting.
#' @param data Data frame including latitude, azimuth and horizon altitude of sites, created
#' using \code{\link{reduct.compass}} or \code{\link{reduct.theodolite}}, plus the uncertainty in
#' azimuthal measurement, which should be added into a column named \emph{Azimuth.Uncertainty}.
#' @param type (Optional) Type of data visualizarion you want to conduct the test on. Current
#' options are 'curv' for \code{\link{curvigram}} and 'hist' for \code{\link{histogram}}. Defaults to 'curv' for curvigram.
#' @param ncores (Optional) Number of processing cores to use for parallelisation. Defaults to the number of
#' available cores minus 1.
#' @param nsims (Optional) Number of simulations to run. The higher this number the slower this process will
#' be, but the lower it is the less power the method has. Defaults to 2000 as a base minimum to test for
#' significance at the p=0.0005 level, but the recommended value is 10,000.
#' @param conf (Optional) Confidence level for p-value calculation and output. Defaults to 0.95,
#' i.e. a 95\% confidence envelope.
#' @param prec (Optional) Smallest possible azimuth for the random sampler, i.e. precision being considered.
#' Defaults to 0.01º.
#' @param range (Optional) Range of declination values to consider.
#' @param verbose (Optional) Boolean to decide whether or not
#' @param ... Other parameters to be passed on to \code{\link{histogram}} or \code{\link{curvigram}}, as appropriate.
#' @export
#' @import parallel foreach doParallel
#' @importFrom rootSolve uniroot.all
#' @seealso \code{\link{reduct.compass}}, \code{\link{reduct.theodolite}}
#' @examples
#' \dontrun{
#' loc <- c(35,-7)
#' mag.az <- c(89.5, 105, 109.5)
#' data <- reduct.compass(loc, mag.az, "2016/04/02", alt=c(1,2,0))
#' data$Azimuth.Uncertainty <- 1  # adds the information on the preision of the azimuthal meaurement
#' sig <- sigTest(data)
#'
#' plot(sig)
#' }
sigTest <- function(data, type='curv', ncores=parallel::detectCores()-1, nsims=2000, conf=.95, prec=.01, range, verbose=T, ...) {

  az <- data$True.Azimuth
  az.unc <- data$Azimuth.Uncertainty
  lat <- data$Latitude
  alt <- data$Altitude

  decs <- az2dec(az, lat, alt)
  unc <- pmax(abs(az2dec(az-az.unc, lat, alt) - az2dec(az, lat, alt)) , abs(az2dec(az, lat, alt) - az2dec(az+az.unc, lat, alt)))

  if (type=='hist') {
    ff <- histogram
    tt <- 'Histogram'
  } else if (type=='curv') {
    ff <- curvigram
    tt <- 'Curvigram'
  } else {
    stop('Type not recognized.')
  }

  ## empirical SPD
  if (verbose) { cat(paste0('Creating Empirical ', tt, '...')) }
  if (!missing(range)) {
    empirical <- ff(decs, unc, ...)
  } else {
    empirical <- ff(decs, unc, ...)
    range <- empirical$metadata$range
  }
  norm <- empirical$metadata$norm
  if (verbose) { cat('Done.\n') }


  ## bootstrapping
  require(foreach)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, library(skyscapeR))
  if (verbose) { cat(paste0('Running ', nsims,' simulations on ', ncores, ' processing cores. This may take a while...')) }

  res <- matrix(NA, nsims, length(empirical$y))
  res <- foreach::foreach(i = 1:nsims, .combine=rbind, .inorder = F) %dopar% {
    simAz <- sample(seq(0, 360, prec), length(az), replace=T)
    simUnc <- sample(unc, replace=F)
    simDecs <- az2dec(simAz, lat, alt)
    # c(ks.test(decs, simDecs)$p.value, ff(simDecs, unc, ...)$density)
    ff(simDecs, unc, norm = norm, range = range)$data$density
  }
  parallel::stopCluster(cl)
  if (verbose) { cat('Done.\n') }

  # KS test p-values
  # ks <- as.numeric(res[,1])
  # res <- res[,-1]

  ## confidence envelope
  zScore.sim <- matrix(0, nrow=nsims, ncol=length(empirical$data$density))
  zMean <- colMeans(res)
  zStd <- apply(res, 2, sd)

  zScore.emp <- (empirical$data$density - zMean)/zStd
  zScore.sim <- apply(res, 2, scale)

  if (verbose) { cat(paste0('Performing a 2-tailed test at the ', conf*100, '% significance level.\n')) }
  lvl.up <- 1-(1-conf)/2; lvl.dn <- (1-conf)/2
  upper <- apply(zScore.sim, 2, quantile, probs = lvl.up, na.rm = T)
  upCI <- zMean + upper*zStd
  upCI[is.na(upCI)] <- 0
  lower <- apply(zScore.sim, 2, quantile, probs = lvl.dn, na.rm = T)
  loCI <- zMean + lower*zStd
  loCI[is.na(loCI)] <- 0

  ## global p-value
  above <- which(zScore.emp > upper); emp.stat <- sum(zScore.emp[above] - upper[above])
  below <- which(zScore.emp < lower); emp.stat <- emp.stat + sum(lower[below] - zScore.emp[below])

  sim.stat <- abs(apply(zScore.sim, 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper))
  sim.stat <- sim.stat + abs(apply(zScore.sim, 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower))
  global.p <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)

  ## local p-value
  ind <- split(above, cumsum(c(1,diff(above) > 1))); ind <- ind[which(lengths(ind) > 1)]
  local <- data.frame(startDec=0, endDec=0, p.value=0, type=NA)
  if (length(ind)>0) {
    for (j in 1:NROW(ind)) {
      emp.stat <- sum(zScore.emp[ind[[j]]] - upper[ind[[j]]])

      sim.stat <- abs(apply(zScore.sim[,ind[[j]]], 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper[ind[[j]]]))
      # sim.stat <- sim.stat + abs(apply(zScore.sim[,ind[[j]]], 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower[ind[[j]]]))

      local[j,]$startDec <- min(empirical$data$dec[ind[[j]]])
      local[j,]$endDec <- max(empirical$data$dec[ind[[j]]])
      local[j,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
      local[j,]$type <- '+'
    }
  }

  ind <- split(below, cumsum(c(1,diff(below) > 1))); ind <- ind[which(lengths(ind) > 1)]
  if (length(ind)>0) {
    for (k in 1:NROW(ind)) {
      emp.stat <- sum(lower[ind[[k]]] - zScore.emp[ind[[k]]])

      # sim.stat <- abs(apply(zScore.sim[,ind[[k]]], 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper[ind[[k]]]))
      sim.stat <- abs(apply(zScore.sim[,ind[[k]]], 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower[ind[[k]]]))

      local[j+k,]$startDec <- min(empirical$data$dec[ind[[k]]])
      local[j+k,]$endDec <- max(empirical$data$dec[ind[[k]]])
      local[j+k,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
      local[j+k,]$type <- '-'
    }
  }

  ## output
  out <- c()
  out$metadata$type <- type
  out$metadata$nsims <- nsims
  out$metadata$conf <- conf
  out$metadata$global.p.value <- global.p
  out$metadata$local <- local
  out$data$CE.mean <- zMean
  out$data$CE.upper <- upCI
  out$data$CE.lower <- loCI
  out$data$empirical <- empirical
  # out$data$ks.p.values <- ks
  class(out) <- 'skyscapeR.sigTest'

  return(out)
}
