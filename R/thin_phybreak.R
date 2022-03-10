#' Remove posterior samples from a phybreak-object.
#' 
#' The function removes (all or some) posterior samples by thinning and/or removing the first part of the chain, 
#'   but keeps the state of variables and parameters intact.
#' 
#' @param x An object of class \code{phybreak}.
#' @param thin The thinning interval to use.
#' @param nkeep the number of most recent samples to keep. If \code{nkeep = Inf} (default), the whole chain is
#'   thinned and returned. Otherwise, samples from the tail are kept.
#' @return The \code{phybreak}-object provided as input, but with fewer posterior samples.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' MCMCstate <- thin.phybreak(MCMCstate, thin = 2)
#' @importFrom coda thin
#' @export
thin.phybreak <- function(x, thin = 1, nkeep = Inf, ...) {
    ### tests
    if(nkeep * thin > length(x$s$logLik) & nkeep < Inf) {
      stop("'nkeep * thin' should not be larger than the current chain length (unless 'nkeep = Inf')")
    }
  
    if (!inherits(x, "phybreak")) 
        stop("object should be of class \"phybreak\"")
    if (nkeep == 0) {
        return(within(x, s <- list(
          inftimes = c(),
          infectors = c(),
          nodetimes = c(), 
          nodehosts = c(), 
          nodeparents = c(), 
          mu = c(), 
          mG = c(), 
          mS = c(), 
          wh.s = c(), 
          wh.e = c(), 
          wh.0 = c(), 
          logLik = c())))
    }
    tokeep <- seq(thin, length(x$s$logLik), thin)
    tokeep <- tail(tokeep, nkeep)
    return(within(x, s <- list(
      inftimes = s$inftimes[, tokeep],
      infectors = s$infectors[, tokeep],
      nodetimes = s$nodetimes[, tokeep], 
      nodehosts = s$nodehosts[, tokeep], 
      nodeparents = s$nodeparents[, tokeep],
      introductions = s$introductions[tokeep],
      mu = s$mu[tokeep], 
      mG = s$mG[tokeep], 
      mS = s$mS[tokeep],
      tG = s$tG[tokeep],
      tS = s$tS[tokeep],
      wh.h = s$wh.h[tokeep],
      wh.s = s$wh.s[tokeep], 
      wh.e = s$wh.e[tokeep], 
      wh.0 = s$wh.0[tokeep],
      dist.e = s$dist.e[tokeep],
      dist.s = s$dist.s[tokeep],
      dist.m = s$dist.m[tokeep],
      logLik = s$logLik[tokeep],
      historyinf = s$historyinf[tokeep],
      heat = s$heat[tokeep],
      chain = s$chain[tokeep],
      hist.mean = s$hist.mean[tokeep])))
}
