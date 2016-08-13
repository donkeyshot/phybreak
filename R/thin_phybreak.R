### thin samples in phybreak.object ###


#' Remove posterior samples from a phybreak-object.
#' 
#' The function removes (all or some) posterior samples by thinning and/or removing the first part of the chain, 
#'   but keeps the state of variables and parameters intact.
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param thin The thinning interval to use.
#' @param nkeep the number of most recent samples to keep. If \code{nkeep = Inf} (default), the whole chain is
#'   thinned and returned. Otherwise, samples from the tail are kept.
#' @return The \code{phybreak}-object provided as input, but with fewer posterior samples.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references Klinkenberg et al, in prep.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 20)
#' MCMCstate <- phybreak(data = simulation)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 50)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 100, thin = 10)
#' 
#' MCMCstate <- thin.phybreak(MCMCstate, thin = 2)
#' @export
thin.phybreak <- function(phybreak.object, thin = 1, nkeep = Inf) {
    ### tests
    if(nkeep * thin > length(phybreak.object$s$logLik & nkeep < Inf)) {
      stop("'nkeep * thin' should not be larger than the current chain length (unless 'nkeep = Inf')")
    }
  
    if (!inherits(phybreak.object, "phybreak")) 
        stop("object should be of class \"phybreak\"")
    if (nkeep == 0) {
        return(within(phybreak.object, s <- list(nodetimes = c(), nodehosts = c(), nodeparents = c(), mu = c(), mG = c(), mS = c(), 
            slope = c(), logLik = c())))
    }
    tokeep <- seq(thin, length(phybreak.object$s$logLik), thin)
    tokeep <- tail(tokeep, nkeep)
    return(within(phybreak.object, s <- list(nodetimes = s$nodetimes[, tokeep], nodehosts = s$nodehosts[, tokeep], nodeparents = s$nodeparents[, 
        tokeep], mu = s$mu[tokeep], mG = s$mG[tokeep], mS = s$mS[tokeep], slope = s$slope[tokeep], logLik = s$logLik[tokeep])))
}
