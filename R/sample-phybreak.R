### run mcmc chain, take samples from the chain, and return the updated phylo-object ###

### phybreak functions called ### .build.pbe .updatehost .updatehost.keepphylo .update.mG .update.mS .update.wh
### .update.mu .destroy.pbe


#' Sampling from a phybreak MCMC-chain.
#' 
#' Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree 
#'   (plus associated parameters), within a \code{phybreak}-object.
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param nsample The number of samples to take.
#' @param thin The thinning to use (values after every \code{thin}'th iteration will be included in the posterior). 
#'   Each iteration does one update of all parameters and tree updates with each host as focal host once.
#' @param keepphylo The proportion of transmission tree updates keeping the phylotree intact.
#' @return The \code{phybreak}-object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' 
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample.phybreak <- function(phybreak.object, nsample, thin, keepphylo = 0.2) {
    ### tests
    if(nsample < 1) stop("nsample should be positive")
    if(thin < 1) stop("thin should be positive")
    if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  
    ### create room in s to add the new posterior samples
    s.post <- list(nodetimes = with(phybreak.object, cbind(s$nodetimes, matrix(NA, nrow = 2 * p$obs - 1, ncol = nsample))), nodehosts = with(phybreak.object, 
        cbind(s$nodehosts, matrix(NA, nrow = 2 * p$obs - 1, ncol = nsample))), nodeparents = with(phybreak.object, cbind(s$nodeparents, 
        matrix(NA, nrow = 3 * p$obs - 1, ncol = nsample))), mu = c(phybreak.object$s$mu, rep(NA, nsample)), mG = c(phybreak.object$s$mG, 
        rep(NA, nsample)), mS = c(phybreak.object$s$mS, rep(NA, nsample)), slope = c(phybreak.object$s$slope, rep(NA, nsample)), 
        logLikseq = c(phybreak.object$s$logLikseq, rep(NA, nsample)))
    
    .build.pbe(phybreak.object)
    
    for (sa in tail(1:length(s.post$mu), nsample)) {
        for (rep in 1:thin) {
            for (i in sample(phybreak.object$p$obs)) {
                if (runif(1) < 1 - keepphylo) 
                  .updatehost(i) else .updatehost.keepphylo(i)
            }
            if (phybreak.object$h$est.mG) 
                .update.mG()
            if (phybreak.object$h$est.mS) 
                .update.mS()
            if (phybreak.object$h$est.wh) 
                .update.wh()
            .update.mu()
        }
        s.post$nodetimes[, sa] <- tail(.pbe0$v$nodetimes, -phybreak.object$p$obs)
        s.post$nodehosts[, sa] <- tail(.pbe0$v$nodehosts, -phybreak.object$p$obs)
        s.post$nodeparents[, sa] <- .pbe0$v$nodeparents
        s.post$mu[sa] <- .pbe0$p$mu
        s.post$mG[sa] <- .pbe0$p$mean.gen
        s.post$mS[sa] <- .pbe0$p$mean.sample
        s.post$slope[sa] <- .pbe0$p$wh.slope
        s.post$logLikseq[sa] <- .pbe0$logLikseq + .pbe0$logLiksam + +.pbe0$logLikgen + .pbe0$logLikcoal
    }
    
    res <- .destroy.pbe(s.post)
    
    
    return(res)
    
}
