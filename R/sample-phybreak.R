### run mcmc chain, take samples from the chain, and return the updated phylo-object ###



#' Sampling from a phybreak MCMC-chain.
#' 
#' Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree 
#'   (plus associated parameters), within a \code{phybreak}-object.
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param nsample The number of samples to take.
#' @param thin The thinning to use (values after every \code{thin}'th iteration will be included in the posterior). 
#'   Each iteration does one update of all parameters and tree updates with each host as focal host once.
#' @param keepphylo The proportion of tree updates keeping the phylotree intact. If there is more than one
#'   sample per host, keepphylo should be 0. If set to NULL (default), this is done automatically, otherwise it is set to 0.2.
#' @param phylotopology_only The proportion of tree updates in which only the within-host minitree topology is sampled, and 
#'   the transmission tree as well as coalescence times are kept unchanged.
#' @return The \code{phybreak}-object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' 
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample.phybreak <- function(phybreak.object, nsample, thin = 1, keepphylo = NULL, phylotopology_only = 0) {
    ### tests
    if(nsample < 1) stop("nsample should be positive")
    if(thin < 1) stop("thin should be positive")
    if(is.null(keepphylo)) {
      if(any(duplicated(phybreak.object$d$hostnames))) {
        keepphylo <- 0
        message("keepphylo = 0")
      } else {
        keepphylo <- 0.2
        message("keepphylo = 0.2")
      }
    }
    if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
    if(phylotopology_only < 0 | phylotopology_only > 1) stop("phylotopology_only should be a fraction")
    if(phylotopology_only + keepphylo > 1) stop("keepphylo + phylotopology_only should be a fraction")
  
    ### create room in s to add the new posterior samples
    s.post <- list(nodetimes = with(phybreak.object, cbind(s$nodetimes, matrix(NA, nrow = d$nsamples + p$obs - 1, ncol = nsample))), 
                   nodehosts = with(phybreak.object, cbind(s$nodehosts, matrix(NA, nrow = d$nsamples + p$obs - 1, ncol = nsample))), 
                   nodeparents = with(phybreak.object, cbind(s$nodeparents, matrix(NA, nrow = 2 * d$nsamples + p$obs - 1, ncol = nsample))), 
                   mu = c(phybreak.object$s$mu, rep(NA, nsample)), 
                   mG = c(phybreak.object$s$mG, rep(NA, nsample)), 
                   mS = c(phybreak.object$s$mS, rep(NA, nsample)), 
                   slope = c(phybreak.object$s$slope, rep(NA, nsample)), 
                   logLik = c(phybreak.object$s$logLik, rep(NA, nsample)))
    
    .build.pbe(phybreak.object)
    
    message(paste0("  sample      logLik         mu  gen.mean  sam.mean parsimony (nSNPs = ", .pbe0$d$nSNPs, ")"))
    
    curtime <- Sys.time()
    
    for (sa in tail(1:length(s.post$mu), nsample)) {
      
      if(Sys.time() - curtime > 10) {
        message(paste0(
          stringr::str_pad(sa, 8),
          stringr::str_pad(round(.pbe0$logLikseq + .pbe0$logLiksam + .pbe0$logLikgen + .pbe0$logLikcoal, 2), 12),
          stringr::str_pad(signif(.pbe0$p$mu, 3), 11),
          stringr::str_pad(signif(.pbe0$p$mean.gen, 3), 10),
          stringr::str_pad(signif(.pbe0$p$mean.sample, 3), 10),
          stringr::str_pad(phangorn::parsimony(
            phybreak2phylo(.pbe0$v), .pbe0$d$sequences), 10)))
        curtime <- Sys.time()
      }
      
      for (rep in 1:thin) {
            for (i in sample(phybreak.object$p$obs)) {
              if (runif(1) < 1 - keepphylo - phylotopology_only) 
                .updatehost(i) else  if (runif(1) < keepphylo/(keepphylo + phylotopology_only)) {
                  .updatehost.keepphylo(i)
                } else .updatehost.topology(i)
            }
            if (phybreak.object$h$est.mG) 
                .update.mG()
            if (phybreak.object$h$est.mS) 
                .update.mS()
            if (phybreak.object$h$est.wh) 
                .update.wh()
            .update.mu()
        }
        s.post$nodetimes[, sa] <- .pbe0$v$nodetimes[.pbe0$v$nodetypes %in% c("c", "t")]
        s.post$nodehosts[, sa] <- .pbe0$v$nodehosts[.pbe0$v$nodetypes %in% c("c", "t")]
        s.post$nodeparents[, sa] <- .pbe0$v$nodeparents
        s.post$mu[sa] <- .pbe0$p$mu
        s.post$mG[sa] <- .pbe0$p$mean.gen
        s.post$mS[sa] <- .pbe0$p$mean.sample
        s.post$slope[sa] <- .pbe0$p$wh.slope
        s.post$logLik[sa] <- .pbe0$logLikseq + .pbe0$logLiksam + +.pbe0$logLikgen + .pbe0$logLikcoal
    }
    
    res <- .destroy.pbe(s.post)
    
    
    return(res)
    
}
