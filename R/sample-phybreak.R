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
#' @param keepphylo The proportion of transmission tree updates keeping the phylotree intact.
#' @return The \code{phybreak}-object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' 
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample.phybreak <- function(phybreak.object, nsample, thin, keepphylo = 0.2, phylotopology_only = 0) {
  ### tests
  if(nsample < 1) stop("nsample should be positive")
  if(thin < 1) stop("thin should be positive")
  if((phybreak.object$h$est.reass | phybreak.object$p$reass.prob > 0) & keepphylo > 0 ) {
    warning("reassortment allowed: keepphylo = 0", immediate. = TRUE)
    keepphylo <- 0
  }
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(phylotopology_only < 0 | phylotopology_only > 1) stop("phylotopology_only should be a fraction")
  if(phylotopology_only + keepphylo > 1) stop("keepphylo + phylotopology_only should be a fraction")
  

  s.post <- with(phybreak.object,
                 list(
                   nodetimes = abind::abind(s$nodetimes, array(NA, dim = c(d$ngenes, 2 * p$obs - 1, nsample))),
                   nodehosts = cbind(s$nodehosts, matrix(NA, nrow = 2 * p$obs - 1, ncol = nsample)),
                   nodeparents = abind::abind(s$nodeparents, array(NA, dim = c(d$ngenes, 3 * p$obs - 1, nsample))),
                   hostreassortment = cbind(s$hostreassortment, matrix(NA, nrow = p$obs, ncol = nsample)),
                   mu = c(s$mu, rep(NA, nsample)),
                   mG = c(s$mG, rep(NA, nsample)),
                   mS = c(s$mS, rep(NA, nsample)),
                   slope = c(s$slope, rep(NA, nsample)),
                   reass = c(s$reass.prob, rep(NA, nsample)),
                   logLik = c(s$logLik, rep(NA, nsample))
                 )
  )

  .build.pbe(phybreak.object)
  
  curtime <- Sys.time()
  
  for (sa in tail(1:length(s.post$mu), nsample)) { 
    # mcmc Diagnostics
    if(Sys.time() - curtime > 10) {
      cat(paste0("sample ", sa, ": logLik = ", 
                 round(.pbe0$logLikseq + .pbe0$logLiksam + .pbe0$logLikgen + .pbe0$logLikcoal, 2),
                 "; mu = ", signif(.pbe0$p$mu, 3), 
                 "; mean.gen = ", signif(.pbe0$p$mean.gen, 3),
                 "; mean.sample = ", signif(.pbe0$p$mean.sample, 3),
                 "; parsimony ="), sapply(1:.pbe0$d$ngenes, 
                                          function(x) {
                                            phangorn::parsimony(phybreak2phylo(.pbe0$v, gene = x), 
                                                                .pbe0$d$sequences[[x]])
                                            }),
                 "(nSNPs =", .pbe0$d$nSNPs, ")\n"
      )
      curtime <- Sys.time()
    }
    
    for (rep in 1:thin) {
      for (i in sample(phybreak.object$p$obs)) {
        if (runif(1) < 1 - keepphylo - phylotopology_only) {
          .updatehost(i) 
        } else if (runif(1) < keepphylo/(keepphylo + phylotopology_only)) {
          .updatehost.keepphylo(i)
        } else {
          .updatehost.topology(i)
        }
      }
      if(phybreak.object$h$est.mG)
        .update.mG()
      if(phybreak.object$h$est.mS)
        .update.mS()
      if(phybreak.object$h$est.wh)
        .update.wh()
      if(phybreak.object$h$est.reass)
        .update.reass()
      .update.mu()
    }
    
    s.post$nodetimes[,, sa] <- t(tail(t(.pbe0$v$nodetimes), -.pbe0$p$obs))
    s.post$nodehosts[, sa] <- tail(.pbe0$v$nodehosts, -.pbe0$p$obs)
    s.post$nodeparents[,, sa] <- .pbe0$v$nodeparents
    s.post$hostreassortment[, sa] <- .pbe0$v$hostreassortment
    s.post$mu[sa] <- .pbe0$p$mu
    s.post$mG[sa] <- .pbe0$p$mean.gen
    s.post$mS[sa] <- .pbe0$p$mean.sample
    s.post$slope[sa] <- .pbe0$p$wh.slope
    s.post$reass.prob[sa] <- .pbe0$p$reass.prob
    s.post$logLik[sa] <- .pbe0$logLikseq + .pbe0$logLiksam  +.pbe0$logLikgen + .pbe0$logLikcoal
  }
  
  res <- .destroy.pbe(s.post)
  return(res)
}
