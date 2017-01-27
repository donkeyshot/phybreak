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
sample.phybreak <- function(phybreak.object, nsample, thin, keepphylo = 0.2,phylotopology_only = 0) {
  ### tests
  if(nsample < 1) stop("nsample should be positive")
  if(thin < 1) stop("thin should be positive")
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(phylotopology_only < 0 | phylotopology_only > 1) stop("phylotopology_only should be a fraction")
  if(phylotopology_only + keepphylo > 1) stop("keepphylo + phylotopology_only should be a fraction")
  if(length(phybreak.object$d$nSNPs)>1 & keepphylo >0 ) {
    warning("Multiple genes detected, keepphylo set to zero!")
    keepphylo <- 0
  }
  
  Ngenes <- dim(phybreak.object$v$nodetimes)[1]
  GeneNames <- paste("gene",1:Ngenes,sep="")
  
  s.post <- list(nodetimes =lapply(1:Ngenes, function(x) with(phybreak.object, cbind(s$nodetimes[[x]], matrix(NA, nrow = 2 * p$obs - 1, ncol = nsample)))),
                 nodehosts = with(phybreak.object, cbind(s$nodehosts, matrix(NA, nrow = 2 * p$obs - 1, ncol = nsample))),
                 nodeparents = lapply(1:Ngenes, function(x) with(phybreak.object, cbind(s$nodeparents[[x]], matrix(NA, nrow = 3 * p$obs - 1, ncol = nsample)))),
                 mu = c(phybreak.object$s$mu, rep(NA, nsample)),
                 mG = c(phybreak.object$s$mG, rep(NA, nsample)),
                 mS = c(phybreak.object$s$mS, rep(NA, nsample)),
                 slope = c(phybreak.object$s$slope, rep(NA, nsample)),
                 logLikseq = with(phybreak.object, cbind(phybreak.object$s$logLikseq, matrix(NA, nrow =Ngenes, ncol = nsample))),
                 rho = c(phybreak.object$s$rho, rep(NA, phybreak.object$h$est.rho*(nsample))),
                 reassortment = with(phybreak.object, cbind(s$reassortment, matrix(rep(rep(NA, h$est.rho*(nsample)), p$obs), nrow = p$obs)))
                )
  names(s.post$nodetimes) <- GeneNames
  names(s.post$nodeparents) <- GeneNames
  
  .build.pbe(phybreak.object)
  
  curtime <- Sys.time()
  
  for (sa in tail(1:length(s.post$mu), nsample)) { 
    # mcmc Diagnostics
    if(Sys.time() - curtime > 10) {
      mcmcdiagnostics(.pbe0,sa,Ngenes)
      curtime <- Sys.time()
    }
    for (rep in 1:thin) {
      if (phybreak.object$h$est.rho) .update.rho()   # Set hosts in which gene reassortment possibly occured.
      for (i in sample(phybreak.object$p$obs)) {
        if (runif(1) < 1 - keepphylo - phylotopology_only) {
          .updatehost(i) 
        } else if (runif(1) < keepphylo/(keepphylo + phylotopology_only)) {
          .updatehost.keepphylo(i)
        } else {
          .updatehost.topology(i)
        }
      }
      if (phybreak.object$h$est.mG)
        .update.mG()
      if (phybreak.object$h$est.mS)
        .update.mS()
      if (phybreak.object$h$est.wh)
        .update.wh()
      .update.mu()
    }
    
    for (gene in 1:Ngenes) {
      s.post$nodetimes[[gene]][, sa]  <- tail(.pbe0$v$nodetimes[gene, ], -phybreak.object$p$obs)
      s.post$nodeparents[[gene]][, sa] <- .pbe0$v$nodeparents[gene, ]
      if (gene == 1){
        s.post$mu[sa] <- .pbe0$p$mu
        s.post$logLikseq[, sa] <- .pbe0$logLikseq + .pbe0$logLiksam  +.pbe0$logLikgen + .pbe0$logLikcoal
        s.post$nodehosts[, sa] <- tail(.pbe0$v$nodehosts, -phybreak.object$p$obs)
        s.post$mG[sa] <- .pbe0$p$mean.gen
        s.post$mS[sa] <- .pbe0$p$mean.sample
        s.post$slope[sa] <- .pbe0$p$wh.slope
        if(Ngenes > 1 & .pbe0$h$est.rho){
          s.post$rho[sa] <- .pbe0$p$rho
          s.post$reassortment[, sa] <- .pbe0$v$reassortment
        }
      }
    }
  }
  
  res <- .destroy.pbe(s.post)
  return(res)
}