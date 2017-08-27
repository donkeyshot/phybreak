### run mcmc chain, take samples from the chain, and return the updated phylo-object ###



#' Sampling from a phybreak MCMC-chain.
#' 
#' Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree 
#'   (plus associated parameters), within a \code{phybreak} object (\code{sample.phybreak} is deprecated).
#' 
#' @param x An object of class \code{phybreak}.
#' @param nsample The number of samples to take.
#' @param thin The thinning to use (values after every \code{thin}'th iteration will be included in the posterior). 
#'   Each iteration does one update of all parameters and tree updates with each host as focal host once.
#' @param keepphylo The proportion of tree updates keeping the phylotree intact, only possible if there is one sample per host
#'   and the \code{wh.model = "linear"}. If \code{NULL} (default), it is set to 0.2 only in that case, otherwise to 0.
#' @param withinhost_only The proportion of tree updates in which only the within-host minitree is sampled, and 
#'   the transmission tree and infection times are kept unchanged.
#' @return The \code{phybreak} object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample_phybreak <- function(x, nsample, thin = 1, keepphylo = NULL, withinhost_only = 0) {
    ### tests
    if(nsample < 1) stop("nsample should be positive")
    if(thin < 1) stop("thin should be positive")
    if(is.null(keepphylo)) {
      if(any(duplicated(x$d$hostnames)) || !(x$p$wh.model %in% c(3, "linear")) ) {
        keepphylo <- 0
        message("keepphylo = 0")
      } else {
        keepphylo <- 0.2
        message("keepphylo = 0.2")
      }
    }
    if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
    if(withinhost_only < 0 | withinhost_only > 1) stop("withinhost_only should be a fraction")
    if(withinhost_only + keepphylo > 1) stop("keepphylo + withinhost_only should be a fraction")
  
    ### create room in s to add the new posterior samples
    s.post <- list(inftimes = with(x, cbind(s$inftimes, matrix(NA, nrow = p$obs, ncol = nsample))),
                   infectors = with(x, cbind(s$infectors, matrix(NA, nrow = p$obs, ncol = nsample))),
                   nodetimes = with(x, cbind(s$nodetimes, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                   nodehosts = with(x, cbind(s$nodehosts, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                   nodeparents = with(x, cbind(s$nodeparents, matrix(NA, nrow = 2 * d$nsamples - 1, ncol = nsample))), 
                   mu = c(x$s$mu, rep(NA, nsample)), 
                   mG = c(x$s$mG, rep(NA, nsample)), 
                   mS = c(x$s$mS, rep(NA, nsample)), 
                   wh.s = c(x$s$wh.s, rep(NA, nsample)), 
                   wh.e = c(x$s$wh.e, rep(NA, nsample)), 
                   wh.0 = c(x$s$wh.0, rep(NA, nsample)), 
                   logLik = c(x$s$logLik, rep(NA, nsample)))
    
    .build.pbe(x)
    
    message(paste0("  sample      logLik         mu  gen.mean  sam.mean parsimony (nSNPs = ", .pbe0$d$nSNPs, ")"))
    
    curtime <- Sys.time()
    
    for (sa in tail(1:length(s.post$mu), nsample)) {
      
      if(Sys.time() - curtime > 10) {
        message(paste0(
          stringr::str_pad(sa, 8),
          stringr::str_pad(round(.pbe0$logLikseq + .pbe0$logLiksam + .pbe0$logLikgen + .pbe0$logLikcoal, 2), 12),
          stringr::str_pad(signif(.pbe0$p$mu, 3), 11),
          stringr::str_pad(signif(.pbe0$p$gen.mean, 3), 10),
          stringr::str_pad(signif(.pbe0$p$sample.mean, 3), 10),
          stringr::str_pad(phangorn::parsimony(
            phybreak2phylo(.pbe0$v), .pbe0$d$sequences), 10)))
        curtime <- Sys.time()
      }
      
      for (rep in 1:thin) {
        for (i in sample(x$p$obs)) {
          if (runif(1) < 1 - keepphylo - withinhost_only) 
            .updatehost(i) else  if (runif(1) < keepphylo/(keepphylo + withinhost_only)) {
              .updatehost.keepphylo(i)
            } else .updatehost.withinhost(i)
        }
        if (x$h$est.mG) 
          .update.mG()
        if (x$h$est.mS) 
          .update.mS()
        if (x$h$est.wh.s) 
          .update.wh.slope()
        if (x$h$est.wh.e) 
          .update.wh.exponent()
        if (x$h$est.wh.0) 
          .update.wh.level()
        .update.mu()
      }
        s.post$inftimes[, sa] <- .pbe0$v$inftimes
        s.post$infectors[, sa] <- .pbe0$v$infectors
        s.post$nodetimes[, sa] <- .pbe0$v$nodetimes[.pbe0$v$nodetypes == "c"]
        s.post$nodehosts[, sa] <- .pbe0$v$nodehosts[.pbe0$v$nodetypes == "c"]
        s.post$nodeparents[, sa] <- .pbe0$v$nodeparents
        s.post$mu[sa] <- .pbe0$p$mu
        s.post$mG[sa] <- .pbe0$p$gen.mean
        s.post$mS[sa] <- .pbe0$p$sample.mean
        s.post$wh.s[sa] <- .pbe0$p$wh.slope
        s.post$wh.e[sa] <- .pbe0$p$wh.exponent
        s.post$wh.0[sa] <- .pbe0$p$wh.level
        s.post$logLik[sa] <- .pbe0$logLikseq + .pbe0$logLiksam + +.pbe0$logLikgen + .pbe0$logLikcoal
    }
    
    res <- .destroy.pbe(s.post)
    
    
    return(res)
    
}



#' @rdname sample_phybreak
#' @export
sample.phybreak <- function(...) {
  .Deprecated("sample_phybreak")
  sample_phybreak(...)
}

