### run mcmc chain and return the updated phylo-object ###

### phybreak functions called ### .build.pbe .updatehost .updatehost.keepphylo .update.mG .update.mS .update.wh
### .update.mu .destroy.pbe

#' MCMC updating of a phybreak-object.
#' 
#' This function allows the MCMC chain to burn in. If used after samples have been taken (with \code{\link{sample.phybreak}}), 
#'   these samples will be returned unchanged in the output.
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param ncycles Number of iterations to be carried out. Each iteration does one update of all parameters and
#'   tree updates with each host as focal host once.
#' @param keepphylo The proportion of tree updates keeping the phylotree intact. If there is more than one
#'   sample per host, keepphylo should be 0. If set to NULL (default), this is done automatically, otherwise it is set to 0.2.
#' @param phylotopology_only The proportion of tree updates in which only the within-host minitree topology is sampled, and 
#'   the transmission tree as well as coalescence times are kept unchanged.
#' @return The \code{phybreak}-object provided as input, with variables and parameters changed due to the updating.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' 
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 50)
#' @export
burnin.phybreak <- function(phybreak.object, ncycles, keepphylo = NULL, phylotopology_only = 0) {
  ### tests
  if(ncycles < 1) stop("ncycles should be positive")
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
  
  .build.pbe(phybreak.object)
  
  message(paste0("   cycle      logLik         mu  gen.mean  sam.mean parsimony (nSNPs = ", .pbe0$d$nSNPs, ")"))
  
  curtime <- Sys.time()
  
  for (rep in 1:ncycles) {
    if(Sys.time() - curtime > 10) {
      message(paste0(
        stringr::str_pad(rep, 8),
        stringr::str_pad(round(.pbe0$logLikseq + .pbe0$logLiksam + .pbe0$logLikgen + .pbe0$logLikcoal, 2), 12),
        stringr::str_pad(signif(.pbe0$p$mu, 3), 11),
        stringr::str_pad(signif(.pbe0$p$mean.gen, 3), 10),
        stringr::str_pad(signif(.pbe0$p$mean.sample, 3), 10),
        stringr::str_pad(phangorn::parsimony(
          phybreak2phylo(.pbe0$v), .pbe0$d$sequences), 10)))
      curtime <- Sys.time()
    }
        

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
  
  res <- .destroy.pbe(phybreak.object$s)
  
  
  return(res)
}
