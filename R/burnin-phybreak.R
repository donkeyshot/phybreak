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
#' @param keepphylo The proportion of transmission tree updates keeping the phylotree intact.
#' @return The \code{phybreak}-object provided as input, with variables and parameters changed due to the updating.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#'
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 50)
#' @export
burnin.phybreak <- function(phybreak.object, ncycles, keepphylo = 0.2, phylotopology_only = 0) {
  ### tests
  if(ncycles < 1) stop("ncycles should be positive")
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(phylotopology_only < 0 | phylotopology_only > 1) stop("phylotopology_only should be a fraction")
  if(phylotopology_only + keepphylo > 1) stop("keepphylo + phylotopology_only should be a fraction")
  if(length(phybreak.object$d$nSNPs) > 1 & keepphylo > 0 ) {
    warning("Multiple genes detected, keepphylo set to zero!")
    keepphylo <- 0
  }
  Ngenes <- dim(phybreak.object$v$nodetimes)[1]
  .build.pbe(phybreak.object)
  curtime <- Sys.time()
  
  for (rep in 1:ncycles) {
    if(Sys.time() - curtime > 10) {
      mcmcdiagnostics(.pbe0,rep,Ngenes)
      curtime <- Sys.time()
    }
    for (i in sample(phybreak.object$p$obs)) {
      if (runif(1) < 1 - keepphylo)                  # If this is true, a new phylogeny + transmission tree is chosen
        .updatehost(i) else .updatehost.keepphylo(i) # Else a new trans tree is chosen while keeping phylo.
    }
    if (phybreak.object$h$est.mG)
      .update.mG()
    if (phybreak.object$h$est.mS)
      .update.mS()
    if (phybreak.object$h$est.wh)
      .update.wh()
    if (phybreak.object$h$est.reass) 
      .update.reass()     
    .update.mu()
  }
  
  res <- .destroy.pbe(phybreak.object$s)
  
  return(res)
}

## Function for displaying mcmc diagnostics given environment .pbe0, current mcmc cycle and the number of genes
# Called by burnin-phybreak.R and sample-phybreak.R
mcmcdiagnostics <- function(envir, cycle, Ngenes){
  cat(paste0("cycle ", cycle, ": logLik = ",
             round(sum(envir$logLikseq) + envir$logLiksam + envir$logLikgen + envir$logLikcoal, 2),
             "; mu = ", signif(envir$p$mu, 3),
             "; mean.gen = ", signif(envir$p$mean.gen, 3),
             "; mean.sample = ", signif(envir$p$mean.sample, 3), "\n"))
  for(Gene in 1:Ngenes){cat(paste0("parsimony gene", Gene),
                            paste0("= ", phangorn::parsimony(phybreak2phylo(envir$v,gene = Gene), envir$d$sequences[[Gene]])),
                            "                 ", paste0("nSNPs gene ", Gene), paste0("= ", envir$d$nSNPs[Gene]),
                            "\n")}; cat(paste0("\n"))
}
