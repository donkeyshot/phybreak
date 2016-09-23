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
#' @param keepphylo The proportion of transmission tree updates keeping the phylotree intact.
#' @return The \code{phybreak}-object provided as input, with variables and parameters changed due to the updating.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' 
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 50)
#' @export
burnin.phybreak <- function(phybreak.object, ncycles, keepphylo = 0.2) {
    ### tests
    if(ncycles < 1) stop("ncycles should be positive")
    if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  
    .build.pbe(phybreak.object)
    
    
    for (rep in 1:ncycles) {
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
    
    res <- .destroy.pbe(phybreak.object$s)
    
    
    return(res)
}
