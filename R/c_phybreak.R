#' Combine mcmc-chains of phybreak-objects.
#' 
#' The function combines the sample slots of its arguments. The variable slot (current state) is taken from the first argument
#' 
#' @param ... Objects of class \code{phybreak} to be concatenated.
#' @return The first \code{phybreak}-object with mcmc samples of all objects.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate1 <- phybreak(dataset = simulation)
#' MCMCstate1 <- burnin_phybreak(MCMCstate1, ncycles = 20)
#' MCMCstate1 <- sample_phybreak(MCMCstate1, nsample = 50, thin = 2)
#' 
#' MCMCstate2 <- phybreak(dataset = simulation)
#' MCMCstate2 <- burnin_phybreak(MCMCstate2, ncycles = 20)
#' MCMCstate2 <- sample_phybreak(MCMCstate2, nsample = 50, thin = 2)
#' 
#' MCMCstate <- c(MCMCstate1, MCMCstate2)
#' @importFrom coda thin
#' @export
c.phybreak <- function(...) {
  obj <- list(...)
  if(!all(sapply(obj, inherits, what = "phybreak"))) {
    stop("all objects should be of class phybreak")
  }
  
  res <- obj[[1]]
  first_with_samples <- 1
  while(first_with_samples <= length(obj) && is.null(res$s$mu)) {
    res$s <- obj[[first_with_samples]]$s
    first_with_samples <- first_with_samples + 1
  }
  
  if(length(obj) <= first_with_samples) return(res)

  for(i in (1 + first_with_samples):length(obj)) {
    for(j in 1:length(res$s)) {
      if(is.null(dim(res$s[[j]]))) {
        res$s[[j]] <- c(res$s[[j]], obj[[i]]$s[[j]])
      } else {
        res$s[[j]] <- cbind(res$s[[j]], obj[[i]]$s[[j]])
      }
    }
  }
  
  return(res)
}