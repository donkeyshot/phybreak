#' Quantiles of sampled parameter values.
#' @param whichpars Which parameters to return. Either a vector with parameter names, or \code{"all"} for all parameters, or 
#'   \code{"posterior"} for parameters for which a posterior is sampled.
#' @param probs Numeric vector of probabilities with values in \eqn{[0,1]}.
#' @param type An integer between 1 and 9 selecting one of the nine quantile algorithms (see \code{\link[stats]{quantile}}).
#' @return A named matrix with desired quantiles. For parameters that are not estimated, all quantiles give the fixed value.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 100)
#' 
#' quantile(MCMCstate)
##' @export
quantile.phybreak <- function(x, whichpars = "posterior", probs = c(0.025, 0.5, 0.975), type = 7, ...) {
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  if(whichpars == "posterior") {
    whichpars <- c("mu", 
                   if(x$h$est.mG) "gen.mean",
                   if(x$h$est.mS) "sample.mean",
                   if(x$h$est.wh.s) "wh.slope",
                   if(x$h$est.wh.e) "wh.exponent",
                   if(x$h$est.wh.0) "wh.level",
                   if(x$h$est.dist.e) "dist.exponent",
                   if(x$h$est.dist.s) "dist.scale",
                   if(x$h$est.dist.m) "dist.mean")
  } else if(whichpars == "all") {
    whichpars <- setdiff(names(x$p), c("obs", "wh.model"))
  } else if(!all(whichpars %in% names(x$p))) {
    warning("some parameters in whichpars don't exist and are left out")
    whichpars <- whichpars[whichpars %in% names(x$p)]
  }
  
  pars <- as.data.frame(x$p)[whichpars][1, ]
  pars <- matrix(pars, nrow = length(whichpars), ncol = length(probs), 
                 dimnames = list(whichpars, paste0("Q", 100 * probs)))
  if("mu" %in% whichpars) pars["mu", ] <- quantile(x$s$mu, probs = probs, type = type)
  if("sample.mean" %in% whichpars) pars["sample.mean", ] <- quantile(x$s$mS, probs = probs, type = type)
  if("gen.mean" %in% whichpars) pars["gen.mean", ] <- quantile(x$s$mG, probs = probs, type = type)
  if("wh.slope" %in% whichpars) pars["wh.slope", ] <- quantile(x$s$wh.s, probs = probs, type = type)
  if("wh.exponent" %in% whichpars) pars["wh.exponent", ] <- quantile(x$s$wh.e, probs = probs, type = type)
  if("wh.level" %in% whichpars) pars["wh.level", ] <- quantile(x$s$wh.0, probs = probs, type = type)
  if("dist.exponent" %in% whichpars) pars["dist.exponent", ] <- quantile(x$s$dist.e, probs = probs, type = type)
  if("dist.scale" %in% whichpars) pars["dist.scale", ] <- quantile(x$s$dist.s, probs = probs, type = type)
  if("dist.mean" %in% whichpars) pars["dist.mean", ] <- quantile(x$s$dist.m, probs = probs, type = type)
  
  return(pars)
}  
