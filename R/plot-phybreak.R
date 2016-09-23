### plot functions for phybreak-object ###

### makes use of plot function in 'simmap', colouring the branches by host in which they reside, but no unique colour per host.
### calls: transtree get.phylo


#' Plotting a phybreak object.
#' 
#' Plots a \code{phybreak}-object as phylogenetic tree with coloured branches indicating hosts. The default 
#'   is to plot the current state, but any posterior sample can be chosen, as well as various consensus trees.
#' 
#' @param x An object of class \code{phybreak}.
#' @param plot.which Either \code{"sample"} to plot the current state or a selected posterior sample, 
#'   \code{"mpc"} or \code{"mtcc"} to plot a consensus transmission tree (see \code{\link{transtree}}) or \code{"mcc"}
#'   to plot the maximum clade credibility tree (see \code{\link{phylotree}}).
#' @param samplenr If \code{plot.which = "sample"}, this indicates which posterior tree should be plotted: 
#'   \code{samplenr = 0} to plot the current state.
#' @param ... Additional options for \code{\link[phytools]{plotSimmap}}.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' plot(MCMCstate, plot.which = "mpc")
#' @export
plot.phybreak <- function(x, plot.which = c("sample", "mpc", "mtcc", "mcc"), samplenr = 0, ...) {
  if(!("phytools" %in% .packages(TRUE))) {
    stop("package 'phytools' should be installed for this function")
  }
  if (!inherits(x, "phybreak")) {
    stop("object x must be of class \"phybreak\"")
  }
  plot.which <- plot.which[which(plot.which %in% c("sample", "mpc", "mtcc", "mcc"))[1]]
  if(is.na(plot.which)) stop("no valid 'plot.which'")
  if(plot.which == "sample" & samplenr > length(x$s$logLik)) {
    warning("requested 'samplenr' not available; current state used")
  }
  
  if (plot.which == "mpc") {
    phytools::plotSimmap(suppressWarnings(transtree(x, "mpc", phylo.class = TRUE)), colors = setNames(nm = c("black", "red", "blue")), ...)
  } else if (plot.which == "mtcc") {
    phytools::plotSimmap(suppressWarnings(transtree(x, "mtcc", phylo.class = TRUE)), colors = setNames(nm = c("black", "red", "blue")), ...)
  } else if (plot.which == "mcc") {
    phytools::plotSimmap(suppressWarnings(phylotree(x, phylo.class = TRUE)), colors = setNames(nm = c("black", "red", "blue")), ...)
  } else {
    phytools::plotSimmap(suppressWarnings(get.phylo(x, samplenr, TRUE)), colors = setNames(nm = c("black", "red", "blue")), ...)
  }
}
