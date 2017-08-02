### plot functions for phybreak-object ###

### makes use of plot function in 'simmap', colouring the branches by host in which they reside, but no unique colour per host.
### calls: transtree get.phylo


#' Plotting a phybreak object.
#' 
#' Plots a \code{phybreak}-object twice: (1) as transmission tree and (2) as phylogenetic tree, using the default graphical parameters
#'   of \code{\link{plotTrans}} and \code{\link{plotPhylo}}. The default 
#'   is to plot the current state, but any posterior sample can be chosen, as well as various consensus trees. Consensus tree "edmonds"
#'   plots only a transmission tree, consensus tree "mcc" only a phylogenetic tree.
#' 
#' @param x An object of class \code{phybreak}.
#' @param plot.which Either \code{"sample"} to plot the current state or a selected posterior sample, 
#'   \code{"mpc"} or \code{"mtcc"} to plot a consensus transmission tree (see \code{\link{transtree}}) or \code{"mcc"}
#'   to plot the maximum clade credibility tree (see \code{\link{phylotree}}).
#' @param samplenr If \code{plot.which = "sample"}, this indicates which posterior tree should be plotted: 
#'   \code{samplenr = 0} to plot the current state.
#' @param ... Some methods for this generic require additional arguments. None are used in this method.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' plot(MCMCstate, plot.which = "mpc")
#' @export
plot.phybreak <- function(x, plot.which = c("sample", "edmonds", "mpc", "mtcc", "mcc"), samplenr = 0, ...) {

  plot.which <- plot.which[which(plot.which %in% c("sample", "edmonds", "mpc", "mtcc", "mcc"))[1]]
  if(is.na(plot.which)) stop("no valid 'plot.which'")
  if(plot.which == "sample" & samplenr > length(x$s$logLik)) {
    warning("requested 'samplenr' not available; current state used")
    samplenr <- 0
  }
  
  devAskNewPage(TRUE)
  on.exit(expr = devAskNewPage(FALSE))
  
  if (plot.which != "mcc") {
    plotTrans(x, plot.which, samplenr)
  } 
  if(plot.which != "edmonds") {
    plotPhylo(x, plot.which, samplenr)
  }
  
}
