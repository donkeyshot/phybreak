### plot functions for phybreakdata-object ###

### makes use of plot function in 'simmap', colouring the branches by host in which they reside, but no unique colour per host.
### calls: transtree get.phylo


#' Plotting a phybreakdata object.
#' 
#' Plots a \code{phybreakdata}-object twice: (1) as transmission tree and (2) as phylogenetic tree, using the default graphical parameters
#'   of \code{\link{plotTrans}} and \code{\link{plotPhylo}}. The default 
#'   is to plot the current state, but any posterior sample can be chosen, as well as various consensus trees. Consensus tree "edmonds"
#'   plots only a transmission tree, consensus tree "mcc" only a phylogenetic tree.
#' 
#' @param x An object of class \code{phybreakdata}.
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
plot.phybreakdata <- function(x, ...) {


  devAskNewPage(TRUE)
  on.exit(expr = devAskNewPage(FALSE))
  plotTrans(x, plot.which = "sample", samplenr = 0)
  plotPhylo(x, plot.which = "sample", samplenr = 0)

}
