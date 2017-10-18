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
#' # First make a phybreakdata-object by simulation.
#' simulation <- sim_phybreak(obsize = 5)
#' 
#' plot(simulation)
#' @export
plot.phybreakdata <- function(x, ...) {


  devAskNewPage(TRUE)
  on.exit(expr = devAskNewPage(FALSE))
  plotTrans(x, plot.which = "sample", samplenr = 0)
  plotPhylo(x, plot.which = "sample", samplenr = 0)

}
