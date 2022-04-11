#' Plotting a phybreak object phylogenetic tree.
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
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' plot(MCMCstate, plot.which = "mpc")
#' @export
plotPhylo <- function(x, plot.which = c("sample", "mpc", "mtcc", "mcc"), samplenr = 0, ...) {
  if(!("phytools" %in% .packages(TRUE))) {
    stop("package 'phytools' should be installed for this function")
  }
  if(!inherits(x, c("phybreak", "phybreakdata"))) {
    stop("object x must be of class \"phybreak\" or \"phybreakdata\"")
  }
  
  if(inherits(x, "phybreakdata")) {
    if(exists("sim.infection.times", x) && exists("sim.infectors", x) && exists("sim.tree", x)) {
      samplenames <- names(x$sample.times)
      x <- transphylo2phybreak(x)
      x$d$names <- samplenames
      plot.which <- "sample"
      samplenr <- 0
    } else {
      stop("object x should contain sim.infection.times, sim.infectors, and sim.tree")
    }
    
  }
  
  plot.which <- plot.which[which(plot.which %in% c("sample", "mpc", "mtcc", "mcc"))[1]]
  if(is.na(plot.which)) stop("no valid 'plot.which'")
  if(plot.which == "sample" & samplenr > length(x$s$logLik)) {
    warning("requested 'samplenr' not available; current state used")
    samplenr <- 0
  }


    
  if (plot.which == "mpc") {
    simmapplot <- suppressWarnings(transtree(x, "mpc", phylo.class = TRUE))
  } else if (plot.which == "mtcc") {
    simmapplot <- suppressWarnings(transtree(x, "mtcc", phylo.class = TRUE))
  } else if (plot.which == "mcc") {
    simmapplot <- suppressWarnings(phylotree(x, phylo.class = TRUE))
  } else if (samplenr == 0) {
    simmapplot <- suppressWarnings(phybreak2phylo(x$v, x$d$names, simmap = TRUE))
  } else {
    x$v <- with(x, 
                 list(
                   inftimes = s$inftimes[, samplenr],
                   infectors = s$infectors[, samplenr],
                   nodetimes = c(v$nodetimes[v$nodetypes %in% c("s", "x")], s$nodetimes[, samplenr]),
                   nodeparents = s$nodeparents[, samplenr],
                   nodehosts = c(v$nodehosts[v$nodetypes %in% c("s", "x")], s$nodehosts[, samplenr]),
                   nodetypes = v$nodetypes
                 ))
    
    simmapplot <- suppressWarnings(phybreak2phylo(x$v, x$d$names, simmap = TRUE))
  }
  phytools::plotSimmap(simmapplot, mar = par("mar"), colors = setNames(nm = unique(as.character(simmapplot$node.state))), ...)
  rootnodetime <- x$d$reference.date + x$v$inftimes[1] + ifelse(sum(x$v$infectors == 0) > 1, min(x$v$nodetimes) - min(x$v$inftimes), simmapplot$root.edge)
  
  labs <- pretty(par("xaxp")[1:2] + rootnodetime)
  axis(1, at = labs - rootnodetime, labels = labs)
}
