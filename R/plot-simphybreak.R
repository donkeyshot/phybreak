### plot functions for simphybreak-object ###

### makes use of plot function in 'simmap', colouring the branches by host in which they reside, but no unique colour per host.
### calls: transtree get.phylo


#' Plotting a phybreak object.
#' 
#' Plots a \code{simphybreak}-object as phylogenetic tree with coloured branches indicating hosts. The default 
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
plot.simphybreak <- function(x, ...) {
  
}




.makephylosimmap.simphybreak <- function(simphybreak.object) {
  ##extract phylotree from obkData
  samtimes <- simphybreak.object$sample.times
  inftimes <- simphybreak.object$infection.times
  infectors <- simphybreak.object$infectors
  infectors <- match(infectors, names(infectors))
  infectors[is.na(infectors)] <- 0
  obs <- length(samtimes)
  phytree <- simphybreak.object$trees[[1]]
  phylobegin <- phytree$edge[,1]
  phyloend <- phytree$edge[,2]
  phylolengths <- phytree$edge.length

  ##initialize nodetimes with infection times, nodeparents with phylotree
  ##initialize edgelengths to help build the tree
  nodetimes <- c(rep(NA,2*obs-1), inftimes)
  nodeparents <- c(head(phylobegin[order(phyloend)], obs),
                   NA,
                   tail(phylobegin[order(phyloend)], obs - 2),
                   rep(NA, obs))
  edgelengths <- c(head(phylolengths[order(phyloend)], obs),
                   NA, tail(phylolengths[order(phyloend)], obs - 2),
                   rep(NA, obs))
  ##link first infection to root node of phylotree
  nodeparents[obs + 1] <- which(infectors == 0) + 2 * obs - 1
  edgelengths[obs + 1] <- samtimes[which(infectors == 0)] - ape::node.depth.edgelength(phytree)[1]
  ##place transmission nodes between sampling and coalescent nodes for hosts without secondary cases
  #a. place coalescent node before transmission node
  nodeparents[(2*obs - 1) + setdiff(1:obs, infectors)] <- nodeparents[setdiff(1:obs, infectors)]
  #b. place transmission node before sampling node
  nodeparents[setdiff(1:obs, infectors)] <- (2*obs - 1) + setdiff(1:obs, infectors)
  
  ##complete nodetimes by adding branch lengths starting from root
  while(any(is.na(nodetimes))) {
    nodetimes[1:(2*obs - 1)] <- nodetimes[nodeparents[1:(2*obs - 1)]] + edgelengths[1:(2*obs - 1)]
  }
  nodetimes <- round(nodetimes, digits = 12) #prevent numerical problems
  
  ##initialize nodehosts, first sampling and transmission nodes...
  nodehosts <- c(1:obs, which(infectors == 0), rep(NA, obs - 2), infectors)   
  ##... then coalescent nodes that are parent of sampling nodes of infectors
  nodehosts[nodeparents[sort(unique(infectors))[-1]]] <- sort(unique(infectors))[-1]
  
  ##complete nodehosts by going backwards in phylotree
  while(any(is.na(nodehosts))) {
    #which are coalescent nodes with known parent & host?; order these by parent
    whichnodes <- tail(which(!is.na(nodehosts) & !is.na(nodeparents)),-obs)
    whichnodes <- whichnodes[order(nodeparents[whichnodes])]
    
    #which of these nodes have the same parent node?
    sameparent.wn <- which(duplicated(nodeparents[whichnodes]))
    #determine the host of these parent nodes
    for(i in sameparent.wn) {
      if(nodehosts[whichnodes[i]] == nodehosts[whichnodes[i-1]]) {
        #...if host of children nodes are the same, this is also host of parent
        nodehosts[nodeparents[whichnodes[i]]] <- nodehosts[whichnodes[i]]
      } else {
        #...if host of children nodes are different, host of parent node is the
        #...infector if one infects the other, or the common infector of both
        posshosts <- nodehosts[whichnodes[i - 1:0]]
        posshosts <- c(posshosts, infectors[posshosts])
        nodehosts[nodeparents[whichnodes[i]]] <- posshosts[duplicated(posshosts)]
      }
    }
  }
  
  ##complete nodeparents: place transmission nodes for hosts with secondary infections
  #a. place coalescent node before transmission node
  nodeparents[nodehosts[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != 
                                nodehosts[(obs + 1):(2 * obs - 1)]) + obs] + 2 * obs - 1] <- 
    nodeparents[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != nodehosts[(obs + 1):(2 * obs - 1)]) + obs]
  #b. place transmission node before sampling node
  nodeparents[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != nodehosts[(obs + 1):(2 * obs - 1)]) + obs] <- 
    nodehosts[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != nodehosts[(obs + 1):(2 * obs - 1)]) + obs] + 2 * obs - 1
  nodeparents[nodehosts == 0] <- 0
  print(nodetimes)
  print(nodeparents)
  print(nodehosts)
  print(names(infectors))
  
  .makephylosimmap.phybreak(nodetimes, nodeparents, nodehosts, names(infectors))
}
