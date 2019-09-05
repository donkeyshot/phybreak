#' Plotting a phybreak object transmission and phylogenetic tree.
#' 
#' Plots a \code{phybreak}-object as a transmission tree with coloured phylogenetic trees for each host. The default 
#'   is to plot the current state, but any posterior sample can be chosen, as well as various consensus trees.
#' 
#' @param x An object of class \code{phybreak}.
#' @param plot.which Either \code{"sample"} to plot the current state or a selected posterior sample, 
#'   \code{"mpc"} or \code{"mtcc"} to plot a consensus transmission tree (see \code{\link{transtree}}) or \code{"mcc"}
#'   to plot the maximum clade credibility tree (see \code{\link{phylotree}}).
#' @param samplenr If \code{plot.which = "sample"}, this indicates which posterior tree should be plotted: 
#'   \code{samplenr = 0} to plot the current state.
#' @param mar Plot margins.
#' @param tree.lwd Line width of phylogenetic trees.
#' @param tree.col Colours of the host phylogenetic trees. Defaults to evenly spaced hues in \code{hcl}, 
#'   with \code{c = 100} and \code{l = 65}.
#' @param hostlabel Whether to print the host names.
#' @param hostlabel.cex Size of host names. Defaults to a value between 0.5 and 1
#'   depending on outbreak size.
#' @param hostlabel.pos Position of host names, either just right from the host's phylogenetic tree,
#'   or aligned at the right hand side of the plot (\code{"right"}).
#' @param hostlabel.col Colours of the host names. Defaults to same as tree colours.
#' @param samplelabel Whether to print the sample names.
#' @param samplelabel.cex Size of sample names. Defaults to a value between 0.4 and 0.8
#'   depending on outbreak size.
#' @param samplelabel.pos Position of sample names, either just right from the host's phylogenetic tree,
#'   or aligned at the right hand side of the plot (\code{"right"}).
#' @param samplelabel.col Colours of the sample names. Defaults to same as tree colours, matching the host.
#' @param label.space Scales the space at the right-hand side to place the host and sample names.
#' @param host.col Colour of shading behind the host phylogenetic trees. Defaults to same as tree colours.
#' @param host.alpha Transparancy of shading behind the host phylogenetic trees.
#' @param cline.lty Line type of vertical lines showing the transmission contacts.
#' @param cline.lwd Line width of vertical lines showing the transmission contacts.
#' @param cline.col Colour of vertical lines showing the transmission contacts. If \code{NULL}, the tree colours are used.
#' @param cpoint.pch Character \code{par("pch")} used for points indicating the transmission contacts, 
#'   both at the tips and roots of the phylogenetic trees. 
#' @param cpoint.cex Size of transmission contact points.
#' @param cpoint.col Colour of points indicating the transmission contacts, associated with the secondary case,
#'   both at the tips and roots of the phylogenetic trees. Defaults to same as tree colours.
#' @param xlab X-axis title.
#' @param axis.cex Size of tick labels.
#' @param title.cex Size of X-axis title.
#' @param ... Further graphical parameters (see details).
#' @section Details: 
#'   Graphical parameters can be added by using names in the format \code{prefix.parameter} for the
#'   different parts of the plot. The \code{parameter} will then be passed on to the appropriate 
#'   graphics function, e.g. \code{tree.lty} to change the line type of the phylogenetic tree. The following 
#'   prefixes can be used: \code{tree} for the phylogenetic trees, \code{hostlabel} for the host names,
#'   \code{samplelabel} for the sample names, \code{host} for the host's shading areas, \code{cline} for 
#'   the vertical contact lines, \code{cpoint} for the transmission contact points, \code{axis} for 
#'   the X-axis, and \code{title} for the X-axis title.
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
plotPhyloTrans <- function(x, plot.which = c("sample", "mpc", "mtcc", "mcc"), samplenr = 0, 
                           mar = 0.1 + c(4, 0, 0, 0), 
                           tree.lwd = 1, tree.col = NULL, 
                           hostlabel = TRUE, hostlabel.cex = NULL, 
                           hostlabel.pos = "adjacent", hostlabel.col = tree.col,
                           samplelabel = FALSE, samplelabel.cex = NULL, 
                           samplelabel.pos = "right", samplelabel.col = tree.col,
                           label.space = 0.15, 
                           host.col = tree.col, host.alpha = 0.2,
                           cline.lty = 3, cline.lwd = 1, cline.col = "black", 
                           cpoint.pch = 20, cpoint.cex = 1, cpoint.col = tree.col, 
                           xlab = "Time", axis.cex = 1, title.cex = 1, ...) {
  ### tests ###
  if(!("phytools" %in% .packages(TRUE))) {
    stop("package 'phytools' should be installed for this function")
  }
  if(!inherits(x, c("phybreak", "phybreakdata"))) {
    stop("object x must be of class \"phybreak\" or \"phybreakdata\"")
  }
  
  ### change phybreakdata into phybreak object ###
  if(inherits(x, "phybreakdata")) {
    if(exists("sim.infection.times", x) && exists("sim.infectors", x) && exists("sim.tree", x)) {
      samplenames <- names(x$sample.times)
      sampletimes <- x$sample.times
      x <- transphylo2phybreak(x)
      x$d$names <- samplenames
      x$d$sample.times <- sampletimes
      plot.which <- "sample"
      samplenr <- 0
    } else {
      stop("object x should contain sim.infection.times, sim.infectors, and sim.tree")
    }
  }
  
  ### get inputdata for plotting the tree ###
  plot.which <- plot.which[which(plot.which %in% c("sample", "mpc", "mtcc", "mcc"))[1]]
  if(is.na(plot.which)) stop("no valid 'plot.which'")
  if(plot.which == "sample" & samplenr > length(x$s$logLik)) {
    warning("requested 'samplenr' not available; current state used")
    samplenr <- 0
  }
  if(plot.which %in% c("mpc", "mtcc", "mcc") && length(x$s$logLik) == 0) {
    warning("no samples available; current state used")
    plot.which <- "sample"
    samplenr <- 0
  }
  
  if (plot.which == "mpc") {
    samplenr <- .mpcinfector(x, length(x$s$logLik), TRUE, FALSE)
  } else if (plot.which == "mtcc") {
    samplenr <- tail(.mtcctree(x$s$infectors, 
                               x$s$inftimes, 
                               c(x$p$obs, length(x$s$logLik))), 1)
  } else if (plot.which == "mcc") {
    samplenr <- tail(.mcctree(x$s$nodeparents, 
                           x$s$nodetimes[1:(x$d$nsamples - 1), ], 
                           c(x$d$nsamples, length(x$s$logLik))), 1)
  } 
  if (samplenr > 0) {
    x$v <- list(
      inftimes = x$s$inftimes[, samplenr],
      infectors = x$s$infectors[, samplenr],
      nodetimes = c(x$v$nodetimes[x$v$nodetypes %in% c("s", "x")], x$s$nodetimes[, samplenr]),
      nodeparents = x$s$nodeparents[, samplenr],
      nodehosts = c(x$v$nodehosts[x$v$nodetypes %in% c("s", "x")], x$s$nodehosts[, samplenr]),
      nodetypes = x$v$nodetypes
    )
  }
  plotinput <- list(d = list(names = x$d$names,
                             hostnames = x$d$hostnames[1:length(x$v$inftimes)],
                             sample.times = x$d$sample.times,
                             reference.date = x$d$reference.date), 
                    v = phybreak2environment(x$v))
  plotinput$v$nodetimes <- plotinput$v$nodetimes
  
  makephylotransplot(plotinput, mar = mar, 
                     tree.lwd = tree.lwd, tree.col = tree.col, 
                     hostlabel = hostlabel, hostlabel.cex = hostlabel.cex, 
                     hostlabel.pos = hostlabel.pos, hostlabel.col = hostlabel.col,
                     samplelabel = samplelabel, samplelabel.cex = samplelabel.cex, 
                     samplelabel.pos = samplelabel.pos, samplelabel.col = samplelabel.col,
                     label.space = label.space, 
                     host.col = host.col, host.alpha = host.alpha,
                     cline.lty = cline.lty, cline.lwd = cline.lwd, cline.col = cline.col, 
                     cpoint.pch = cpoint.pch, cpoint.cex = cpoint.cex, cpoint.col = cpoint.col, 
                     xlab = xlab, axis.cex = axis.cex, title.cex = title.cex, ...)
}


makephylotransplot <- function(plotinput, mar = 0.1 + c(4, 0, 0, 0), 
                               tree.lwd = 1, tree.col = NULL, 
                               hostlabel = TRUE, hostlabel.cex = NULL, 
                               hostlabel.pos = "adjacent", hostlabel.col = tree.col,
                               samplelabel = FALSE, samplelabel.cex = NULL, 
                               samplelabel.pos = "right", samplelabel.col = tree.col,
                               label.space = 0.15, 
                               host.col = tree.col, host.alpha = 0.2,
                               cline.lty = 3, cline.lwd = 1, cline.col = "black", 
                               cpoint.pch = 20, cpoint.cex = 1, cpoint.col = tree.col, 
                               xlab = "Time", axis.cex = 1, title.cex = 1, ...) {
  oldmar <- par("mar")
  par(mar = mar)
  on.exit(par(mar = oldmar))
  obs <- length(plotinput$v$inftimes)
  
  ### determine order for plotting hosts ###
  inputhosts <- plotinput$d$hostnames
  inputinfectors <- c("index", plotinput$d$hostnames)[1 + plotinput$v$infectors]
  names(inputinfectors) <- inputhosts
  hostorder <- order(plotinput$v$inftimes)
  inputhosts <- inputhosts[hostorder]
  inputinfectors <- inputinfectors[hostorder]
  plotrank <- rankhostsforplot(inputhosts, inputinfectors)
  hostorder[plotrank] <- hostorder
  
  #########################################
  ### build the complete phylotranstree ###
  #########################################
  
  ### start with all tips as separate trees ###
  treelist <- list()
  # move along all nodes
  for(i in 1:length(plotinput$v$nodeparents)) {
    # only tips inside hosts
    if(plotinput$v$nodetypes[i] != "c" &&
       plotinput$v$nodehosts[i] > 0) {
      # parentnode of tip
      parentnode <- plotinput$v$nodeparents[i]
      # treat index case and pre-index as one host (split later)  
      if(plotinput$v$nodehosts[parentnode] == 0) {
        parentnode <- plotinput$v$nodeparents[parentnode]
      }
      treelist <- c(
        treelist,
        list(
          list(
            # line coordinates of tree
            x1vec = if(parentnode > 0) {
              plotinput$v$nodetimes[parentnode]
            } else min(plotinput$v$inftimes),
            x2vec = plotinput$v$nodetimes[i],
            y1vec = 0,
            y2vec = 0,
            # nodes of tree
            nodevec = i,
            rootnode = parentnode,
            xroot = if(parentnode > 0) {
              plotinput$v$nodetimes[parentnode]
            } else min(plotinput$v$inftimes),
            # give score for 'verticality' of tree, depending on tips
            # transmitting to lower or higher infectees
            yscore = if(plotinput$v$nodetypes[i] %in% c("t", "b")) {
              if(which(hostorder == plotinput$v$nodehosts[i]) >
                 which(hostorder == plotinput$v$nodehosts[
                   which(plotinput$v$nodeparents == i)
                   ])) {
                -1
              } else 1
            } else 0,
            host = plotinput$v$nodehosts[i],
            # is the within-host tree complete?
            completeYN = if(parentnode == 0 ||
                            plotinput$v$nodetypes[parentnode] %in% c("t", "b")) {
              TRUE
            } else FALSE
          )
        )
      )
    }
  }
  
  ### step by step coalesce all trees within hosts ###
  ### and collect them by host ###
  hosttrees <- list()
  while(length(treelist) > 0) {
    # next trees to coalesce, by going back in time
    xroots <- sapply(treelist, function(xx) xx$xroot)
    tojoin <- which(xroots == max(xroots))
    trees2join <- treelist[tojoin]
    treelist <- treelist[-tojoin]
    
    # order trees by yscore
    treeorder <- order(sapply(trees2join, function(xx) xx$yscore))
    trees2join <- trees2join[treeorder]
    
    # if trees are complete: collect them and add to hosttrees
    if(trees2join[[1]]$completeYN) {
      # finalize host
      hosttrees <- c(hosttrees, list(trees2join))
      # pass yscores on to transmission nodes in infector
      trnodes <- unlist(sapply(trees2join, 
                               function(xx) 
                                 which(xx$rootnode == 
                                         sapply(treelist, function(xxx) xxx$nodevec[1]))) )
      for(i in seq_along(trnodes)) {
        treelist[[trnodes[i]]]$yscore <- 
          treelist[[trnodes[i]]]$yscore + trees2join[[i]]$yscore/10
      }
    } else {
      # trees are not complete: coalesce within host
      
      # create room vertically by moving higher tree upwards
      yshift <- 1 + max(trees2join[[1]]$y1vec)
      trees2join[[2]]$y1vec <- trees2join[[2]]$y1vec + yshift
      trees2join[[2]]$y2vec <- trees2join[[2]]$y2vec + yshift
      
      # find new rootnode: treat index case and pre-index as one host (split later)  
      newrootnode <- plotinput$v$nodeparents[trees2join[[1]]$rootnode]
      if(newrootnode > 0 &&
         plotinput$v$nodehosts[newrootnode] == 0 &&
         plotinput$v$nodetypes[newrootnode] %in% c("t", "b")) {
        newrootnode <- plotinput$v$nodeparents[newrootnode]
      }
      
      # x-coordinate of rootnode: infection time of index case if
      # root of complete phylotree inside index case (root branch)
      newxroot <- if(newrootnode > 0) {
        plotinput$v$nodetimes[newrootnode]
      } else min(c(trees2join[[1]]$xroot, plotinput$v$inftimes))
      
      newyscore <- trees2join[[1]]$yscore + trees2join[[2]]$yscore
      newcompleteYN <- 
        if(newrootnode == 0 ||
           plotinput$v$nodetypes[newrootnode] %in% c("t", "b")) {
          TRUE
        } else FALSE
      
      newtree <- list(
        # join trees, add root and vertical line
        x1vec = c(newxroot, 
                  trees2join[[1]]$xroot, 
                  unlist(sapply(trees2join, function(xx) xx$x1vec))),
        x2vec = c(trees2join[[1]]$xroot, 
                  trees2join[[1]]$xroot, 
                  unlist(sapply(trees2join, function(xx) xx$x2vec))),
        y1vec = c(mean(c(trees2join[[1]]$y1vec[1], trees2join[[2]]$y1vec[1])),
                  trees2join[[1]]$y1vec[1],
                  unlist(sapply(trees2join, function(xx) xx$y1vec))),
        y2vec = c(mean(c(trees2join[[1]]$y1vec[1], trees2join[[2]]$y1vec[1])),
                  trees2join[[2]]$y1vec[1],
                  unlist(sapply(trees2join, function(xx) xx$y2vec))),
        # nodes of tree
        nodevec = c(trees2join[[1]]$rootnode, trees2join[[1]]$rootnode,
                    unlist(sapply(trees2join, function(xx) xx$nodevec))),
        rootnode = newrootnode,
        xroot = newxroot,
        # new score for tree 'verticality'
        yscore = newyscore,
        host = trees2join[[1]]$host,
        completeYN = newcompleteYN
      )
      treelist <- c(treelist, list(newtree))
    }
  }
  
  ### join the trees per host ###
  # final hosttree is index case and is already one tree
  for(i in (length(hosttrees)-1):1) {
    # order trees as transmission nodes in infector
    infector <- plotinput$v$infectors[hosttrees[[i]][[1]]$host]
    infectorhosttree <- 
      which(sapply(hosttrees, function(x) x[[1]]$host) == infector)
    reorderhosttree <-
      order(match(sapply(hosttrees[[i]], function(xx) xx$rootnode),
                  hosttrees[[infectorhosttree]][[1]]$nodevec))
    hosttrees[[i]] <- hosttrees[[i]][reorderhosttree]
    
    # create room vertically by moving higher trees upwards
    yshifts <- cumsum(c(0, 
                        sapply(hosttrees[[i]], function(xx) 1 + max(xx$y1vec))))
    for(j in 1:length(hosttrees[[i]])) {
      hosttrees[[i]][[j]]$y1vec <- hosttrees[[i]][[j]]$y1vec + yshifts[j]
      hosttrees[[i]][[j]]$y2vec <- hosttrees[[i]][[j]]$y2vec + yshifts[j]
    }
    
    # combine all trees into a single list
    hosttrees[[i]] <- list(list(
      x1vec = unlist(sapply(hosttrees[[i]], function(xx) xx$x1vec)),
      x2vec = unlist(sapply(hosttrees[[i]], function(xx) xx$x2vec)),
      y1vec = unlist(sapply(hosttrees[[i]], function(xx) xx$y1vec)),
      y2vec = unlist(sapply(hosttrees[[i]], function(xx) xx$y2vec)),
      nodevec = unlist(sapply(hosttrees[[i]], function(xx) xx$nodevec)),
      host = hosttrees[[i]][[1]]$host
    ))
  }
  
  ### finalize complete tree ###
  #order trees correctly
  hosttrees <- unlist(hosttrees, recursive = F)
  hosttrees <- hosttrees[match(hostorder, sapply(hosttrees, function(xx) xx$host))]
  
  # create room vertically by moving higher trees upwards
  yshifts <- cumsum(c(0, 
                      sapply(hosttrees, function(xx) 1 + max(xx$y1vec))))
  for(j in 1:length(hosttrees)) {
    hosttrees[[j]]$y1vec <- hosttrees[[j]]$y1vec + yshifts[j]
    hosttrees[[j]]$y2vec <- hosttrees[[j]]$y2vec + yshifts[j]
  }
  
  # combine all trees into a single list
  completetree <- data.frame(
    x1vec = unlist(sapply(hosttrees, function(xx) xx$x1vec)),
    x2vec = unlist(sapply(hosttrees, function(xx) xx$x2vec)),
    y1vec = unlist(sapply(hosttrees, function(xx) xx$y1vec)),
    y2vec = unlist(sapply(hosttrees, function(xx) xx$y2vec)),
    nodevec = unlist(sapply(hosttrees, function(xx) xx$nodevec)),
    hostvec = rep(hostorder, sapply(hosttrees, function(xx) length(xx$nodevec)))
  )
  
  # extract pre-index part from index case
  preindexrows <- completetree$x1vec <= min(plotinput$v$inftimes)
  introtree <- completetree[preindexrows, ]
  introtree$x2vec <- pmin(introtree$x2vec, min(plotinput$v$inftimes))
  introtree$hostvec <- 0
  nodevec2replace <- introtree$x2vec == min(plotinput$v$inftimes)
  introtree$nodevec[nodevec2replace] <-
    plotinput$v$nodeparents[introtree$nodevec[nodevec2replace]]
  # remove pre-index part from index case
  postindexrows <- completetree$x2vec > min(plotinput$v$inftimes)
  resttree <- completetree[postindexrows, ] 
  resttree$x1vec <- pmax(resttree$x1vec, min(plotinput$v$inftimes))
  # add pre-index part as separate subtree
  completetree <- rbind(introtree, resttree)
  completetree$hostvec <- factor(completetree$hostvec)
  
  ###########################################
  ### prepare vectors needed for plotting ###
  ###########################################
  
  ### input for phylotrees ###
  xphylo1 <- completetree$x1vec + plotinput$d$reference.date
  xphylo2 <- completetree$x2vec + plotinput$d$reference.date
  yphylo1 <- completetree$y1vec 
  yphylo2 <- completetree$y2vec 
  colphylo <- as.numeric(completetree$hostvec)
  
  ### input for transmission links ###
  xtrans <- plotinput$v$inftimes[plotinput$v$infectors > 0]
  ytrans <- sapply(xtrans,
                   function(xx) c(
                     min(completetree$y1vec[completetree$x1vec == xx]),
                     max(completetree$y1vec[completetree$x1vec == xx]),
                     min(completetree$y2vec[completetree$x2vec == xx]),
                     max(completetree$y2vec[completetree$x2vec == xx])))
  xtrans <- xtrans + plotinput$d$reference.date
  ytrans1 <- apply(ytrans, 2, min)
  ytrans2 <- apply(ytrans, 2, max)
  
  ### input for transmission nodes ###
  xnodes <- c(
    completetree$x2vec[completetree$nodevec %in%
                         which(plotinput$v$nodetypes %in% c("t", "b"))],
    unlist(sapply(hosttrees, 
                  function(xx) xx$x1vec[xx$x1vec == min(xx$x1vec) & xx$x1vec >= min(plotinput$v$inftimes)]))
  ) + plotinput$d$reference.date
  ynodes <- c(
    completetree$y1vec[completetree$nodevec %in%
                         which(plotinput$v$nodetypes %in% c("t", "b"))],
    unlist(sapply(hosttrees, 
                  function(xx) xx$y1vec[xx$x1vec == min(xx$x1vec) & xx$x1vec >= min(plotinput$v$inftimes)]))
  )
  nodeIDs <- completetree$nodevec[completetree$nodevec %in%
                                    which(plotinput$v$nodetypes %in% c("t", "b"))]
  infecteenodes <- match(nodeIDs, plotinput$v$nodeparents)
  infecteehosts <- c(
    plotinput$v$nodehosts[infecteenodes], 
    unlist(sapply(hosttrees, function(xx) rep(xx$host, sum(xx$x1vec == min(xx$x1vec) & xx$x1vec >= min(plotinput$v$inftimes))))))
  
  ### input for hosts ###
  xhost1 <- plotinput$v$inftimes + plotinput$d$reference.date
  xhost2 <- sapply(1:obs, function(xx) max(completetree$x2vec[xx == completetree$hostvec])) + plotinput$d$reference.date
  yhost1 <- sapply(1:obs, function(xx) min(completetree$y1vec[xx == completetree$hostvec])) - 0.25
  yhost2 <- sapply(1:obs, function(xx) max(completetree$y1vec[xx == completetree$hostvec])) + 0.25
  
  ### input for host labels ###
  if(hostlabel.pos != "right") {
    xhostname <- xhost2
  } else {
    xhostname <- max(xhost2)
  }
  yhostname <- (yhost1 + yhost2)/2
  labelhostname <- paste0(" ", plotinput$d$hostnames)

  ### input for sample labels ###
  samplenodes <- which(completetree$nodevec %in% which(plotinput$v$nodetypes %in% c("s", "x")))
  if(samplelabel.pos != "right") {
    xsamplename <- xhost2[colphylo[samplenodes] - 1]
  } else {
    xsamplename <- max(xhost2)
  }
  ysamplename <- completetree$y1vec[samplenodes]
  hostsamplename <- colphylo[samplenodes] - 1
  labelsamplename <- paste0(" ", plotinput$d$names[completetree$nodevec[samplenodes]])
    
  ### colours ###
  sqrtobs <- floor(sqrt(obs))
  if(is.null(tree.col)) {
    treecolours <- c("#888888", 
                     hcl(unlist(sapply(1:sqrtobs - 1, 
                                       function(xx) seq(xx, obs - 1, sqrtobs))) * 360/obs, 
                         c = 100, l = 65))
  } else {
    treecolours <- c("#888888",
                     rep_len(tree.col, length.out = obs))
  }
  if(is.null(hostlabel.col)) {
    hostlabelcolours <- tail(treecolours, -1)
  } else {
    hostlabelcolours <- rep_len(hostlabel.col, length.out = obs)
  }
  if(is.null(samplelabel.col)) {
    samplelabelcolours <- tail(treecolours, -1)
  } else {
    samplelabelcolours <- rep_len(samplelabel.col, length.out = obs)
  }
  if(is.null(host.col)) {
    hostcolours <- tail(treecolours, -1)
  } else {
    hostcolours <- rep_len(host.col, length.out = obs)
  }
  if(is.null(cline.col)) {
    clinecolours <- tail(treecolours, -1)[plotinput$v$infectors > 0]
  } else {
    clinecolours <- rep_len(cline.col, length.out = obs)[plotinput$v$infectors > 0]
  }
  if(is.null(cpoint.col)) {
    cpointcolours <- tail(treecolours, -1)
  } else {
    cpointcolours <- rep_len(cpoint.col, length.out = obs)
  }
  
  ### some smart graphical parameters
  if(is.null(hostlabel.cex)) hostlabel.cex <- max(0.5, min(1, 30/obs))
  if(is.null(samplelabel.cex)) samplelabel.cex <- 0.8 * hostlabel.cex
  tmin <- min(xphylo1)
  tmax <- max(xphylo2)
  
  
  
  ### initialize plot
  plot.new()
  par(cex = 1)
  if(hostlabel || samplelabel) {
    plot.window(xlim = c(tmin, tmax + label.space * hostlabel.cex * (tmax - tmin)), 
                ylim = c(0, max(yphylo2)))
  } else {
    plot.window(xlim = c(tmin, tmax), 
                ylim = c(0, max(yphylo2)))
  }
  
  ### X-axis
  do.call(Axis,
          c(list(x = c(round(tmin), round(tmax)),
                 side = 1,
                 cex.axis = axis.cex),
            graphicalparameters("axis", 1, ...)))

  ### Axis title
  do.call(title,
          c(list(xlab = xlab,
                 cex.lab = title.cex,
                 line = par("mar")[1]*5/8),
            graphicalparameters("title", 1, ...)
          )
  )
  
  ### hosts
  do.call(rect,
          c(list(xleft = xhost1, 
                 xright = xhost2,
                 ybottom = yhost1, 
                 ytop = yhost2,
                 col = adjustcolor(hostcolours, alpha = host.alpha), 
                 border = NA),
            graphicalparameters("host", 1:obs, ...)))

  ### tree
  do.call(segments,
          c(list(x0 = xphylo1, 
                 x1 = xphylo2,
                 y0 = yphylo1, 
                 y1 = yphylo2,
                 col = treecolours[colphylo], 
                 lwd = tree.lwd),
            graphicalparameters("tree", 1 + (colphylo + obs - 1) %% (1 + obs), ...)))
  
  ### contact lines
  do.call(segments,
          c(list(x0 = xtrans, 
                 y0 = ytrans1, 
                 y1 = ytrans2,
                 lty = cline.lty,
                 lwd = cline.lwd,
                 col = clinecolours),
            graphicalparameters("cline", which(plotinput$v$infectors > 0, ...))))
  
  ### contact points
  do.call(points,
          c(list(x = xnodes, 
                 y = ynodes, 
                 pch = cpoint.pch,
                 cex = cpoint.cex,
                 col = cpointcolours[infecteehosts]),
            graphicalparameters("cpoint", infecteehosts, ...)))
  
  ### host labels
  if(hostlabel) {
    do.call(text,
            c(list(x = xhostname, 
                   y = yhostname, 
                   label = labelhostname,
                   col = hostlabelcolours, 
                   cex = hostlabel.cex,
                   adj = 0),
              graphicalparameters("hostlabel", 1:obs, ...)))
  }
  
  ### sample labels
  if(samplelabel) {
    do.call(text,
            c(list(x = xsamplename, 
                   y = ysamplename, 
                   label = labelsamplename,
                   col = samplelabelcolours[hostsamplename],
                   cex = samplelabel.cex,
                   adj = 0),
              graphicalparameters("samplelabel", 1:obs, ...)))
  }
  
}
