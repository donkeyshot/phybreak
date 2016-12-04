phybreak2phylo <- function(vars, samplenames = c(), simmap = FALSE) {
  ### extract variables
  nodetimes <- vars$nodetimes
  nodeparents <- vars$nodeparents
  nodehosts <- vars$nodehosts
  nodetypes <- vars$nodetypes
  Nhosts <- sum(nodetypes == "t")
  Nsamples <- sum(nodetypes == "s")
  if(is.null(samplenames)) samplenames <- 1:Nsamples
  if(length(samplenames) != Nsamples) {
    warning("length of samplenames does not match number of samples; samplenames not used")
    samplenames <- 1:Nsamples
  }
  
  ### topology for phylo
  indext <- (1:length(nodetimes))[nodeparents == 0]
  indexc <- (1:length(nodetimes))[nodeparents == indext]
  edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
  edgeend <- (1:length(nodetimes))[nodeparents != 0 & nodeparents != indext]
  edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]
  edgestarttypes <- nodetypes[edgestart]
  edgeendtypes <- nodetypes[edgeend]
  
  toremove <- edgeendtypes == "t"
  toreplace <- edgestarttypes == "t" 
  edgelengths[toreplace] <- edgelengths[toreplace] + edgelengths[toremove][match(edgestart[toreplace], edgeend[toremove])]
  edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace], edgeend[toremove])]
  edgestart <- edgestart[!toremove]
  edgeend <- edgeend[!toremove]
  edgelengths <- edgelengths[!toremove]
  
  edges <- matrix(c(edgestart, edgeend), ncol = 2)
  
  rootedge <- min(nodetimes[nodetypes == "c"]) - min(nodetimes[nodetypes == "t"])
  
  ### items for simmap (phytools)
  if(simmap) {
    whichcolor <- function(iors, path) {
      if (path[1] == 0) {
        return(length(path) %% 3)
      } else {
        return(whichcolor(iors, c(iors[path[1]], path)))
      }
    }
    
    hoststates <- 1 + sapply(1:Nhosts, whichcolor, iors = nodehosts[nodetypes == "t"])
    tipstates <- c("red", "blue", "black")[hoststates[nodehosts[nodetypes == "s"]]]
    nodestates <- matrix(tipstates[nodehosts[edges]], ncol = 2)
    edgemaps <- list()
    for (i in 1:nrow(edges)) {
      if (nodestates[i, 1] == nodestates[i, 2]) {
        edgemaps[[i]] <- edgelengths[i]
        names(edgemaps[[i]]) <- nodestates[i, 1]
      } else {
        nodes <- edges[i, ]
        ntimes <- nodetimes[c(nodes[1], nodeparents[nodes[2]], nodes[2])]
        edgemaps[[i]] <- ntimes[2:3] - ntimes[1:2]
        names(edgemaps[[i]]) <- nodestates[i, ]
      }
    }
    
    mappededge <- matrix(nrow = nrow(edges), ncol = 3)
    for (i in 1:nrow(edges)) {
      mappededge[i, ] <- edgemaps[[i]][c("blue", "black", "red")]
    }
    mappededge[is.na(mappededge)] <- 0
    colnames(mappededge) <- c("blue", "black", "red")
    rownames(mappededge) <- paste0(edges[, 1], ",", edges[, 2])
  }
  
  
  ### give first coalescent node number Nsamples+1
  if (indexc != Nsamples + 1) {
    edgestart[edgestart == indexc] <- 0
    edgeend[edgeend == indexc] <- 0
    edgestart[edgestart == Nsamples + 1] <- indexc
    edgeend[edgeend == Nsamples + 1] <- indexc
    edgestart[edgestart == 0] <- Nsamples + 1
    edgeend[edgeend == 0] <- Nsamples + 1
  }
  edges <- matrix(c(edgestart, edgeend), ncol = 2)
  
  ### make preliminary result for ladderizing
  respreladder <- list(edge = edges, edge.length = edgelengths, Nnode = Nsamples - 1, tip.label = samplenames)
  class(respreladder) <- c("phylo")
  respreladder <- ape::reorder.phylo(respreladder)
  
  ### ladderize to determine reordering of edges
  edgeladder <- ape::ladderize(respreladder)$edge[, 2]
  newedgeorder <- match(edgeladder, edgeend)
  newtiporder <- match(respreladder$tip.label, samplenames)
  
  if(simmap) {
    reorderedmaps <- list()
    for (i in 1:nrow(edges)) {
      reorderedmaps[[i]] <- edgemaps[[newedgeorder[i]]]
    }

    res <- list(edge = edges[newedgeorder, ], edge.length = edgelengths[newedgeorder], Nnode = Nsamples - 1, 
                tip.label = samplenames[newtiporder], root.edge = rootedge,
                node.state = nodestates[newedgeorder, ], states = tipstates[newtiporder], maps = reorderedmaps, mapped.edge = mappededge[newedgeorder, 
                                                                                                                                         ])
    class(res) <- c("simmap", "phylo")
  } else {
    res <- list(edge = edges[newedgeorder, ], edge.length = edgelengths[newedgeorder], Nnode = Nsamples - 1, 
                tip.label = samplenames[newtiporder], root.edge = rootedge)
    class(res) <- c("phylo")
  }
  
  
  attr(res, "order") <- "cladewise"
  
  return(res)
  
}

phybreak2trans <- function(vars, hostnames = c(), reference.date = 0) {
  ### extract variables
  nodetimes <- vars$nodetimes
  nodeparents <- vars$nodeparents
  nodehosts <- vars$nodehosts
  nodetypes <- vars$nodetypes
  Nhosts <- sum(nodetypes == "t")
  Nsamples <- sum(nodetypes == "s")
  if(is.null(hostnames)) hostnames <- paste0("host.",1:Nhosts)
  if(length(hostnames) != Nhosts) {
    warning("length of hostnames does not match number of hosts; hostnames not used")
    hostnames <- paste0("host.",1:Nhosts)
  }
  
  ### make new variables
  samtimes <- nodetimes[nodetypes == "s"] + reference.date
  names(samtimes) <- hostnames[nodehosts[nodetypes == "s"]]
  inftimes <- nodetimes[nodetypes == "t"] + reference.date
  names(inftimes) <- hostnames
  infectors <- c("index",hostnames)[1 + nodehosts[nodetypes == "t"]]
  names(infectors) <- hostnames
  
  ### return result
  return(list(
    sample.times = samtimes,
    sim.infection.times = inftimes,
    sim.infectors = infectors
  ))
}


### vars should contain $sample.times
### vars may contain $sim.infection.times, $sim.infectors, and $sim.tree
### vars may contain $sample.hosts, but this is not yet used
### if resample = TRUE, then resamplepars should contain
###     $mean.sample, $shape.sample, $mean.gen, $shape.gen, $wh.model, $wh.slope
transphylo2phybreak <- function(vars, resample = FALSE, resamplepars = NULL) {

  ### extract variables
  refdate <- min(vars$sample.times)
  samtimes <- as.numeric(vars$sample.times - refdate)
  Nsamples <- length(samtimes)
  hostnames <- names(vars$sample.times)
  ##### NB: adjustment needed for .rinftimes to deal with multiple samples #####
  if(is.null(vars$sim.infection.times) | is.null(vars$sim.infectors) | is.null(vars$sim.tree) | resample == TRUE) {
    resample <- TRUE
    inftimes <- .rinftimes(samtimes, resamplepars$mean.sample, resamplepars$shape.sample)
    infectors <- .rinfectors(inftimes, resamplepars$mean.gen, resamplepars$shape.gen)
  } else {
    inftimes <- as.numeric(vars$sim.infection.times - refdate)
    infectors <- match(vars$sim.infectors, hostnames)
    infectors[is.na(infectors)] <- 0
  }
  Nhosts <- length(inftimes)
  if(is.null(hostnames)) hostnames <- 1:Nhosts

  ##### NB: to be implemented later
  if(is.null(vars$sample.hosts)) {
    if(Nhosts != Nsamples) {
      stop("samples cannot be linked to hosts: add sample.hosts data")
    }
    samhosts <- 1:Nsamples
  } else {
    if(length(vars$sample.hosts) != Nsamples) {
      if(Nhosts != Nsamples) {
        stop("samples cannot be linked to hosts: adjust length of sample.hosts data")
      }
      warning("sample.hosts has incorrect length and is ignored")
      samhosts <- hostnames
    } else {
      samhosts <- match(vars$sample.hosts, hostnames)
    }
  }
  

  nodetypes <- c(rep("s", Nsamples), rep("c", Nsamples - 1),
                 rep("t", Nhosts))  #node type (sampling, coalescent, transmission)
  
  if(resample) {
    ##nodehosts, coalescent nodes in hosts with secondary cases or multiple samples
    coalnodehosts <- sort(c(infectors, samhosts))
    coalnodehosts <- coalnodehosts[duplicated(coalnodehosts)]
    nodehosts <- c(samhosts, coalnodehosts, infectors)
    
    ##initialize nodetimes with sampling and infection times
    nodetimes <- c(samtimes, rep(NA, Nsamples - 1), inftimes)
    ##sample coalescent times
    for(i in 1:Nhosts) {
      nodetimes[nodehosts == i & nodetypes == "c"] <-   #change the times of the coalescence nodes in host i...
        inftimes[i] +                      #...to the infection time +
        .samplecoaltimes(nodetimes[nodehosts == i & nodetypes != "c"] - inftimes[i],
                         resamplepars$wh.model, resamplepars$wh.slope)  #...sampled coalescent times
    }
    
    ##initialize nodeparents 
    nodeparents <- 0*nodetimes
    ##sample nodeparents
    for(i in 1:Nhosts) {
      nodeparents[nodehosts == i] <-     #change the parent nodes of all nodes in host i...
        .sampletopology(which(nodehosts == i), nodetimes[nodehosts == i], 
                        nodetypes[nodehosts == i], i + 2*Nsamples - 1, resamplepars$wh.model)
      #...to a correct topology, randomized where possible
    }
    
  } else {
    ##extract phylotree
    phytree <- vars$sim.tree
    phylobegin <- phytree$edge[, 1]
    phyloend <- phytree$edge[, 2]
    phylolengths <- phytree$edge.length

    ##determine hosts with coalescence nodes
    coalnodehosts <- sort(c(infectors, samhosts))
    coalnodehosts <- coalnodehosts[duplicated(coalnodehosts)]
    
    ##initialize nodetimes with infection times, nodeparents with phylotree
    ##initialize edgelengths to help build the tree
    nodetimes <- c(rep(NA, 2*Nsamples - 1), inftimes)
    nodeparents <- c(head(phylobegin[order(phyloend)], Nsamples),
                     NA,
                     tail(phylobegin[order(phyloend)], Nsamples - 2),
                     rep(NA, Nhosts))
    edgelengths <- c(head(phylolengths[order(phyloend)], Nsamples),
                     NA, tail(phylolengths[order(phyloend)], Nsamples - 2),
                     rep(NA, Nhosts))
    ##link first infection to root node of phylotree
    nodeparents[Nsamples + 1] <- which(infectors == 0) + 2 * Nsamples - 1
    edgelengths[Nsamples + 1] <- -inftimes[which(infectors == 0)] - ape::node.depth.edgelength(phytree)[1]
    
    ##complete nodetimes by adding branch lengths starting from root
    while(any(is.na(nodetimes))) {
      nodetimes[1:(2*Nsamples - 1)] <- nodetimes[nodeparents[1:(2*Nsamples - 1)]] + edgelengths[1:(2*Nsamples - 1)]
    }

    ##place transmission nodes between sampling and coalescent nodes for hosts without secondary cases and with one sample
    #a. place coalescent node before transmission node
    nodeparents[(2*Nsamples - 1) + setdiff(1:Nhosts, coalnodehosts)] <- nodeparents[samhosts[setdiff(1:Nhosts, coalnodehosts)]]
    #b. place transmission node before sampling node
    nodeparents[samhosts[setdiff(1:Nhosts, infectors)]] <- (2*Nsamples - 1) + setdiff(1:Nhosts, infectors)
    
    ##initialize nodehosts, first sampling and transmission nodes...
    nodehosts <- c(samhosts, which(infectors == 0), rep(NA, Nsamples - 2), infectors)   
    ##... then coalescent nodes that are parent of sampling nodes of infectors
    nodehosts[nodeparents[sort(unique(infectors))[-1]]] <- sort(unique(infectors))[-1]
    
    ##complete nodehosts by going backwards in phylotree
    while(any(is.na(nodehosts))) {
      #which are coalescent nodes with known parent & host?; order these by parent
      whichnodes <- which(!is.na(nodehosts) & !is.na(nodeparents) & nodetypes != "s")
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
    edgestosplit <- (nodehosts[nodeparents] != nodehosts & 
                       nodetypes == "c" & 
                       nodehosts[nodeparents] != 0)
    tnodestoadd <- nodehosts[edgestosplit] + 2*Nsamples - 1
    nodeparents[tnodestoadd] <- nodeparents[edgestosplit]
    nodeparents[edgestosplit] <- tnodestoadd
    nodeparents[nodehosts == 0] <- 0
  }
  
  return(list(
    d = list(hostnames = hostnames,
             samplenames = hostnames,
             referencedate = refdate),
    v = list(nodetimes = round(nodetimes, digits = 12),
             nodehosts = nodehosts,
             nodeparents = nodeparents,
             nodetypes = nodetypes)
  ))
}


obkData2phybreak <- function(data, resample = FALSE, resamplepars = NULL) {
  ### extract variables
  samtimes <- as.numeric(OutbreakTools::get.dates(data, "dna"))
  names(samtimes) <- OutbreakTools::get.individuals(data)
  if(!resample) {
    inftimes <- as.numeric(OutbreakTools::get.dates(data, "individuals"))
    infectors <- OutbreakTools::get.data(simulobk, "infector")
    tree <- OutbreakTools::get.trees(data)[[1]]  
  } else {
    inftimes <- NULL
    infectors <- NULL
    tree <- NULL
  }
  
  varslist <- list(
    sample.times = samtimes,
    sim.infection.times = inftimes,
    sim.infectors = infectors,
    sim.tree = tree
  )
  
  return(transphylo2phybreak(vars = varslist, resample = resample, resamplepars = resamplepars))
}
  
 