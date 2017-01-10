phybreak2phylo <- function(vars, samplenames = c(), simmap = FALSE) {
  ### extract variables
  nodetimes <- vars$nodetimes
  nodeparents <- as.integer(vars$nodeparents)
  nodehosts <- as.integer(vars$nodehosts)
  nodetypes <- vars$nodetypes
  Nhosts <- sum(nodetypes == "t")
  Nsamples <- sum(nodetypes == "s")
  if(is.null(samplenames)) samplenames <- 1:Nsamples
  if(length(samplenames) != Nsamples) {
    warning("length of samplenames does not match number of samples; samplenames not used")
    samplenames <- 1:Nsamples
  }
  
  ### give first coalescent node number Nsamples + 1 for compatibility with phylo
  currentrootnode <- which(nodeparents == which(nodeparents == 0))
  if(currentrootnode != Nhosts + 1) {
    #swap hosts and times
    nodehosts[c(Nhosts + 1, currentrootnode)] <- nodehosts[c(currentrootnode, Nhosts + 1)]
    nodetimes[c(Nhosts + 1, currentrootnode)] <- nodetimes[c(currentrootnode, Nhosts + 1)]
    
    #swap parents
    nodeparents[c(Nhosts + 1, currentrootnode)] <- nodeparents[c(currentrootnode, Nhosts + 1)]

    #change parents of all nodes
    tooldroot <- nodeparents == currentrootnode
    tonewroot <- nodeparents == Nhosts + 1
    nodeparents[tooldroot] <- Nhosts + 1
    nodeparents[tonewroot] <- currentrootnode
  }
  
  ### topology for phylo
  edges <- removetransnodes(nodeparents, nodetypes)
 
  edgelengths <- apply(edges, 1, function(x) nodetimes[x[2]] - nodetimes[x[1]])
  
  rootedge <- min(nodetimes[nodetypes == "c"]) - min(nodetimes[nodetypes == "t"])
  
  ### items for simmap (phytools)
  if(simmap) {
    
    hostcolors <- c("red", "blue", "black")[1 + sapply(1:Nhosts, whichgeneration, infectors = nodehosts[nodetypes == "t"]) %% 3]
    tipstates <- hostcolors[nodehosts[nodetypes == "s"]]
    nodestates <- apply(edges, 1:2, function(x) tipstates[nodehosts[x]])
    
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
  
  

  ### make preliminary result for ladderizing
  res_preladder <- list(edge = edges, 
                       edge.length = edgelengths, 
                       Nnode = Nsamples - 1, 
                       tip.label = samplenames)
  class(res_preladder) <- c("phylo")
  res_preladder <- ape::reorder.phylo(res_preladder)
  
  ### ladderize to determine reordering of edges
  edgeladder <- ape::ladderize(res_preladder)$edge[, 2]
  newedgeorder <- match(edgeladder, edges[, 2])
  newtiporder <- match(res_preladder$tip.label, samplenames)
  
  if(simmap) {
    
    reorderedmaps <- list()
    for (i in 1:nrow(edges)) {
      reorderedmaps[[i]] <- edgemaps[[newedgeorder[i]]]
    }

    res <- list(edge = edges[newedgeorder, ], 
                edge.length = edgelengths[newedgeorder], 
                Nnode = Nsamples - 1L, 
                tip.label = samplenames[newtiporder], 
                root.edge = rootedge,
                node.state = nodestates[newedgeorder, ], 
                states = tipstates[newtiporder], 
                maps = reorderedmaps, 
                mapped.edge = mappededge[newedgeorder, ])
    class(res) <- c("simmap", "phylo")
    
  } else {
    
    res <- list(edge = edges[newedgeorder, ], 
                edge.length = edgelengths[newedgeorder], 
                Nnode = Nsamples - 1L, 
                tip.label = samplenames[newtiporder], 
                root.edge = rootedge)
    class(res) <- c("phylo")
    
  }
  
  attr(res, "order") <- "cladewise"
  
  return(res)
  
}

phybreak2trans <- function(vars, hostnames = c(), reference.date = 0) {
  ### extract variables
  nodetimes <- vars$nodetimes
  nodehosts <- vars$nodehosts
  nodetypes <- vars$nodetypes
  Nhosts <- sum(nodetypes == "t")
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
  if(is.null(vars$sim.infection.times) | is.null(vars$sim.infectors) | resample) {
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
  
  if(is.null(vars$sim.tree) | resample) {
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
    phytree <- vars$sim.tree
    nodetimes <- c(ape::node.depth.edgelength(phytree) - min(ape::node.depth.edgelength(phytree)[1:Nsamples]),
                   inftimes)
    nodehosts_parents <- make_nodehostsparents(phytree$edge, samhosts, infectors)
    nodehosts <- nodehosts_parents[, 1]
    nodeparents <- nodehosts_parents[, 2]
    }
  
  return(list(
    d = list(hostnames = hostnames,
             samplenames = hostnames,
             reference.date = refdate),
    v = list(nodetimes = round(nodetimes, digits = 12),
             nodehosts = nodehosts,
             nodeparents = nodeparents,
             nodetypes = nodetypes)
  ))
}


obkData2phybreak <- function(data, resample = FALSE, resamplepars = NULL) {
  ### extract variables
  samtimes <- OutbreakTools::get.dates(data, "dna")
  names(samtimes) <- OutbreakTools::get.individuals(data)
  if(!resample) {
    inftimes <- OutbreakTools::get.dates(data, "individuals")
    infectors <- OutbreakTools::get.data(data, "infector")
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
  
 

### remove transmission nodes from nodeparents object, and return 'edge' item in phylo-object
removetransnodes <- function(nodeparents, nodetypes) {
  # index node
  indext <- which(nodeparents == 0)

  # all edges apart from index edge
  edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
  edgeend <- (1:length(nodetypes))[nodeparents != 0 & nodeparents != indext]
  
  # remove edges ending in transmission node
  # replace parent transmission nodes by their parent coalescent node
  toremove <- nodetypes[edgeend] == "t"
  toreplace <- nodetypes[edgestart] == "t" 
  edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace], edgeend[toremove])]
  edgestart <- edgestart[!toremove]
  edgeend <- edgeend[!toremove]
  
  return(matrix(c(edgestart, edgeend), ncol = 2))
}

### to which generation does hostID belong in the transmission tree?
whichgeneration <- function(infectors, hostID) {
  if (hostID[1] == 0) {
    return(length(hostID) - 2L)
  } else {
    return(whichgeneration(infectors, c(infectors[hostID[1]], hostID)))
  }
}


### random infection times given sampling times and sampling interval distribution
.rinftimes <- function(st, meanS, shapeS) {
  ### tests
  if(class(st) != "numeric" && class(st) != "integer") {
    stop(".rinftimes called with non-numeric sampling times")
  }
  if(meanS <= 0) stop(".rinftimes called with non-positive mean sampling interval")
  if(shapeS <= 0) stop(".rinftimes called with non-positive shape parameter")
  
  ### function body
  st - rgamma(length(st), shape = shapeS, scale = meanS/shapeS)
}

### random infectors times given infection times and generation interval distribution
.rinfectors <- function(it, meanG, shapeG) {
  ### tests
  if(class(it) != "numeric" && class(it) != "integer") {
    stop(".rinfectors called with non-numeric infection times")
  }
  if(sum(it == min(it)) > 1) stop("rinfectors with >1 index case")
  if(meanG <= 0) stop(".rinfectors called with non-positive mean generation interval")
  if(shapeG <= 0) stop(".rinfectors called with non-positive shape parameter")
  
  ### function body
  res <- rep(0,length(it))
  for(i in 1:length(it)) {
    if(it[i] > min(it)) {
      dist <- dgamma(it[i] - it, shape = shapeG, scale = meanG/shapeG)
      dist[i] <- 0
      res[i] <- sample(length(it), 1, prob = dist)
    }
  }
  return(res)
}

###turns phylogenetic tree & infectors into nodehosts & nodeparents (wrapper function)
make_nodehostsparents <- function(edges, samplehosts, infectors) {
  rootnode <- which(!(1:(1 + nrow(edges))) %in% edges[, 2])
  indextnode <- which(infectors == 0) + 2 * length(samplehosts) - 1
  nodeparents <- c(c(indextnode, edges[, 1])[order(c(rootnode, edges[, 2]))], rep(NA, length(infectors)))
  nodehosts <- c(samplehosts, rep(NA, length(samplehosts) - 1), infectors)
  
  nodehp <- cbind(nodehosts, nodeparents)
  res <- makenodehp(nodehp, rootnode, infectors, length(samplehosts))
  res[indextnode, 2] <- 0
  return(res)
}

###turns phylogenetic tree & infectors into nodehosts & nodeparents (worker function)
makenodehp <- function(nodehp, nodeID, infectors, Nsamples) {
  # identify children in phylotree
  childnodes <- which(nodehp[, 2] == nodeID)
  
  # For each child's host that is unknown, call this function with child as nodeID
  if(is.na(nodehp[childnodes[1], 1])) {
    nodehp <- makenodehp(nodehp, childnodes[1], infectors, Nsamples)
  }
  if(is.na(nodehp[childnodes[2], 1])) {
    nodehp <- makenodehp(nodehp, childnodes[2], infectors, Nsamples)
  }

  # Determine node's host and place transmission nodes in downstream edges
  if(nodehp[childnodes[1], 1] == nodehp[childnodes[2], 1]) {
    nodehp[nodeID, 1] <- nodehp[childnodes[1], 1]
  } else {
    if(nodehp[childnodes[1], 1] != infectors[nodehp[childnodes[2], 1]]) {
      nodehp[nodeID, 1] <- nodehp[infectors[nodehp[childnodes[1], 1]], 1]
      transmissionnode <- nodehp[childnodes[1], 1] + 2 * Nsamples - 1
      nodehp[transmissionnode, 2] <- nodeID
      nodehp[childnodes[1], 2] <- transmissionnode
    }
    if(nodehp[childnodes[2], 1] != infectors[nodehp[childnodes[1], 1]]) {
      nodehp[nodeID, 1] <- nodehp[infectors[nodehp[childnodes[2], 1]], 1]
      transmissionnode <- nodehp[childnodes[2], 1] + 2 * Nsamples - 1
      nodehp[transmissionnode, 2] <- nodeID
      nodehp[childnodes[2], 2] <- transmissionnode
    }
  }

  # return tree with resolved node + downstream nodes
  return(nodehp)
}

