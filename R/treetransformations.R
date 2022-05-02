phybreak2phylo <- function(vars, samplenames = c(), simmap = FALSE) {
  ### extract variables
  nodetimes <- vars$nodetimes
  nodeparents <- as.integer(vars$nodeparents)
  nodehosts <- as.integer(vars$nodehosts) + 1L
  nodetypes <- vars$nodetypes
  inftimes <- c(-Inf, vars$inftimes)
  infectors <- c(0, vars$infectors + 1L)
  nhosts <- length(inftimes)
  nsamples <- sum(nodetypes %in% c("s", "x"))
  if(is.null(samplenames)) samplenames <- 1:nsamples
  if(length(samplenames) != nsamples) {
    warning("length of samplenames does not match number of samples; samplenames not used")
    samplenames <- 1:nsamples
  }
  
  ### give first coalescent node number nsamples + 1 for compatibility with phylo
  currentrootnode <- which(nodeparents == 0)
  if(currentrootnode != nsamples + 1) {
    #swap hosts and times
    nodehosts[c(nsamples + 1, currentrootnode)] <- nodehosts[c(currentrootnode, nsamples + 1)]
    nodetimes[c(nsamples + 1, currentrootnode)] <- nodetimes[c(currentrootnode, nsamples + 1)]
    
    #swap parents
    nodeparents[c(nsamples + 1, currentrootnode)] <- nodeparents[c(currentrootnode, nsamples + 1)]
    
    #change parents of all nodes
    tooldroot <- nodeparents == currentrootnode
    tonewroot <- nodeparents == nsamples + 1
    nodeparents[tooldroot] <- nsamples + 1
    nodeparents[tonewroot] <- currentrootnode
  }
  
  ### topology for phylo
  
  edges <- cbind(nodeparents, 1:(2*nsamples - 1))[-(nsamples + 1), ]
  
  edgelengths <- apply(edges, 1, function(x) nodetimes[x[2]] - nodetimes[x[1]])
  
  rootedge <- max(0, min(nodetimes) - min(inftimes[-1]))
  
  ### items for simmap (phytools)
  if(simmap) {
    
    hostcolors <- c("red", "blue", "black")[1 + sapply(1:nhosts, whichgeneration, infectors = infectors) %% 3]
    tipstates <- hostcolors[nodehosts[nodetypes %in% c("s", "x")]]
    nodestates <- apply(edges, 1:2, function(x) hostcolors[nodehosts[x]])
    
    edgemaps <- list()
    for(i in 1:nrow(edges)) {
      starthost <- nodehosts[edges[i, 1]]
      endhost <- nodehosts[edges[i, 2]]
      edgehosts <- rev(c(setdiff(.ptr(infectors, endhost), .ptr(infectors, starthost)), starthost))
      ntimesstart <- c(nodetimes[edges[i, 1]], inftimes[edgehosts[-1]])
      ntimesend <- c(inftimes[edgehosts[-1]], nodetimes[edges[i, 2]])
      edgemaps[[i]] <- ntimesend - ntimesstart
      names(edgemaps[[i]]) <- hostcolors[edgehosts]
    }
    
    
    mappededge <- matrix(nrow = nrow(edges), ncol = 3)
    for (i in 1:nrow(edges)) {
      mappededge[i, ] <- sapply(c("blue", "black", "red"), function(x) sum(edgemaps[[i]][names(edgemaps[[i]]) == x]))
    }
    mappededge[is.na(mappededge)] <- 0
    colnames(mappededge) <- c("blue", "black", "red")
    rownames(mappededge) <- paste0(edges[, 1], ",", edges[, 2])
  }
  
  
  
  ### make preliminary result for ladderizing
  res_preladder <- list(edge = edges, 
                        edge.length = edgelengths, 
                        Nnode = nsamples - 1, 
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
                Nnode = nsamples - 1L, 
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
                Nnode = nsamples - 1L, 
                tip.label = samplenames[newtiporder], 
                root.edge = rootedge)
    class(res) <- c("phylo")
    
  }
  
  attr(res, "order") <- "cladewise"
  
  return(res)
  
}

phybreak2trans <- function(vars, hostnames = c(), reference.date = 0,
                           culling.times = NULL) {
  ### extract variables
  nodetimes <- vars$nodetimes
  nodehosts <- vars$nodehosts
  nodetypes <- vars$nodetypes
  nhosts <- length(vars$infectors)
  if(is.null(hostnames)) hostnames <- paste0("host.",1:nhosts)
  hostnames <- unique(hostnames)
  if(length(hostnames) < nhosts) {
    pre_hostnames <- paste0("host.",1:nhosts)
    pre_hostnames[1:length(hostnames)] <- hostnames
    hostnames <- pre_hostnames
  } else if(length(hostnames) > nhosts) {
    hostnames <- hostnames[1:nhosts]
  }
  
  ### make new variables
  samtimes <- nodetimes[nodetypes %in% c("s", "x")] + reference.date
  names(samtimes) <- hostnames[nodehosts[nodetypes %in% c("s", "x")]]
  inftimes <- vars$inftimes + reference.date
  names(inftimes) <- hostnames
  infectors <- c("index", hostnames)[1 + vars$infectors]
  names(infectors) <- hostnames
  
  ### return result
  return(list(
    sample.times = samtimes,
    culling.times = culling.times,
    sim.infection.times = inftimes,
    sim.infectors = infectors
  ))
}


transphylo2phybreak <- function(vars, resample = FALSE, resamplepars = NULL, 
                                introductions = 1) {
  
  ### extract and order samples
  refdate <- min(vars$sample.times)
  samtimes <- vars$sample.times - refdate
  if(exists("culling.times", vars))
    cultimes <- vars$culling.times - refdate
  else 
    cultimes <- NULL
  nsamples <- length(samtimes)
  if(exists("sample.hosts", vars)) {
    samhosts <- vars$sample.hosts
    
    if(length(samtimes) != length(samhosts)) {
      stop("lengths of sample.times and sample.hosts not equal")
    }
    
    if(!is.null(names(samtimes))) {
      if(!is.null(names(samhosts))) {
        if(!(all(names(samtimes) %in% names(samhosts)))) {
          stop("sample.times and sample.hosts don't have same names")
        }
        samhosts <- samhosts[names(samtimes)]
      } else {
        names(samhosts) <- names(samtimes)
      }
    }
    
    # samhosts <- samhosts[order(samtimes)]
    if(is.null(names(samhosts))) {
      names(samhosts) <- 1:nsamples
    }
  } else {
    if(!is.null(names(samtimes))) {
      samhosts <- names(samtimes)
      names(samhosts) <- names(samtimes)
      # samhosts <- samhosts[order(samtimes)]
    } else {
      samhosts <- 1:nsamples
      names(samhosts) <- 1:nsamples
    }
  }
  samtimes <- as.numeric(samtimes)
  cultimes <- as.numeric(cultimes)
  samplenames <- names(samhosts)
  hostnames <- unique(samhosts)
  samhosts <- match(samhosts, hostnames)
  nhosts <- length(hostnames)

  ### infection times and infectors
  if(is.null(vars$sim.infection.times) | is.null(vars$sim.infectors) | resample) {
    resample <- TRUE
    inftimes <- .rinftimes(samtimes[1:nhosts], resamplepars$sample.mean, resamplepars$sample.shape)
    infectors <- .rinfectors(inftimes, introductions, p = resamplepars, 
                             v = list(cultimes=cultimes, nodetimes = samtimes))
  } else {
    inftimes <- as.numeric(vars$sim.infection.times - refdate)
    infectors <- match(vars$sim.infectors, hostnames)
    infectors[is.na(infectors)] <- 0
    if(length(hostnames) != length(infectors) | length(hostnames) != length(inftimes)) {
      stop("numbers of infection times and infectors should match the number of sampled hosts")
    }
  }
  
  if(is.null(vars$sim.tree) | resample) {
    res <- list(
      inftimes = inftimes,
      infectors = infectors,
      cultimes = cultimes,
      nodetypes = c(rep("s", nhosts), rep("x", nsamples - nhosts), rep("c", nsamples - 1), rep("t", nhosts)),  
      #node type (primary sampling, extra sampling, coalescent)
      nodetimes = c(samtimes, samtimes[-1], inftimes),
      nodehosts = c(samhosts, rep(-1, length(samhosts) - 1), infectors),
      nodeparents = rep(-1, 2 * nsamples + nhosts - 1)
    )
    list2env(list(v = res, p = resamplepars, d = list(nsamples = nsamples)), pbe1)

    if(resamplepars$wh.bottleneck == "wide") {
      invisible(sapply(1:nhosts, rewire_pullnodes_wide))
    } else {
      invisible(sapply(0:nhosts, rewire_pullnodes_complete))
    }
    res <- environment2phybreak(pbe1$v)
  } else {
    res <- list(
      inftimes = inftimes,
      infectors = infectors,
      nodetypes = c(rep("s", nhosts), rep("x", nsamples - nhosts), rep("c", nsamples - 1)),  
      #node type (primary sampling, extra sampling, coalescent)
      nodetimes = c(samtimes, samtimes[-1]),
      nodehosts = c(samhosts, rep(-1, length(samhosts) - 1)),
      nodeparents = rep(-1, 2 * nsamples - 1)
    )
    phytree <- vars$sim.tree
    res$nodetimes <- c(ape::node.depth.edgelength(phytree) - min(ape::node.depth.edgelength(phytree)[1:nsamples]))
    res$nodeparents <- c(0, phytree$edge[,1])[order(c(nsamples + 1, phytree$edge[, 2]))]
    list2env(list(v = res, p = list(wh.model = "single")), pbe1)
    make_nodehosts(pbe1)
    res <- pbe1$v
  }
  
  return(list(
    d = list(names = samplenames,
             hostnames = hostnames,
             reference.date = refdate),
    v = res
  ))
}


phybreak2environment <- function(vars) {
  obs <- sum(vars$nodetypes == "s")
  nsamples <- sum(vars$nodetypes %in% c("s", "x"))
  
  if (obs < length(vars$infectors)){
    newnode <- 2*nsamples
    vars$nodetypes[newnode] <- "t"
    vars$nodetimes[newnode] <- vars$inftimes[1]
    vars$nodehosts[newnode] <- 0
    vars$nodeparents[which(vars$nodeparents == 0)] <- newnode
    vars$nodeparents[newnode] <- 0
    
  }
  
  for(edge in 1:(obs)) {
    curedge <- edge
    parentnode <- vars$nodeparents[curedge]
    endhost <- vars$nodehosts[curedge]
    if (parentnode > 0){
      parenthost <- vars$nodehosts[parentnode]
      while(endhost == parenthost) {
        curedge <- parentnode
        parentnode <- vars$nodeparents[curedge]
        endhost <- vars$nodehosts[curedge]
        parenthost <- max(vars$nodehosts[parentnode], 0)
      }
    }
    newnode <- length(vars$nodehosts) + 1
    vars$nodetypes[newnode] <- "t"
    vars$nodetimes[newnode] <- vars$inftimes[endhost]
    vars$nodehosts[newnode] <- vars$infectors[endhost]
    vars$nodeparents[newnode] <- parentnode
    vars$nodeparents[curedge] <- newnode
  }
  
  xthextra <- length(vars$nodetypes)
  for(edge in (obs + 1):xthextra) {
    curedge <- edge
    parentnode <- max(vars$nodeparents[curedge], 1)
    if(vars$nodetypes[parentnode] == "c") {
      endhost <- vars$nodehosts[curedge]
      parenthost <- vars$nodehosts[parentnode]
      while(parenthost != endhost) {
        newnode <- xthextra <- xthextra + 1L
        vars$nodetypes[newnode] <- "b"
        vars$nodetimes[newnode] <- vars$inftimes[endhost]
        vars$nodehosts[newnode] <- vars$infectors[endhost]
        vars$nodeparents[newnode] <- parentnode
        vars$nodeparents[curedge] <- newnode
        curedge <- newnode
        endhost <- vars$nodehosts[curedge]
      }
    }
  }
  
  vars$tree <- sapply(1:obs, function(i) tail(.ptr(vars$infectors, i), 1))
  
  return(vars)
} 

environment2phybreak <- function(varsenv) {
  varlength <- sum(varsenv$nodetypes %in% c("s", "x", "c"))
  
  np <- c(1, 1 + varsenv$nodeparents)
  np[np == 0] <- which(np == 0)
  while(any(np[1:(varlength + 1)] > 1 + varlength)) {
    np[np > 1 + varlength] <- np[np][np > 1 + varlength]
  }
  
  varsenv$nodeparents <- (np - 1)[-1][1:varlength]
  varsenv$nodehosts <- varsenv$nodehosts[1:varlength]
  varsenv$nodetimes <- varsenv$nodetimes[1:varlength]
  varsenv$nodetypes <- varsenv$nodetypes[1:varlength]
  
  return(varsenv)
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

### random infectors given infection times and generation interval distribution
.rinfectors <- function(it, introductions, p, v) {
  ### tests
  if(class(it) != "numeric" && class(it) != "integer") {
    stop(".rinfectors called with non-numeric infection times")
  }
  if(sum(it == min(it)) > introductions) stop("rinfectors with too many index cases")
  if(p$trans.model ==  "gamma"){
    if(p$gen.mean <= 0) stop(".rinfectors called with non-positive mean generation interval")
    if(p$gen.shape <= 0) stop(".rinfectors called with non-positive shape parameter")
  }
  
  ### function body
  res <- rep(0, length(it))
  for(i in 1:length(it)) {
    if(it[i] > min(it)) {
      dist <- infect_distribution(it[i], it, 
                                  cultimes = v$cultimes,
                                  nodetimes = v$nodetimes[1:length(it)], 
                                  p = p)
      dist[i] <- 0
      res[i] <- sample(length(it), 1, prob = dist)
    }
  }
  # add additional introductions
  res[sample((1:length(it))[it > min(it)], introductions - 1)] <- 0
  return(res)
}


make_nodehosts <- function(phybreakenv) {
  whmodel2 <- all(phybreakenv$v$nodetimes[phybreakenv$v$nodeparents[phybreakenv$v$nodetypes == "s"]] <= 
                    phybreakenv$v$inftimes) + 1
  for(i in 1:length(phybreakenv$v$nodetypes %in% c("s", "x"))) {
    nodepath <- .ptr(as.integer(phybreakenv$v$nodeparents), i)
    hostpath <- c(.ptr(phybreakenv$v$infectors, phybreakenv$v$nodehosts[i]), 0)
    hostpathtimes <- c(phybreakenv$v$inftimes[hostpath], -Inf) + 
      switch(whmodel2,
             1e-10,
             -1e-10)
    phybreakenv$v$nodehosts[nodepath] <- 
      invisible(sapply(nodepath, function(x) hostpath[hostpathtimes < phybreakenv$v$nodetimes[x]][1]))
  }
}

# pullnodes_complete_multiple_introductions <- function(currentID) {
#   nodes <- which(pbe1$v$nodehosts == currentID)
#   transnodes <- which(pbe1$v$nodetypes == "t")
#   transnode.currentID <- transnodes[currentID]
#   if (any(pbe1$v$infectors == currentID) | length(nodes) > 1){
#     innodes <- pbe1$v$infectors == currentID
#     edgesin <- nodes
#     edgeintimes <- pbe1$v$nodetimes[nodes]
#     transnode.inftime <- pbe1$v$nodetimes[transnodes[currentID]]
#     if (currentID == 1){
#       coaltimes <- sapply(edgeintimes, function(t) runif(1, transnode.inftime, t))
#       coaltimes <- coaltimes[-which(coaltimes == max(coaltimes))]
#     } else {
#       coaltimes <- sample_coaltimes(edgeintimes, transnode.inftime, pbe1$p)
#     }
#     
#     if(length(coaltimes) >= 1){
#       coalnodes <- which(pbe1$v$nodetypes == "c" &
#                            pbe1$v$nodehosts == -1)[1:length(coaltimes)]
#       pbe1$v$nodehosts[coalnodes] <- currentID
#       pbe1$v$nodetimes[coalnodes] <- coaltimes
#     } else {
#       coalnodes <- c()
#     }
#     
#     nodeorder <- order(c(coaltimes, edgeintimes))
#     edgeend <- c(coalnodes, edgesin)[nodeorder]
#     edgeendtimes <- c(coaltimes, edgeintimes)[nodeorder]
#     
#     # sample topology of minitree within hostID
#     edgestart <- sample_topology(edgeend,
#                                  edgeendtimes,
#                                  c(rep("c", length(coaltimes)),
#                                    rep("x", length(edgeintimes)))[nodeorder],
#                                  transnode.currentID)
#     
#     pbe1$v$nodeparents[edgeend] <- edgestart
#     if (pbe1$v$nodeparents[transnode.currentID] == -1 &
#         pbe1$v$infectors[currentID] == 0)
#       pbe1$v$nodeparents[transnode.currentID] <- 0
#   } else {
#     pbe1$v$nodeparents[currentID-1] <- transnodes[currentID]
#     if (pbe1$v$nodeparents[transnode.currentID] == -1 &
#         pbe1$v$infectors[currentID] == 0)
#       pbe1$v$nodeparents[transnode.currentID] <- 0
#   }
# }
