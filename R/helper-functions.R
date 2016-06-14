### helper functions called from multiple phybreak functions ###


### sample coalescent times in a host, given the tip times since infection,
### ...the wh-model and the slope (used if WHmodel = 3)
### called from:
# phybreak
# .sim.phylotree
### calls:
# .sctwh3  ## C++ function
.samplecoaltimes <- function(tleaves, WHmodel = 3, slope = 1) {
  ### tests
  if(min(tleaves) < 0) stop(".samplecoaltimes with negative tip times")
  if(!any(WHmodel == 1:3)) stop(paste0(".samplecoaltimes called with WHmodel = ",WHmodel))
  if(WHmodel == 3 && slope < 0) stop(".samplecoaltimes called with negative slope")
  
  ### function body
  if(length(tleaves) < 2) return(c())
  
  switch(
    WHmodel,
    #coalescence at transmission
    return(head(sort(tleaves),-1)),
    #coalescence at infection
    return(0*tleaves[-1]),
    {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(log(tleaves)/(slope), decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- sort(exp(slope*tnodetrans))
      # make sure that all branches will have positive length
      res <- apply(cbind(res,
                         min(10^-5,tleaves/length(tleaves))*(1:length(res))),1,max)
      
      return(res)
    }
    
  )
}


### sample tree topology in a host, given the node IDs, ntimes, and types,
### ...the root node and the WHmodel
### called from:
# phybreak
# .sim.phylotree
.sampletopology <- function(nIDs, ntimes, ntypes, rootnode, WHmodel = 3) {
  ### tests
  if(!any(WHmodel == 1:3)) stop(paste0(".sampletopology called with WHmodel = ",WHmodel))

  ### function body
  if(length(nIDs) == 1) return(rootnode)
  switch(
    WHmodel,
    #coalescence at transmission
    {
      cnodes <- nIDs[ntypes=="c"]
      cnodeparents <- c(rootnode,head(cnodes,-1))
      leafparents <- c(cnodes,tail(cnodes,1))
      leafparents <- leafparents[rank(ntimes[ntypes != "c"],ties.method="first")]
      res <- c(head(leafparents,sum(ntypes=="s")),
               cnodeparents,
               tail(leafparents,sum(ntypes=="t")))
      return(res)
    },
    #coalescence at infection
    {
      cnodes <- sort(nIDs[ntypes=="c"],decreasing = TRUE)
      res <- c(rep(NA,sum(ntypes=="s")),
               rootnode,tail(-nIDs[ntypes=="c"],-1),
               rep(NA,sum(ntypes=="t")))
      for(i in cnodes) {
        res[sample(which(is.na(res)),2)] <- i
        res[res == -i] <- NA
      }
      return(res)
    }
  )
  IDs <- nIDs[order(ntimes,ntypes)]
  tys <- ntypes[order(ntimes,ntypes)]
  if(tys[1] != "c") {
    print(c(nIDs, ntimes, ntypes, rootnode))
    stop("host topology does not start with coalescence node")
  }
  res <- rep(rootnode, length(nIDs))
  tochoose <- rep(IDs[1], 2)
  for(i in 2:length(nIDs)) {
    res[i] <- tochoose[1]
    if(tys[i] == "c") {
      tochoose <- sample(c(tochoose[-1], IDs[i], IDs[i]))
    } else {
      tochoose <- tochoose[-1]
    }
  }
  return(res[order(IDs)])
}


### make a class 'phylo' object from class 'phybreak' input
### called by:
# get.phylo
# sim.phybreak
### required packages:
# phangorn
.makephylo.phybreak <- function(nodetimes, nodeparents, nodenames) {
  ###topology
  Nhosts <- (1+length(nodetimes))/3
  indext <- (1:length(nodetimes))[nodeparents == 0]
  indexc <- (1:length(nodetimes))[nodeparents == indext]
  edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
  edgeend <- (1:length(nodetimes))[nodeparents != 0 & nodeparents != indext]
  edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]
  
  toremove <- (edgeend >= 2*Nhosts)
  toreplace <- (edgestart >= 2*Nhosts)
  edgelengths[toreplace] <- edgelengths[toreplace] + edgelengths[toremove][match(edgestart[toreplace],edgeend[toremove])]
  edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace],edgeend[toremove])]
  edgestart <- edgestart[!toremove]
  edgeend <- edgeend[!toremove]
  edgelengths <- edgelengths[!toremove]
  
  ###give first coalescent node number Nhosts+1
  if(indexc != Nhosts + 1) {
    edgestart[edgestart == indexc] <- 0
    edgeend[edgeend == indexc] <- 0
    edgestart[edgestart == Nhosts + 1] <- indexc
    edgeend[edgeend == Nhosts + 1] <- indexc
    edgestart[edgestart == 0] <- Nhosts + 1
    edgeend[edgeend == 0] <- Nhosts + 1
  }
  
  edges <- matrix(c(edgestart,edgeend),ncol=2)
  
  toname <- nodenames
  
  
  res <- list(
    edge = edges,
    edge.length = edgelengths,
    Nnode = Nhosts - 1,
    tip.label = toname
  )
  class(res) <- "phylo"
  res <- reorder(res)
  res <- ladderize(res)
  return(res)
  
}


### make a class 'phylo' object from class 'phybreak' input, including Simmap items
### called by:
# get.phylo
### required packages:
# phytools
# phangorn
.makephylosimmap.phybreak <- function(nodetimes, nodeparents, nodehosts, nodenames) {
  ###topology for phylo
  Nhosts <- (1+length(nodetimes))/3
  indext <- (1:length(nodetimes))[nodeparents == 0]
  indexc <- (1:length(nodetimes))[nodeparents == indext]
  edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
  edgeend <- (1:length(nodetimes))[nodeparents != 0 & nodeparents != indext]
  edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]
  
  toremove <- (edgeend >= 2*Nhosts)
  toreplace <- (edgestart >= 2*Nhosts)
  edgelengths[toreplace] <- edgelengths[toreplace] + edgelengths[toremove][match(edgestart[toreplace],edgeend[toremove])]
  edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace],edgeend[toremove])]
  edgestart <- edgestart[!toremove]
  edgeend <- edgeend[!toremove]
  edgelengths <- edgelengths[!toremove]
  
  edges <- matrix(c(edgestart,edgeend),ncol=2)
  
  ###items for simmap (phytools)  
  whichcolor <- function(iors, path) {
    if(path[1] == 0) {
      return(length(path) %% 3)
    } else {
      return(whichcolor(iors, c(iors[path[1]], path)))
    }
  }
  
  tipstates <- c("red","blue","black")[1+sapply(1:Nhosts, whichcolor, 
                                                iors=tail(nodehosts,Nhosts))]
  nodestates <- matrix(tipstates[nodehosts[edges]], ncol = 2)
  edgemaps <- list()
  for(i in 1:(2 * (Nhosts - 1))) {
    if(nodestates[i,1] == nodestates[i,2]) {
      edgemaps[[i]] <- edgelengths[i]
      names(edgemaps[[i]]) <- nodestates[i,1]
    } else {
      nodes <- edges[i,]
      ntimes <- nodetimes[c(nodes[1],nodeparents[nodes[2]],nodes[2])]
      edgemaps[[i]] <- ntimes[2:3] - ntimes[1:2]
      names(edgemaps[[i]]) <- nodestates[i,]
    }
  }
  
  mappededge <- matrix(nrow = 2 * (Nhosts - 1), ncol = 3)
  for(i in 1:(2 * (Nhosts - 1))) {
    mappededge[i,] <- edgemaps[[i]][c("blue","black","red")]
  }
  mappededge[is.na(mappededge)] <- 0
  colnames(mappededge) <- c("blue","black","red")
  rownames(mappededge) <- paste0(edges[,1],",",edges[,2])
  
  
  ###give first coalescent node number Nhosts+1
  if(indexc != Nhosts + 1) {
    edgestart[edgestart == indexc] <- 0
    edgeend[edgeend == indexc] <- 0
    edgestart[edgestart == Nhosts + 1] <- indexc
    edgeend[edgeend == Nhosts + 1] <- indexc
    edgestart[edgestart == 0] <- Nhosts + 1
    edgeend[edgeend == 0] <- Nhosts + 1
  }
  edges <- matrix(c(edgestart,edgeend),ncol=2)
  
  ###make preliminary result for ladderizing
  respreladder <- list(
    edge = edges,
    edge.length = edgelengths,
    Nnode = Nhosts - 1,
    tip.label = nodenames
  )
  class(respreladder) <- c("phylo")
  respreladder <- reorder(respreladder)
  
  ###ladderize to determine reordering of edges
  edgeladder <- ladderize(respreladder)$edge[,2]
  newedgeorder <- match(edgeladder, edgeend)
  newtiporder <- match(respreladder$tip.label, nodenames)
  
  reorderedmaps <- list()
  for(i in 1:(2*(Nhosts-1))) {
    reorderedmaps[[i]] <- edgemaps[[newedgeorder[i]]]
  }
  
  
  res <- list(
    edge = edges[newedgeorder,],
    edge.length = edgelengths[newedgeorder],
    Nnode = Nhosts - 1,
    tip.label = nodenames[newtiporder],
    node.state = nodestates[newedgeorder,],
    states = tipstates[newtiporder],
    maps = reorderedmaps,
    mapped.edge = mappededge[newedgeorder,]
  )
  class(res) <- c("simmap","phylo")
  attr(res,"order") <- "cladewise"
  
  return(res)
  
}

### calculate the log-likelihood of sampling intervals
### called from:
# logLik.phybreak
# .build.phybreakenv
# .propose.phybreakenv
.lik.gentimes <- function(obs, shapeG, meanG, nodetimes, nodehosts, nodetypes) {
  sum(dgamma(nodetimes[nodetypes == "t" & nodehosts > 0] -
               nodetimes[nodehosts[nodetypes == "t" & nodehosts > 0] + 2*obs -1],
             shape = shapeG, scale = meanG/shapeG, log=TRUE))
}

### calculate the log-likelihood of generation intervals
### called from:
# logLik.phybreak
# .build.phybreakenv
# .propose.phybreakenv
.lik.sampletimes <- function(shapeS, meanS, nodetimes, nodetypes) {
  sum(dgamma(nodetimes[nodetypes == "s"] -
               nodetimes[nodetypes == "t"],
             shape = shapeS, scale = meanS/shapeS, log=TRUE))
}

### calculate the log-likelihood of coalescent intervals
### called from:
# logLik.phybreak
# .build.phybreakenv
# .propose.phybreakenv
.lik.coaltimes <- function(obs, wh.model, slope, nodetimes, nodehosts, nodetypes) {
  if(wh.model == 1 || wh.model == 2) return(0)
  
  coalnodes <- nodetypes == "c"
  orderednodes <- order(nodehosts, nodetimes)
  orderedtouse <- orderednodes[c(duplicated(nodehosts[orderednodes])[-1], FALSE)]
  #only use hosts with secondary infections
  
  ##make vectors with information on intervals between nodes
  coalno <- c(FALSE, head(coalnodes[orderedtouse],-1)) #interval starts with coalescence
  nodeho <- nodehosts[orderedtouse] #host in which interval resides
  coalmultipliers <- choose(2 + cumsum(2*coalno - 1),2) #coalescence coefficient
 
  ##from t to tau (time since infection)
  whtimes <- nodetimes - c(0,tail(nodetimes,obs))[1+nodehosts]
  
  noderates <- 1/(slope*whtimes[orderedtouse])  
  #coalescence rate (per pair of lineages)
  nodeescrates <- log(whtimes[orderedtouse])/(slope)
  #cumulative coalescence rate since infection of host (per pair of lineages)
  

  escratediffs <- nodeescrates - c(0, head(nodeescrates,-1))
  escratediffs[!duplicated(nodeho)] <- nodeescrates[!duplicated(nodeho)]
  #cumulative coalescence rate within interval (per pair of lineages)
  

  ##First: coalescence rates at coalescence nodes
  logcoalrates <- log(noderates[c(coalno[-1],FALSE)])
  
  #Second: probability to escape coalescence in all intervals
  logescapes <- -escratediffs*coalmultipliers
  
  
  return(sum(logcoalrates) + sum(logescapes))
  
}
