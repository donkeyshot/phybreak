### functions to simulate mini-trees ###

sample_coaltimes <- function(tiptimes, parameters) {
  ### tests
  if(min(tiptimes) < 0) stop("sample_coaltimes with negative tip times")

  ### function body
  if(length(tiptimes) < 2) return(c())
  
  switch(
    parameters$wh.model,
    #coalescence at transmission
    single = head(sort(tiptimes), -1),
    #coalescence at infection
    infinite = 0*tiptimes[-1],
    #linear increase
    linear = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(log(tiptimes)/(parameters$wh.slope), decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- sort(exp(parameters$wh.slope * tnodetrans))
      # make sure that all branches will have positive length
      res <- apply(cbind(res,
                         min(10^-5,tiptimes/length(tiptimes))*(1:length(res))), 1, max)
      
      return(res)
    },
    #exponential increase
    exponential = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(-exp(-parameters$wh.exponent * tiptimes)/(parameters$wh.exponent * parameters$wh.level), 
                     decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- sort(-log(-parameters$wh.exponent * parameters$wh.level * tnodetrans)/(parameters$wh.exponent))

      return(res)
    },
    #constant level
    constant = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(tiptimes/parameters$wh.level, decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- sort(parameters$wh.level * tnodetrans)
      
      return(res)
    }
  )
}

sample_topology <- function(nodeIDs, nodetimes, nodetypes, infectornodes) {
  res <- rep(-1, length(nodeIDs))
  tochoose <- infectornodes
  for(i in 1:length(nodeIDs)) {
    res[i] <- tochoose[1]
    if(nodetypes[i] == "c") {
      tochoose <- sample(c(tochoose[-1], nodeIDs[i], nodeIDs[i]))
    } else {
      tochoose <- tochoose[-1]
    }
  }
  return(res)
}

sample_singlecoaltime <- function(oldtiptimes, oldcoaltimes, newtiptime, parameters) {
  ### return -Inf if there is no existing tree
  if(length(oldtiptimes) + length((oldcoaltimes)) == 0) return(-Inf)
  
  ### add 0s to prevent problems with empty vectors
  oldtiptimes <- c(0, oldtiptimes)
  oldcoaltimes <- c(0, oldcoaltimes)
  
  ### function body
  switch(
    parameters$wh.model,
    #coalescence at transmission
    single = return(min(newtiptime, max(oldtiptimes))),
    #coalescence at infection
    infinite = return(0),
    linear = {
      # transform times so that fixed rate 1 can be used
      transtiptimes <- log(oldtiptimes)/parameters$wh.slope
      transcoaltimes <- log(oldcoaltimes)/parameters$wh.slope
      transcurtime <- log(newtiptime)/parameters$wh.slope
    },
    exponential = {
      transtiptimes <- -exp(-parameters$wh.exponent * oldtiptimes)/(parameters$wh.exponent * parameters$wh.level)
      transcoaltimes <- -exp(-parameters$wh.exponent * oldcoaltimes)/(parameters$wh.exponent * parameters$wh.level)
      transcurtime <- -exp(-parameters$wh.exponent * newtiptime)/(parameters$wh.exponent * parameters$wh.level)
    },
    constant = {
      transtiptimes <- oldtiptimes/parameters$wh.level
      transcoaltimes <- oldcoaltimes/parameters$wh.level
      transcurtime <- newtiptime/parameters$wh.level
    }
  )

  # start with number of edges at current time and time of next node (backwards next)
  curnedge <- sum(transtiptimes >= transcurtime) - sum(transcoaltimes >= transcurtime)
  transnexttime <- max(c(-Inf, transtiptimes[transtiptimes < transcurtime],  
                         transcoaltimes[transcoaltimes < transcurtime]))
  
  # traverse minitree node by node, subtracting coalescence rate untile cumulative rate reaches 0
  frailty <- stats::rexp(1)
  while(curnedge * (transcurtime - transnexttime) < frailty) {
    transcurtime <- transnexttime
    curnedge <- sum(transtiptimes >= transcurtime) - sum(transcoaltimes >= transcurtime)
    transnexttime <- max(c(-Inf, transtiptimes[transtiptimes < transcurtime],  
                           transcoaltimes[transcoaltimes < transcurtime]))
    frailty <- stats::rexp(1)
  }
  
  # calculate transformed node time
  transreturntime <- transcurtime - frailty / curnedge
  
  # transform to real time
  switch(
    parameters$wh.model, single =, infinite = ,
    linear = exp(parameters$wh.slope * transreturntime),
    exponential = -log(-parameters$wh.exponent * parameters$wh.level * transreturntime)/(parameters$wh.exponent),
    constant = parameters$wh.level * transreturntime
  )
}

sample_singlechildnode <- function(nodeIDs, nodeparents, nodetimes, newnodetime) {
  candidatenodes <- nodeIDs[nodetimes >= newnodetime & 
                           c(-Inf, nodetimes)[1 + match(nodeparents, nodeIDs, nomatch = 0)] < newnodetime]
  if(length(candidatenodes) == 0) candidatenodes <- nodeIDs[nodetimes == max(nodetimes)]
  return(sample(rep(candidatenodes, 2), 1))
}


### sample coalescent times in a host, given the tip times since infection,
### ...the wh-model and the slope (used if WHmodel = 3)
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
.sampletopology <- function(nIDs, ntimes, ntypes, rootnode, WHmodel = 3) {
  ### tests
  if(!any(WHmodel == 1:3)) stop(paste0(".sampletopology called with WHmodel = ",WHmodel))
  
  ### function body
  if(length(nIDs) == 1) return(rootnode)
  switch(
    WHmodel,
    #coalescence at transmission
    {
      cnodes <- nIDs[ntypes == "c"]
      cnodeparents <- c(rootnode, head(cnodes, -1))
      leafparents <- c(cnodes,tail(cnodes,1))
      leafparents <- leafparents[rank(ntimes[ntypes != "c"],ties.method="first")]
      res <- c(head(leafparents,sum(ntypes %in% c("s", "x"))),
               cnodeparents,
               tail(leafparents,sum(ntypes=="t")))
      return(res)
    },
    #coalescence at infection
    {
      cnodes <- sort(nIDs[ntypes=="c"],decreasing = TRUE)
      res <- c(rep(NA,sum(ntypes %in% c("s", "x"))),
               rootnode,tail(-nIDs[ntypes=="c"],-1),
               rep(NA,sum(ntypes=="t")))
      for(i in cnodes) {
        res[sample(which(is.na(res)),2)] <- i
        res[res == -i] <- NA
      }
      return(res)
    }
  )
  IDs <- nIDs[order(ntimes, ntypes)]
  tys <- ntypes[order(ntimes, ntypes)]
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



.sampleextracoaltime <- function(ntimes, ntypes, newtiptime, WHmodel = 3, slope = 1) {
  ### tests
  if(min(ntimes) < 0) stop(".sampleextracoaltime with negative node times")
  if(!any(WHmodel == 1:3)) stop(paste0(".sampleextracoaltime called with WHmodel = ",WHmodel))
  if(WHmodel == 3 && slope < 0) stop(".sampleextracoaltime called with negative slope")
  
  ### function body
  switch(
    WHmodel,
    #coalescence at transmission
    return(min(newtiptime, max(ntimes))),
    #coalescence at infection
    return(0),
    {
      # transform times so that fixed rate 1 can be used
      transntimes <- log(ntimes)/slope
      transcurtime <- log(newtiptime)/slope
      
      # random cumulative coalescence rate
      rcumcoalrate <- stats::rexp(1)
      
      # start with number of edges at current time and time of next node (backwards next)
      curnedge <- sum(transntimes - transcurtime >= 0 & ntypes != "c") - sum(transntimes - transcurtime >= 0 & ntypes == "c")
      transnexttime <- max(c(-Inf, transntimes[transntimes < transcurtime]))
      
      # traverse minitree node by node, subtracting coalescence rate untile cumulative rate reaches 0
      while(curnedge * (transcurtime - transnexttime) < rcumcoalrate) {
        rcumcoalrate <- rcumcoalrate - curnedge * (transcurtime - transnexttime)
        transcurtime <- transnexttime
        curnedge <- sum(transntimes - transcurtime >= 0 & ntypes != "c") - sum(transntimes - transcurtime >= 0 & ntypes == "c")
        transnexttime <- max(c(-Inf, transntimes[transntimes < transcurtime]))
      }
      
      # calculate transformed node time, and transform to real time
      transreturntime <- transcurtime - rcumcoalrate / curnedge
      returntime <- exp(slope * transreturntime)
      
      return(returntime)
    }
  )
}

.sampleextraupnode <- function(nIDs, nparents, ntimes, newnodetime) {
  candidatenodes <- nIDs[ntimes > newnodetime & 
                           c(-Inf,ntimes)[1 + match(nparents, nIDs, nomatch = 0)] < newnodetime]
  if(length(candidatenodes) == 0) candidatenodes <- nIDs[ntimes == max(ntimes)]
  return(sample(rep(candidatenodes, 2), 1))
}

