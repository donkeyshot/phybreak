### functions to simulate mini-trees ###

sample_coaltimes <- function(tiptimes, inftime, parameters) {
  ### tests
  if(min(tiptimes) < inftime) stop("sample_coaltimes with negative tip times")

  ### function body
  if(length(tiptimes) < 2) return(c())
  
  switch(
    parameters$wh.model,
    #coalescence at transmission
    single = head(sort(tiptimes), -1),
    no_coalescent = head(sort(tiptimes), -1),
    #coalescence at infection
    infinite = inftime + 0*tiptimes[-1],
    #linear increase
    linear = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(log(parameters$wh.level/parameters$wh.slope + tiptimes - inftime)/(parameters$wh.slope), decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- sort(exp(parameters$wh.slope * tnodetrans))
      res <- apply(cbind(res,
                         min(10^-5, (parameters$wh.level/parameters$wh.slope + tiptimes - inftime)/length(tiptimes))),
                   1, max)
      res <- res + inftime - parameters$wh.level/parameters$wh.slope
      
      return(res)
    },
    #exponential increase
    exponential = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(-exp(-parameters$wh.exponent * (tiptimes - inftime))/(parameters$wh.exponent * parameters$wh.level), 
                     decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- inftime + sort(-log(-parameters$wh.exponent * parameters$wh.level * tnodetrans)/(parameters$wh.exponent))

      return(res)
    },
    #constant level
    constant = {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort((tiptimes - inftime)/parameters$wh.level, decreasing = TRUE)
      tnodetrans <- .sctwh3(ttrans)
      
      res <- inftime + sort(parameters$wh.level * tnodetrans)
      
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

sample_singlecoaltime <- function(oldtiptimes, oldcoaltimes, newtiptime, inftime, parameters) {
  ### return -Inf if there is no existing tree
  if(length(oldtiptimes) + length((oldcoaltimes)) == 0) return(-Inf)
  
  ### add 0s to prevent problems with empty vectors
  oldtiptimes <- c(inftime, oldtiptimes)
  oldcoaltimes <- c(inftime, oldcoaltimes)
  
  ### function body
  switch(
    parameters$wh.model,
    #coalescence at transmission
    single = return(min(newtiptime, max(oldtiptimes))),
    no_coalescent = return(min(newtiptime, max(oldtiptimes))),
    #coalescence at infection
    infinite = return(inftime),
    linear = {
      # transform times so that fixed rate 1 can be used
      transtiptimes <- log(parameters$wh.level/parameters$wh.slope + oldtiptimes - inftime)/parameters$wh.slope
      transcoaltimes <- log(parameters$wh.level/parameters$wh.slope + oldcoaltimes - inftime)/parameters$wh.slope
      transcurtime <- log(parameters$wh.level/parameters$wh.slope + newtiptime - inftime)/parameters$wh.slope
    },
    exponential = {
      transtiptimes <- -exp(-parameters$wh.exponent * (oldtiptimes - inftime))/(parameters$wh.exponent * parameters$wh.level)
      transcoaltimes <- -exp(-parameters$wh.exponent * (oldcoaltimes - inftime))/(parameters$wh.exponent * parameters$wh.level)
      transcurtime <- -exp(-parameters$wh.exponent * (newtiptime - inftime))/(parameters$wh.exponent * parameters$wh.level)
    },
    constant = {
      transtiptimes <- (oldtiptimes - inftime)/parameters$wh.level
      transcoaltimes <- (oldcoaltimes - inftime)/parameters$wh.level
      transcurtime <- (newtiptime - inftime)/parameters$wh.level
    }
  )

  # start with number of edges at current time and time of next node (backwards next)
  curnedge <- sum(transtiptimes >= transcurtime) - sum(transcoaltimes >= transcurtime)
  transnexttime <- max(c(-Inf, transtiptimes[transtiptimes < transcurtime],  
                         transcoaltimes[transcoaltimes < transcurtime]))
  
  # traverse minitree node by node, subtracting coalescence rate until cumulative rate reaches 0
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
    parameters$wh.model, single =, no_coalescent = , infinite = ,
    linear = {
      res <- min(c(10^-5, (oldcoaltimes[-1] - inftime)/2))
      res <- max(res, exp(parameters$wh.slope * transreturntime))
      res <- inftime + res - parameters$wh.level/parameters$wh.slope
      return(res)
      },
    exponential = inftime - log(-parameters$wh.exponent * parameters$wh.level * transreturntime)/(parameters$wh.exponent),
    constant = inftime + parameters$wh.level * transreturntime
  )
}

sample_singlechildnode <- function(nodeIDs, nodeparents, nodetimes, newnodetime) {
  candidatenodes <- nodeIDs[nodetimes >= newnodetime & 
                           c(-Inf, nodetimes)[1 + match(nodeparents, nodeIDs, nomatch = 0)] < newnodetime]
  if(length(candidatenodes) == 0) candidatenodes <- nodeIDs[nodetimes == max(nodetimes)]
  return(sample(rep(candidatenodes, 2), 1))
}


### Should I ever consider weighted topology sampling, the topology should be sampled backwards in time instead of forwards.
### This function does exactly that and is slightly less efficient than the now-used forward sampling
# sample_topology_backwards <- function(nodeIDs, nodetypes, infectornodes) {
#   res <- rep(-1, length(nodeIDs))
#   tochoose <- c()
#   for(i in length(nodeIDs):1) {
#     if(nodetypes[i] == "c") {
#       tojoin <- sample(length(tochoose), 2)
#       res[tochoose[tojoin]] <- nodeIDs[i]
#       tochoose <- c(tochoose[-tojoin], i)
#     } else {
#       tochoose <- c(tochoose, i)
#     }
#   }
#   res[res == -1] <- infectornodes
#   return(res)
# }
