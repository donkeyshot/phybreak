### functions to simulate mini-trees ###


### sample coalescent times in a host, given the tip times since infection,
### ...the wh-model and the slope (used if WHmodel = 3)
.samplecoaltimes <- function(tleaves, WHmodel = 3, slope = 1, ngenes = 1, reassortment = FALSE) {
  ### tests
  if(min(tleaves) < 0) stop(".samplecoaltimes with negative tip times")
  if(!any(WHmodel == 1:3)) stop(paste0(".samplecoaltimes called with WHmodel = ",WHmodel))
  if(WHmodel == 3 && slope < 0) stop(".samplecoaltimes called with negative slope")
  
  ### function body
  if(length(tleaves) < 2) return(c())
  
  switch(
    WHmodel,
    #coalescence at transmission
    return(t(replicate(ngenes, head(sort(tleaves),-1)))),
    #coalescence at infection
    return(matrix(0, nrow = ngenes, ncol = length(tleaves) - 1)),
    {
      # transform times so that fixed rate 1 can be used
      ttrans <- sort(log(tleaves)/(slope), decreasing = TRUE)
      if(reassortment) {
        tnodetrans <- t(replicate(ngenes, .sctwh3(ttrans)))
      } else {
        tnodetrans <- matrix(rep(.sctwh3(ttrans), each = ngenes), nrow = ngenes)
      }
      
      res <- t(matrix(apply(exp(slope*tnodetrans), 1, sort), ncol = ngenes))
      # make sure that all branches will have positive length
      res <- pmax(res, matrix(rep(min(10^-6,tleaves/length(tleaves)) * (1:ncol(res)), each = ngenes), nrow = ngenes))
      
      return(res)
    }
    
  )
}


### sample tree topology in a host, given the node IDs, ntimes, and types,
### ...the root node and the WHmodel
.sampletopology <- function(nIDs, ntimes, ntypes, rootnode, WHmodel = 3, reassortment = FALSE) {
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
      leafparents <- leafparents[rank(ntimes[1, ntypes != "c"], ties.method="first")]
      res <- c(head(leafparents,sum(ntypes == "s")),
               cnodeparents,
               tail(leafparents,sum(ntypes == "t")))
      return(t(replicate(nrow(ntimes), res)))
    },
    #coalescence at infection
    {
      cnodes <- sort(nIDs[ntypes=="c"], decreasing = TRUE)
      res <- c(rep(NA,sum(ntypes=="s")),
               rootnode,tail(-nIDs[ntypes=="c"],-1),
               rep(NA,sum(ntypes=="t")))
      for(i in cnodes) {
        res[sample(which(is.na(res)),2)] <- i
        res[res == -i] <- NA
      }
      return(t(replicate(nrow(ntimes), res)))
    }
  )
  if(reassortment) {
    IDs <- t(apply(ntimes, 1, function(x) nIDs[order(x, ntypes)]))
    tys <- t(apply(ntimes, 1, function(x) ntypes[order(x, ntypes)]))
    if(any(tys[, 1] != "c")) {
      print(c(nIDs, ntimes, ntypes, rootnode))
      stop("host topology does not start with coalescence node")
    }
    res <- matrix(rootnode, nrow = nrow(ntimes), ncol = length(nIDs))
    for(gene in 1:nrow(res)) {
      tochoose <- rep(IDs[gene, 1], 2)
      for(i in 2:length(nIDs)) {
        res[gene, i] <- tochoose[1]
        if(tys[gene, i] == "c") {
          tochoose <- sample(c(tochoose[-1], IDs[gene, i], IDs[gene, i]))
        } else {
          tochoose <- tochoose[-1]
        }
      }
    }
    res <- (t(sapply(1:nrow(ntimes), function(x) res[x, order(IDs[x, ])])))
  } else {
    IDs <- nIDs[order(ntimes[1, ], ntypes)]
    tys <- ntypes[order(ntimes[1, ], ntypes)]
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
    res <- t(replicate(nrow(ntimes), res[order(IDs)]))
  }
  return(res)
}
