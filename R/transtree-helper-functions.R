### function to obtain the most likely infector for each host, possibly containing cycles and multiple roots,
### plus (if requested) means and standard deviations of infection times
### called by:
# transtree
### calls:
# .postinfector
.transtreecount <- function(phybreak.object, samplesize, includetimes) {
  ### initialize some constants
  chainlength <- length(phybreak.object$s$mu)
  obsize <- phybreak.object$p$obs
  samplerange <- (chainlength-samplesize+1):chainlength
  
  ### get result
  res <- t(matrix(with(phybreak.object,
                       apply(s$nodehosts[obsize:(2*obsize-1),samplerange],
                             1,.postinfector,support = TRUE)),nrow = 2))
  
  ### get time summaries
  if(includetimes) {
  timesums <- with(phybreak.object,
                   rowSums(s$nodetimes[obsize:(2*obsize-1),samplerange]*(
                     s$nodehosts[obsize:(2*obsize-1),samplerange] == res[,1]
                   )))  
  timesumsqs <- with(phybreak.object,
                     rowSums((s$nodetimes[obsize:(2*obsize-1),samplerange]^2)*(
                       s$nodehosts[obsize:(2*obsize-1),samplerange] == res[,1]
                     )))  
  res <- cbind(res, timesums/res[,2])
  res <- cbind(res, sqrt((timesumsqs - res[,3]^2*res[,2])/(res[,2]-1)))
  }
  
  ### return result
  return(res)
}


### function to obtain a transmission tree based on most likely infectors, using Edmonds's algorithm 
### to remove cycles, applied with multiple root candidates, 
### plus (if requested) means and standard deviations of infection times
### called by
# transtree
### calls:
# .edmondsiterative
.transtreeedmonds <- function(phybreak.object, samplesize, includetimes) {
  ### initialize some constants
  chainlength <- length(phybreak.object$s$mu)
  obsize <- phybreak.object$p$obs
  samplerange <- (chainlength-samplesize+1):chainlength
  
  ### obtaining the result in steps
  
  # matrix with support for each infector (row) per host (column), with 0 as maximum, and a column for the index
  supportmatrix <- cbind(
    c(0, rep(-1,obsize)),
    apply(1 + phybreak.object$s$nodehosts[obsize:(2*obsize-1), samplerange],
          1, tabulate, nbins = obsize + 1))
  supportmatrix <- supportmatrix - rep(apply(supportmatrix, 2, max), each = obsize + 1)
  
  # vector with hosts, ordered by support to be index, and vector with these supports
  candidateindex <- order(apply(
    phybreak.object$s$nodehosts[obsize:(2*obsize-1), samplerange] == 0,
    1, sum), decreasing = TRUE)
  indexsupports <- apply(
    phybreak.object$s$nodehosts[obsize:(2*obsize-1), samplerange] == 0,
    1, sum)[candidateindex]
  
  ## index host candidate by candidate, make a tree with edmonds's algorithm
  ## after each tree, test for each remaining candidate if their support-to-be-index
  ## is smaller than the support for their infector in the last tree.
  ## If so, stop making new trees, as no better tree will come out. Then select the best tree.
  # first initialize some variables
  bestYN <- FALSE
  nextcandidate <- 1
  alltrees <- c()
  allsupports <- c()
  # then make trees as long as bestYN == FALSE
  while(!bestYN) {
    # make copy of supportmatrix, maximally supporting the next candidate index
    suppmat <- supportmatrix
    suppmat[1,] <- -samplesize
    suppmat[1,candidateindex[nextcandidate] + 1] <- 0
    
    # make the tree
    thistree <- .edmondsiterative(suppmat, samplesize, obsize)[-1] - 1
    alltrees <- c(alltrees, thistree)
    thissupport <- rowSums(phybreak.object$s$nodehosts[obsize:(2*obsize-1), 
                                                                  samplerange] == thistree)              
    allsupports <- c(allsupports, thissupport)

    # test if any unused candidate index has higher index support than infector support in last tree
    bestYN <- TRUE                 
    for(i in tail(candidateindex,-nextcandidate)) {
      if(thissupport[i] < indexsupports[candidateindex==i]) {
        bestYN <- FALSE
      }
    }
    nextcandidate <- nextcandidate + 1
  }
  
  # find the tree with maximum support
  dim(alltrees) <- c(obsize,length(alltrees)/obsize)
  dim(allsupports) <- dim(alltrees)
  besttree <- alltrees[,which(colSums(allsupports)==max(colSums(allsupports)))]
  treesupport <- allsupports[,which(colSums(allsupports)==max(colSums(allsupports)))]
  res <- matrix(c(besttree,
                  treesupport),
                ncol=2)
  
  ### get time summaries
  if(includetimes) {
    timesums <- with(phybreak.object,
                     rowSums(s$nodetimes[obsize:(2*obsize-1),samplerange]*(
                       s$nodehosts[obsize:(2*obsize-1),samplerange] == res[,1]
                     )))  
    timesumsqs <- with(phybreak.object,
                       rowSums((s$nodetimes[obsize:(2*obsize-1),samplerange]^2)*(
                         s$nodehosts[obsize:(2*obsize-1),samplerange] == res[,1]
                       )))  
    res <- cbind(res, timesums/res[,2])
    res <- cbind(res, sqrt((timesumsqs - res[,3]^2*res[,2])/(res[,2]-1)))
  }
  
  ### return result
  return(res)
}


### function to obtain a transmission tree based on most likely infectors, selecting the tree 
### among posterior trees with highest support,
### plus (if requested) means and standard deviations of infection times
### called by
# transtree
### calls:
#return 'maximum parent credibility' tree with posterior support and means and standard deviations of 
#infection times, and infection times of the mpc tree itself
.mpcinfector <- function(phybreak.object, samplesize, phylo.class = FALSE, includetimes = FALSE) {
  ### initialize some constants
  chainlength <- length(phybreak.object$s$mu)
  obsize <- phybreak.object$p$obs
  samplerange <- (chainlength-samplesize+1):chainlength
  posteriorsamples <- phybreak.object$s$nodehosts[obsize:(2*obsize-1), samplerange]

  ### obtaining the result in steps

  # matrix containing posterior support for each infector in each tree
  allsupports <- matrix(0, nrow = obsize, ncol=samplesize)
  for(hostID in 1:obsize) {
    for(infector in 0:obsize) {
      allsupports[hostID, posteriorsamples[hostID,] == infector] <- 
        sum(posteriorsamples[hostID,] == infector)
    }
  }
  
  # 
  besttree <- which.max(colSums(log(allsupports)))
  if(phylo.class) return(besttree + samplerange[1] - 1)
  res <- matrix(
    c(
      posteriorsamples[, besttree],
      allsupports[, besttree]
    ),ncol = 2)
  
  if(includetimes) {
    timesums <- with(phybreak.object,
                     rowSums(s$nodetimes[obsize:(2*obsize-1),samplerange]*(
                       s$nodehosts[obsize:(2*obsize-1),samplerange] == res[,1]
                     )))  
    timesumsqs <- with(phybreak.object,
                       rowSums((s$nodetimes[obsize:(2*obsize-1),samplerange]^2)*(
                         posteriorsamples == res[,1]
                       )))  
    res <- cbind(res, timesums/res[,2])
    res <- cbind(res, sqrt((timesumsqs - res[,3]^2*res[,2])/(res[,2]-1)))
    res <- cbind(res, phybreak.object$s$nodetimes[obsize:(2*obsize-1),
                                                  bestpars + samplerange[1] - 1])
  }
  
  
  return(res)
}




### function returning the ranknr'th most frequent entry in inf.chain, plus the frequency if support = TRUE
### called by:
# .transtreecount
# .inflist (function 'infectorsets')
.postinfector <- function(inf.chain, ranknr = 1, support = FALSE) {
  chainlength <- length(inf.chain)
  ordinf <- rle(sort(inf.chain))
  if(length(ordinf$lengths) < ranknr) {
    if(support) return(c(0,0)) else return(0)
  }
  ord <- order(ordinf$lengths,decreasing = TRUE)
  if(support) {
    return(c(ordinf[[2]][ord[ranknr]],ordinf[[1]][ord[ranknr]]))
  } else {
    return(ordinf[[2]][ord[ranknr]])
  }
}


### function iteratively removing cycles by edmond's algorithm, as long as there are cycles
### called by:
# .transtreeedmonds
# .edmondsiterative
### calls:
# .cycleYN
# .pathtocycle
# .edmondsiterative
.edmondsiterative <- function(suppmat, samplesize, obs) {
  # get most likely infectors, determine which hosts are in cycles,
  # return ML infectors if there are no cycles
  parentset <- apply(suppmat, 2, which.max)
  cycleIDs <- which(sapply(1:obs, 
                           .cycleYN, 
                           parset = parentset))
  if(length(cycleIDs) == 0) {
    return(parentset)
  }
  
  # take one cycle, and infector for each member
  cycletonode <- unique(.pathtocycle(parentset, cycleIDs[1]))
  cycleinfector <- parentset[cycletonode]
  
  # determine for each host which cycle member is their most likely infector,
  # and for each host which cycle member would be their infectee
  whichincoming <- cycletonode[apply(suppmat[cycletonode,],2,which.max)]
  whichoutgoing <- cycletonode[apply(suppmat[,cycletonode],1,which.max)]

  ## turn the cycle into a single host: in the suppmat, use the position of host #1 for this cycle
  ## move all infector support for any cycle member to host #1
  ## remove all infector support for other hosts 
  # let all infectee-candidates be infected by the cycle (entered into position of host #1)
  suppmat[cycletonode[1],] <- 
    apply(suppmat[cycletonode,], 2, max)
  # let all infectee-candidates not be infected by the other positions (removing support)
  suppmat[tail(cycletonode,-1),] <- -samplesize
  # remove mutual support for cycle members
  suppmat[cycletonode,cycletonode] <- -samplesize
  # adjust the support for infectors of the cycle: the maximum of supports among the cycle members
  suppmat[,cycletonode[1]] <- 
    apply(suppmat[,cycletonode], 1, max) - max(suppmat[,cycletonode])
  # let the other positions all take the index as infector (removing their role in the tree)
  suppmat[,tail(cycletonode,-1)] <- -samplesize
  suppmat[1,tail(cycletonode,-1)] <- 0
  
  # use the adjusted support matrix to get a proper transmission tree
  treeres <- .edmondsiterative(suppmat, samplesize, obs)
  
  # determine which hosts are infected by the cycle,
  # and the infector of the cycle
  incoming <- which(treeres == cycletonode[1])
  outgoing <- treeres[cycletonode[1]]
  
  # let the cycle infectees be infected by their most likely infector
  treeres[incoming] <- whichincoming[incoming]
  # resolve the cycle by choosing the best infectee for the selected infector
  cycleinfector[cycletonode == whichoutgoing[outgoing]] <- outgoing
  treeres[cycletonode] <- cycleinfector
  
  
  return(treeres)
}

### returns TRUE if ID is in a cycle, FALSE if not, based on vector of infectors
### called by:
# .edmondsiterative
### calls:
# .pathtocycle
.cycleYN <- function(parset, ID) {
  .pathtocycle(parset, ID)[1] == ID && ID != 1
}

### return path to the root from IDs[1], iteratively calling this function
### stop if a cycle is encountered
### called by:
# .edmondsiterative
# .cycleYN
# .pathtocycle
### calls:
# .pathtocycle
.pathtocycle <- function(parset, IDs) {
  if(IDs[1] == 1 | length(IDs) > length(unique(IDs))) {
    return(IDs)
  } else {
    return(.pathtocycle(parset, c(parset[IDs[1]],IDs)))
  }
}
