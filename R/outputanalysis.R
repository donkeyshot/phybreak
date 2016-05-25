###obtain current state variables
get.infectors <- function(phybreak.object) {
  return(tail(phybreak.object$v$nodehosts,phybreak.object$p$obs))
}

get.parameters <- function(phybreak.object) {
  return(unlist(phybreak.object$p))
}

###MCMC analysis
#include 'coda'

make.coda <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  res <- with(phybreak.object,
              cbind(t(s$nodetimes[p$obs:(2*p$obs-1),(chainlength-samplesize+1):chainlength]),
                    t(s$nodehosts[p$obs:(2*p$obs-1),(chainlength-samplesize+1):chainlength])))
  parnames <- c(paste0("tinf.",phybreak.object$d$names),
                     paste0("parent.",phybreak.object$d$names))
  if(phybreak.object$h$est.mS) {
    res <- cbind(phybreak.object$s$mS[(chainlength-samplesize+1):chainlength], res)
    parnames <- c("mS",parnames)
  }
  if(phybreak.object$h$est.mG) {
    res <- cbind(phybreak.object$s$mG[(chainlength-samplesize+1):chainlength], res)
    parnames <- c("mG",parnames)
  }
  if(phybreak.object$h$est.wh) {
    if(phybreak.object$p$wh.model == 3) {
      res <- cbind(phybreak.object$s$r0[(chainlength-samplesize+1):chainlength], res)
      parnames <- c("r0",parnames)
    } else {
      res <- cbind(phybreak.object$s$slope[(chainlength-samplesize+1):chainlength], res)
      parnames <- c("slope",parnames)
    }
  }
  res <- cbind(phybreak.object$s$mu[(chainlength-samplesize+1):chainlength], res,
               phybreak.object$s$logLik[(chainlength-samplesize+1):chainlength])
  parnames <- c("mu",parnames,"logLik")
  colnames(res) <- parnames
  return(mcmc(res))
}

###phangorn functionality
#include 'phangorn'
.makephylo <- function(nodetimes, nodeparents, nodenames) {
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
.makephylosimmap <- function(nodetimes, nodeparents, nodehosts, nodenames) {
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

.makephyDat <- function(SNPs, SNPfreqs, nodenames) {
  Nhosts <- nrow(SNPs)

  ###sequences
  dnadata <- t(matrix(rep(t(SNPs),rep(SNPfreqs,Nhosts)),ncol=Nhosts))
  rownames(dnadata) <- nodenames

  return(phyDat(dnadata))
}


make.phylo.phybreak <- function(phybreak.object, samplenr = 0, simmap = FALSE) {
  if(samplenr > length(phybreak.object$s$mu)) {
    warning("requested 'samplenr' not available; current state used")
    samplenr <- 0
  }

  if(samplenr == 0) {
    nodetimes <- phybreak.object$v$nodetimes
    nodeparents <- phybreak.object$v$nodeparents
    nodehosts <- phybreak.object$v$nodehosts
  } else {
    nodetimes <- with(phybreak.object,
                      c(v$nodetimes[1:p$obs],
                        s$nodetimes[,samplenr]))
    nodeparents <- phybreak.object$s$nodeparents[,samplenr]
    nodehosts <- with(phybreak.object,
                      c(v$nodehosts[1:p$obs],
                        s$nodehosts[,samplenr]))
  }
  nodenames <- phybreak.object$d$names

  if(simmap) {
    return(.makephylosimmap(nodetimes, nodeparents, nodehosts, nodenames))
  } else {
    return(.makephylo(nodetimes, nodeparents, nodenames))
  }
}

treeplot.phybreak <- function(phybreak.object, plot.mpc = FALSE, samplenr = 0, ...) {
  if(plot.mpc) {
    plotSimmap(MLtrans(phybreak.object, "mpc", phylo.class = TRUE),
               colors = setNames(nm = c("black","red","blue")), ...)
  } else {
    plotSimmap(make.phylo.phybreak(phybreak.object, samplenr, TRUE),
               colors = setNames(nm = c("black","red","blue")), ...)
  }
}

make.multiPhylo <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  treelist <- vector('list', samplesize)
  for(i in 1:samplesize){
    treelist[[i]] <-
      make.phylo.phybreak(phybreak.object,
                 i + chainlength - samplesize)
  }
  class(treelist) <- "multiPhylo"
  return(treelist)
}

make.phyDat <- function(phybreak.object) {
  nodenames <- phybreak.object$d$names
  return(with(phybreak.object$d,
              .makephyDat(SNP, SNPfr, nodenames)))
}


###infector analysis
#return 'ranknr'th most likely infector from a posterior chain
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



#return ordered most likely infectors, with posterior support
.inflist <- function(inf.chain, nhosts) {
  sapply(1:nhosts, .postinfector, inf.chain=inf.chain, support = TRUE)
}

#return ordered most likely infectors for multiple hosts, with posterior support
.infarray <- function(inf.matrix) {
  res <- t(apply(inf.matrix,MARGIN = 1, FUN = .inflist, nhosts = nrow(inf.matrix)))
  dim(res) <- c(nrow(inf.matrix), 2, nrow(inf.matrix))
  return(res)
}

#return distribution of infectors for hostID, including 0-counts, starting with root 0
.infectordist <- function(hostID, posterior, obs) {
  res <- tabulate(posterior[hostID,] + 1,nbins = obs+1)
  return(res - max(res))
}

#return support for infector of hostID
.infectorsupport <- function(hostID, infectors, posterior) {
  sum(posterior[hostID,]==infectors[hostID])
}

#return most frequent number in vector
.MLinfector <- function(infdistcolumn) {
  which(infdistcolumn == max(infdistcolumn))[1]
}

#return path to the root from IDs[1], iteratively calling this function
#stop if a cycle is encountered
.pathtoroot <- function(parset, IDs) {
  if(IDs[1] == 1 | length(IDs) > length(unique(IDs))) {
    return(IDs)
  } else {
    return(.pathtoroot(parset, c(parset[IDs[1]],IDs)))
  }
}

#return if ID is in a cycle
.cycleYN <- function(parset, ID) {
  .pathtoroot(parset, ID)[1] == ID && ID != 1
}

#iterative call to edmond's algorithm, as long as there are cycles
.edmondsiterative <- function(infdistmat, samplesize, obs) {
  parentset <- apply(infdistmat, 2, .MLinfector)
  cycleIDs <- which(sapply(1:obs, 
                           .cycleYN, 
                           parset = parentset))
  if(length(cycleIDs) == 0) {
    return(parentset)
  }
  cycletonode <- unique(.pathtoroot(parentset, cycleIDs[1]))
  cycleinfector <- parentset[cycletonode]
  
  whichincoming <- cycletonode[apply(infdistmat[cycletonode,],2,which.max)]
  whichoutgoing <- cycletonode[apply(infdistmat[,cycletonode],1,which.max)]
  
  infdistmat[cycletonode[1],] <- 
    apply(infdistmat[cycletonode,], 2, max)
  
  infdistmat[tail(cycletonode,-1),] <- -samplesize
  
  infdistmat[cycletonode,cycletonode] <- -samplesize
  
  infdistmat[,cycletonode[1]] <- 
    apply(infdistmat[,cycletonode], 1, max) - max(infdistmat[,cycletonode])
  
  infdistmat[,tail(cycletonode,-1)] <- -samplesize
  infdistmat[1,tail(cycletonode,-1)] <- 0
  
  treeres <- .edmondsiterative(infdistmat, samplesize, obs)
  
  
  incoming <- which(treeres == cycletonode[1])
  outgoing <- treeres[cycletonode[1]]
  
  treeres[incoming] <- whichincoming[incoming]
  cycleinfector[cycletonode == whichoutgoing[outgoing]] <- outgoing
  treeres[cycletonode] <- cycleinfector
  
  
  return(treeres)
  
}


#worker function to obtain a transmission tree based on most likely infector, possibly containing cycles and multiple roots,
#plus means and standard deviations of infection times
.transtreecount <- function(phybreak.object, samplesize) {
  chainlength <- length(phybreak.object$s$mu)
  obsize <- phybreak.object$p$obs
  samplerange <- (chainlength-samplesize+1):chainlength
  
  res <- t(matrix(with(phybreak.object,
                       apply(s$nodehosts[obsize:(2*obsize-1),samplerange],
                             1,.postinfector,support = TRUE)),nrow = 2))
  
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
  return(res)
}


#worker function to obtain a transmission tree based on most likely infectors, using Edmonds's algorithm 
#to remove cycles, applied with multiple root candidates, plus means and standard deviations of infection times
.transtreeedmonds <- function(phybreak.object, samplesize) {
  chainlength <- length(phybreak.object$s$mu)
  obsize <- phybreak.object$p$obs
  samplerange <- (chainlength-samplesize+1):chainlength
  
  infdistmatrix <- cbind(
    c(0, rep(-1,obsize)),
    sapply(1:obsize,.infectordist,
           posterior=phybreak.object$s$nodehosts[obsize:(2*obsize-1), samplerange],
           obs=obsize))
  candidateroots <- as.numeric(
    names( sort( table( 
      apply(phybreak.object$s$nodehosts[obsize:(2*obsize-1),]==0,2,which)
      ),decreasing = TRUE)))
  rootsupports <- as.numeric(
    sort( table(
      apply(phybreak.object$s$nodehosts[obsize:(2*obsize-1),]==0,2,which)
      ),decreasing = TRUE))
  bestYN <- FALSE
  nextcandidate <- 1
  alltrees <- c()
  allsupports <- c()
  while(!bestYN) {
    infdistmat <- infdistmatrix
    infdistmat[1,] <- -samplesize
    infdistmat[1,candidateroots[nextcandidate] + 1] <- 0
    
    thistree <- .edmondsiterative(infdistmat, samplesize, obsize)[-1] - 1
    alltrees <- c(alltrees, thistree)
    thissupport <- sapply(1:obsize, .infectorsupport, 
                          infectors = thistree,
                          posterior = phybreak.object$s$nodehosts[obsize:(2*obsize-1), 
                                                                  samplerange])              
    allsupports <- c(allsupports, thissupport)
    
    bestYN <- TRUE                 
    for(i in tail(candidateroots,-nextcandidate)) {
      if(thissupport[i] < rootsupports[candidateroots==i]) {
        bestYN <- FALSE
      }
    }
    nextcandidate <- nextcandidate + 1
  }
  
  dim(alltrees) <- c(obsize,length(alltrees)/obsize)
  dim(allsupports) <- dim(alltrees)
  
  besttree <- alltrees[,which(colSums(allsupports)==max(colSums(allsupports)))]

  treesupport <- sapply(1:obsize, .infectorsupport, 
                        infectors = besttree,
                        posterior = phybreak.object$s$nodehosts[obsize:(2*obsize-1), 
                                                                samplerange])
  res <- matrix(c(besttree,
                  treesupport),
                ncol=2)
  
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
  
  return(res)
  
}


#return 'maximum parent credibility' tree with posterior support and means and standard deviations of 
#infection times, and infection times of the mpc tree itself
.mpcinfector <- function(phybreak.object, samplesize, phylo.class = FALSE) {
  chainlength <- length(phybreak.object$s$mu)
  obsize <- phybreak.object$p$obs
  samplerange <- (chainlength-samplesize+1):chainlength
  
  parcounts <- apply(
    phybreak.object$s$nodehosts[obsize:(2*obsize-1), samplerange] + 1,
    1, tabulate, nbins=obsize + 1)
  parcred <- function(hostID, samplenr) {
    parcounts[
      phybreak.object$s$nodehosts[hostID + obsize - 1, samplenr] + 1,
      hostID]
  }
  parcreds <- colSums(
    log(
      matrix(
        mapply(
          parcred,
          rep(1:obsize,samplesize),
          rep(samplerange,each=obsize)
        ),nrow = obsize
      )
    )
  )
  bestpars <- which(parcreds == max(parcreds))[1]
  if(phylo.class) return(bestpars + samplerange[1] - 1)
  res <- matrix(
    c(
      phybreak.object$s$nodehosts[obsize:(2*obsize-1),
                                  bestpars + samplerange[1] - 1],
      mapply(parcred,
             1:obsize,
             rep(bestpars + samplerange[1] - 1,obsize)
      )
    ),ncol = 2
  )
  
  
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
  res <- cbind(res, phybreak.object$s$nodetimes[obsize:(2*obsize-1),
                                                bestpars + samplerange[1] - 1])
  
  
  return(res)
}

#### WRAPPER to get posterior mean transmission trees
MLtrans <- function(phybreak.object,
                    method = c("count", "edmonds", "mpc", "mtcc", "cc.construct"),
                    samplesize = Inf, infector.name = TRUE, support.ntrees = FALSE,
                    phylo.class = FALSE) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs
  
  res <- c()
  
  if(method[1] == "count") {
    res <- .transtreecount(phybreak.object, samplesize)
  }
  if(method[1] == "edmonds") {
    res <- .transtreeedmonds(phybreak.object, samplesize)
  }
  if(method[1] == "mpc") {
    res <- .mpcinfector(phybreak.object, samplesize, phylo.class)
    if(phylo.class) return(make.phylo.phybreak(phybreak.object, res, TRUE))
  }
  if(method[1] == "mtcc") {
    res <- matrix(.CCtranstree2(
      phybreak.object$s$nodehosts[obs:(2*obs-1),
                                  (1:samplesize) + chainlength - samplesize],
      phybreak.object$s$nodetimes[obs:(2*obs-1),
                                  (1:samplesize) + chainlength - samplesize],
      c(obs, samplesize)
    ), ncol = 5)
  }
  if(method[1] == "cc.construct") {
    res <- matrix(.CCtranstreeconstruct(
      phybreak.object$s$nodehosts[obs:(2*obs-1),
                                  (1:samplesize) + chainlength - samplesize],
      phybreak.object$s$nodetimes[obs:(2*obs-1),
                                  (1:samplesize) + chainlength - samplesize],
      c(obs, samplesize)
    ), ncol = 4)
  }
  
  if(length(res) == 0) {
    stop("incorrect method provided, choose \"count\", \"edmonds\",
\"mpc\", \"mtcc\", or \"cc.construct\"")
  }
  
  if(infector.name) {
    infectors.out <- matrix(c("index",phybreak.object$d$names)[1+res[,1]],
                            ncol = 1,
                            dimnames=list(phybreak.object$d$names,"infector"))
  } else {
    infectors.out <- matrix(res[,1],ncol = 1,
                            dimnames=list(phybreak.object$d$names,"infector"))
  }
  if(support.ntrees) {
    support.out <- res[,2]
  } else {
    support.out <- res[,2] / samplesize
  }
  if(method[1] == "mpc" || method[1] == "mtcc") {
    return(
      data.frame(
        infectors = infectors.out,
        support = support.out,
        inftime.mean = res[,3],
        inftime.sd = res[,4],
        inftime.mc = res[,5]
      )
    )
  } else {
    return(
      data.frame(
        infectors = infectors.out,
        support = support.out,
        inftime.mean = res[,3],
        inftime.sd = res[,4]
      )
    )
  }
  
  
}

##sets of possible infectors per host, based on 'count' method
infectorsets <- function(phybreak.object, percentile = 0.95, 
                         minsupport = 0, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs
  
  inffreqs <- .infarray(phybreak.object$s$nodehosts[obs:(2*obs-1),])
  
  includemat <- inffreqs[,2,] > minsupport * samplesize
  cumfreqincl <- t(rbind(rep(TRUE, obs),
                         apply(inffreqs[,2,],1,cumsum) < percentile * samplesize)[-(obs+1),])
  includemat <- includemat & cumfreqincl
  
  res <- list()
  for(i in 1:obs) {
    res[[i]] <- inffreqs[i,1,][includemat[i,]]
  }
  names(res) <- phybreak.object$d$names
  return(res)
}


#count equal elements in two vectors (such as infectors)
.equal.infectors <- function(parset1, parset2) {
  sum(parset1 == parset2)
}

#return distance matrix (class dist) between trees, based on infector count
treedists.phybreak <- function(phybreak.object, samplesize = 100, thin = 10) {
  obsize <- phybreak.object$p$obs
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize*thin > chainlength & samplesize < Inf) {
    warning("desired 'samplesize*thin' larger than number of available samples")
  }
  samplesize <- min(samplesize, floor(chainlength/thin))

  res <- matrix(NA, nrow = samplesize, ncol = samplesize)
  for(i in 1:samplesize) {
    for(j in 1:i) {
      res[i,j] <- .equal.infectors(phybreak.object$s$nodehosts[obsize:(2*obsize-1),(chainlength - thin*samplesize) + thin*i],
                                   phybreak.object$s$nodehosts[obsize:(2*obsize-1),(chainlength - thin*samplesize) + thin*j])
      res[j,i] <- res[i,j]
    }
  }
  return(as.dist(obsize-res))
}


.makephyloparset <- function(parentset) {
  obs <- (1 + length(parentset))/3
  res <- parentset
  while(max(res) >= 2*obs) {
    res[res >= 2*obs] <- parentset[res][res >= 2*obs]
  }
  return(res[1 : (2*obs-1)])
}

.makephyloparsets <- function(nodeparentsets) {
  apply(nodeparentsets,
        MARGIN = 2, .makephyloparset)
}

.makephylo2 <- function(nodeparents, nodetimes, nodenames) {
  ###topology
  Nhosts <- (1+length(nodeparents))/2
  indexc <- (1:length(nodeparents))[nodeparents == 0]
  edgestart <- nodeparents[nodeparents != 0]
  edgeend <- (1:length(nodeparents))[nodeparents != 0]
  edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]
  
  if(indexc != Nhosts + 1) {
    edgestart[edgestart == indexc] <- 0
    edgeend[edgeend == indexc] <- 0
    edgestart[edgestart == Nhosts + 1] <- indexc
    edgeend[edgeend == Nhosts + 1] <- indexc
    edgestart[edgestart == 0] <- Nhosts + 1
    edgeend[edgeend == 0] <- Nhosts + 1
  }
  
  edges <- matrix(c(edgestart,edgeend),ncol=2)
  
  
  
  res <- list(
    edge = edges,
    edge.length = edgelengths,
    Nnode = Nhosts - 1,
    tip.label = nodenames
  )
  class(res) <- "phylo"
  res <- reorder(res)
  res <- ladderize(res)
  return(res)
  
}


MLphylo <- function(phybreak.object,
                    method = c("mcc", "cc.construct"),
                    phylo.class = TRUE, mc.times = TRUE, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs
  
  res <- c()
  
  if(method[1] == "mcc") {
    res <- matrix(.CCphylotree(
      .makephyloparsets(phybreak.object$s$nodeparents[,
                                                     (1:samplesize) + chainlength - samplesize]),
      phybreak.object$s$nodetimes[1:(obs - 1),
                                    (1:samplesize) + chainlength - samplesize],
      c(obs, samplesize)
    ), ncol = 5)
    res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs]
    res[1:obs, 5] <- phybreak.object$v$nodetimes[1:obs]
  }
  if(method[1] == "cc.construct") {
    res <- matrix(.CCphylotreeconstruct(
      .makephyloparsets(phybreak.object$s$nodeparents[,
                                                     (1:samplesize) + chainlength - samplesize]),
      phybreak.object$s$nodetimes[1:(obs - 1),
                                  (1:samplesize) + chainlength - samplesize],
      c(obs, samplesize)
    ), ncol = 4)
    res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs]
  }
  
  if(length(res) == 0) {
    stop("incorrect method provided, choose \"mcc\" or \"cc.construct\"")
  }

  parents.out <- matrix(res[,1],ncol = 1,
                          dimnames=list(1:(2*obs-1),"parent"))
  if(method[1] == "mcc") {
    if(phylo.class) {
      if(mc.times) {
        return(.makephylo2(res[,1],res[,5],phybreak.object$d$names))
      } else {
        return(.makephylo2(res[,1],res[,3],phybreak.object$d$names))
      }
    } else {
      return(
        data.frame(
          parents = parents.out,
          support = res[,2],
          nodetime.mean = res[,3],
          nodetime.sd = res[,4],
          nodetime.mc = res[,5]
        )
      )
    }
  } else {
    if(phylo.class) {
      return(.makephylo2(res[,1],res[,3],phybreak.object$d$names))
    } else
    return(
      data.frame(
        parents = parents.out,
        support = res[,2],
        nodetime.mean = res[,3],
        nodetime.sd = res[,4]
      )
    )
  }
  
  
}


comparephybreak.infectors <- function(phybreak.object, sim.object) {
  res <- rep(NA,length(phybreak.object$s$logLik))
  for(i in 1:length(res)) {
    res[i] <- .equal.infectors(tail(phybreak.object$s$nodehosts[,i],phybreak.object$p$obs),
                              sim.object@individuals$infector)
  }
  return(res)
}

# eqinfch <- equal.infectors.chain(curstate,simulatie)
#
# equal.infectors(curstate,simulatie)
#
#
#
# dist.two.hosts <- function(parset, h1, h2) {
#   ptr1 <- pathtoroot(parset, h1)
#   ptr2 <- pathtoroot(parset, h2)
#   length(union(ptr1,ptr2)) - length(intersect(ptr1,ptr2))
# }
#
# dist.all.parents <- function(parset1, parset2) {
#   distmat <- matrix(0,ncol=2, nrow=length(parset1))
#   paths1 <- matrix(0,ncol=length(parset1), nrow=length(parset1))
#   paths2 <- matrix(0,ncol=length(parset1), nrow=length(parset1))
#
#   for(i in 1:length(parset1)) {
#     path1 <- pathtoroot(parset1, i)
#     path2 <- pathtoroot(parset2, i)
#     paths1[i,1:length(path1)] <- path1
#     paths2[i,1:length(path2)] <- path2
#   }
#
#   for(i in 1:length(parset1)) {
#     distmat[i,1] <- length(union(paths1[parset1[i],],paths1[parset2[i],])) -
#       length(intersect(paths1[parset1[i],],paths1[parset2[i],]))
#     distmat[i,2] <- length(union(paths2[parset1[i],],paths2[parset2[i],])) -
#       length(intersect(paths2[parset1[i],],paths2[parset2[i],]))
#   }
#   return(distmat)
# }
#
# dist.post.to.sim <- function(parset.post, parset.sim) {
#   distvec <- rep(0, length(parset.sim))
#   paths <- matrix(0,ncol=length(parset.sim), nrow=length(parset.sim))
#
#   for(i in 1:length(parset.sim)) {
#     path <- pathtoroot(parset.sim, i)
#     paths[i,1:length(path)] <- path
#   }
#
#   for(i in 1:length(parset.sim)) {
#     distvec[i] <- length(union(paths[parset.post[i],],paths[parset.sim[i],])) -
#       length(intersect(paths[parset.post[i],],paths[parset.sim[i],]))
#   }
#   return(distvec)
# }
#
#
# MLinfector <- function(inf.chain) {
#   as.numeric(names(sort(table(inf.chain),decreasing=TRUE))[1])
# }
#
# MLinfectors <- function(phybreak.object) {
#   with(phybreak.object,apply(s$nodehosts[p$obs:(2*p$obs-1),],1,MLinfector))
# }
#
# ps1 <- MLinfectors(curstate)
# ps2 <- simulatie$trueoutbreak$nodehosts[400:599]
#
#
# make.transmission.clades <- function(parset) {
#   res <- matrix(FALSE,nrow=length(parset),ncol=length(parset))
#   for(i in 1:length(parset)) {
#     res[cbind(pathtoroot(parset,i),i)] <- TRUE
#   }
#   return(res)
# }
#
# trclades <- make.transmission.clades(ps2)
# plot(curstate$s$nodehosts[200,])
# postclades <- array(FALSE,dim=c(1000,200,200))
# for(i in 1:1000) {
#   postclades[i,,] <- make.transmission.clades(curstate$s$nodehosts[200:399,i])
# }
# cladescores <- array(0,dim=c(1000,200))
# for(i in 1:1000) {
#   for(j in 1:200) {
#     if(cladescores[i,j] == 0) {
#       cladeidenticalQ <- array(TRUE,dim=c(1000,200))
#       for(k in 1:200) {
#         cladeidenticalQ[postclades[i,j,k] != postclades[,,k]] <- FALSE
#       }
#       cladescores[cladeidenticalQ] <- sum(cladeidenticalQ)
#     }
#   }
# }
# treescores <- rowSums(cladescores)
# which(treescores == max(treescores))
# ps3 <- curstate$s$nodehosts[200:399,770]
# sum(ps2==ps3)
# idscores <- rep(0,1000)
# for(i in 1:1000) {
#   idscores[i] <- sum(ps2 == curstate$s$nodehosts[200:399,i])
# }

