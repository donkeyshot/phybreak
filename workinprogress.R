setwd("S:/R/klinkend/Phylomodellen/ResultsJan16/Resultsplus")

res <- readRDS("res1__wh1")

setwd("S:/R/klinkend/Phylomodellen/ResultsJan16/Simulations")

load("sims_wh1")

library(phybreak)
simulatie <- sim.phybreak.gentime(20)
curstate <- make.phybreak.obkData(simulatie)
curstate <- burnin.phybreak(curstate, 500)
curstate <- sample.phybreak(curstate, 100, 10)

for(i in 1:2786) {
  for(j in (i+1):2787) {
    if(all(testres2[(i*20-19):(i*20)] == testres2[(j*20-19):(j*20)])) {
      print(c(i,j))
    }
  }
}


makephyloparset <- function(parentset) {
  obs <- (1 + length(parentset))/3
  res <- parentset
  while(max(res) >= 2*obs) {
    res[res >= 2*obs] <- parentset[res][res >= 2*obs]
  }
  return(res[1 : (2*obs-1)])
}


makephyloparsets <- function(nodeparentsets) {
  apply(nodeparentsets,
        MARGIN = 2, makephyloparset)
}

CCphyloscores.phybreak <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs


  return(.MLphylotree_MCC(
    makephyloparsets(phybreak.object$s$nodeparents[,
                                (1:samplesize) + chainlength - samplesize]),
    c(obs, samplesize)
  ))

}


CCphylotreeconstruct.phybreak <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs


  return(.CCphylotreeconstruct(
    makephyloparsets(phybreak.object$s$nodeparents[,
                                                   (1:samplesize) + chainlength - samplesize]),
    c(obs, samplesize)
  ))

}

.makephylo2 <- function(nodeparents) {
  ###topology
  Nhosts <- (1+length(nodeparents))/2
  indexc <- (1:length(nodeparents))[nodeparents == 0]
  edgestart <- nodeparents[nodeparents != 0]
  edgeend <- (1:length(nodeparents))[nodeparents != 0]
  edgelengths <- rep(1,length(nodeparents)-1)

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
    tip.label = 1:Nhosts
  )
  class(res) <- "phylo"
  res <- reorder(res)
  return(res)

}


CCtranstree.phybreak <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs

  return(.CCtranstree(
    phybreak.object$s$nodehosts[obs:(2*obs-1),
                                (1:samplesize) + chainlength - samplesize],
    phybreak.object$s$nodetimes[obs:(2*obs-1),
                                (1:samplesize) + chainlength - samplesize],
    c(obs, samplesize)
  ))

}

CCtransconstruct.phybreak <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs

  return(.CCtranstreeconstruct(
    phybreak.object$s$nodehosts[obs:(2*obs-1),
                                (1:samplesize) + chainlength - samplesize],
    c(obs, samplesize)
  ))

}



CCphylotree.phybreak <- function(phybreak.object, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs


  return(.CCphylotree(
    makephyloparsets(phybreak.object$s$nodeparents[,
                                                   (1:samplesize) + chainlength - samplesize]),
    c(obs, samplesize)
  ))

}


.postinfector <- function(inf.chain, ranknr = 1, support = FALSE) {
  chainlength <- length(inf.chain)
  ordinf <- rle(sort(inf.chain))
  if(length(ordinf$lengths) < ranknr) {
    if(support) return(c(0,0)) else return(0)
  }
  ord <- order(ordinf$lengths,decreasing = TRUE)
  if(support) {
    return(c(ordinf[[2]][ord[ranknr]],ordinf[[1]][ord[ranknr]]/chainlength))
  } else {
    return(ordinf[[2]][ord[ranknr]])
  }
}

pcmatrix <- matrix(0, nrow = 50, ncol = 1000)
parmatrix <- res[[2]]$s$nodehosts[50:99,]+1
parcounts <- apply(parmatrix,1,tabulate,nbins=51)
parcred <- function(hostID, samplenr, pcounts, pmat) {
  pcounts[pmat[hostID, samplenr], hostID]
}
result <- matrix(mapply(parcred, rep(1:50,1000),rep(1:1000,each=50),MoreArgs = list(pcounts = parcounts,pmat = parmatrix)),nrow=50)
parcreds <- colSums(log(result))
which(parcreds == max(parcreds))

MPC.tree <- function(phybreak.object, support = FALSE, infectornames = TRUE, samplesize = Inf) {
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  samplesize <- min(samplesize, chainlength)
  res <- matrix(with(phybreak.object,
                     apply(s$nodehosts[p$obs:(2*p$obs-1),(chainlength-samplesize+1):chainlength],
                           1,.postinfector,support = support)),nrow = 1 + support)
  if(infectornames) res[1,] <- c("index",phybreak.object$d$names)[1+res[1,]]
  if(!support) {
    return(as.data.frame(t(res),ncol=1,dimnames=list(phybreak.object$d$names,"infectors")))
  } else {
    rownames(res) <- c("infectors","support")
    colnames(res) <- phybreak.object$d$names
    return(as.data.frame(t(res)))
  }
}




posttrees <- c(make.multiPhylo(res[[3]])[500:1000],
               make.multiPhylo(res[[3]])[1:500])
posttrees <- c(posttrees[502:1001],posttrees[1:501])
truestate <- make.phybreak.obkData(simulaties[[3]],use.tree = TRUE)
posttrees[[1001]] <- make.phylo.phybreak(truestate)
mccdistmat <- rspr.matrix(posttrees[991:1001],"full")
mccdistmat
ccphydists <- CCphyloscores.phybreak(res[[3]])
which(ccphydists == max(ccphydists))
rspr(posttrees[[1001]],posttrees[[557]])

testpars <- res[[1]]$s$nodeparents[,1]
testpars[testpars >= 100] <- testpars[testpars][testpars >= 100]

saveRDS(simulaties, "sims_wh1RDS.rds")
simulaties2 <- readRDS("sims_wh1RDS.rds")

correctconstruct <- rep(0,10)
for(i in 1:10){
  correctconstruct[i] <- sum(MLtrans(res[[i]],"cc.construct") == simulaties[[i]]@individuals$infector)
}
correctMCC <- rep(0,10)
for(i in 1:10){
  correctMCC[i] <- sum(MLtrans(res[[i]],"mcc") == simulaties[[i]]@individuals$infector)
}
correctMLtreeinf <- rep(0,10)
for(i in 1:10) {
  correctMLtreeinf[i] <- sum(MLtrans(res[[i]],"edmunds") == simulaties[[i]]@individuals$infector)
}
correctMLinf <- rep(0,10)
for(i in 1:10) {
  correctMLinf[i] <- sum(MLtrans(res[[i]],"count") == simulaties[[i]]@individuals$infector)
}



whichtree <- which(CCscores.phybreak(res[[3]]) == max(CCscores.phybreak(res[[3]])))
sum(res[[3]]$s$nodehosts[50:99,whichtree[2]] == res[[3]]$s$nodehosts[50:99,whichtree[1]])
