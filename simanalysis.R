setwd("S:/R/klinkend/Phylomodellen/ResultsJan16/Resultsplus")

res <- readRDS("res2__wh2")

setwd("S:/R/klinkend/Phylomodellen/ResultsJan16/Simulations")

load("sims_wh2")

library(phybreak)
setwd("S:/R/klinkend/Phylomodellen/ResultsJan16/Resultsplus")

resfiles <- list.files()
simfiles <- list.files("../Simulations")

#define required output evaluation functions
{
infsetcoverage <- function(phybreak.object, sim.obkData.object, percentile = 0.95, 
                           minsupport = 0, samplesize = Inf) {
  infsets <- infectorsets(phybreak.object, percentile, minsupport, samplesize)
  res <- rep(FALSE, phybreak.object$p$obs)
  for(i in 1:length(res)) {
    res[i] <- sim.obkData.object@individuals$infector[i] %in% infsets[[i]]
  }
  return(sum(res))
}

makephyloparset <- function(parentset) {
  obs <- (1 + length(parentset))/3
  res <- parentset
  while(max(res) >= 2*obs) {
    res[res >= 2*obs] <- parentset[res][res >= 2*obs]
  }
  return(res[1 : (2*obs-1)])
}


pathtoroot <- function(parset, IDs) {
  if(IDs[1] == 0 | length(IDs) > length(unique(IDs))) {
    return(IDs)
  } else {
    return(pathtoroot(parset, c(parset[IDs[1]],IDs)))
  }
}

#
dist.two.hosts <- function(parset, h1, h2) {
  ptr1 <- pathtoroot(parset, h1)
  ptr2 <- pathtoroot(parset, h2)
  length(union(ptr1,ptr2)) - length(intersect(ptr1,ptr2))
}
#which posterior clades are correct?
equal.clades.trans <- function(parset.post, parset.sim) {
  parset.post <- unlist(parset.post)
  pathspost <- matrix(FALSE,ncol=length(parset.post), nrow=length(parset.post))
  pathssim <- matrix(FALSE,ncol=length(parset.post), nrow=length(parset.post))
  
  for(i in 1:length(parset.post)) {
    path1 <- pathtoroot(parset.post, i)[-1]
    path2 <- pathtoroot(parset.sim, i)[-1]
    pathspost[i,path1] <- TRUE
    pathssim[i,path2] <- TRUE
  }
  res <- rep(FALSE, length(parset.post))
  
  for(i in 1:length(parset.post)) {
    res[i] <- any(apply(pathspost[,i] == pathssim, 2, all))
  }
  return(res)
}

equal.clades.phylo <- function(parset.post, parset.sim) {
  parset.post <- unlist(parset.post)
  obs <- (length(parset.post) + 1)/2
  pathspost <- matrix(FALSE,ncol=length(parset.post), nrow=length(parset.post))
  pathssim <- matrix(FALSE,ncol=length(parset.post), nrow=length(parset.post))
  
  for(i in 1:length(parset.post)) {
    path1 <- pathtoroot(parset.post, i)[-1]
    path2 <- pathtoroot(parset.sim, i)[-1]
    pathspost[i,path1] <- TRUE
    pathssim[i,path2] <- TRUE
  }
  
  res <- rep(FALSE, length(parset.post))
  
  for(i in obs:(2*obs-1)) {
    res[i] <- any(apply(pathspost[1:obs,i] == pathssim[1:obs,], 2, all))
  }
  return(res[obs:(2*obs-1)])
}


#how far in the posterior tree is the true infector?
dist.post.to.sim <- function(parset.post, parset.sim) {
  parset.post <- unlist(parset.post)
  distvec <- rep(0, length(parset.post))
  paths <- matrix(0,ncol=length(parset.post), nrow=length(parset.post))
  
  for(i in 1:length(parset.post)) {
    path <- pathtoroot(parset.post, i)
    paths[i,1:length(path)] <- path
  }
  
  for(i in 1:length(parset.post)) {
    distvec[i] <- length(union(paths[parset.sim[i],],paths[parset.post[i],])) -
      length(intersect(paths[parset.sim[i],],paths[parset.post[i],]))
  }
  return(distvec)
}
}
#define function 'rspr(tree1, tree2)'
source("S:/R/klinkend/Phylomodellen/rsprCO/R/rspR_rev.r")

nextsims <- 8
load(paste0("../Simulations/sims_",
            substring(resfiles[nextsims], 7)))
res <- readRDS(resfiles[nextsims])


alloutputs <- array(rep(NA,12*74*208),dim=c(208,74,12), dimnames = list(
  resfiles,c("xburnin", "samplesize", "tr.all.tim", "tr.count.cov50", "tr.count.cov80",
    "tr.count.ior.0", "tr.count.ior.50", "tr.count.ior.80", "tr.count.iorset",
    "tr.count.tim.0", "tr.count.tim.50", "tr.count.tim.80",
    "tr.edmunds.cov50", "tr.edmunds.cov80",
    "tr.edmunds.ior.0", "tr.edmunds.ior.50", "tr.edmunds.ior.80",
    "tr.edmunds.tim.0", "tr.edmunds.tim.50", "tr.edmunds.tim.80",
    "tr.edmunds.cla.0", "tr.edmunds.cla.50", "tr.edmunds.cla.80",
    "tr.edmunds.2dg.0", "tr.edmunds.2dg.50", "tr.edmunds.2dg.80",
    "tr.mpc.cov50", "tr.mpc.cov80",
    "tr.mpc.ior.0", "tr.mpc.ior.50", "tr.mpc.ior.80",
    "tr.mpc.tim.0", "tr.mpc.tim.50", "tr.mpc.tim.80",
    "tr.mpc.cla.0", "tr.mpc.cla.50", "tr.mpc.cla.80",
    "tr.mpc.2dg.0", "tr.mpc.2dg.50", "tr.mpc.2dg.80",
    "tr.mtcc.cov50", "tr.mtcc.cov80",
    "tr.mtcc.ior.0", "tr.mtcc.ior.50", "tr.mtcc.ior.80",
    "tr.mtcc.tim.0", "tr.mtcc.tim.50", "tr.mtcc.tim.80",
    "tr.mtcc.cla.0", "tr.mtcc.cla.50", "tr.mtcc.cla.80",
    "tr.mtcc.2dg.0", "tr.mtcc.2dg.50", "tr.mtcc.2dg.80",
    "tr.construct.cov50", "tr.construct.cov80",
    "tr.construct.ior.0", "tr.construct.ior.50", "tr.construct.ior.80",
    "tr.construct.tim.0", "tr.construct.tim.50", "tr.construct.tim.80",
    "tr.construct.cla.0", "tr.construct.cla.50", "tr.construct.cla.80",
    "tr.construct.2dg.0", "tr.construct.2dg.50", "tr.construct.2dg.80",
    "ph.mcc.cov50", "ph.mcc.cov80",
    "ph.mcc.cla.0", "ph.mcc.cla.50", "ph.mcc.cla.80",
    "ph.mcc.rspr"#,
#     "ph.construct.cov50", "ph.construct.cov80",
#     "ph.construct.cla.0", "ph.construct.cla.50", "ph.construct.cla.80",
#     "ph.construct.rspr"
  ),
  c("sim1","sim2","sim3","sim4","sim5","sim6",
    "sim7","sim8","sim9","sim10","total","mean")
))

for(nextsims in 64:208) {
  load(paste0("../Simulations/sims_",
              substring(resfiles[nextsims], 7)))
  res <- readRDS(resfiles[nextsims])
  
  for(i in 1:10) {
    print(i)
    obs <- res[[i]]$p$obs
    parsimonies <- sapply(make.multiPhylo(res[[i]]), 
                          parsimony, 
                          data = make.phyDat(res[[i]]))
    
    ##additional burnin steps
    xburnin <- 100*ceiling((which(parsimonies == min(parsimonies))[1]-1)/100)
    alloutputs[nextsims,"xburnin",i] <- xburnin
    
    ##required sample size
    samplesize <- 1000
    codares <- make.coda(res[[i]])
    if(length(parsimonies) > 1000 && samplesize + xburnin < length(parsimonies)) {
      effsize <- min(effectiveSize(codares[(xburnin) + 1:samplesize,c("mu","logLik")]))
      while(effsize < 200 && samplesize + xburnin < length(parsimonies)) {
        samplesize <- samplesize + 100
        effsize <- min(effectiveSize(codares[(xburnin) + 1:samplesize,c("mu","logLik")]))
      }
    }
    alloutputs[nextsims,"samplesize",i] <- samplesize
    
    ##infection time coverage by posterior percentiles
    codares <- make.coda(res[[i]],samplesize)
    truetimes <- as.numeric(simulaties[[i]]@individuals$date - min(simulaties[[i]]@dna@meta$date))
    tr.all.tim <- sum((summary(codares[,1+(1:obs)],quantiles = c(.025,.975))[[2]][,1]<truetimes &
                         summary(codares[,1+(1:obs)],quantiles = c(.025,.975))[[2]][,2]>truetimes))
    alloutputs[nextsims,"tr.all.tim",i] <- tr.all.tim
    
    ##infector and infection time coverage for "count"
    truetree <- simulaties[[i]]@individuals$infector
    posttree <- MLtrans(res[[i]], samplesize = samplesize)
    posttree[is.nan(posttree[,"inftime.sd"]),"inftime.sd"] <- 
      posttree[is.nan(posttree[,"inftime.sd"]),"inftime.mean"]
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.count.cov50",i] <- sum(postcov50)
    alloutputs[nextsims,"tr.count.cov80",i] <- sum(postcov80)
    alloutputs[nextsims,"tr.count.ior.0",i] <- sum((posttree["infector"] == truetree))
    alloutputs[nextsims,"tr.count.ior.50",i] <- sum((posttree["infector"] == truetree) & postcov50)
    alloutputs[nextsims,"tr.count.ior.80",i] <- sum((posttree["infector"] == truetree) & postcov80)
    alloutputs[nextsims,"tr.count.iorset",i] <- infsetcoverage(res[[i]], simulaties[[i]], samplesize = samplesize)
    alloutputs[nextsims,"tr.count.tim.0",i] <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
    alloutputs[nextsims,"tr.count.tim.50",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
    alloutputs[nextsims,"tr.count.tim.80",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
    
    ##infector and infection time coverage for "edmunds"
    posttree <- MLtrans(res[[i]], method = "edmunds", samplesize = samplesize)
    posttree[is.nan(posttree[,"inftime.sd"]),"inftime.sd"] <- 
      posttree[is.nan(posttree[,"inftime.sd"]),"inftime.mean"]
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.edmunds.cov50",i] <- sum(postcov50)
    alloutputs[nextsims,"tr.edmunds.cov80",i] <- sum(postcov80)
    alloutputs[nextsims,"tr.edmunds.ior.0",i] <- sum((posttree["infector"] == truetree))
    alloutputs[nextsims,"tr.edmunds.ior.50",i] <- sum((posttree["infector"] == truetree) & postcov50)
    alloutputs[nextsims,"tr.edmunds.ior.80",i] <- sum((posttree["infector"] == truetree) & postcov80)
    alloutputs[nextsims,"tr.edmunds.tim.0",i] <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
    alloutputs[nextsims,"tr.edmunds.tim.50",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
    alloutputs[nextsims,"tr.edmunds.tim.80",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
    alloutputs[nextsims,"tr.edmunds.2dg.0",i] <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
    alloutputs[nextsims,"tr.edmunds.2dg.50",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
    alloutputs[nextsims,"tr.edmunds.2dg.80",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
    alloutputs[nextsims,"tr.edmunds.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.edmunds.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.edmunds.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##infector and infection time coverage for "mpc"
    posttree <- MLtrans(res[[i]], method = "mpc", samplesize = samplesize)
    posttree[is.nan(posttree[,"inftime.sd"]),"inftime.sd"] <- 
      posttree[is.nan(posttree[,"inftime.sd"]),"inftime.mean"]
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.mpc.cov50",i] <- sum(postcov50)
    alloutputs[nextsims,"tr.mpc.cov80",i] <- sum(postcov80)
    alloutputs[nextsims,"tr.mpc.ior.0",i] <- sum((posttree["infector"] == truetree))
    alloutputs[nextsims,"tr.mpc.ior.50",i] <- sum((posttree["infector"] == truetree) & postcov50)
    alloutputs[nextsims,"tr.mpc.ior.80",i] <- sum((posttree["infector"] == truetree) & postcov80)
    alloutputs[nextsims,"tr.mpc.tim.0",i] <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
    alloutputs[nextsims,"tr.mpc.tim.50",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
    alloutputs[nextsims,"tr.mpc.tim.80",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
    alloutputs[nextsims,"tr.mpc.2dg.0",i] <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
    alloutputs[nextsims,"tr.mpc.2dg.50",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
    alloutputs[nextsims,"tr.mpc.2dg.80",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
    alloutputs[nextsims,"tr.mpc.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.mpc.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.mpc.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##infector and infection time coverage for "mtcc"
    posttree <- MLtrans(res[[i]], method = "mtcc", samplesize = samplesize)
    posttree[is.nan(posttree[,"inftime.sd"]),"inftime.sd"] <- 
      posttree[is.nan(posttree[,"inftime.sd"]),"inftime.mean"]
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.mtcc.cov50",i] <- sum(postcov50)
    alloutputs[nextsims,"tr.mtcc.cov80",i] <- sum(postcov80)
    alloutputs[nextsims,"tr.mtcc.ior.0",i] <- sum((posttree["infector"] == truetree))
    alloutputs[nextsims,"tr.mtcc.ior.50",i] <- sum((posttree["infector"] == truetree) & postcov50)
    alloutputs[nextsims,"tr.mtcc.ior.80",i] <- sum((posttree["infector"] == truetree) & postcov80)
    alloutputs[nextsims,"tr.mtcc.tim.0",i] <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
    alloutputs[nextsims,"tr.mtcc.tim.50",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
    alloutputs[nextsims,"tr.mtcc.tim.80",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
    alloutputs[nextsims,"tr.mtcc.2dg.0",i] <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
    alloutputs[nextsims,"tr.mtcc.2dg.50",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
    alloutputs[nextsims,"tr.mtcc.2dg.80",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
    alloutputs[nextsims,"tr.mtcc.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.mtcc.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.mtcc.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##infector and infection time coverage for "cc.construct"
    posttree <- MLtrans(res[[i]], method = "cc.construct", samplesize = samplesize)
    posttree[is.nan(posttree[,"inftime.sd"]),"inftime.sd"] <- 
      posttree[is.nan(posttree[,"inftime.sd"]),"inftime.mean"]
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.construct.cov50",i] <- sum(postcov50)
    alloutputs[nextsims,"tr.construct.cov80",i] <- sum(postcov80)
    alloutputs[nextsims,"tr.construct.ior.0",i] <- sum((posttree["infector"] == truetree))
    alloutputs[nextsims,"tr.construct.ior.50",i] <- sum((posttree["infector"] == truetree) & postcov50)
    alloutputs[nextsims,"tr.construct.ior.80",i] <- sum((posttree["infector"] == truetree) & postcov80)
    alloutputs[nextsims,"tr.construct.tim.0",i] <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
    alloutputs[nextsims,"tr.construct.tim.50",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
    alloutputs[nextsims,"tr.construct.tim.80",i] <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
    alloutputs[nextsims,"tr.construct.2dg.0",i] <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
    alloutputs[nextsims,"tr.construct.2dg.50",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
    alloutputs[nextsims,"tr.construct.2dg.80",i] <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
    alloutputs[nextsims,"tr.construct.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.construct.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.construct.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##clade coverage for "mcc"
    truetree <- makephyloparset(
      make.phybreak.obkData(
        simulaties[[i]],use.tree = TRUE
      )$v$nodeparents)
    posttree <- MLphylo(res[[i]], phylo.class = FALSE, samplesize = samplesize)
    postcov50 <- (posttree["support"] >= .5*samplesize)[obs:(2*obs-1)]
    postcov80 <- (posttree["support"] >= .8*samplesize)[obs:(2*obs-1)]
    alloutputs[nextsims,"ph.mcc.cov50",i] <- sum(postcov50)
    alloutputs[nextsims,"ph.mcc.cov80",i] <- sum(postcov80)
    alloutputs[nextsims,"ph.mcc.cla.0",i] <- sum(equal.clades.phylo(posttree["parent"],truetree))
    alloutputs[nextsims,"ph.mcc.cla.50",i] <- sum(equal.clades.phylo(posttree["parent"],truetree) & postcov50)
    alloutputs[nextsims,"ph.mcc.cla.80",i] <- sum(equal.clades.phylo(posttree["parent"],truetree) & postcov80)
    alloutputs[nextsims,"ph.mcc.rspr",i] <- rspr(simulaties[[i]]@trees[[1]],MLphylo(res[[i]]))
    
    ##clade coverage for "cc.construct"
    #   posttree <- MLphylo(res[[i]], "cc.construct", phylo.class = FALSE, samplesize = samplesize)
    #   postcov50 <- (posttree["support"] >= .5*samplesize)[50:99]
    #   postcov80 <- (posttree["support"] >= .8*samplesize)[50:99]
    #   alloutputs[nextsims,"ph.construct.cov50",i] <- sum(postcov50)
    #   alloutputs[nextsims,"ph.construct.cov80",i] <- sum(postcov80)
    #   alloutputs[nextsims,"ph.construct.cla.0",i] <- sum(equal.clades.phylo(posttree["parent"],truetree))
    #   alloutputs[nextsims,"ph.construct.cla.50",i] <- sum(equal.clades.phylo(posttree["parent"],truetree) & postcov50)
    #   alloutputs[nextsims,"ph.construct.cla.80",i] <- sum(equal.clades.phylo(posttree["parent"],truetree) & postcov80)
    #   alloutputs[nextsims,"ph.construct.rspr",i] <- rspr(simulaties[[i]]@trees[[1]],MLphylo(res[[i]],"cc.construct"))
    
  }
  
  alloutputs[nextsims,,11] <- rowSums(alloutputs[nextsims,,1:10])
  alloutputs[nextsims,,12] <- alloutputs[nextsims,,11]/
    c(10,10,10*obs,10*obs,
      rep(c(10*obs,10*obs,
            alloutputs[nextsims,"tr.count.cov50",11],
            alloutputs[nextsims,"tr.count.cov80",11]), 2),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.edmunds.cov50",11],
            alloutputs[nextsims,"tr.edmunds.cov80",11]), 4),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.mpc.cov50",11],
            alloutputs[nextsims,"tr.mpc.cov80",11]), 4),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.mtcc.cov50",11],
            alloutputs[nextsims,"tr.mtcc.cov80",11]), 4),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.construct.cov50",11],
            alloutputs[nextsims,"tr.construct.cov80",11]), 4),
      10*obs, 10*obs, 
      10*obs, alloutputs[nextsims,"ph.mcc.cov50", 11],
      alloutputs[nextsims,"ph.mcc.cov80", 11], 10#,
      #     10*obs, 10*obs, 
      #     10*obs, alloutputs[nextsims,"ph.construct.cov50", 11],
      #     alloutputs[nextsims,"ph.construct.cov80", 11], 10
    )
  
  saveRDS(alloutputs, file = "resultsplussummary")
  
}


for(nextsims in 1:208) {
  load(paste0("../Simulations/sims_",
              substring(resfiles[nextsims], 7)))
  res <- readRDS(resfiles[nextsims])
  print(resfiles[nextsims])
  
  for(i in 1:10) {
    print(i)
    obs <- res[[i]]$p$obs

    ##additional burnin steps
    xburnin <- alloutputs[nextsims,"xburnin",i]
    
    ##required sample size
    samplesize <- alloutputs[nextsims,"samplesize",i]
    
    ##infector and infection time coverage for "count"
    truetree <- simulaties[[i]]@individuals$infector

    ##infector and infection time coverage for "edmunds"
    posttree <- MLtrans(res[[i]], method = "edmunds", samplesize = samplesize)
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.edmunds.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.edmunds.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.edmunds.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##infector and infection time coverage for "mpc"
    posttree <- MLtrans(res[[i]], method = "mpc", samplesize = samplesize)
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.mpc.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.mpc.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.mpc.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##infector and infection time coverage for "mtcc"
    posttree <- MLtrans(res[[i]], method = "mtcc", samplesize = samplesize)
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.mtcc.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.mtcc.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.mtcc.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    
    ##infector and infection time coverage for "cc.construct"
    posttree <- MLtrans(res[[i]], method = "cc.construct", samplesize = samplesize)
    postcov50 <- posttree["support"] >= .5*samplesize
    postcov80 <- posttree["support"] >= .8*samplesize
    alloutputs[nextsims,"tr.construct.cla.0",i] <- sum(equal.clades.trans(posttree["infector"],truetree))
    alloutputs[nextsims,"tr.construct.cla.50",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov50)
    alloutputs[nextsims,"tr.construct.cla.80",i] <- sum(equal.clades.trans(posttree["infector"],truetree) & postcov80)
    

  }
  
  alloutputs[nextsims,,11] <- rowSums(alloutputs[nextsims,,1:10])
  alloutputs[nextsims,,12] <- alloutputs[nextsims,,11]/
    c(10,10,10*obs,10*obs,
      rep(c(10*obs,10*obs,
            alloutputs[nextsims,"tr.count.cov50",11],
            alloutputs[nextsims,"tr.count.cov80",11]), 2),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.edmunds.cov50",11],
            alloutputs[nextsims,"tr.edmunds.cov80",11]), 4),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.mpc.cov50",11],
            alloutputs[nextsims,"tr.mpc.cov80",11]), 4),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.mtcc.cov50",11],
            alloutputs[nextsims,"tr.mtcc.cov80",11]), 4),
      10*obs, 10*obs, 
      rep(c(10*obs,alloutputs[nextsims,"tr.construct.cov50",11],
            alloutputs[nextsims,"tr.construct.cov80",11]), 4),
      10*obs, 10*obs, 
      10*obs, alloutputs[nextsims,"ph.mcc.cov50", 11],
      alloutputs[nextsims,"ph.mcc.cov80", 11], 10#,
      #     10*obs, 10*obs, 
      #     10*obs, alloutputs[nextsims,"ph.construct.cov50", 11],
      #     alloutputs[nextsims,"ph.construct.cov80", 11], 10
    )
  
  saveRDS(alloutputs, file = "resultsplussummary")
  
}



alloutputs <- readRDS("resultsplussummary")


allres <- array(alloutputs, dim = c(13,4,4,74,12), dimnames = list(
  episim = c("default", substring(dimnames(alloutputs)[[1]][2:13], 11)),
  wh.sim = paste0("sim_",substring(dimnames(alloutputs)[[1]][seq(1,40,13)], 7)),
  wh.res = paste0("res_",substring(dimnames(alloutputs)[[1]][seq(1,40,13)], 7)),
  output = dimnames(alloutputs)[[2]],
  which = dimnames(alloutputs)[[3]]
))


comparebybest <- function(episim, output) {
  f <- function(wh.sim) {
    rowSums(
    rep(apply(
      allres[episim, wh.sim,, output,1:10],2,max
    ),each=4) == 
      allres[episim, wh.sim,, output,1:10]
  )
  }
  sapply(dimnames(allres)[[2]], f)
} 
summarybybest <- function(episim, output) {
  mat <- comparebybest(episim, output)
  return(array(c(rowSums(mat),sum(diag(mat))), 
               dimnames = list(model = c(dimnames(allres)[[3]], "res_own"))))
}

rowSums(sapply(dimnames(allres)[[1]], summarybybest, output = "tr.edmunds.cla.0"))

summarybycount <- function(episim, output) {
  mat <- allres[episim,,,output,"mean"]/4
  return(array(c(colSums(mat),sum(diag(mat))), 
               dimnames = list(model = c(dimnames(allres)[[3]], "res_own"))))
}

rowSums(sapply(dimnames(allres)[[1]], summarybycount, output = "tr.edmunds.cla.0"))
sapply(dimnames(allres)[[1]], summarybycount, output = "tr.count.ior.0")

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
