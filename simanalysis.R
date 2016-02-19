comparephybreak.infectors <- function(phybreak.object, sim.object) {
  res <- rep(NA,length(phybreak.object$s$logLik))
  for(i in 1:length(res)) {
    res[i] <- .equal.infectors(tail(phybreak.object$s$nodehosts[,i],phybreak.object$p$obs),
                               sim.object@individuals$infector)
  }
  return(res)
}



infsetcoverage <- function(phybreak.object, sim.obkData.object, percentile = 0.95, 
                           minsupport = 0, samplesize = Inf) {
  infsets <- infectorsets(phybreak.object, percentile, minsupport, samplesize)
  res <- rep(FALSE, phybreak.object$p$obs)
  for(i in 1:length(res)) {
    res[i] <- sim.obkData.object@individuals$infector[i] %in% infsets[[i]]
  }
  return(sum(res))
}


parsimonies <- sapply(make.multiPhylo(res[[1]]), 
                      parsimony, 
                      data = make.phyDat(res[[1]]))
##additional burnin steps
xburnin <- 100*ceiling((which(parsimonies == min(parsimonies))[1]-1)/100)

##required sample size
samplesize <- 1000
effsize <- min(effectiveSize(codares[(xburnin) + 1:samplesize,c("mu","logLik")]))
while(effsize < 200) {
  samplesize <- samplesize + 100
  effsize <- min(effectiveSize(codares[(1+xburnin) + 1:samplesize,c("mu","logLik")]))
}

##infection time coverage by posterior percentiles
codares <- make.coda(res[[1]],samplesize)
truetimes <- as.numeric(simulaties[[1]]@individuals$date - min(simulaties[[1]]@dna@meta$date))
inftimeperc <- sum((summary(codares[,2:51],quantiles = c(.025,.975))[[2]][,1]<truetimes &
  summary(codares[,2:51],quantiles = c(.025,.975))[[2]][,2]>truetimes))/res[[1]]$p$obs

##infector and infection time coverage for "count"
truetree <- simulaties[[1]]@individuals$infector
posttree <- MLtrans(res[[1]], samplesize = samplesize)
postcov50 <- posttree["support"] >= .5*samplesize
postcov80 <- posttree["support"] >= .8*samplesize
postcoverage.count.50 <- sum(postcov50)
postcoverage.count.80 <- sum(postcov80)
postinfector.count.0 <- sum((posttree["infector"] == truetree))
postinfector.count.50 <- sum((posttree["infector"] == truetree) & postcov50)
postinfector.count.80 <- sum((posttree["infector"] == truetree) & postcov80)
postinfectorset.count <- infsetcoverage(res[[1]], simulaties[[1]], minsupport = 0.05, samplesize = samplesize)
postinftimes.count.0 <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
postinftimes.count.50 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
postinftimes.count.80 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])

##infector and infection time coverage for "edmunds"
posttree <- MLtrans(res[[1]], method = "edmunds", samplesize = samplesize)
postcov50 <- posttree["support"] >= .5*samplesize
postcov80 <- posttree["support"] >= .8*samplesize
postcoverage.edmunds.50 <- sum(postcov50)
postcoverage.edmunds.80 <- sum(postcov80)
postinfector.edmunds.0 <- sum((posttree["infector"] == truetree))
postinfector.edmunds.50 <- sum((posttree["infector"] == truetree) & postcov50)
postinfector.edmunds.80 <- sum((posttree["infector"] == truetree) & postcov80)
postinftimes.edmunds.0 <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
postinftimes.edmunds.50 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
postinftimes.edmunds.80 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
post2nddegree.edmunds.0 <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
post2nddegree.edmunds.50 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
post2nddegree.edmunds.80 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
postclades.edmunds.0 <- sum(equal.clades(posttree["infector"],truetree))
postclades.edmunds.50 <- sum(equal.clades(posttree["infector"],truetree) & postcov50)
postclades.edmunds.80 <- sum(equal.clades(posttree["infector"],truetree) & postcov80)

##infector and infection time coverage for "mpc"
posttree <- MLtrans(res[[1]], method = "mpc", samplesize = samplesize)
postcov50 <- posttree["support"] >= .5*samplesize
postcov80 <- posttree["support"] >= .8*samplesize
postcoverage.mpc.50 <- sum(postcov50)
postcoverage.mpc.80 <- sum(postcov80)
postinfector.mpc.0 <- sum((posttree["infector"] == truetree))
postinfector.mpc.50 <- sum((posttree["infector"] == truetree) & postcov50)
postinfector.mpc.80 <- sum((posttree["infector"] == truetree) & postcov80)
postinftimes.mpc.0 <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
postinftimes.mpc.50 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
postinftimes.mpc.80 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
post2nddegree.mpc.0 <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
post2nddegree.mpc.50 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
post2nddegree.mpc.80 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
postclades.mpc.0 <- sum(equal.clades(posttree["infector"],truetree))
postclades.mpc.50 <- sum(equal.clades(posttree["infector"],truetree) & postcov50)
postclades.mpc.80 <- sum(equal.clades(posttree["infector"],truetree) & postcov80)

##infector and infection time coverage for "mtcc"
posttree <- MLtrans(res[[1]], method = "mtcc", samplesize = samplesize)
postcov50 <- posttree["support"] >= .5*samplesize
postcov80 <- posttree["support"] >= .8*samplesize
postcoverage.mtcc.50 <- sum(postcov50)
postcoverage.mtcc.80 <- sum(postcov80)
postinfector.mtcc.0 <- sum((posttree["infector"] == truetree))
postinfector.mtcc.50 <- sum((posttree["infector"] == truetree) & postcov50)
postinfector.mtcc.80 <- sum((posttree["infector"] == truetree) & postcov80)
postinftimes.mtcc.0 <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
postinftimes.mtcc.50 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
postinftimes.mtcc.80 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
post2nddegree.mtcc.0 <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
post2nddegree.mtcc.50 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
post2nddegree.mtcc.80 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
postclades.mtcc.0 <- sum(equal.clades(posttree["infector"],truetree))
postclades.mtcc.50 <- sum(equal.clades(posttree["infector"],truetree) & postcov50)
postclades.mtcc.80 <- sum(equal.clades(posttree["infector"],truetree) & postcov80)

##infector and infection time coverage for "cc.construct"
posttree <- MLtrans(res[[1]], method = "cc.construct", samplesize = samplesize)
postcov50 <- posttree["support"] >= .5*samplesize
postcov80 <- posttree["support"] >= .8*samplesize
postcoverage.construct.50 <- sum(postcov50)
postcoverage.construct.80 <- sum(postcov80)
postinfector.construct.0 <- sum((posttree["infector"] == truetree))
postinfector.construct.50 <- sum((posttree["infector"] == truetree) & postcov50)
postinfector.construct.80 <- sum((posttree["infector"] == truetree) & postcov80)
postinftimes.construct.0 <- sum(abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)
postinftimes.construct.50 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov50])
postinftimes.construct.80 <- sum((abs(posttree["inftime.mean"]-truetimes) - 2*posttree["inftime.sd"] < 0)[postcov80])
post2nddegree.construct.0 <- sum(dist.post.to.sim(posttree["infector"],truetree)<3)
post2nddegree.construct.50 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov50)
post2nddegree.construct.80 <- sum((dist.post.to.sim(posttree["infector"],truetree)<3) & postcov80)
postclades.construct.0 <- sum(equal.clades(posttree["infector"],truetree))
postclades.construct.50 <- sum(equal.clades(posttree["infector"],truetree) & postcov50)
postclades.construct.80 <- sum(equal.clades(posttree["infector"],truetree) & postcov80)

posttree <- MLphylo(res[[1]], phylo.class = FALSE, samplesize = samplesize)




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
equal.clades <- function(parset.post, parset.sim) {
  parset.post <- unlist(parset.post)
  pathspost <- matrix(FALSE,ncol=length(parset.post), nrow=length(parset.post))
  pathssim <- matrix(FALSE,ncol=length(parset.post), nrow=length(parset.post))

  for(i in 1:length(parset.post)) {
    path1 <- pathtoroot(parset.post, i)
    path2 <- pathtoroot(parset.sim, i)
    pathspost[i,1:length(path1)] <- path1
    pathssim[i,1:length(path2)] <- path2
  }
  res <- rep(FALSE, length(parset.post))
  
  for(i in 1:length(parset.post)) {
    res[i] <- any(apply(pathspost[,i] == pathssim, 2, all))
  }
  return(res)
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
