### calculate the log-likelihood ###

### phybreak functions called ###
# .likseq  ##C++
# .lik.gentimes
# .lik.sampletimes
# .lik.coaltimes


logLik.phybreak <- function(phybreak.object, genetic = TRUE, withinhost = TRUE,
                            sampling = TRUE, generation = TRUE) {
  res <- 0
  if(genetic) {
    res <- res + with(phybreak.object, .likseq(t(d$SNP), d$SNPfr,
                                            v$nodeparents, v$nodetimes, p$mu,p$obs))
  } 
  if(generation) {
    res <- res + with(phybreak.object, .lik.gentimes(p$obs, p$shape.gen, p$mean.gen,
                                                  v$nodetimes, v$nodehosts, v$nodetypes))
  } 
  if(sampling) {
    res <- res + with(phybreak.object, .lik.sampletimes(p$shape.sample, p$mean.sample,
                                                     v$nodetimes, v$nodetypes))
  } 
  if(withinhost) {
    res <- res + with(phybreak.object, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, 
                                                   v$nodetimes, v$nodehosts, v$nodetypes))
  }
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
