logLik.phybreak <- function(object, genetic = TRUE, withinhost = TRUE, sampling = TRUE, generation = TRUE, ...) {
  res <- 0
  if (genetic) {
    resSeq <- rep(0, dim(object$v$nodetimes)[1])
    for (ii in 1:dim(object$v$nodetimes)[1]) {
      object2 <- extractTree(object,ii)              # Extract tree for gene ii 
      res <- res + with(object2, .likseq(matrix(unlist(d$sequences), ncol = p$obs),
                                         attr(d$sequences, "weight"),
                                         v$nodeparents, v$nodetimes, p$mu, p$obs))
    }
  }
  
  if (generation) {res <- res + with(object, .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes[1, ], v$nodehosts, v$nodetypes)) }
  if (sampling) {res <- res + with(object, .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes[1, ], v$nodetypes))   }
  if (withinhost) { res <- res+ with(object, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes,v$reassortment)) }
  
  attributes(res) <- list(
    nobs = object$p$obs,
    df = 1 + object$h$est.mG + object$h$est.mS + object$h$est.wh,
    genetic = genetic, withinhost = withinhost, sampling = sampling, generation = generation
  )
  class(res) <- "logLik"
  return(res)
}


### calculate the log-likelihood of sampling intervals called from: logLik.phybreak .build.phybreakenv .propose.phybreakenv
.lik.gentimes <- function(obs, shapeG, meanG, nodetimes, nodehosts, nodetypes) {
  sum(dgamma(nodetimes[nodetypes == "t" & nodehosts > 0] - nodetimes[nodehosts[nodetypes == "t" & nodehosts > 0] + 2 * obs -
                                                                       1], shape = shapeG, scale = meanG/shapeG, log = TRUE))
}

### calculate the log-likelihood of generation intervals called from: logLik.phybreak .build.phybreakenv .propose.phybreakenv
.lik.sampletimes <- function(shapeS, meanS, nodetimes, nodetypes) {
  sum(dgamma(nodetimes[nodetypes == "s"] - nodetimes[nodetypes == "t"], shape = shapeS, scale = meanS/shapeS, log = TRUE))
}

# Specify function inputs
### calculate the log-likelihood of coalescent intervals called from: logLik.phybreak .build.phybreakenv .propose.phybreakenv
.lik.coaltimes <- function(obs, wh.model, slope, nodetimes, nodehosts, nodetypes, reassortment = c()) {
  if (wh.model == 1 || wh.model == 2)
    return(0)
  Ngenes <- dim(nodetimes)[1]
   
  coalnodes <- nodetypes == "c"
  hostR1 <- which(reassortment==1)
  orderednodes <- t(sapply(1:Ngenes, function(x) order(nodehosts, nodetimes[x, ])))
  orderedtouseAll <- t(sapply(1:Ngenes, function(x) orderednodes[x,c(duplicated(nodehosts[orderednodes[x,]])[-1], FALSE)]))
  
  # Separate orderedtouseAll by the absence/occurence of reassortment
  orderedtouseidxR1 <- matrix(orderedtouseAll[,c(which( (nodehosts[orderedtouseAll[1, ]] %in% hostR1) == 1))], nrow = Ngenes) # Reassortment occured
  orderedtouseidxR0 <- matrix(orderedtouseAll[,c(which( (nodehosts[orderedtouseAll[1, ]] %in% hostR1) == 0))], nrow = Ngenes)[1, ] # Absence of reassortment
  res <- 0
  for (gene in 1:Ngenes){
    orderedtouse <- orderedtouseidxR1[gene, ]
    ## make vectors with information on intervals between nodes
    coalno <- c(FALSE, head(coalnodes[orderedtouse], -1))    # interval starts with coalescence
    nodeho <- nodehosts[orderedtouse]                        # host in which interval resides
    coalmultipliers <- choose(2 + cumsum(2 * coalno - 1), 2) # coalescence coefficient
    
    ## from t to tau (time since infection)
    whtimes <- nodetimes[gene, ] - c(0, tail(nodetimes[gene, ], obs))[1 + nodehosts]
    noderates <- 1/(slope * whtimes[orderedtouse])
    
    # coalescence rate (per pair of lineages)
    nodeescrates <- log(whtimes[orderedtouse])/(slope)
    # cumulative coalescence rate since infection of host (per pair of lineages)
    
    escratediffs <- nodeescrates - c(0, head(nodeescrates, -1))
    escratediffs[!duplicated(nodeho)] <- nodeescrates[!duplicated(nodeho)]
    # cumulative coalescence rate within interval (per pair of lineages)
    
    ## First: coalescence rates at coalescence nodes
    logcoalrates <- log(noderates[c(coalno[-1], FALSE)])
    
    # Second: probability to escape coalescence in all intervals
    logescapes <- -escratediffs * coalmultipliers
    
    res <- res + sum(logcoalrates) + sum(logescapes)
  }
  
  ## Repeat previous part for hosts in which reassortment did not occur
  orderedtouse <- orderedtouseidxR0
  coalno <- c(FALSE, head(coalnodes[orderedtouse], -1))
  nodeho <- nodehosts[orderedtouse]
  coalmultipliers <- choose(2 + cumsum(2 * coalno - 1), 2)
  
  ## from t to tau (time since infection)
  whtimes <- nodetimes[1, ] - c(0, tail(nodetimes[1, ], obs))[1 + nodehosts]
  noderates <- 1/(slope * whtimes[orderedtouse])
  
  # coalescence rate (per pair of lineages)
  nodeescrates <- log(whtimes[orderedtouse])/(slope)
  
  # cumulative coalescence rate since infection of host (per pair of lineages)
  escratediffs <- nodeescrates - c(0, head(nodeescrates, -1))
  escratediffs[!duplicated(nodeho)] <- nodeescrates[!duplicated(nodeho)]
  # cumulative coalescence rate within interval (per pair of lineages)
  
  ## First: coalescence rates at coalescence nodes
  logcoalrates <- log(noderates[c(coalno[-1], FALSE)])
  
  # Second: probability to escape coalescence in all intervals
  logescapes <- -escratediffs * coalmultipliers
  
  return( res + sum(logcoalrates) + sum(logescapes))
}
