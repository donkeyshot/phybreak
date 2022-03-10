add_history <- function(d, vars, pars, h, s, build = FALSE, hist.inf = histtime){
  ### Store a back-up of all variables
  if (build){
    v <- vars
    p <- pars
    copy2pbe2("d", environment())
    copy2pbe2("v", environment())
    copy2pbe2("p", environment())
    copy2pbe2("h", environment())
  }
  
  ### Check if history for 1 introduction should be added
  introductions <- sum(vars$infectors==0)
  if (introductions == 1 & build) buildoneintro <- TRUE
  else buildoneintro <- FALSE
  
  h.dist.median <- median(1/h$dist[-1,-1])
  h$dist <- rbind(1/h.dist.median, cbind(1/h.dist.median, h$dist))
  
  ### Extract variables 
  nodetimes <- vars$nodetimes
  nodeparents <- as.integer(vars$nodeparents)
  nodehosts <- as.integer(vars$nodehosts) + 1
  nodetypes <- vars$nodetypes
  inftimes <- c(hist.inf, vars$inftimes)
  infectors <- c(0, vars$infectors + 1)
  nhosts <- length(inftimes)
  nsamples <- sum(nodetypes %in% c("s", "x"))
  ncoalnodes <- sum(nodetypes == "c")
  ntransnodes <- sum(nodetypes == "t")
  
  names(nodetimes) <- c(rep("s", nsamples), 
                        rep("c", ncoalnodes),
                        rep("t", ntransnodes))
  names(nodeparents) <- c(rep("s", nsamples), 
                        rep("c", ncoalnodes),
                        rep("t", ntransnodes))
  
  ### Create artificial host
  nodehosts <- c(nodehosts[1:nsamples], 
                 rep(1, introductions -1), nodehosts[nodetypes == "c"],
                 0, nodehosts[nodetypes == "t"])
  
  ### Add artificial types
  nodetypes <- c(rep("s", nhosts-1),
                 rep("x", nsamples-nhosts+1),
                 rep("c", nsamples - 1),
                 rep("t", ntransnodes + 1))
    
  
  ### Add artificial times
  edgesin <- which(nodehosts == 1 & nodetypes != "c")
  edgeintimes <- inftimes[infectors == 1]
  
  newcoaltimes <- sample_coaltimes(edgeintimes, inftimes[1], pars)
  if(length(newcoaltimes)>0)
    newcoaltimes[1] <- sample(rnorm(1, mean = h$mrca.av, sd = h$mrca.sd))
  
  nodetimes <- c(nodetimes[1:nsamples], 
                 newcoaltimes, 
                 nodetimes[names(nodetimes) == "c"],
                 inftimes[1],
                 nodetimes[names(nodetimes) == "t"])

  ### Add artificial parents
  coalnodes <- which(nodetypes == "c")[1:length(newcoaltimes)]
  transnode <- which(nodetypes == "t")[1]
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]
  
  ### Sample topology of history host
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)
  
  
  nodeparents <- sapply(nodeparents, function(xx){
    if (xx == 0) xx <- 1
    else if (xx <= sum(nsamples + ncoalnodes)) xx <- xx + length(newcoaltimes)
    else xx <- xx + length(newcoaltimes) + 1
  })
  nodeparents <- c(nodeparents[names(nodeparents) == "s"],
                        rep(NA, length(newcoaltimes)),
                        nodeparents[names(nodeparents) == "c"],
                        0,
                        nodeparents[names(nodeparents) == "t"]) 
  
  if (buildoneintro) nodeparents[transnode+which(infectors == 1) -1] <- transnode
  else nodeparents[edgeend] <- edgestart
  
  ### Make element v
  v <- list()
  v$inftimes <- inftimes
  v$infectors <- infectors
  v$nodetypes <- nodetypes
  v$nodetimes <- nodetimes
  v$nodehosts <- nodehosts
  v$nodeparents <- nodeparents
  v$cultimes <- vars$cultimes
  
  if (sum(v$infectors==1) == sum(pbe0$hs$infectors==1)) {
    if (all(which(v$infectors==1) == which(pbe0$hs$infectors==1))) {
      if (all(v$inftimes[v$infectors==1] == pbe0$hs$inftimes[v$infectors==1])) {
        v$inftimes[1] <- pbe0$hs$inftimes[1]
        v$nodetimes[pars$obs+1:length(newcoaltimes)] <- pbe0$hs$nodetimes[pars$obs+1:length(newcoaltimes)]
        v$nodetimes[sum(v$nodetypes!="t")+1] <- v$inftimes[1]
        
        for(n in names(pbe1))
          copy2pbe2(n, pbe1)
        
        hs <- v
        copy2pbe1("v", environment())
        copy2pbe1("hs", environment())
        return()
      }
    }
  }
  
  return(list("v" = v, "h" = h))
}

remove_history <- function(x = NULL, keepenv = FALSE){
  ### Remove history host
  if (is.null(x)) {
    v <- pbe0$v
    p <- pbe0$p
    d <- pbe0$d
  } else {
    v <- x$v
    p <- x$p
    d <- x$d
  }
  
  v$inftimes <- tail(v$inftimes, p$obs)
  v$infectors <- tail(v$infectors, p$obs) - 1
  intros <- sum(v$infectors == 0)
  
  keephosts <- setdiff(which(v$nodehosts != 0), 
                       intersect(which(v$nodehosts == 1), which(v$nodetypes == "c")))
  v$nodehosts <- v$nodehosts[keephosts] - 1
  v$nodetypes <- v$nodetypes[keephosts]
  v$nodetimes <- v$nodetimes[keephosts]
  
  remnodes <- setdiff(1:length(v$nodeparents), keephosts)
  transnode <- 2*d$nsamples - length(remnodes) + 1
  nodeparents <- v$nodeparents
  nodeparents[nodeparents %in% remnodes] <- 0
  for (i in seq_along(remnodes)){
    remnode <- remnodes[i] - i + 1
    nodeparents <- sapply(nodeparents, function(xx) {
      if (xx > remnode) xx - 1 
      else xx
    })
    nodeparents <- nodeparents[-remnode]
  }
  v$nodeparents <- nodeparents
  
  if (is.null(x)){
    if (keepenv)
      copy2pbe2("v", environment())
    else {
      pbe0$h$dist <- pbe0$h$dist[,2:ncol(pbe0$h$dist)][2:nrow(pbe0$h$dist),]
      copy2pbe0("v", environment())
    }
  } else {
    return(v)
  }
}
