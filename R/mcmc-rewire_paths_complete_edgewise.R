rewire_pathA_complete_edgewise <- function() {
  ### First, deconnect hostID and change minitree in hostID
  # calculate old topology likelihood
  logLiktopo_old <- lik_topology_host(pbe1, pbe1$hostID)

  # coalescent nodes and time order
  coalnodes_hostID <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
  coalnodes_rank <- rank(pbe1$v$nodetimes[coalnodes_hostID])

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
  
  # new coalescent times
  coaltimes_new <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"],
                                    pbe1$tinf.prop, pbe1$p)[coalnodes_rank]
  
  # move inftime and nodetimes
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop
  pbe1$v$nodetimes[coalnodes_hostID] <- coaltimes_new

  # check and report consistency
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == pbe1$hostID]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  
  # calculate new topology likelihood
  logLiktoporatio <- lik_topology_host(pbe1, pbe1$hostID) - logLiktopo_old
  copy2pbe1("logLiktoporatio", environment())
  
  ### Second, reconnect hostID
}

rewire_pathB_complete_edgewise <- function() {
  ### Identify new index
  newindex <- which(pbe1$v$inftimes == sort(pbe1$v$inftimes)[2])
  
  ### First, deconnect hostID and change minitree in hostID
  # transmission node of new index
  transnode_ni <- 2*pbe1$d$nsamples - 1 + newindex
  
  # take coalescent node from hostID
  coalnode_tomove <- take_cnode(transnode_ni)

  # calculate old topology likelihood
  logLiktopo_old <- lik_topology_host(pbe1, pbe1$hostID)

  # coalescent nodes and time order
  coalnodes_hostID <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
  coalnodes_rank <- rank(pbe1$v$nodetimes[coalnodes_hostID])

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # new coalescent times
  coaltimes_new <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"],
                                    pbe1$tinf.prop, pbe1$p)[coalnodes_rank]
  
  # move inftime and nodetimes
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop
  pbe1$v$nodetimes[coalnodes_hostID] <- coaltimes_new
  
  # check and report consistency
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == pbe1$hostID]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  
  # calculate new topology likelihood
  logLiktoporatio <- lik_topology_host(pbe1, pbe1$hostID) - logLiktopo_old
  copy2pbe1("logLiktoporatio", environment())
  
  ### Second, change tree topology and reconnect hostID
  pbe1$v$infectors[pbe1$hostID] <- pbe1$infector.proposed.ID
  pbe1$v$infectors[newindex] <- 0L
  pbe1$v$nodeparents[c(transnode, transnode_ni)] <- pbe1$v$nodeparents[c(transnode_ni, transnode)]
  pbe1$v$nodehosts[transnode_ni] <- 0L
  pbe1$v$nodehosts[transnode] <- pbe1$infector.proposed.ID
  rewire_pullnodes_complete(pbe1$infector.proposed.ID)
}

rewire_pathCF_complete_edgewise <- function() {
  ### Identify new infector and old infector
  newinfector <- which(pbe1$v$infectors == pbe1$hostID)
  newinfector <- newinfector[which(pbe1$v$inftimes[newinfector] == min(pbe1$v$inftimes[newinfector]))]
  oldinfector <- pbe1$v$infectors[pbe1$hostID]
  
  ### First, deconnect hostID and new infector and change minitrees in hostID and new infector
  # transmission node of new infector
  transnode_ni <- 2*pbe1$d$nsamples - 1 + newinfector
  
  # take coalescent node from hostID
  coalnode_tomove <- take_cnode(transnode_ni)
  
  # calculate old topology likelihood
  logLiktopo_old <- lik_topology_host(pbe1, pbe1$hostID) + lik_topology_host(pbe1, newinfector)

  # coalescent nodes and time orders
  coalnodes_hostID <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
  coalnodes_rank_hostID <- rank(pbe1$v$nodetimes[coalnodes_hostID])
  coalnodes_ni <- which(pbe1$v$nodehosts == newinfector & pbe1$v$nodetypes == "c")
  coalnodes_rank_ni <- rank(pbe1$v$nodetimes[coalnodes_ni])

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
  
  # new coalescent times
  coaltimes_new_hostID <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"],
                                           pbe1$v$inftimes[newinfector], pbe1$p)[coalnodes_rank_hostID]
  coaltimes_new_ni <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == newinfector & pbe1$v$nodetypes != "c"],
                                                pbe1$v$inftimes[pbe1$hostID], pbe1$p)[coalnodes_rank_ni]
  
  # move inftime and nodetimes
  pbe1$v$inftimes[c(newinfector, pbe1$hostID)] <- pbe1$v$inftimes[c(pbe1$hostID, newinfector)]
  pbe1$v$nodetimes[c(transnode, transnode_ni)] <- pbe1$v$nodetimes[c(transnode_ni, transnode)]
  pbe1$v$nodetimes[coalnodes_hostID] <- coaltimes_new_hostID
  pbe1$v$nodetimes[coalnodes_ni] <- coaltimes_new_ni
  
  # check and report consistency
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == pbe1$hostID]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == newinfector] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == newinfector]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  
  # calculate new topology likelihood
  logLiktoporatio <- lik_topology_host(pbe1, pbe1$hostID) + lik_topology_host(pbe1, newinfector) - logLiktopo_old
  copy2pbe1("logLiktoporatio", environment())

  ### Second, change tree topology and reconnect hostID
  pbe1$v$infectors[newinfector] <- oldinfector
  pbe1$v$infectors[pbe1$hostID] <- newinfector
  pbe1$v$nodeparents[c(transnode, transnode_ni)] <- pbe1$v$nodeparents[c(transnode_ni, transnode)]
  pbe1$v$nodehosts[transnode_ni] <- oldinfector
  pbe1$v$nodehosts[transnode] <- newinfector
  rewire_pullnodes_complete(newinfector)
}

rewire_pathD_complete_edgewise <- function() {
  ### Identify old index
  oldindex <- which(pbe1$v$infectors == 0)
  
  ### First, deconnect hostID and change minitree in hostID
  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
  
  # take coalescent node from infector
  coalnode_tomove <- take_cnode(transnode)
  
  # calculate old topology likelihood
  logLiktopo_old <- lik_topology_host(pbe1, pbe1$hostID)

  # coalescent nodes and time order
  coalnodes_hostID <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
  coalnodes_rank <- rank(pbe1$v$nodetimes[coalnodes_hostID])

  # transmission node of oldindex
  transnode_oi <- 2*pbe1$d$nsamples - 1 + oldindex
  
  # new coalescent times
  coaltimes_new <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"],
                                    pbe1$tinf.prop, pbe1$p)[coalnodes_rank]

  # move inftime and nodetimes
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop
  pbe1$v$nodetimes[coalnodes_hostID] <- coaltimes_new
  
  # check and report consistency
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == pbe1$hostID]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  
  # calculate new topology likelihood
  logLiktoporatio <- lik_topology_host(pbe1, pbe1$hostID) - logLiktopo_old
  copy2pbe1("logLiktoporatio", environment())
  
  ### Second, change tree topology and reconnect hostID
  pbe1$v$infectors[oldindex] <- pbe1$hostID
  pbe1$v$infectors[pbe1$hostID] <- 0L
  pbe1$v$nodeparents[c(transnode, transnode_oi)] <- pbe1$v$nodeparents[c(transnode_oi, transnode)]
  pbe1$v$nodehosts[transnode] <- 0L
  pbe1$v$nodehosts[transnode_oi] <- pbe1$hostID
  rewire_pullnodes_complete(pbe1$hostID)
}

rewire_pathE_complete_edgewise <- function() {
  ### First, deconnect hostID and change minitree in hostID
  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
  
  # take coalescent node from infector
  coalnode_tomove <- take_cnode(transnode)
  
  # calculate old topology likelihood
  logLiktopo_old <- lik_topology_host(pbe1, pbe1$hostID)

  # coalescent nodes and time order
  coalnodes_hostID <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
  coalnodes_rank <- rank(pbe1$v$nodetimes[coalnodes_hostID])

  # new coalescent times
  coaltimes_new <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"],
                                    pbe1$tinf.prop, pbe1$p)[coalnodes_rank]

  # move inftime and nodetimes
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop
  pbe1$v$nodetimes[coalnodes_hostID] <- coaltimes_new
  
  # check and report consistency
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == pbe1$hostID] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == pbe1$hostID]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  
  # calculate new topology likelihood
  logLiktoporatio <- lik_topology_host(pbe1, pbe1$hostID) - logLiktopo_old
  copy2pbe1("logLiktoporatio", environment())
  
  ### Second, change tree topology and reconnect hostID
  pbe1$v$infectors[pbe1$hostID] <- pbe1$infector.proposed.ID
  pbe1$v$nodehosts[transnode] <- pbe1$infector.proposed.ID
  rewire_pullnodes_complete(pbe1$infector.proposed.ID)
}



