rewire_change_infector_complete <- function(ID, newinfector) {
  btnodes <- 2 * pbe0$d$nsamples + ID - 1
  
  if(pbe1$v$infectors[ID] == 0) {
    coalnodes <- c()
    pbe1$v$nodeparents[btnodes] <- -1L
  } else {
    coalnodes <- take_cnode(btnodes)
  }
  
  pbe1$v$nodehosts[btnodes] <- newinfector
  pbe1$v$infectors[ID] <- newinfector
  return(coalnodes)
}

rewire_within_complete_edgewise <- function(ID, tinf) {
  coalnodes <- which(pbe1$v$nodehosts == ID & pbe1$v$nodetypes == "c")
  coaltimes_old <- pbe1$v$nodetimes[coalnodes]
  
  coaltimes_new <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == ID & pbe1$v$nodetypes != "c"],
                                    tinf, pbe1$p)
  coaltimes_new <- coaltimes_new[coaltimes_new > tinf]
  
  pbe1$logLiktoporatio <- pbe1$logLiktoporatio - lik_topology_host(pbe1, ID)
  
  coaltimes_new <- coaltimes_new[rank(coaltimes_old)]
  
  pbe1$v$nodetimes[coalnodes] <- coaltimes_new
  btnodes <- 2 * pbe0$d$nsamples + ID - 1
  pbe1$v$nodetimes[btnodes] <- tinf
  pbe1$v$inftimes[ID] <- tinf
  
  if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == ID] - 
         pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == ID]] < 0)) {
    pbe1$logLiktoporatio <- -Inf
    return()
  } 
  
  pbe1$logLiktoporatio <- pbe1$logLiktoporatio + lik_topology_host(pbe1, ID)
}

rewire_pathA_complete_edgewise <- function() {
  rewire_change_infector_complete(pbe1$hostID, 0L)
  rewire_within_complete_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_pullnodes_complete(0L)
  }
}


rewire_pathB_complete_edgewise <- function() {
  ### Identify new index
  newindex <- which(pbe1$v$inftimes == sort(pbe1$v$inftimes)[2])
  
  rewire_change_infector_complete(newindex, 0L)
  rewire_change_infector_complete(pbe1$hostID, pbe1$infector.proposed.ID)
  rewire_within_complete_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_pullnodes_complete(0)
    rewire_pullnodes_complete(pbe1$infector.proposed.ID)
  }
}

rewire_pathCF_complete_edgewise <- function() {
  ### Identify new infector and old infector
  newinfector <- which(pbe1$v$infectors == pbe1$hostID)
  newinfector <- newinfector[which(pbe1$v$inftimes[newinfector] == min(pbe1$v$inftimes[newinfector]))]
  oldinfector <- pbe1$v$infectors[pbe1$hostID]
  
  firstinfectiontime <- pbe1$v$inftimes[pbe1$hostID]
  secondinfectiontime <- pbe1$v$inftimes[newinfector]
  
  rewire_change_infector_complete(newinfector, oldinfector)
  rewire_within_complete_edgewise(newinfector, firstinfectiontime)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_change_infector_complete(pbe1$hostID, newinfector)
    rewire_within_complete_edgewise(pbe1$hostID, secondinfectiontime)
  }
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_pullnodes_complete(newinfector)
    rewire_pullnodes_complete(oldinfector)
  }
}

rewire_pathD_complete_edgewise <- function() {
  currentindex <- which(pbe1$v$infectors == 0)
  
  rewire_change_infector_complete(pbe1$hostID, 0L)
  rewire_within_complete_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_change_infector_complete(currentindex, pbe1$hostID)
    rewire_pullnodes_complete(pbe1$hostID)
    rewire_pullnodes_complete(0L)
  }
}

rewire_pathE_complete_edgewise <- function() {
  rewire_change_infector_complete(pbe1$hostID, pbe1$infector.proposed.ID)
  rewire_within_complete_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) rewire_pullnodes_complete(pbe1$infector.proposed.ID)
}
