swap_heats <- function(heats, likelihood){#, heats, pbe_j, pbe_k){
  # proposed_swap <- sample(n, 2)
  # pbe_j <- envirs[[proposed_swap[1]]] #get(sprintf("pbe0_%s",proposed_swap[1]))
  # logLik_j <- sum(c(pbe_j$logLikcoal, pbe_j$logLikdist,  pbe_j$logLikgen, 
  #                   pbe_j$logLiksam, pbe_j$logLikseq))
  # heat_j <- pbe_j$heat
  # 
  # pbe_k <- envirs[[proposed_swap[2]]] #get(sprintf("pbe0_%s",proposed_swap[2]))
  # logLik_k <- sum(c(pbe_k$logLikcoal, pbe_k$logLikdist,  pbe_k$logLikgen, 
  #                   pbe_k$logLiksam, pbe_k$logLikseq))
  # heat_k <- pbe_k$heat
  # 
  # if (heat_j >= heat_k)
  #   logacceptanceprob <- (heat_j - heat_k)*(logLik_k - logLik_j)
  # else 
  #   logacceptanceprob <- (heat_k - heat_j)*(logLik_j - logLik_k)
  # 
  # if (runif(1) < exp(logacceptanceprob)) {
  #   heat <- heat_j
  #   assign("heat", get("heat", environment()), envirs[[proposed_swap[2]]])
  #   heat <- heat_k
  #   assign("heat", get("heat", environment()), envirs[[proposed_swap[1]]])
  # }
  # 
  # return(envirs)
  
  proposed_swap <- sample(length(heats), 2)
  Lj <- likelihood[proposed_swap[1]]
  Lk <- likelihood[proposed_swap[2]]
  heatj <- heats[proposed_swap[1]]
  heatk <- heats[proposed_swap[2]]
  
  if (heatj >= heatk)
    logacceptanceprob <- (heatj - heatk)*(Lk - Lj)
  else
    logacceptanceprob <- (heatk - heatj)*(Lj - Lk)

  if (runif(1) < exp(logacceptanceprob)) {
    heats[proposed_swap] <- heats[rev(proposed_swap)]
  }
  
  return(heats)

}
