lik_genetic <- function(infectors, sequences, samplehosts, mu, bn) {
  pr_tree_general(
    transmissiontree = infectors,
    haplotypes = matrix(
      attr(sequences, "allLevels")[do.call(rbind, sequences)],
      nrow = length(sequences)),
    haplotypefrequencies = attr(sequences, "weight"),
    statecount = 4L,
    statecodes = matrix(attr(sequences, "contrast"), 
                        ncol = 4,
                        dimnames = list(
                          attr(sequences, "allLevels"),
                          attr(sequences, "levels")
                        )),
    samplehosts = samplehosts, mu = 1 - exp(-mu), bn = 1 - exp(-bn)
  )
}

pr_tree_general <- function(transmissiontree, haplotypes,
                            haplotypefrequencies, statecount,
                            statecodes, samplehosts, mu, bn) {
  allstates <- make_allstates(statecount)
  
  transitionmatrix <- make_transitionmatrix(allstates, mu, bn)
  
  observationprobabilities <- make_observationprobabilities(
    statecodes, haplotypes, obs = length(transmissiontree),
    allstates, samplehosts)
  
  conditionalprobabilities <-
    pr_subtree(which(transmissiontree == 0), transmissiontree,
               transitionmatrix, observationprobabilities)
  priorprobabilities <- 
    rowSums(transitionmatrix[, rowSums(allstates) == 1]) / ncol(allstates)
  
  toreturn <- log(colSums(conditionalprobabilities * priorprobabilities))
  toreturn <- sum(toreturn * haplotypefrequencies)
  
  return(toreturn)
}




### all possible start and finish states in a host and observation
make_allstates <- function(statecount) {
  toreturn <- c()
  for(i in (statecount - 1):0) {
    toreturn <- cbind(
      toreturn, rep(c(F, T), each = 2 ^ i)
    )
  }
  toreturn <- toreturn[-1,]
  return(toreturn)
}

### transition matrix between all states from previous to current host
### makes use of transition probabilities from 
# - number of variants in previous host
# - number of variants in previous & current host (identical)
# - number of variants in current host
make_transitionmatrix <- function(allstates, mu, bn) {
  statecount <- ncol(allstates)
  changearray <- array(0, dim = c(statecount, statecount, statecount))
  
  for(i in 1:statecount) {
    for(j in 1:i) {
      for(k in j:(statecount + j - i)) {
        cstart <- c(rep(1, i), rep(0, statecount - i))
        cfinish <- c(rep(1, j), rep(0, statecount - k), rep(1, k - j))
        changearray[i, j, k] <-
          pr_change(cstart, cfinish, mu, bn)
      }
    }
  }
  
  toreturn <- apply(allstates, 1, 
                    function(x) apply(allstates, 1, 
                                      function(y) pr_change_fromarray(x, y, changearray)))
  
  return(toreturn)
}
pr_change <- function(c_start, c_finish, mu, bn) {
  sum(
    sapply(1:sum(c_start),
           function(x) {
             bn ^ (x - 1) *
               (1 - bn) ^ (sum(c_start) > x) *
               (choose(sum(c_start * c_finish), x) /
                  choose(sum(c_start), x)) *
               mu ^ (sum(c_finish) - x) *
               (1 - mu) ^ (length(c_start) - sum(c_finish))
           })
  )
}
pr_change_fromarray <- function(c_start, c_finish, change_array) {
  if(sum(c_start * c_finish)) {
    change_array[sum(c_start), sum(c_start * c_finish), sum(c_finish)]
  } else 0
}

### for each host the conditional probabilities to observe the data
make_observationprobabilities <- function(statecodes, haplotypes,
                                          obs, allstates, samplehosts) {
  dataarray <- statecodes[haplotypes, ]
  dim(dataarray) <- c(dim(haplotypes), ncol(statecodes))
  
  toreturn <- array(0, dim = c(ncol(haplotypes), 
                               obs, 
                               nrow(allstates)))
  for(haplo in 1:ncol(haplotypes)) {
    for(host in 1:obs) {
      toreturn[haplo, host, ] <-
        apply(
          apply(dataarray[samplehosts == host, haplo, , drop = FALSE], 1,
                function(x) apply(allstates, 1, pr_datapoint, c_data = x)
          ), 1, prod)
    }
  }
  
  return(toreturn)
  
}
pr_datapoint <- function(c_host, c_data) {
  sum(c_host * c_data)/sum(c_host)
}

### the worker function: pruning the transmission tree
pr_subtree <- function(host, trtree, trmatrix, obsprobs) {
  p_leaves <- obsprobs[, host,]
  
  if(host %in% trtree) {
    infectees <- which(host == trtree)
    p_leaves <- p_leaves * t(apply(
      array(sapply(infectees, pr_subtree, 
                   trtree = trtree, trmatrix = trmatrix, 
                   obsprobs = obsprobs),
            dim = c(dim(p_leaves)[2:1], length(infectees))), 
      1:2, prod))
  }
  
  toreturn <- apply(p_leaves, 1, 
                    function(x) colSums(trmatrix * x))
  return(toreturn)
  
}
