### consensus phylogenetic tree from phybreak-object ###


### maximum clade credibility tree
### calls:
# .mcctree  ##C++
# .makephyloparsets
# .makephylo2

phylotree <- function(phybreak.object,
                    samplesize = Inf, support = "proportion", clade.times = TRUE,
                    time.quantiles = c(.025, 0.5, 0.975), phylo.class = FALSE) {
  chainlength <- length(phybreak.object$s$mu)
  
  ### tests
  if(chainlength == 0) stop("no sampled trees available")
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  if(support != "proportion" && support != "count") {
    warning("support is given as proportion")
  }
  if(class(time.quantiles) != "numeric") {
    stop("time.quantiles should be numeric")
  }
  if(max(time.quantiles) > 1 || min(time.quantiles) < 0) {
    stop("time.quantiles should be between 0 and 1")
  }
  
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs
  
  res <- c()
  
#  if(method[1] == "mcc") {
    res <- .mcctree(
      .makephyloparsets(phybreak.object$s$nodeparents[,
                                                      (1:samplesize) + chainlength - samplesize]),
      phybreak.object$s$nodetimes[1:(obs - 1),
                                  (1:samplesize) + chainlength - samplesize],
      c(obs, samplesize))
    if(phylo.class) return(get.phylo(phybreak.object, tail(res,1) + chainlength - samplesize, TRUE))
    res <- matrix(head(res,-1), ncol = 5)
    res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs]
    res[1:obs, 5] <- phybreak.object$v$nodetimes[1:obs]
#  }
#   if(method[1] == "cc.construct") {
#     res <- matrix(.CCphylotreeconstruct(
#       .makephyloparsets(phybreak.object$s$nodeparents[,
#                                                       (1:samplesize) + chainlength - samplesize]),
#       phybreak.object$s$nodetimes[1:(obs - 1),
#                                   (1:samplesize) + chainlength - samplesize],
#       c(obs, samplesize)
#     ), ncol = 4)
#     res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs]
#   }
#   
#   if(length(res) == 0) {
#     stop("incorrect method provided, choose \"mcc\" or \"cc.construct\"")
#   }
  
  parents.out <- matrix(res[,1],ncol = 1,
                        dimnames=list(1:(2*obs-1),"parent"))
#   if(method[1] == "mcc") {
  return(
    data.frame(
      parents = parents.out,
      support = res[,2],
      node.times.mean = res[,3],
      node.times.sd = res[,4],
      node.times.mc.tree = res[,5]
    )
  )
  #   } else {
#       return(
#         data.frame(
#           parents = parents.out,
#           support = res[,2],
#           nodetime.mean = res[,3],
#           nodetime.sd = res[,4]
#         )
#       )
#   }
  
  
}
