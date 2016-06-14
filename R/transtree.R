### consensus transmission tree from phybreak-object ###


### several methods available, default is "count" = counting infectors without resolving
### cycles or multiple/no index cases. Other options are "edmunds" = counting infectors
### with resolving cycles and multiple/no index cases; "mpc" = maximum parent credibility,
### selecting the sampled posterior tree with best support (summed count of infectors per host);
### "mtcc" = maximum transmission cluster credibility, selecting the sampled posterior
### tree with best support (summed count of clusters, where cluster = host + all progeny)
### calls:
# .transtreecount
# .transtreeedmonds
# .mpcinfector
# .mtcctree  ##C++
# get.phylo

transtree <- function(phybreak.object,
                    method = c("count", "edmonds", "mpc", "mtcc"),
                    samplesize = Inf, infector.name = TRUE, support = "proportion",
                    infection.times = "all", time.quantiles = c(.025, 0.5, 0.975),
                    phylo.class = FALSE) {
  chainlength <- length(phybreak.object$s$mu)
  
  ### tests
  if(chainlength == 0) stop("no sampled trees available")
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  if(support != "proportion" && support != "count") {
    warning("support is given as proportion")
  }
  if(infection.times != "all" && infection.times != "infector" && infection.times != "infector.sd") {
    warning("infection time summaries based on all samples")
  }
  if(infection.times == "all" || infection.times == "infector") {
    if(class(time.quantiles) != "numeric") {
      stop("time.quantiles should be numeric")
    }
    if(max(time.quantiles) > 1 || min(time.quantiles) < 0) {
      stop("time.quantiles should be between 0 and 1")
    }
  }
  
  ### initialization
  samplesize <- min(samplesize, chainlength)
  obs <- phybreak.object$p$obs
  res <- c()
  
  ### decision tree to use correct model
  if(method[1] == "count") {
    res <- .transtreecount(phybreak.object, samplesize, infection.times == "infector.sd")
  }
  if(method[1] == "edmonds") {
    res <- .transtreeedmonds(phybreak.object, samplesize, infection.times == "infector.sd")
  }
  if(method[1] == "mpc") {
    res <- .mpcinfector(phybreak.object, samplesize, phylo.class, infection.times == "infector.sd")
    if(phylo.class) return(get.phylo(phybreak.object, res, TRUE))
  }
  if(method[1] == "mtcc") {
    res <- .mtcctree(
      phybreak.object$s$nodehosts[obs:(2*obs-1),
                                  (1:samplesize) + chainlength - samplesize],
      phybreak.object$s$nodetimes[obs:(2*obs-1),
                                  (1:samplesize) + chainlength - samplesize],
      c(obs, samplesize)
    )
    if(phylo.class) return(get.phylo(phybreak.object, tail(res,1) + chainlength - samplesize, TRUE))
    res <- matrix(head(res,-1), ncol = 5)
  }
#   if(method[1] == "cc.construct") {
#     res <- matrix(.CCtranstreeconstruct(
#       phybreak.object$s$nodehosts[obs:(2*obs-1),
#                                   (1:samplesize) + chainlength - samplesize],
#       phybreak.object$s$nodetimes[obs:(2*obs-1),
#                                   (1:samplesize) + chainlength - samplesize],
#       c(obs, samplesize)
#     ), ncol = 4)
#   }
  ### test: no correct method provided
  if(length(res) == 0) {
    stop("incorrect method provided, choose \"count\", \"edmonds\",
         \"mpc\", or \"mtcc\"")
  }

  ### build output
  # infectors
  if(infector.name) {
    infectors.out <- matrix(c("index",phybreak.object$d$names)[1+res[,1]],
                            ncol = 1,
                            dimnames=list(phybreak.object$d$names,"infector"))
  } else {
    infectors.out <- matrix(res[,1],ncol = 1,
                            dimnames=list(phybreak.object$d$names,"infector"))
  }
  # support
  if(support == "count") {
    support.out <- res[,2]
  } else {
    support.out <- res[,2] / samplesize
  } 
  # times
  if(infection.times == "infector") {
    posttimes <- phybreak.object$s$nodetimes[obs:(2*obs-1),
                                             (1:samplesize) + chainlength - samplesize]
    posttimes[res[,1] != phybreak.object$s$nodehosts[obs:(2*obs-1),
                                                                  (1:samplesize) + chainlength - samplesize]] <- NA
    time.out <- t(apply(posttimes,1,quantile,probs=time.quantiles,na.rm=TRUE))
  } else if(infection.times == "infector.sd") {
    if(method[1] == "mpc" || method[1] == "mtcc") {
      time.out <- res[,3:5]
      colnames(time.out) <- c("mean","sd","mc.tree")
    } else {
      time.out <- res[,3:4]
      colnames(time.out) <- c("mean","sd")
    }
  } else {
    posttimes <- phybreak.object$s$nodetimes[obs:(2*obs-1),
                                             (1:samplesize) + chainlength - samplesize]
    time.out <- t(apply(posttimes,1,quantile,probs=time.quantiles))
  } 
  
  ### return the result
  return(
    data.frame(
      infectors = infectors.out,
      support = support.out,
      inf.times = time.out
    )
  )
}
