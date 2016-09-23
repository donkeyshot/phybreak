### phybreak class constructor ###

### phybreak functions called ###
# .rinftimes(sampling times, sample.mean, sample.shape)
# .rinfectors(infection times, generation.mean, generation.shape)
# .samplecoaltimes(tip times, wh.model, wh.slope)
# .sampletopology(nodes, nodetimes, nodetypes, rootnode, wh.model)
# .distmatrix(SNP, SNPfrequencies)

#' Create a phybreak-object from data and prior distributions.
#' 
#' phybreak takes as data either an \code{'obkData'}-object or a \code{matrix} with sequences 
#' (individuals in rows, nucleotides in columns). If only sequences are used, a vector with 
#' sampling times is also needed. Parameter values are used as initial values in the MCMC-chain 
#' or kept fixed. All variables are initialized by random samples from the prior distribution, 
#' unless the tree in the obkData is used (\code{use.tree = TRUE}).
#' 
#' @param data The sequence data. These can be given as objects of class \code{'DNAbin'} or \code{'phyDat'}, or in a \code{matrix}
#'  with nucleotides, each row a host, each column a nucleotide. All nucleotides that are not \code{'a'}, \code{'c'}, \code{'g'}, 
#'  or \code{'t'}, will be turned into \code{'n'}. If the sequences are named, these names will be used. It is also possible
#'  to provide the data in object of class \code{'obkData'}, containing sequences and sampling times as metadata with 
#'  these sequences.
#' @param times Vector of sampling times (not needed if the data are of class \code{'obkData'}). If the vector is named,
#'   these names will be used to identify the hosts.
#' @param mu Initial value for mutation rate (defined per site per unit of time). 
#'   NOTE: mutation is defined as assignment of a random nucleotide at a particular site; this could be the 
#'   nucleotide that was there before the mutation event. Therefore, the actual rate of change of nucleotides 
#'   is \code{0.75*mu}.
#' @param gen.shape Shape parameter of the generation interval distribution (not estimated).
#' @param gen.mean Initial value for the mean generation interval, i.e. the interval between infection of a secondary
#'   case by a primary case.
#' @param sample.shape Shape parameter of the sampling interval distribution (not estimated), i.e. the interval between
#'   infection and sampling of a host.
#' @param sample.mean Initial value for the mean sampling interval.
#' @param wh.model The model for within-host pathogen dynamics (effective pathogen population size = 
#'   N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Options are:
#'   \enumerate{
#'     \item Effective size = 0, so coalescence occurs 'just before' transmission in the infector
#'     \item Effective size = Inf, so coalescence occurs 'just after' transmission in the infectee
#'     \item Effective size at time t after infection = \code{wh.slope * t}
#'   }
#' @param wh.slope Initial value for the within-host slope, used if \code{wh.model = 3}.
#' @param est.gen.mean Whether to estimate the mean generation interval or keep it fixed. 
#' @param prior.mean.gen.mean Mean of the (gamma) prior distribution of mean generation interval \code{mG} 
#'   (only if \code{est.gen.mean = TRUE}).
#' @param prior.mean.gen.sd Standard deviation of the (gamma) prior distribution of mean generation interval \code{mG} 
#'   (only if \code{est.gen.mean = TRUE}).
#' @param est.sample.mean Whether to estimate the mean sampling interval or keep it fixed. 
#' @param prior.mean.sample.mean Mean of the (gamma) prior distribution of mean sampling interval \code{mS} 
#'   (only if \code{est.sample.mean = TRUE}).
#' @param prior.mean.sample.sd Standard deviation of the (gamma) prior distribution of mean sampling interval \code{mS} 
#'   (only if \code{est.sample.mean = TRUE}).
#' @param est.wh.slope Whether to estimate the within-host slope or keep it fixed. 
#' @param prior.wh.shape Shape parameter of the (gamma) prior distribution of \code{slope} 
#'   (only if \code{est.wh.slope = TRUE}).
#' @param prior.wh.mean Mean of the (gamma) prior distribution of \code{slope} 
#'   (only if \code{est.wh.slope = TRUE}).
#' @param use.tree Whether to use the transmission and phylogenetic tree given in data of class \code{'obkData'}, 
#'   to create a \code{phybreak}-object with an exact copy of the outbreak. This requires more data in \code{data}: 
#'   the slot \code{individuals} with vectors \code{infector} and \code{date}, and the slot \code{trees} with at least 
#'   one phylogenetic tree. Such data can be simulated with \code{\link{sim.phybreak}}.
#' @return An object of class \code{phybreak} with the following elements
#'   \describe{
#'     \item{d}{a \code{list} with data, i.e. names, sequences, sampling times, and total number of SNPs.}
#'     \item{v}{a \code{list} with current state of all nodes in the tree: times, hosts in which they reside,
#'       parent nodes, node types (sampling, coalescent, or transmission)}
#'     \item{p}{a \code{list} with the parameter values}
#'     \item{h}{a \code{list} with helper information for the MCMC-method: \code{si.mu} and \code{si.wh} for 
#'       efficiently proposing \code{mu} and \code{slope}, matrix \code{dist} with weights for infector sampling 
#'       based on sequence distances, \code{logical}s \code{est.mG}, \code{est.mS}, and \code{est.wh.slope} whether to estimate 
#'       mean generation interval \code{mG}, mean sampling interval \code{mS}, and within-host \code{slope}, 
#'       and parameters for the priors of \code{mG}, \code{mS}, and \code{slope}.
#'     }
#'     \item{s}{an empty \code{list} that will contain vector and matrices with the posterior samples; 
#'       in matrices, the rows are nodes in the phylogenetic tree, the columns are the samples}
#'   }
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' simulation <- sim.phybreak()
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' 
#' \dontrun{
#'  ### only if 'OutbreakTools' is installed
#'  simulation <- sim.phybreak(output-class = "obkData")
#'  MCMCstate <- phybreak(data = simulation, use.tree = TRUE)
#' }
#' 
#' sampletimedata <- c(0,2,2,4,4)
#' sampleSNPdata <- matrix(c("a","a","a","a","a",
#'                           "a","c","c","c","c",
#'                           "t","t","t","g","g"), nrow = 5)
#' MCMCstate <- phybreak(data = sampleSNPdata, times = sampletimedata)
#' @export
phybreak <- function(data, times = NULL, 
         mu = .0001, gen.shape = 3, gen.mean = 1,
         sample.shape = 3, sample.mean = 1, 
         wh.model = 3, wh.slope = 1,
         est.gen.mean = TRUE, prior.mean.gen.mean = 1, prior.mean.gen.sd = Inf,
         est.sample.mean = TRUE, prior.mean.sample.mean = 1, prior.mean.sample.sd = Inf,
         est.wh.slope = TRUE, prior.wh.shape = 1, prior.wh.mean = 1,
         use.tree = FALSE) {
  
  ###########################################
  ### check for correct classes and sizes ###
  ###########################################
  if(!any(class(data) %in% c("DNAbin", "phyDat", "matrix", "obkData"))) {
    stop("data should be of class \"DNAbin\", \"phyDat\", \"character\" or \"obkData\"")
  }
  if(inherits(data, "matrix") && !inherits(data[1], "character")) {
    stop("data matrix should contain \"character\" elements")
  }
  if(inherits(data, "DNAbin")) {
    data <- as.character(data)
  }
  if(inherits(data, "phyDat")) {
    data <- as.character(data)
  }
  if(!inherits(data, "matrix") & use.tree) {
    warning("tree can only be used if provided in obkData-format; random tree will be generated")
  }
  if(inherits(data, "obkData")) {
    if(!("OutbreakTools" %in% .packages(TRUE))) {
      stop("package 'OutbreakTools' not installed, while data-class is \"obkData\"")
    }
    if(!("OutbreakTools" %in% .packages(FALSE))) {
      warning("package 'OutbreakTools' is not attached")
    }
    if(is.null(data@dna)) {
      stop("no sequence data in data@dna")
    }
    obs <- nrow(data@dna@dna[[1]])
    if(length(data@dna@meta$date) != obs) {
      stop("length of data@dna@meta$date does not match number of sequences")
    }
    if(use.tree & length(data@individuals$infector) != obs) {
      stop("use.tree = TRUE, but length of data@individuals$infector does not match number of sequences")
    }
    if(use.tree & length(data@individuals$date) != obs) {
      stop("use.tree = TRUE, but length of data@individuals$date does not match number of sequences")
    }
  } else {
    if(!(class(times) %in% c("Date", "numeric", "integer"))) {
      stop("times should be numeric or of class \"Date\"")
    }
    if(length(times) != nrow(data)) {
      stop("data and times have different number of observations")
    }
    if((!is.null(names(times)) && !is.null(row.names(data))) && !all(names(times) %in% row.names(data))) {
      stop("names in sequence data and times are not equal")
    }
  }
  numFALSE <- 
    unlist(
      lapply(
        list(mu, gen.shape, gen.mean, sample.shape, sample.mean,
             wh.model, wh.slope, prior.mean.gen.mean, prior.mean.gen.sd,
             prior.mean.sample.mean, prior.mean.sample.sd,
             prior.wh.shape, prior.wh.mean),
        class
      )
    ) != "numeric"
  if(any(numFALSE)) {
    stop(paste0("parameters ",
                c("mu", "gen.shape", "gen.mean", "sample.shape", "sample.mean",
                  "wh.model", "wh.slope", "prior.mean.gen.mean", "prior.mean.gen.sd",
                  "prior.mean.shape.mean", "prior.mean.shape.sd",
                  "prior.wh.shape", "prior.wh.mean")[numFALSE],
                " should be numeric"))
  }
  numNEGATIVE <- 
    c(mu, gen.shape, gen.mean, sample.shape, sample.mean,
             wh.model, wh.slope, prior.mean.gen.mean, prior.mean.gen.sd,
             prior.mean.sample.mean, prior.mean.sample.sd,
             prior.wh.shape, prior.wh.mean) <= 0
  if(any(numNEGATIVE)) {
    stop(paste0("parameters ",
                c("mu", "gen.shape", "gen.mean", "sample.shape", "sample.mean",
                  "wh.model", "wh.slope", "prior.mean.gen.mean", "prior.mean.gen.sd",
                  "prior.mean.shape.mean", "prior.mean.shape.sd",
                  "prior.wh.shape", "prior.wh.mean")[numNEGATIVE],
                " should be positive"))
  }
  logFALSE <- 
    unlist(
      lapply(
        list(est.gen.mean, est.sample.mean, est.wh.slope, use.tree),
        class
      )
    ) != "logical"
  if(any(logFALSE)) {
    stop(paste0("parameters ",
                c("est.gen.mean", "est.sample.mean", "est.wh.slope", "use.tree")[logFALSE],
                " should be logical"))
  }
  
  ### outbreak parameters ###
  #outbreak size
  if(class(data) == "obkData") {
    obs <- nrow(data@dna@dna[[1]])
  } else {
    obs <- length(times)  
  }
  
  #################################
  ### first slot: sequence data ###
  #################################
  #names
  if(class(data) == "obkData") {
    hostnames <- rownames(data@dna@dna[[1]])
  } else {
    if(is.null(names(times))) {
      if(is.null(row.names(data))) {
        hostnames <- paste0("host.",1:obs)
      } else hostnames <- row.names(data)
    } else {
      hostnames <- names(times)
      if(!is.null(row.names(data))) {
        data <- data[match(hostnames, row.names(data)), ]
      }
    }
  }
  
  #sequences (SNP)
  SNP.sample <- c()
  if(class(data) == "obkData") {
    allgenes <- OutbreakTools::get.dna(data)
    for(i in 1:(OutbreakTools::get.nlocus(data))) {
      SNP.sample <- cbind(SNP.sample, as.character(allgenes[[i]]))
    }
    if(length(setdiff(SNP.sample, c("a","c","g","t")))) {
      warning("all nucleotides other than actg are turned into n")
    }
    SNP.sample[SNP.sample != "a" & SNP.sample != "c" & SNP.sample != "g" & SNP.sample != "t"] <- "n"
  } else {
    SNP.sample <- data
    if(length(setdiff(SNP.sample,c("a","c","g","t")))) warning("all nucleotides other than actg are turned into n")
    SNP.sample[SNP.sample != "a" & SNP.sample != "c" & SNP.sample != "g" & SNP.sample != "t"] <- "n"
  }
  seqdata <- phangorn::as.phyDat(SNP.sample)
  #SNP count
  nsnps <- sum( (apply(SNP.sample=="a", 2, any) +
                   apply(SNP.sample=="c", 2, any) +
                   apply(SNP.sample=="g", 2, any) +
                   apply(SNP.sample=="t", 2, any) +
                   apply(SNP.sample=="n", 2, all) - 1))
  
  
  ############################
  ### initialize variables ###
  ############################
  #auxiliary variables
  if(class(data) == "obkData") {
    refdate <- min(OutbreakTools::get.dates(data, "dna"))
  } else refdate <- min(times)
  #transmission tree
  if(class(data) == "obkData") {
    if(use.tree) {
      inftimes <- as.numeric(OutbreakTools::get.dates(data, "individuals") - refdate)
      infectors <- data@individuals$infector
    } else {
      samtimes <- as.numeric(OutbreakTools::get.dates(data, "dna") - refdate)
      inftimes <- .rinftimes(samtimes, sample.mean, sample.shape)
      infectors <- .rinfectors(inftimes, gen.mean, gen.shape)
    }
  } else {
    samtimes <- as.numeric(times - refdate)
    inftimes <- .rinftimes(samtimes, sample.mean, sample.shape)
    infectors <- .rinfectors(inftimes, gen.mean, gen.shape)
  }
  #phylogenetic tree
  nodetypes <- c(rep("s",obs), rep("c",obs - 1),
                 rep("t",obs))  #node type (sampling, coalescent, transmission)
  if(class(data) == "obkData" & use.tree) {
    ##extract phylotree from obkData
    phytree <- OutbreakTools::get.trees(data)[[1]]
    phylobegin <- phytree$edge[,1]
    phyloend <- phytree$edge[,2]
    phylolengths <- phytree$edge.length
    
    ##initialize nodetimes with infection times, nodeparents with phylotree
    ##initialize edgelengths to help build the tree
    nodetimes <- c(rep(NA,2*obs-1), inftimes)
    nodeparents <- c(head(phylobegin[order(phyloend)], obs),
                     NA,
                     tail(phylobegin[order(phyloend)], obs - 2),
                     rep(NA, obs))
    edgelengths <- c(head(phylolengths[order(phyloend)], obs),
                     NA, tail(phylolengths[order(phyloend)], obs - 2),
                     rep(NA, obs))
    ##link first infection to root node of phylotree
    nodeparents[obs + 1] <- which(infectors == 0)+ 2 * obs - 1
    edgelengths[obs + 1] <- -inftimes[which(infectors == 0)] - ape::node.depth.edgelength(phytree)[1]
    ##place transmission nodes between sampling and coalescent nodes for hosts without secondary cases
    #a. place coalescent node before transmission node
    nodeparents[(2*obs - 1) + setdiff(1:obs, infectors)] <- nodeparents[setdiff(1:obs, infectors)]
    #b. place transmission node before sampling node
    nodeparents[setdiff(1:obs,infectors)] <- (2*obs - 1) + setdiff(1:obs, infectors)
    
    ##complete nodetimes by adding branch lengths starting from root
    while(any(is.na(nodetimes))) {
      nodetimes[1:(2*obs - 1)] <- nodetimes[nodeparents[1:(2*obs - 1)]] + edgelengths[1:(2*obs - 1)]
    }
    nodetimes <- round(nodetimes, digits = 12) #prevent numerical problems
    
    ##initialize nodehosts, first sampling and transmission nodes...
    nodehosts <- c(1:obs, which(infectors == 0), rep(NA, obs - 2), infectors)   
    ##... then coalescent nodes that are parent of sampling nodes of infectors
    nodehosts[nodeparents[sort(unique(infectors))[-1]]] <- sort(unique(infectors))[-1]
    
    ##complete nodehosts by going backwards in phylotree
    while(any(is.na(nodehosts))) {
      #which are coalescent nodes with known parent & host?; order these by parent
      whichnodes <- tail(which(!is.na(nodehosts) & !is.na(nodeparents)),-obs)
      whichnodes <- whichnodes[order(nodeparents[whichnodes])]

      #which of these nodes have the same parent node?
      sameparent.wn <- which(duplicated(nodeparents[whichnodes]))
      #determine the host of these parent nodes
      for(i in sameparent.wn) {
        if(nodehosts[whichnodes[i]] == nodehosts[whichnodes[i-1]]) {
          #...if host of children nodes are the same, this is also host of parent
          nodehosts[nodeparents[whichnodes[i]]] <- nodehosts[whichnodes[i]]
        } else {
          #...if host of children nodes are different, host of parent node is the
          #...infector if one infects the other, or the common infector of both
          posshosts <- nodehosts[whichnodes[i - 1:0]]
          posshosts <- c(posshosts, infectors[posshosts])
          nodehosts[nodeparents[whichnodes[i]]] <- posshosts[duplicated(posshosts)]
        }
      }
    }
    
    ##complete nodeparents: place transmission nodes for hosts with secondary infections
    #a. place coalescent node before transmission node
    nodeparents[nodehosts[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != 
                                  nodehosts[(obs + 1):(2 * obs - 1)]) + obs] + 2 * obs - 1] <- 
      nodeparents[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != nodehosts[(obs + 1):(2 * obs - 1)]) + obs]
    #b. place transmission node before sampling node
    nodeparents[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != nodehosts[(obs + 1):(2 * obs - 1)]) + obs] <- 
      nodehosts[which(nodehosts[nodeparents[(obs + 1):(2 * obs - 1)]] != nodehosts[(obs + 1):(2 * obs - 1)]) + obs] + 2 * obs - 1
    nodeparents[nodehosts == 0] <- 0
  } else {
    ##nodehosts, coalescent nodes in hosts with secondary cases
    nodehosts <- c(1:obs, sort(infectors)[-1], infectors)
    
    ##initialize nodetimes with sampling and infection times
    nodetimes <- c(samtimes, rep(NA, obs - 1), inftimes)
    ##sample coalescent times
    for(i in 1:obs) {
      nodetimes[nodehosts == i & nodetypes == "c"] <-   #change the times of the coalescence nodes in host i...
        inftimes[i] +                      #...to the infection time +
        .samplecoaltimes(nodetimes[nodehosts == i & nodetypes != "c"] - inftimes[i],
                         wh.model, wh.slope)  #...sampled coalescent times
    }
    
    ##initialize nodeparents 
    nodeparents <- 0*nodetimes
    ##sample nodeparents
    for(i in 1:obs) {
      nodeparents[nodehosts == i] <-     #change the parent nodes of all nodes in host i...
        .sampletopology(which(nodehosts == i), nodetimes[nodehosts == i], nodetypes[nodehosts == i], i + 2*obs - 1, wh.model)
      #...to a correct topology, randomized where possible
    }
  }
  
  ################################
  ### make the phybreak object ###
  ################################
  res <- list(
    d = list(
      names = hostnames,
      sequences = seqdata,
      sample.times = times,
      nSNPs = nsnps
    ),
    v = list(
      nodetimes = nodetimes,
      nodehosts = nodehosts,
      nodeparents = nodeparents,
      nodetypes = nodetypes
    ),
    p = list(
      obs = obs,
      mu = mu,
      mean.sample = sample.mean,
      mean.gen = gen.mean,
      shape.sample = sample.shape,
      shape.gen = gen.shape,
      wh.model = wh.model,
      wh.slope = wh.slope
    ),
    h = list(si.mu = if(nsnps == 0) 0 else 2.38*sqrt(trigamma(nsnps)),
             si.wh = 2.38*sqrt(trigamma(obs - 1)),
             dist = .distmatrix(SNP.sample),
             est.mG = est.gen.mean,
             est.mS = est.sample.mean,
             est.wh = est.wh.slope,
             mG.av = prior.mean.gen.mean,
             mG.sd = prior.mean.gen.sd,
             mS.av = prior.mean.sample.mean,
             mS.sd = prior.mean.sample.sd,
             wh.sh = prior.wh.shape,
             wh.av = prior.wh.mean),
    s = list(
      nodetimes = c(),
      nodehosts = c(),
      nodeparents = c(),
      mu = c(),
      mG = c(),
      mS = c(),
      slope = c(),
      logLik = c()
    )
    )
  
  class(res) <- c("phybreak", "list")
  
  return(res)
}


### random infection times given sampling times and sampling interval distribution
### called from: 
# phybreak
.rinftimes <- function(st, meanS, shapeS) {
  ### tests
  if(class(st) != "numeric" && class(st) != "integer") {
    stop(".rinftimes called with non-numeric sampling times")
  }
  if(meanS <= 0) stop(".rinftimes called with non-positive mean sampling interval")
  if(shapeS <= 0) stop(".rinftimes called with non-positive shape parameter")
  
  ### function body
  st - rgamma(length(st), shape = shapeS, scale = meanS/shapeS)
}

### random infectors times given infection times and generation interval distribution
### called from:
# phybreak
.rinfectors <- function(it, meanG, shapeG) {
  ### tests
  if(class(it) != "numeric" && class(it) != "integer") {
    stop(".rinfectors called with non-numeric infection times")
  }
  if(sum(it == min(it)) > 1) stop("rinfectors with >1 index case")
  if(meanG <= 0) stop(".rinfectors called with non-positive mean generation interval")
  if(shapeG <= 0) stop(".rinfectors called with non-positive shape parameter")
  
  ### function body
  res <- rep(0,length(it))
  for(i in 1:length(it)) {
    if(it[i] > min(it)) {
      dist <- dgamma(it[i] - it, shape = shapeG, scale = meanG/shapeG)
      dist[i] <- 0
      res[i] <- sample(length(it), 1, prob = dist)
    }
  }
  return(res)
}

### distance matrix between sequences given SNP data
### called from:
# phybreak
.distmatrix <- function(seqmatrix) {

  # count SNPs excluding "n"
  res <- as.matrix(ape::dist.dna(ape::as.DNAbin(seqmatrix), model = "N", pairwise.deletion = TRUE))

  # prob of SNP per nucleotide in most distant entry
  nscore <- max(res)/ncol(seqmatrix)

  # add nscore for each missing nucleotide
  nmatrix <- seqmatrix == "n"
  countns <- function(s1, s2) sum(nmatrix[s1,] | nmatrix[s2,])  
  for(i in 1:nrow(res)) {
    res[i, ] <- res[i, ] + sapply(1:nrow(res), countns, s2 = i) * nscore
  }
  
  #add 1 to avoid division by 0, and make distances proportional
  return((res + 1) / max(res + 1))
}



