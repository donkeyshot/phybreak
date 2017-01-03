### phybreak class constructor ###


#' Create a phybreak-object from data and prior distributions.
#' 
#' phybreak takes as data either an \code{'obkData'}-object or a \code{\link{phybreakdata}}-object with sequences 
#' (individuals in rows, nucleotides in columns). Both \code{'obkData'} and \code{\link{phybreakdata}} 
#' contain at least sequences and sampling times, and potentially more. Parameter values are used as initial values in the MCMC-chain 
#' or kept fixed. All variables are initialized by random samples from the prior distribution, 
#' unless a complete tree is given in the data and should be used (\code{use.tree = TRUE}).
#' It is also possible to provide only sequences as data, and sampling times separately.
#' 
#' @param data An object with sequences plus additional data.
#'  (class \code{'obkData'} or \code{'phybreakdata'}). All nucleotides that are not \code{'a'}, \code{'c'}, \code{'g'}, 
#'  or \code{'t'}, will be turned into \code{'n'}.  
#'  
#'  If the data are provided as an object of class \code{'obkData'}, these should contain sequences and sampling times 
#'  as metadata with these sequences. The object may also contain \code{infector} and infection \code{date} vectors in the
#'  \code{individuals} slot, plus (at least) one tree in the \code{'trees'} slot (class \code{'multiPhylo'}).
#'  
#'  Data provided as an object of class \code{'phybreakdata'} contain \code{sequences} and \code{sampling.times},
#'  and potentially \code{infection.times}, \code{infectors}, and \code{trees}. Prepare your data in this format
#'  by \code{\link{phybreakdata}} or by simulation with \code{\link{sim.phybreak}}.
#'  
#'  It is also possible to provide only sequences as data, (class \code{'DNAbin'}, \code{'phyDat'}, or a \code{matrix} with nucleotides, 
#'  each row a host, each column a nucleotide), and corresponding sampling times in the separate \code{times} argument.
#' @param times Vector of sampling times, needed if the data consist of only sequences. If the vector is named,
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
#' simulation <- sim.phybreak(obsize = 10)
#' MCMCstate <- phybreak(data = simulation)
#' 
#' simulation <- sim.phybreak(obsize = 10)
#' MCMCstate <- phybreak(data = simulation, use.tree = TRUE)
#' 
#' 
#' sampletimedata <- c(0,2,2,4,4)
#' sampleSNPdata <- matrix(c("a","a","a","a","a",
#'                           "a","c","c","c","c",
#'                           "t","t","t","g","g"), nrow = 5)
#' dataset <- phybreakdata(sequences = sampleSNPdata, sample.times = sampletimedata)
#' MCMCstate <- phybreak(data = dataset)
#' 
#' ### also possible without 'phybreakdata' as intermediate, 
#' ### but not with additional data (future implementation)
#' MCMCstate <- phybreak(data = sampleSNPdata, times = sampletimedata)
#' @export
phybreak <- function(data, times = NULL,
         mu = NULL, gen.shape = 3, gen.mean = 1,
         sample.shape = 3, sample.mean = 1, 
         wh.model = 3, wh.slope = 1,
         est.gen.mean = TRUE, prior.mean.gen.mean = 1, prior.mean.gen.sd = Inf,
         est.sample.mean = TRUE, prior.mean.sample.mean = 1, prior.mean.sample.sd = Inf,
         est.wh.slope = TRUE, prior.wh.shape = 2, prior.wh.mean = 1,
         use.tree = FALSE) {
  
  ###########################################
  ### check for correct classes and sizes ###
  ###########################################
  data <- testdataclass_phybreak(data, times)
  if(use.tree) testfortree_phybreak(data)
  testargumentsclass_phybreak(environment())
  
  ### outbreak parameters ###
  
  ########################
  ### first slot: data ###
  ########################
  dataslot <- list()
  
  #names
  if(inherits(data, "obkData")) {
    dataslot$names <- OutbreakTools::get.individuals(data)
  } else {
    dataslot$names <- names(data$sample.times)
  }
  
  #sequences (SNP)
  SNP.sample <- c()
  if(inherits(data, "obkData")) {
    allgenes <- OutbreakTools::get.dna(data)
    for(i in 1:(OutbreakTools::get.nlocus(data))) {
      SNP.sample <- cbind(SNP.sample, as.character(allgenes[[i]]))
    }
    if(length(setdiff(SNP.sample, c("a","c","g","t")))) {
      warning("all nucleotides other than actg are turned into n")
    }
    SNP.sample[SNP.sample != "a" & SNP.sample != "c" & SNP.sample != "g" & SNP.sample != "t"] <- "n"
    dataslot$sequences <- phangorn::as.phyDat(SNP.sample)
  } else {
    dataslot$sequences <- data$sequences
  }
  
  #sample times
  if(inherits(data, "obkData")) {
    dataslot$sample.times <- OutbreakTools::get.dates(data, "dna")
  } else {
    dataslot$sample.times <- data$sample.times
  }
  names(dataslot$sample.times) <- dataslot$names
  
  #SNP count
  seqmat <- as.character(dataslot$sequences)
  dataslot$nSNPs <- as.integer(
    sum(apply(seqmat=="a", 2, any) +
          apply(seqmat=="c", 2, any) +
          apply(seqmat=="g", 2, any) +
          apply(seqmat=="t", 2, any) +
          apply(seqmat=="n", 2, all) - 1)
  )
  
  ##############################
  ### third slot: parameters ###
  ##############################
  parameterslot <- list(
    obs = length(dataslot$sample.times),
    mu = NULL,
    mean.sample = sample.mean,
    mean.gen = gen.mean,
    shape.sample = sample.shape,
    shape.gen = gen.shape,
    wh.model = wh.model,
    wh.slope = wh.slope
  )
  
  ##############################
  ### second slot: variables ###
  ##############################
  if(inherits(data, "phybreakdata")) {
    phybreakvariables <- transphylo2phybreak(data, resample = !use.tree, resamplepars = parameterslot)
  } else {
    phybreakvariables <- obkData2phybreak(data, resample = !use.tree, resamplepars = parameterslot)
  }
  variableslot <- phybreakvariables$v
  dataslot$reference.date <- phybreakvariables$d$reference.date
  
  #################
  # parameters$mu #
  #################
  if(is.null(mu)) {
    treelength <- with(variableslot, sum(nodetimes[nodeparents != 0] - nodetimes[nodeparents]))
    parameterslot$mu <- (dataslot$nSNPs / ncol(seqmat)) / treelength
  } else {
    parameterslot$mu <- mu
  }
  
  #################################
  ### fourth slot: helper input ###
  #################################
  helperslot <- list(si.mu = if(dataslot$nSNPs == 0) 0 else 2.38*sqrt(trigamma(dataslot$nSNPs)),
                     si.wh = 2.38*sqrt(trigamma(parameterslot$obs - 1)),
                     dist = distmatrix_phybreak(seqmat),
                     est.mG = est.gen.mean,
                     est.mS = est.sample.mean,
                     est.wh = est.wh.slope,
                     mG.av = prior.mean.gen.mean,
                     mG.sd = prior.mean.gen.sd,
                     mS.av = prior.mean.sample.mean,
                     mS.sd = prior.mean.sample.sd,
                     wh.sh = prior.wh.shape,
                     wh.av = prior.wh.mean)
  
  ###########################
  ### fifth slot: samples ###
  ###########################
  sampleslot <- list(
    nodetimes = c(),
    nodehosts = c(),
    nodeparents = c(),
    mu = c(),
    mG = c(),
    mS = c(),
    slope = c(),
    logLik = c()
  )

  ################################
  ### make the phybreak object ###
  ################################
  res <- list(
    d = dataslot,
    v = variableslot,
    p = parameterslot,
    h = helperslot,
    s = sampleslot
    )
  
  class(res) <- c("phybreak", "list")
  
  return(res)
}


### Test data class
testdataclass_phybreak <- function(data, times) {
  if(inherits(data, c("DNAbin", "phyDat", "matrix"))) {
    data <- phybreakdata(sequences = data, sample.times = times)
  }
  
  if(!inherits(data, c("obkData", "phybreakdata"))) {
    stop("data should be of class \"obkData\" or \"phybreakdata\"")
  }
  
  if(inherits(data, "obkData")) {
    if(!("OutbreakTools" %in% .packages(TRUE))) {
      stop("package 'OutbreakTools' not installed, while data-class is \"obkData\"")
    }
    if(!("OutbreakTools" %in% .packages(FALSE))) {
      warning("package 'OutbreakTools' is not attached")
    }
  }
  
  return(data)
}

### Test for presence of tree
testfortree_phybreak <- function(data) {
  if(inherits(data, "phybreakdata")) {
    if(is.null(data$sim.infection.times) | is.null(data$sim.infectors)) {
      warning("transmission tree can only be used if provided in data; random tree will be generated")
    }
    if(is.null(data$sim.tree)) {
      warning("phylogenetic tree can only be used if provided in data; random tree will be generated")
    }
  } else {
    if(is.null(OutbreakTools::get.dates(data, "individuals")) | is.null(OutbreakTools::get.data(data, "infector"))) {
      warning("transmission tree can only be used if provided in data; random tree will be generated")
    }
    if(is.null(OutbreakTools::get.trees(data))) {
      warning("phylogenetic tree can only be used if provided in data; random tree will be generated")
    }
  }
}

### Test arguments classes
testargumentsclass_phybreak <- function(env) {
  with(env, {
    if(is.null(mu)) mutest <- 1 else mutest <- mu
    numFALSE <- 
      unlist(
        lapply(
          list(mutest, gen.shape, gen.mean, sample.shape, sample.mean,
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
      c(mutest, gen.shape, gen.mean, sample.shape, sample.mean,
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
  })
}



### pseudo-distance matrix between sequences given SNP data
distmatrix_phybreak <- function(seqmatrix) {

  # count SNPs excluding "n"
  res <- as.matrix(ape::dist.dna(ape::as.DNAbin(seqmatrix), model = "N", pairwise.deletion = TRUE))

  # prob of SNP per nucleotide in most distant entry
  nscore <- max(res)/ncol(seqmatrix)

  # add nscore for each missing nucleotide
  nmatrix <- seqmatrix == "n"
  for(i in 1:nrow(res)) {
    res[i, ] <- res[i, ] + sapply(1:nrow(res), 
                                  function(s1, s2) sum(nmatrix[s1,] | nmatrix[s2,]),
                                  s2 = i) * nscore
  }
  
  #add 1 to avoid division by 0, and make distances proportional
  return((res + 1) / max(res + 1))
}


