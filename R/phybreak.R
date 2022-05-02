#' Create a phybreak-object from data and prior distributions.
#' 
#' phybreak takes as data a \code{\link{phybreakdata}}-object with sequences (individuals in rows, 
#' nucleotides in columns) and sampling times, and potentially more. Parameter values are used as 
#' initial values in the MCMC-chain or kept fixed. All variables are initialized by random samples 
#' from the prior distribution, unless a complete tree is given in the data and should be used 
#' (\code{use.tree = TRUE}). It is also possible to provide only sequences as data, and sampling times separately.
#' 
#' @param dataset An object with sequences plus additional data (class \code{'phybreakdata'}). 
#'  Data contain \code{sequences} and \code{sampling.times}, and potentially \code{sim.infection.times}, 
#'  \code{sim.infectors}, and \code{sim.tree}. Prepare your data in this format
#'  by \code{\link{phybreakdata}} or by simulation with \code{\link{sim.phybreak}}.
#'  
#'  It is also possible to provide only sequences as data, (class \code{'DNAbin'}, \code{'phyDat'}, or a \code{matrix} with nucleotides, 
#'  each row a host, each column a nucleotide), and corresponding sampling times in the separate \code{times} argument.
#' @param times Vector of sampling times, needed if the data consist of only sequences. If the vector is named,
#'   these names will be used to identify the hosts.
#' @param mu Initial value for mutation rate (defined per site per unit of time). If \code{NULL} (default), then an initial
#'   value is calculated by dividing the number of SNPs by the product: 0.75 times 'total sequence length' times 'sum of
#'   edge lengths in the initial phylogenetic tree'.
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
#'   N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Names and numbers are allowed.
#'   Options are:
#'   \enumerate{
#'     \item "single": effective size = 0, so coalescence occurs 'just before' transmission in the infector (complete bottleneck)
#'     \item "infinite": effective size = Inf, with complete bottleneck, so coalescence occurs 'just after' transmission in the infectee
#'     \item "linear": effective size at time t after infection = \code{wh.level + wh.slope * t} (complete or wide bottleneck; if complete, \code{wh.level = 0})
#'     \item "exponential": effective size at time t after infection = \code{wh.level * exp(wh.exponent * t)} (wide bottleneck)
#'     \item "constant": effective size = wh.level (wide bottleneck)
#'   }
#' @param wh.bottleneck Whether the bottleneck should be complete or wide, which is only an option if \code{wh.model = "linear"} 
#'   (in that case, \code{"auto"} defaults to \code{"complete"}).
#' @param wh.slope Initial value for the within-host slope, used if \code{wh.model = "linear"}.
#' @param wh.exponent Initial value for the within-host exponent, used if \code{wh.model = "exponential"}
#' @param wh.level Initial value for the within-host effective pathogen size at transmission, used if \code{wh.bottleneck = "wide"}
#'   (if \code{wh.model = "exponential"} or \code{"constant"}, and optional if \code{wh.model = "linear"})
#' @param est.gen.mean Whether to estimate the mean generation interval or keep it fixed. 
#' @param prior.gen.mean.mean Mean of the (gamma) prior distribution of mean generation interval \code{mG} 
#'   (only if \code{est.gen.mean = TRUE}).
#' @param prior.gen.mean.sd Standard deviation of the (gamma) prior distribution of mean generation interval \code{mG} 
#'   (only if \code{est.gen.mean = TRUE}).
#' @param est.sample.mean Whether to estimate the mean sampling interval or keep it fixed. 
#' @param prior.sample.mean.mean Mean of the (gamma) prior distribution of mean sampling interval \code{mS} 
#'   (only if \code{est.sample.mean = TRUE}).
#' @param prior.sample.mean.sd Standard deviation of the (gamma) prior distribution of mean sampling interval \code{mS} 
#'   (only if \code{est.sample.mean = TRUE}).
#' @param est.wh.slope Whether to estimate the within-host slope or keep it fixed. 
#' @param prior.wh.slope.shape Shape parameter of the (gamma) prior distribution of \code{wh.slope} 
#'   (only if \code{est.wh.slope = TRUE}).
#' @param prior.wh.slope.mean Mean of the (gamma) prior distribution of \code{wh.slope} 
#'   (only if \code{est.wh.slope = TRUE}).
#' @param est.wh.exponent Whether to estimate the within-host exponent or keep it fixed. 
#' @param prior.wh.exponent.shape Shape parameter of the (gamma) prior distribution of \code{wh.exponent} 
#'   (only if \code{est.wh.exponent = TRUE}).
#' @param prior.wh.exponent.mean Mean of the (gamma) prior distribution of \code{wh.exponent} 
#'   (only if \code{est.wh.exponent = TRUE}).
#' @param est.wh.level Whether to estimate the within-host level at \code{t = 0} or keep it fixed. 
#' @param prior.wh.level.shape Shape parameter of the (gamma) prior distribution of \code{wh.level} 
#'   (only if \code{est.wh.level = TRUE}).
#' @param prior.wh.level.mean Mean of the (gamma) prior distribution of \code{wh.level} 
#'   (only if \code{est.wh.level = TRUE}).
#' @param use.tree Whether to use the transmission and phylogenetic tree given in data of class \code{'obkData'}, 
#'   to create a \code{phybreak}-object with an exact copy of the outbreak. This requires more data in \code{data}: 
#'   the slot \code{individuals} with vectors \code{infector} and \code{date}, and the slot \code{trees} with at least 
#'   one phylogenetic tree. Such data can be simulated with \code{\link{sim.phybreak}}.
#' @param ... If arguments from previous versions of this function are used, they may be interpreted correctly through 
#'   this argument, but it is better to provide the correct argument names.
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
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' simulation <- sim_phybreak(obsize = 10)
#' MCMCstate <- phybreak(dataset = simulation)
#' 
#' simulation <- sim_phybreak(obsize = 10)
#' MCMCstate <- phybreak(dataset = simulation, use.tree = TRUE)
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
#' ### but only with single samples per host and not with additional data (future implementation)
#' MCMCstate <- phybreak(data = sampleSNPdata, times = sampletimedata)
#' @export
phybreak <- function(dataset, times = NULL, 
         mu = NULL, gen.shape = 3, gen.mean = 1, trans.model = "gamma",
         trans.init = 1e-4, trans.growth = 1, trans.sample = 0.1, trans.culling = 5,
         sample.shape = 3, sample.mean = 1, 
         multiple.introductions = FALSE, introductions = 1, intro.rate = 1,
         wh.model = "linear", wh.bottleneck = "auto", wh.history = 1, wh.slope = 1, wh.exponent = 1, wh.level = 0.1,
         dist.model = "power", dist.exponent = 2, dist.scale = 1, dist.mean = 1,
         est.mu = TRUE, prior.mu.mean = 0, prior.mu.sd = 100,
         est.gen.mean = TRUE, prior.gen.mean.mean = 1, prior.gen.mean.sd = Inf,
         est.sample.mean = TRUE, prior.sample.mean.mean = 1, prior.sample.mean.sd = Inf,
         est.intro.rate = TRUE, prior.intro.rate.mean = 1, prior.intro.rate.sd = Inf,
         est.trans.growth = TRUE, est.trans.sample = TRUE, 
         est.wh.slope = TRUE, prior.wh.slope.shape = 3, prior.wh.slope.mean = 1,
         est.wh.exponent = TRUE, prior.wh.exponent.shape = 1, prior.wh.exponent.mean = 1,
         est.wh.level = TRUE, prior.wh.level.shape = 1, prior.wh.level.mean = 0.1,
         est.wh.history = TRUE, prior.wh.history.shape = 1, prior.wh.history.mean = 100,
         est.dist.exponent = TRUE, prior.dist.exponent.shape = 1, prior.dist.exponent.mean = 1,
         est.dist.scale = TRUE, prior.dist.scale.shape = 1, prior.dist.scale.mean = 1,
         est.dist.mean = TRUE, prior.dist.mean.shape = 1, prior.dist.mean.mean = 1,
         use.tree = FALSE, ...) {
  ########################################################
  ### parameter name compatibility with older versions ###
  ########################################################
  old_arguments <- list(...)
  if(exists("oldarguments$prior.mean.gen.mean")) prior.gen.mean.mean <- oldarguments$prior.mean.gen.mean
  if(exists("oldarguments$prior.mean.gen.sd")) prior.gen.mean.sd <- oldarguments$prior.mean.gen.sd
  if(exists("oldarguments$prior.mean.sample.mean")) prior.sample.mean.mean <- oldarguments$prior.mean.sample.mean
  if(exists("oldarguments$prior.mean.sample.sd")) prior.sample.mean.sd <- oldarguments$prior.mean.sample.sd
  if(exists("oldarguments$prior.wh.shape")) prior.wh.slope.shape <- oldarguments$prior.wh.shape
  if(exists("oldarguments$prior.wh.mean")) prior.wh.slope.mean <- oldarguments$prior.wh.mean
  
  
  ###########################################
  ### check for correct classes and sizes ###
  ###########################################
  dataset <- testdataclass_phybreak(dataset, times, ...)
  rownames(as.character(dataset$sequences))
  if(use.tree) introductions <- testfortree_phybreak(dataset, introductions)
  if(introductions > 1) multiple.introductions <- TRUE
  if(!multiple.introductions) {
    intro.rate <- 1
    est.intro.rate <- FALSE
    }
  testargumentsclass_phybreak(environment())
  wh.model <- choose_whmodel(wh.model)
  wh.bottleneck <- choose_whbottleneck(wh.bottleneck, wh.model)
  dist.model <- choose_distmodel(dist.model, dataset$distances)

  ### outbreak parameters ###
  
  ########################
  ### first slot: data ###
  ########################
  dataslot <- list()
  
  #names, sequences, sample times, distances 
  dataslot$names <- names(dataset$sample.times)
  dataslot$hostnames <- dataset$sample.hosts
  dataslot$sequences <- dataset$sequences
  dataslot$sample.times <- dataset$sample.times
  dataslot$culling.times <- dataset$culling.times
  dataslot$locations <- dataset$locations
  dataslot$distances <- dataset$distances
  
  #SNP count
  SNPpatterns <- do.call(rbind, dataslot$sequences)
  dataslot$nSNPs <- as.integer(
    sum(apply(SNPpatterns, 2, 
              function(x) {
                max(0, length(unique(x[x < 5])) - 1)
              }
              ) * attr(dataslot$sequences, "weight")
        )
    )

  #Sample size
  dataslot$nsamples <- length(dataslot$names)
  
  ##############################
  ### third slot: parameters ###
  ##############################
  parameterslot <- list(
    obs = length(unique(dataslot$hostnames)),
    mu = NULL,
    sample.mean = sample.mean,
    gen.mean = gen.mean,
    sample.shape = sample.shape,
    gen.shape = gen.shape,
    trans.init = trans.init,
    trans.growth = trans.growth,
    trans.sample = trans.sample,
    trans.culling = trans.culling,
    intro.rate = intro.rate,
    wh.model = wh.model,
    wh.bottleneck = wh.bottleneck,
    wh.slope = wh.slope,
    wh.exponent = wh.exponent,
    wh.level = wh.level * (wh.bottleneck == "wide"),
    wh.history = wh.history,
    dist.model = dist.model,
    dist.exponent = dist.exponent,
    dist.scale = dist.scale,
    dist.mean = dist.mean,
    trans.model = trans.model,
    mult.intro = multiple.introductions
  )
  
  ##############################
  ### second slot: variables ###
  ##############################
  phybreakvariables <- transphylo2phybreak(dataset, resample = !use.tree, resamplepars = parameterslot,
                                           introductions = introductions)
  variableslot <- phybreakvariables$v
  dataslot$reference.date <- phybreakvariables$d$reference.date
  
  #################
  # parameters$mu #
  #################
  if(is.null(mu)) {
    treelength <- with(variableslot, sum(nodetimes[nodeparents != 0] - nodetimes[nodeparents]))
    curparsimony <- phangorn::parsimony(phybreak2phylo(variableslot), dataslot$sequences)
    sequencelength <- sum(attr(dataslot$sequences, "weight"))
    parameterslot$mu <- (curparsimony / sequencelength) / treelength / 0.75
  } else {
    parameterslot$mu <- mu
  }
  
  #################################
  ### fourth slot: helper input ###
  #################################
  helperslot <- list(si.mu = if(dataslot$nSNPs == 0) 0 else 2.38*sqrt(trigamma(dataslot$nSNPs)),
                     si.ir = 2.38*sqrt(trigamma(introductions)),
                     si.wh = 2.38*sqrt(trigamma(dataslot$nsamples - 1)),
                     si.dist = 2.38*sqrt(trigamma(parameterslot$obs - 1)),
                     dist = distmatrix_phybreak(subset(dataslot$sequences, subset = 1:parameterslot$obs)),
                     est.mu = est.mu,
                     est.mG = est.gen.mean,
                     est.mS = est.sample.mean,
                     est.tG = est.trans.growth,
                     est.tS = est.trans.sample,
                     est.ir = est.intro.rate,
                     est.wh.s = est.wh.slope && wh.model == "linear",
                     est.wh.e = est.wh.exponent && wh.model == "exponential",
                     est.wh.0 = est.wh.level && wh.bottleneck == "wide",
                     est.wh.h = est.wh.history && multiple.introductions,
                     est.dist.e = est.dist.exponent && dist.model %in% c("power", "exponential"),
                     est.dist.s = est.dist.scale && dist.model == "power",
                     est.dist.m = est.dist.mean && dist.model == "poisson",
                     mu.av = prior.mu.mean,
                     mu.sd = prior.mu.sd,
                     mG.av = prior.gen.mean.mean,
                     mG.sd = prior.gen.mean.sd,
                     mS.av = prior.sample.mean.mean,
                     ir.av = prior.intro.rate.mean,
                     ir.sd = prior.intro.rate.sd,
                     mS.sd = prior.sample.mean.sd,
                     wh.s.sh = prior.wh.slope.shape,
                     wh.s.av = prior.wh.slope.mean,
                     wh.e.sh = prior.wh.exponent.shape,
                     wh.e.av = prior.wh.exponent.mean,
                     wh.0.sh = prior.wh.level.shape,
                     wh.0.av = prior.wh.level.mean,
                     wh.h.sh = prior.wh.history.shape,
                     wh.h.av = prior.wh.history.mean,
                     dist.e.sh = prior.dist.exponent.shape,
                     dist.e.av = prior.dist.exponent.mean,
                     dist.s.sh = prior.dist.scale.shape,
                     dist.s.av = prior.dist.scale.mean,
                     dist.m.sh = prior.dist.mean.shape,
                     dist.m.av = prior.dist.mean.mean)
  
  ###########################
  ### fifth slot: samples ###
  ###########################
  sampleslot <- list(
    inftimes = c(),
    infectors = c(),
    nodetimes = c(),
    nodehosts = c(),
    nodeparents = c(),
    introductions = c(),
    mu = c(),
    mG = c(),
    mS = c(),
    tG = c(),
    tS = c(),
    ir = c(),
    wh.s = c(),
    wh.e = c(),
    wh.0 = c(),
    wh.h = c(),
    dist.e = c(),
    dist.s = c(),
    dist.m = c(),
    logLik = c(),
    heat = c()
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


### Test dataset class
testdataclass_phybreak <- function(dataset, times, ...) {
  if(inherits(dataset, c("DNAbin", "phyDat", "matrix"))) {
    dataset <- phybreakdata(sequences = dataset, sample.times = times, ...)
  }
  
  if(!inherits(dataset, "phybreakdata")) {
    stop("dataset should be of class \"phybreakdata\"")
  }
  
  return(dataset)
}

### Test for presence of tree
testfortree_phybreak <- function(dataset, introductions) {
  if(is.null(dataset$sim.infection.times) | is.null(dataset$sim.infectors)) {
    warning("transmission tree can only be used if provided in dataset; random tree will be generated")
  } else {
    simulatedintroductions <- sum(dataset$sim.infectors == "index") * 1
    if(introductions != simulatedintroductions) {
      warning(paste0("Provided transmission tree contains ", simulatedintroductions, " introduction(s) instead of ",
                     introductions,". Transmission tree will be used"))
      introductions <- simulatedintroductions
    }
  }
  if(is.null(dataset$sim.tree)) {
    warning("phylogenetic tree can only be used if provided in dataset; random tree will be generated")
  }
  return(introductions)
}

### Test arguments classes
testargumentsclass_phybreak <- function(env) {
  with(env, {
    if(is.null(mu)) mutest <- 1 else mutest <- mu
    numFALSE <- 
      unlist(
        lapply(
          list(mutest, gen.shape, gen.mean, sample.shape, sample.mean,
               introductions, intro.rate,
               wh.slope, wh.exponent, wh.level, prior.gen.mean.mean, prior.gen.mean.sd,
               prior.sample.mean.mean, prior.sample.mean.sd,
               prior.wh.slope.shape, prior.wh.slope.mean,
               prior.wh.exponent.shape, prior.wh.exponent.mean,
               prior.wh.level.shape, prior.wh.level.mean),
          class
        )
      ) != "numeric"
    if(any(numFALSE)) {
      stop(paste0("parameters ",
                  c("mu", "gen.shape", "gen.mean", "sample.shape", "sample.mean",
                    "introductions", "intro.rate",
                    "wh.slope", "wh.exponent", "wh.level", "prior.gen.mean.mean", "prior.gen.mean.sd",
                    "prior.shape.mean.mean", "prior.shape.mean.sd",
                    "prior.wh.slope.shape", "prior.wh.slope.mean",
                    "prior.wh.exponent.shape", "prior.wh.exponent.mean",
                    "prior.wh.level.shape", "prior.wh.level.mean")[numFALSE],
                  " should be numeric"))
    }
    numNEGATIVE <- 
      c(mutest, gen.shape, gen.mean, sample.shape, sample.mean,
        introductions, intro.rate,
        wh.slope, wh.exponent, wh.level, prior.gen.mean.mean, prior.gen.mean.sd,
        prior.sample.mean.mean, prior.sample.mean.sd,
        prior.wh.slope.shape, prior.wh.slope.mean,
        prior.wh.exponent.shape, prior.wh.exponent.mean,
        prior.wh.level.shape, prior.wh.level.mean
  ) <= 0
    if(any(numNEGATIVE)) {
      stop(paste0("parameters ",
                  c("mu", "gen.shape", "gen.mean", "sample.shape", "sample.mean",
                    "introductions", "intro.rate",
                    "wh.slope", "wh.exponent", "wh.level", "prior.gen.mean.mean", "prior.gen.mean.sd",
                    "prior.shape.mean.mean", "prior.shape.mean.sd",
                    "prior.wh.slope.shape", "prior.wh.slope.mean",
                    "prior.wh.exponent.shape", "prior.wh.exponent.mean",
                    "prior.wh.level.shape", "prior.wh.level.mean")[numNEGATIVE],
                  " should be positive"))
    }
    logFALSE <- 
      unlist(
        lapply(
          list(multiple.introductions, est.gen.mean, est.sample.mean,
               est.wh.slope, est.wh.exponent, est.wh.level, use.tree),
          class
        )
      ) != "logical"
    if(any(logFALSE)) {
      stop(paste0("parameters ",
                  c("multiple.introductions", "est.gen.mean", "est.sample.mean", "est.wh.slope", "est.wh.exponent", 
                    "est.wh.level", "use.tree")[logFALSE],
                  " should be logical"))
    }
  })
}

### set within-host model from possible input values
choose_whmodel <- function(x) {
  whoptions <- c("single", "infinite", "linear", "exponential", "constant")
  if(is.numeric(x)) {
    if(floor(x) %in% 1:5) {
      return(whoptions[x])
    } else {
      stop("pick one of five within-host models")
    }
  } else {
    return(match.arg(x, whoptions))
  }
}

### set within-host bottleneck from possible input values
choose_whbottleneck <- function(x, wh.model) {
  x <- match.arg(x, c("auto", "complete", "wide"))
  if(wh.model %in% c("single", "infinite")) {
    if(x == "wide") message(paste0("wh.model = ", wh.model, " only possible with complete bottleneck"))
    return("complete")
  } else if(wh.model %in% c("exponential", "constant")) {
    if(x == "complete") message(paste0("wh.model = ", wh.model, " only possible with wide bottleneck"))
    return("wide")
  } else {
    if(x == "auto") {
      return("complete")
    } else return(x)
  }
}

### set distance model from possible input values
choose_distmodel <- function(x, distances) {
  x <- match.arg(x, c("power", "exponential", "poisson", "none"))
  if(is.null(distances)) {
    x <- "none"
  }
  return(x)
}


### pseudo-distance matrix between sequences given SNP data
distmatrix_phybreak <- function(sequences) {
  
  # count SNPs excluding "n"
  res <- as.matrix(phangorn::dist.hamming(sequences, exclude = "pairwise", ratio = FALSE))
  
  # prob of SNP per nucleotide in most distant entry
  nscore <- max(res)/sum(attr(sequences, "weight"))
  
  # add nscore for each missing nucleotide
  seqs_n <- do.call(rbind, sequences) == 16
  res <- res + outer(X = 1:length(sequences), Y = 1:length(sequences), 
                     FUN = Vectorize(
                       function(x, y) sum((seqs_n[x,] | seqs_n[y,]) * attr(sequences, "weight"))
                     )) * nscore
  

  #add 1 to avoid division by 0, and make distances proportional
  res <- (res + 1) / max(res + 1)
  return(res / min(res))
}

