#' Log-likelihood of a phybreak-object.
#' 
#' The likelihood of a \code{phybreak}-object is calculated, with the option to include or exclude parts of the 
#' likelihood for genetic data, phylogenetic tree (within-host model), sampling times and generation times.
#' 
#' The sequence likelihood is calculated by Felsenstein's pruning algorithm, assuming a prior probability of 0.25 
#' for each nucleotide. The within-host likelihood is the likelihood of coalescence times given the within-host model 
#' and slope. The generation interval and sampling interval likelihood are log-densities of the gamma distributions 
#' for these variables.
#' 
#' @param object An object of class \code{phybreak}.
#' @param genetic Whether to include the likelihood of the mutation model.
#' @param withinhost Whether to include the likelihood of within-host (coalescent) model.
#' @param sampling Whether to include the likelihood of the sampling model (sampling intervals).
#' @param generation Whether to include the likelihood of the transmission model (generation intervals).
#' @param ... Some methods for this generic require additional arguments. None are used in this method.
#' @return The log-likelihood as an object of class logLik.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' logLik(MCMCstate)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' logLik(MCMCstate)
#' 
#' tree0 <- get_phylo(MCMCstate)
#' seqdata <- get_data(MCMCstate)$sequences
#' phangorn::pml(tree0, seqdata, rate = 0.75*get_parameters(MCMCstate, "mu") 
#' logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
#'        sampling = FALSE, generation = FALSE) #should give the same result as 'pml'
#' @export
logLik.phybreak <- function(object, genetic = TRUE, withinhost = TRUE, sampling = TRUE, generation = TRUE, ...) {
  res <- 0
  if (genetic) {
    res <- res + with(object, .likseq(matrix(unlist(d$sequences), ncol = d$nsamples), 
                                      attr(d$sequences, "weight"), 
                                      v$nodeparents, v$nodetimes, p$mu, d$nsamples))
  }
  if (generation) {
    res <- res + with(object, lik_gentimes(p$gen.shape, p$gen.mean, v$inftimes, v$infectors))
  }
  if (sampling) {
    res <- res + with(object, lik_sampletimes(p$obs, p$sample.shape, p$sample.mean, v$nodetimes, v$inftimes))
  }
  if (withinhost) {
    objectenv <- object
    objectenv$v <- phybreak2environment(objectenv$v)
    res <- res + with(object, lik_coaltimes(objectenv))
  }
  attributes(res) <- list(
    nobs = object$p$obs,
    df = 1 + object$h$est.mG + object$h$est.mS + object$h$est.wh.s + object$h$est.wh.e + object$h$est.wh.0,
    genetic = genetic, withinhost = withinhost, sampling = sampling, generation = generation
  )
  class(res) <- "logLik"
  return(res)
}


### calculate the log-likelihood of sampling intervals 
lik_gentimes <- function(shapeG, meanG, inftimes, infectors) {
  sum(dgamma(inftimes[infectors > 0] - 
               inftimes[infectors[infectors > 0]], 
             shape = shapeG, scale = meanG/shapeG, log = TRUE))
}

### calculate the log-likelihood of generation intervals 
lik_sampletimes <- function(obs, shapeS, meanS, nodetimes, inftimes) {
  sum(dgamma(nodetimes[1:obs] - inftimes, shape = shapeS, scale = meanS/shapeS, log = TRUE))
}


### calculate the log-likelihood of coalescent intervals 
lik_coaltimes <- function(phybreakenv) {
  return(sum(sapply(0:phybreakenv$p$obs, lik_coaltimes_host, phybreakenv = phybreakenv)))
}

### calculate the log-likelihood of coalescent intervals 
lik_coaltimes_host <- function(phybreakenv, hostID) {
  if (phybreakenv$p$wh.model %in% c(1, 2, "single", "infinite")) return(0)
  
  if(hostID > 0) {
    bnodes_from_infector <- which(phybreakenv$v$nodehosts == phybreakenv$v$infectors[hostID] & phybreakenv$v$nodetypes == "b")
    bnodes_in_infectees <- which(phybreakenv$v$nodeparents %in% bnodes_from_infector)
    bottleneck <- 1 + length(phybreakenv$v$nodehosts[bnodes_in_infectees] == hostID)
    starttime <- phybreakenv$v$inftimes[hostID]
    inftime <- 0
  } else {
    bottleneck <- 1
    starttime <- min(phybreakenv$v$inftimes) - phybreakenv$p$sample.mean
    inftime <- min(phybreakenv$v$nodetimes) - 1 - starttime
  }
  
  ### minustimes
  minustimes <- phybreakenv$v$nodetimes[phybreakenv$v$nodehosts == hostID & phybreakenv$v$nodetypes %in% c("s", "x", "t", "b")] - starttime
  
  plustimes <- phybreakenv$v$nodetimes[phybreakenv$v$nodehosts == hostID & phybreakenv$v$nodetypes == "c"] - starttime
  
  eventslist <- c(bottleneck, rep(-1, length(minustimes)), rep(1, length(plustimes)))
  timeslist <- c(inftime, minustimes, plustimes)
  eventslist <- eventslist[order(timeslist)]
  timeslist <- sort(timeslist)
  
  lineagecounts <- cumsum(eventslist)
  coalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                      linear = 1/(phybreakenv$p$wh.slope * plustimes),
                      exponential = 1/(phybreakenv$p$wh.level * exp(phybreakenv$p$wh.exponent * plustimes)),
                      constant = rep(1, length(plustimes))/phybreakenv$p$wh.level)
  logcoalrates <- log(coalrates)
  cumcoalrates <- switch(phybreakenv$p$wh.model, single =, infinite =,
                         linear = c(0, log(timeslist[-1]) / phybreakenv$p$wh.slope),
                         exponential = -1/(phybreakenv$p$wh.level * phybreakenv$p$wh.exponent * 
                                             exp(phybreakenv$p$wh.exponent * timeslist)),
                         constant = timeslist/phybreakenv$p$wh.level)
  logcoalescapes <- -(tail(cumcoalrates, -1) - head(cumcoalrates, -1)) * head(choose(lineagecounts, 2), -1)
  
  return(sum(logcoalrates) + sum(logcoalescapes))
  
}




