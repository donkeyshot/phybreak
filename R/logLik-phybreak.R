### calculate the log-likelihood ###



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
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' logLik(MCMCstate)
#' 
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' logLik(MCMCstate)
#' 
#' tree0 <- get.phylo(MCMCstate)
#' seqdata <- get.seqdata(MCMCstate)
#' pml(tree0, seqdata, rate = 0.75*get.parameters(MCMCstate)["mu"]) 
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
        res <- res + with(object, .lik.gentimes(p$gen.shape, p$gen.mean, v$inftimes, v$infectors))
    }
    if (sampling) {
        res <- res + with(object, .lik.sampletimes(p$obs, p$sample.shape, p$sample.mean, v$nodetimes, v$inftimes))
    }
    if (withinhost) {
        res <- res + with(object, .lik.coaltimes(object))
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
.lik.gentimes <- function(shapeG, meanG, inftimes, infectors) {
    sum(dgamma(inftimes[infectors > 0] - 
                 inftimes[infectors[infectors > 0]], 
               shape = shapeG, scale = meanG/shapeG, log = TRUE))
}

### calculate the log-likelihood of generation intervals 
.lik.sampletimes <- function(obs, shapeS, meanS, nodetimes, inftimes) {
    sum(dgamma(nodetimes[1:obs] - inftimes, shape = shapeS, scale = meanS/shapeS, log = TRUE))
}


### calculate the log-likelihood of coalescent intervals 
.lik.coaltimes <- function(phybreakenv) {
  return(sum(sapply(0:phybreakenv$p$obs, .lik.coaltimes.host, phybreakenv = phybreakenv)))
}

### calculate the log-likelihood of coalescent intervals 
.lik.coaltimes.host <- function(phybreakenv, hostID) {
  if (phybreakenv$p$wh.model %in% c(1, 2)) return(0)
  
  if(hostID > 0) {
    bottleneck <- length(firstnodes(phybreakenv, hostID))
    starttime <- phybreakenv$v$inftimes[hostID]
    inftime <- 0
  } else {
    bottleneck <- 1
    starttime <- min(phybreakenv$v$inftimes) - phybreakenv$p$sample.mean
    inftime <- min(phybreakenv$v$nodetimes) - 1 - starttime
  }
  
  ### minustimes
  minusnodes <- infectee_nodes(phybreakenv, hostID)
  minushosts <- unlist(sapply(phybreakenv$v$nodehosts[minusnodes], secondhostinchain, 
                       infectors = phybreakenv$v$infectors, firsthostID = hostID))
  minustimes <- c(phybreakenv$v$inftimes[minushosts], 
                  phybreakenv$v$nodetimes[phybreakenv$v$nodehosts == hostID & phybreakenv$v$nodetypes %in% c("s", "x")]) - starttime
  
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




