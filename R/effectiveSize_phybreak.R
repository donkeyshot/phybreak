#' Effective sample size
#' @name ESS
NULL

### Calculate ESS, using coda::effectiveSize for all continuous parameters,
### and effectiveSize_infector for all infectors
#' @describeIn ESS Effective sample size of phybreak posterior.
#' @param x An object of class \code{phybreak}.
#' @return Effective sample sizes.
#' @details For all parameters and continuous variables (infection times), \code{effectiveSize_phybreak} calculates
#'   the effective sample size (ESS) with the \code{\link[coda]{effectiveSize}} function in \pkg{coda}. For the infectors, 
#'   a method is used that is similar to the method for the approximate ESS for phylogenetic trees, described in 
#'   Lanfaer et al (2016):
#'   \enumerate{
#'     \item Define as distance measure between two sampled infectors \eqn{D(i,j) = 0} if \eqn{i = j}, 
#'         and \eqn{D(i,j) = 1} if \eqn{i \neq= j}{i <> j}
#'     \item Calculate the mean squared distance f(k) between sampled infectors at intervals k = 1,2,... in the mcmc chain.
#'         The distance will increase with increasing interval k.
#'     \item Use the rate at which f(k) approaches the asymptote to calculate the ESS (see Lanfaer et al, 2016)
#'   }
#'   This method can also be directly called with \code{effectiveSize_factor} for a single vector of categorical variables.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1093/gbe/evw171}{Lanfaer et al. (2016)} Estimating 
#'   the effective sample size of tree topologies from Bayesian phylogenetic analyses. 
#'   \emph{Genome Biol Evol}, \strong{8}(8): 2319-2332.
#' @examples 
#' #First create a phybreak object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' effectiveSize_phybreak(MCMCstate)
#' @export
effectiveSize_phybreak <- function(x) {
  ### tests
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  # Only with at least 2 posterior samples
  if(length(x$s$mu) < 2) stop("not possible with less than 2 posterior samples")
  
  # Get mcmc object and identify the columns with infectors
  mcmcobj <- suppressWarnings(get_mcmc(x))
  whichinfectors <- grepl("infector", colnames(mcmcobj))
  
  # Calculate the ESS per chain with the correct method
  res <- sapply(1:ncol(mcmcobj), function(i) {
    if(whichinfectors[i]) {
      effectiveSize_categorical(mcmcobj[, i])
    } else {
      coda::effectiveSize(mcmcobj[, i])
    }
  })
  
  # Add names and return
  names(res) <- colnames(mcmcobj)
  return(res)
}

#' @describeIn ESS Effective sample size of a categorical variable.
#' @param chain Vector of class \code{factor}, or a numeric vector with integers indicating categories.
#' @export
effectiveSize_factor <- function(chain) {
  # length must be 2 at least
  if(length(chain) < 2) return(NA_real_)  
  
  # maximum interval k to calculate mean distance
  max_int <- min(floor(length(chain)/2), 250L)
  
  # make numerical vector
  chain <- as.numeric(chain)
  
  # no ESS if all samples are identical
  if(length(unique(chain)) == 1) return(NaN)
  
  ### get ballpark m by looking at complete chain (intervals up to half chain length)
  # step 1: squared distance as function of interval
  intervalstart <- floor(length(chain)/(2 * max_int))
  squared_distances <- rep(0, max_int)
  for(i in 1:max_int) {
    samplinginterval <- 1 + (i - 1) * intervalstart
    
    # distance is 0 if identical, 1 if different
    squared_distances[i] <- mean(chain != c(tail(chain, -samplinginterval), head(chain, samplinginterval)))
  }
  
  # step 2: fit to exponential semivariogram to find asymptote
  Dpar <- if(squared_distances[1] >= mean(squared_distances)) {
    mean(squared_distances)
  } else {
    min(max(squared_distances),
        optim(c(.1, .1),
              function(x)
                sum(
                  (x[1] * (1 - exp(-x[2] * (1 + intervalstart*(0:(max_int - 1))))) -
                     squared_distances) ^ 2
                ))$par[1])
  }
  
  # step 3: no ESS if no asymptote has been reached
  if(all(tail(squared_distances, -1) - head(squared_distances, -1) > 0)) return(0)
  
  # step 4: calculate ballpark m
  mpar <- intervalstart * (which(squared_distances >= Dpar)[1] - 1) + 1
  
  
  ### get final m by calculating f(k) for k <= 1.25*mpar
  # step 1: squared distance as function of interval
  intervalend <- min(intervalstart, ceiling(mpar/200))
  sqdists <- rep(0, max_int)
  for(i in 1:max_int) {
    samplinginterval <- i * intervalend
    sqdists[i] <- mean(chain != c(tail(chain, -samplinginterval), head(chain, samplinginterval)))
  }
  
  # step 2: calculate m as lowest k where f is within 5% of asymptote
  mpar <- intervalend * (which(sqdists > 0.95 * Dpar)[1])
  
  ### calculate other parameters of equation
  Dpar <- max(Dpar, sqdists)
  Npar <- length(chain)
  
  ### calculate right-hand side of equation
  if(mpar == 0) {
    return(Npar)
  } else {
    right_hand_side <- intervalend * sum(
      sapply(1:(mpar/intervalend),
             function(x) sum(sqdists[x] * (Npar - x))
      )
    )
    right_hand_side <- right_hand_side + Dpar*((Npar - mpar) * (Npar - mpar + 1)/2)
    right_hand_side <- min(right_hand_side/(2 * Npar ^ 2), Dpar/4)
    
    ### final result up to a maximum of the chain length
    return(min(Npar, 1/(1 - 4 * right_hand_side / Dpar)))
  }
}

