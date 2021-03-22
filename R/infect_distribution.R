#' Infectiousness distribution.
#' 
#' This function calculates the probability of transmission given the infectiousness 
#' distribution for the hosts.
#' 
#' @param time The proposed infection time of the current host
#' @param vars The list of variables describing the nodes in the tree. This list should 
#'   at least include the infection times of the hosts. For \code{p$trans.model = "sample+culling"}
#'   also nodetimes of the sample nodes and the culling times of the hosts are needed as input.
#' @param p The list with parameters for the model.  
#' 
#' @return The transmission probabilities of all hosts infected before the current host. For each host
#'   with an infection time before the proposed infection time, the probability of being an infector
#'   of the current host is calculated according to the given infectiousness function.
#' 
#' @export
infect_distribution <- function(time, inftimes, p, 
                                nodetimes = NULL, cultimes = NULL,  
                                host = NULL, log = FALSE){
  if(is.null(p$trans.model))
    p$trans.model <- "gamma"
  if(!inherits(p$trans.model, "character"))
    stop("trans.model must be a character")
  if(!(p$trans.model %in% c("gamma", "unif", "sampling+culling", "user-defined")))
    stop("transmission model is unknown")
  
  ### Gamma distributed ####
  if(p$trans.model == "gamma") {
    if(log)
      return(dgamma(time - inftimes, 
                    shape = p$gen.shape, 
                    scale = p$gen.mean/p$gen.shape,
                    log = TRUE))
    else 
      return(dgamma(time - inftimes, 
                    shape = p$gen.shape, 
                    scale = p$gen.mean/p$gen.shape,
                    log = FALSE))
  
  ### Uniform distributed ###
  } else if(p$trans.model == "unif") {
    return(dunif(time - inftimes, min = 0, max = max(time-vars$inftimes)))
  
  ### Sampling and culling of hosts ###
  } else if(p$trans.model =="sampling+culling") {
    
    if(is.null(nodetimes) | is.null(cultimes))
      stop("times of first sample and culling times of hosts must be provided")
    # if(length(nodetimes) != length(inftimes)-1)
    #   stop("length of nodetimes must be equal to number of hosts")
    # if(length(cultimes) != length(inftimes)-1)
    #   stop("length of cultimes must be equal to number of hosts")
 
    if(is.null(p$trans.init))
      stop("initial fraction infected is missing")
    if(is.null(p$trans.growth))
      stop("growth factor of infectiousness is missing")
    if(is.null(p$trans.sample))
      stop("reduction factor after first positive sample is missing")
    if(is.null(p$trans.culling))
      stop("decay factor after culling is missing")
    
    a <- (1-p$trans.init)/p$trans.init
    r <- p$trans.growth
    S <- p$trans.sample
    C <- p$trans.culling
    
    if(length(nodetimes) == length(inftimes)-1)
      samtimes <- nodetimes - inftimes[-1]
    else 
      samtimes <- nodetimes - inftimes
    
    if(length(cultimes) == length(inftimes)-1)
      cultimes <- cultimes - inftimes[-1]
    else 
      cultimes <- as.numeric(cultimes - inftimes)
    
    hosttimes <- as.numeric(time - inftimes)
    
    if(is.null(host)){
      if(length(hosttimes) != length(nodetimes)){
        probs <- 0.1
        j <- 1
      } else {
        probs <- c()
        j <- 0
      }
      for (i in 1:length(samtimes)){
        if(hosttimes[i+j] < 0)
          probs <- c(probs, 0)
        else if(hosttimes[i+j] < samtimes[i])
          #probs <- c(probs, 1/(1+exp(-1/2*(hosttimes[i+j]-p$gen.mean))))
          probs <- c(probs, 1/(1+a*exp(-r*hosttimes[i+j])))
        else if(hosttimes[i+j] >= samtimes[i] & hosttimes[i+j] < cultimes[i])
          #probs <- c(probs, p$gen.sample.scale*1/(1+exp(-(hosttimes[i+j]-p$gen.mean))))
          probs <- c(probs, S/(1+a*exp(-r*hosttimes[i+j])))
        else 
          # probs <- c(probs, p$gen.sample.scale*1/(1+exp(-(cultimes[i+j]-p$gen.mean))) * 
          #              exp(-p$gen.culling.scale*(hosttimes[i+j]-cultimes[i])))
          probs <- c(probs, S/(1+a*exp(-r*cultimes[i])) * exp(-C*(hosttimes[i+j]-cultimes[i])))
      }
    } else {
      probs <- c()
      for (t in hosttimes){
        if (t < samtimes[host]){
          #probs <- c(probs, 1/(1+exp(-1/2*(t-p$gen.mean))))
          probs <- c(probs, 1/(1+a*exp(-r*t)))
        } else if (t >= samtimes[host] & t < cultimes[host]){
          #probs <- c(probs, p$gen.sample.scale*1/(1+exp(-(t-p$gen.mean))))
          probs <- c(probs, S*1/(1+a*exp(-r*t)))
        } else if (t >= cultimes[host]){
          # probs <- c(probs, p$gen.sample.scale*1/(1+exp(-(cultimes[host]-p$gen.mean))) * 
          #          exp(-p$gen.culling.scale*(t-cultimes[host])))
          probs <- c(probs, S*1/(1+a*exp(-r*cultimes[host])) * exp(-C*(t-cultimes[host])))
        }
      }
    }
    
    if(log)
      return(log(probs))
    else
      return(probs)
  }
}
