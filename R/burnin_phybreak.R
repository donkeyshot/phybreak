#' MCMC updating of a phybreak-object.
#' 
#' This function allows the MCMC chain to burn in. If used after samples have been taken (with \code{\link{sample_phybreak}}), 
#'   these samples will be returned unchanged in the output (\code{burnin.phybreak} is deprecated).
#' 
#' @param x An object of class \code{phybreak}.
#' @param ncycles Number of iterations to be carried out. Each iteration does one update of all parameters and
#'   tree updates with each host as focal host once.
#' @param classic The proportion of tree updates with the classic protocol (reference see below), in which within-host
#'   minitrees are proposed by simulating coalescent times and tree topology. In the current default protocol only 
#'   coalescent times are proposed with the minitree topology kept intact. This is followed by removing and reconnecting
#'   the sampling tips one by one. This results in better mixing of the mcmc-chain if there is much 
#'   genetic information (many SNPs) and/or if there are many possible within-host minitree topologies 
#'   (e.g. many samples per host). The classic protocol is faster in terms of updates/second and can thus be more efficient
#'   with little genetic information.
#' @param keepphylo The proportion of tree updates keeping the phylogenetic tree intact, only possible if there is one 
#'   sample per host and the \code{wh.model = "linear"} with complete bottleneck.
#' @param withinhost_only The proportion of tree updates in which only the within-host minitree is sampled, and 
#'   the transmission tree and infection times are kept unchanged.
#' @param parameter_frequency The relative frequency by which the model parameters are updated relative to updating each host.
#' @param status_interval The number of seconds between each on-screen print of the progress of the mcmc-chain.
#' @return The \code{phybreak} object provided as input, with variables and parameters changed due to the updating.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 50)
#' @export
burnin_phybreak <- function(x, ncycles, classic = 0, keepphylo = 0, withinhost_only = 0, 
                            parameter_frequency = 1, status_interval = 10, 
                            historydist = 0.5,
                            nchains = 1, heats = NULL, swap = 1) {
  ### tests
  if(ncycles < 1) stop("ncycles should be positive")
  if(is.null(x$p$wh.bottleneck)) {
    x$p$wh.bottleneck <- choose_whbottleneck("auto", x$p$wh.model)
  }
  if(classic < 0 | classic > 1) stop("classic should be a fraction")
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(withinhost_only < 0 | withinhost_only > 1) stop("withinhost_only should be a fraction")
  if(withinhost_only + keepphylo + classic > 1) stop("classic + keepphylo + withinhost_only should be a fraction")
  if(keepphylo > 0) {
    if(any(duplicated(x$d$hostnames)) || !(x$p$wh.model %in% c(3, "linear")) || x$p$wh.bottleneck == "wide") {
      keepphylo <- 0
      warning("model incompatible with keepphylo-updates: they will not be used", immediate. = TRUE)
    } 
  }
  
  ### add distance model if not present
  if(is.null(x$p$dist.model)) {
    x$p$dist.model <- "none"
    x$p$dist.exponent <- 2
    x$p$dist.scale <- 1
    x$p$dist.mean <- 1
    x$s$dist.e <- c()
    x$s$dist.s <- c()
    x$s$dist.m <- c()
    x$h$est.dist.s <- FALSE
    x$h$est.dist.e <- FALSE
    x$h$est.dist.m <- FALSE
  }
  
  protocoldistribution <- c(1 - classic - keepphylo - withinhost_only, classic, keepphylo, withinhost_only)
  historydistribution <- c(historydist, 1 - historydist)
  
  build_pbe(x)
  
  if (is.null(heats))
    heats <- 1/(1+1*(1:nchains-1))
  else if (inherits(heats, "numeric") & length(heats) != nchains)
    stop("length of heats is not the same as number of chains")
  else if (!inherits(heats, "numeric"))
    stop("heats is not a numeric vector")
  
  envirs <- list()
  
  for (n in 1:nchains){
    heat <- heats[n]
    copy2pbe0("heat", environment())
    chain <- n
    copy2pbe0("chain", environment())
    envirs[[n]] <- as.environment(as.list(pbe0, all.names = TRUE))
  }
  
  message(paste0("   cycle      logLik  introductions       mu  gen.mean  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
  #message(paste0("   cycle      logLik  introductions       mu  t.half  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
  print_screen_log(0)
  
  curtime <- Sys.time()
  
  swap_thin <- 0
  shared_heats <- heats
  
  for (rep in 1:ncycles) {
    
    if(Sys.time() - curtime > status_interval) {
      for (i in ls(envir=envirs[[1]])) copy2pbe0(i, envirs[[1]])
      print_screen_log(rep)
      curtime <- Sys.time()
    }
    
    for (i in 1:nchains){
      envirs[[i]]$heat <- shared_heats[i]
    }
    
    envirs <- lapply(envirs, function(e) {
      for (i in ls(envir = e)) copy2pbe0(i,e)
        
      for(i in sample(c(rep(-(1:13), parameter_frequency), 0:x$p$obs))) {
        if(i >= 0) {
          which_protocol <- sample(c("edgewise", "classic", "keepphylo", "withinhost"),
                                   1,
                                   prob = protocoldistribution)
          history <- sample(c(TRUE, FALSE), 1, prob = historydistribution)
          if (i > 0) {
            update_host(i, which_protocol, history)
          } else if (i == 0 & history) {
            update_host(i, which_protocol, TRUE)
          }
        }
      
        if (i == -1 && x$h$est.mu)  update_mu()
        if (i == -2 && x$h$est.mG)  update_mG()
        if (i == -3 && x$h$est.mS)  update_mS()
        if (i == -11 && x$h$est.tG) update_tG()
        if (i == -12 && x$h$est.tS) update_tS()
        if (i == -10 && x$h$est.hist.m) update_hist_mean()
        if (i == -4 && x$h$est.wh.s)  update_wh_slope()
        if (i == -5 && x$h$est.wh.e)  update_wh_exponent()
        if (i == -6 && x$h$est.wh.0)  update_wh_level()
        if (i == -7 && x$h$est.dist.e)  update_dist_exponent()
        if (i == -8 && x$h$est.dist.s)  update_dist_scale()
        if (i == -9 && x$h$est.dist.m)  update_dist_mean()
      }
      
      as.environment(as.list(pbe0, all.names = TRUE))
    })
    
    if(nchains > 1 & rep %% swap == 0){
      shared_lik <- do.call(cbind, lapply(envirs, function(xx){
        sum(xx$logLikcoal, xx$logLikgen, xx$logLiksam, xx$logLikdist, 
            xx$logLikseq)
      }))
      
      shared_heats <- swap_heats(shared_heats, shared_lik)
    } 
    swap_thin <- swap_thin + 1
  }
  
  heats <- do.call(c, lapply(envirs, function(e){
    for (i in ls(envir=e)) copy2pbe0(i, e)
    return(pbe0$heat)
  }))
  
  for (i in ls(envir=envirs[[which(heats == 1)]])) copy2pbe0(i, envirs[[which(heats == 1)]])
  
  res <- list(d = pbe0$d, v = environment2phybreak(pbe0$v), p = pbe0$p, h = pbe0$h, s = x$s,
                     hist = pbe0$v$inftimes[1])
  class(res) <- c("phybreak", "list")
  rm(list = ls(pbe0), envir = pbe0)
  rm(list = ls(pbe1), envir = pbe1)
  
  return(res)
}


#' @rdname burnin_phybreak
#' @export
burnin.phybreak <- function(...) {
  .Deprecated("burnin_phybreak")
  burnin_phybreak(...)
}

print_screen_log <- function(iteration) {
  if(pbe0$p$trans.model == "gamma"){
    message(paste0(
      stringr::str_pad(iteration, 8),
      stringr::str_pad(round(sum(pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikcoal + pbe0$logLikdist), 2), 12),
      stringr::str_pad(signif(sum(pbe0$v$infectors==0), 1), 15),
      stringr::str_pad(signif(pbe0$p$mu, 3), 9),
      stringr::str_pad(signif(pbe0$p$gen.mean, 3), 10),
      stringr::str_pad(signif(pbe0$p$sample.mean, 3), 10),
      stringr::str_pad(phangorn::parsimony(
        phybreak2phylo(environment2phybreak(pbe0$v)), pbe0$d$sequences), 10)))
  } else if (pbe0$p$trans.model == "sampling+culling"){
    message(paste0(
      stringr::str_pad(iteration, 8),
      stringr::str_pad(round(pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikcoal + pbe0$logLikdist, 2), 12),
      stringr::str_pad(signif(sum(pbe0$v$infectors==0), 1), 15),
      stringr::str_pad(signif(pbe0$p$mu, 3), 9),
      stringr::str_pad(signif(log((1-pbe0$p$trans.init)/pbe0$p$trans.init)/pbe0$p$trans.growth, 3), 10),
      stringr::str_pad(signif(pbe0$p$sample.mean, 3), 10),
      stringr::str_pad(phangorn::parsimony(
        phybreak2phylo(environment2phybreak(pbe0$v)), pbe0$d$sequences), 10)))
  }
}
