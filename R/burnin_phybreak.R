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
                            parameter_frequency = 1, status_interval = 10) {
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
  
  protocoldistribution <- c(1 - classic - keepphylo - withinhost_only, classic, keepphylo, withinhost_only)
  
  build_pbe(x)

  message(paste0("   cycle      logLik         mu  gen.mean  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
  print_screen_log(0)
  
  curtime <- Sys.time()

  for (rep in 1:ncycles) {
    if(Sys.time() - curtime > status_interval) {
      print_screen_log(rep)
      curtime <- Sys.time()
    }
    for(i in sample(c(rep(-(1:6), parameter_frequency), 1:x$p$obs))) {
      if(i > 0) {
        which_protocol <- sample(c("edgewise", "classic", "keepphylo", "withinhost"),
                                 1,
                                 prob = protocoldistribution)
        update_host(i, which_protocol)
      }
    
      if (i == -1)  update_mu()
      if (i == -2 && x$h$est.mG)  update_mG()
      if (i == -3 && x$h$est.mS)  update_mS()
      if (i == -4 && x$h$est.wh.s)  update_wh_slope()
      if (i == -5 && x$h$est.wh.e)  update_wh_exponent()
      if (i == -6 && x$h$est.wh.0)  update_wh_level()
    }
  }
  
  res <- destroy_pbe(x$s)
  
  return(res)
}


#' @rdname burnin_phybreak
#' @export
burnin.phybreak <- function(...) {
  .Deprecated("burnin_phybreak")
  burnin_phybreak(...)
}

print_screen_log <- function(iteration) {
  message(paste0(
    stringr::str_pad(iteration, 8),
    stringr::str_pad(round(pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikcoal, 2), 12),
    stringr::str_pad(signif(pbe0$p$mu, 3), 11),
    stringr::str_pad(signif(pbe0$p$gen.mean, 3), 10),
    stringr::str_pad(signif(pbe0$p$sample.mean, 3), 10),
    stringr::str_pad(phangorn::parsimony(
      phybreak2phylo(environment2phybreak(pbe0$v)), pbe0$d$sequences), 10)))
}
