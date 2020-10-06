#' Sampling from a phybreak MCMC-chain.
#' 
#' Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree 
#'   (plus associated parameters), within a \code{phybreak} object (\code{sample.phybreak} is deprecated).
#' 
#' @param x An object of class \code{phybreak}.
#' @param nsample The number of samples to take.
#' @param thin The thinning to use (values after every \code{thin}'th iteration will be included in the posterior). 
#'   Each iteration does one update of all parameters and tree updates with each host as focal host once.
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
#' @return The \code{phybreak} object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample_phybreak <- function(x, nsample, thin = 1, classic = 0, keepphylo = 0, withinhost_only = 0, 
                            parameter_frequency = 1, status_interval = 10) {
    ### tests
    if(nsample < 1) stop("nsample should be positive")
    if(thin < 1) stop("thin should be positive")
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
  
    ### create room in s to add the new posterior samples
    s.post <- list(inftimes = with(x, cbind(s$inftimes, matrix(NA, nrow = p$obs, ncol = nsample))),
                   infectors = with(x, cbind(s$infectors, matrix(NA, nrow = p$obs, ncol = nsample))),
                   nodetimes = with(x, cbind(s$nodetimes, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                   nodehosts = with(x, cbind(s$nodehosts, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                   nodeparents = with(x, cbind(s$nodeparents, matrix(NA, nrow = 2 * d$nsamples - 1, ncol = nsample))),
                   introductions = c(sum(x$s$infectors==0), rep(NA, nsample)),
                   mu = c(x$s$mu, rep(NA, nsample)), 
                   mG = c(x$s$mG, rep(NA, nsample)), 
                   mS = c(x$s$mS, rep(NA, nsample)), 
                   wh.h = c(x$s$wh.h, rep(NA, nsample)), 
                   wh.s = c(x$s$wh.s, rep(NA, nsample)), 
                   wh.e = c(x$s$wh.e, rep(NA, nsample)), 
                   wh.0 = c(x$s$wh.0, rep(NA, nsample)), 
                   dist.e = c(x$s$dist.e, rep(NA, nsample)), 
                   dist.s = c(x$s$dist.s, rep(NA, nsample)), 
                   dist.m = c(x$s$dist.m, rep(NA, nsample)), 
                   logLik = c(x$s$logLik, rep(NA, nsample)))
    
    build_pbe(x)

    message(paste0("  sample      logLik         mu  gen.mean  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
    print_screen_log(length(x$s$mu))
    
    curtime <- Sys.time()
    
    for (sa in tail(1:length(s.post$mu), nsample)) {
      
      for (rep in 1:thin) {
        if(Sys.time() - curtime > status_interval) {
          print_screen_log(sa)
          curtime <- Sys.time()
        }
        for(i in  sample(c(rep(-(1:10), parameter_frequency), 2:x$p$obs))) {
          if(i > 0) {
            which_protocol <- sample(c("edgewise", "classic", "keepphylo", "withinhost"),
                                     1,
                                     prob = protocoldistribution)
            update_host(i, which_protocol)
          }
        
          if (i == -1)  update_mu()
          if (i == -2 && x$h$est.mG)  update_mG()
          if (i == -3 && x$h$est.mS)  update_mS()
          if (i == -10 && x$h$est.wh.h) update_wh_history()
          if (i == -4 && x$h$est.wh.s)  update_wh_slope()
          if (i == -5 && x$h$est.wh.e)  update_wh_exponent()
          if (i == -6 && x$h$est.wh.0)  update_wh_level()
          if (i == -7 && x$h$est.dist.e)  update_dist_exponent()
          if (i == -8 && x$h$est.dist.s)  update_dist_scale()
          if (i == -9 && x$h$est.dist.m)  update_dist_mean()
        }
      }
      remove_history(keepenv = TRUE)
      vars_to_log <- environment2phybreak(pbe0_2$v)
      s.post$inftimes[, sa] <- vars_to_log$inftimes
      s.post$infectors[, sa] <- vars_to_log$infectors
      s.post$nodetimes[, sa] <- c(vars_to_log$nodetimes[vars_to_log$nodetypes == "c"],
                                  rep(NA,x$d$nsamples-1-length(which(vars_to_log$nodetypes == "c"))))
      s.post$nodehosts[, sa] <- c(vars_to_log$nodehosts[vars_to_log$nodetypes == "c"],
                                  rep(NA,x$d$nsamples-1-length(which(vars_to_log$nodetypes == "c"))))
      s.post$nodeparents[, sa] <- c(vars_to_log$nodeparents,
                                    rep(NA,x$d$nsamples-1-length(which(vars_to_log$nodetypes == "c"))))
      s.post$introductions[sa] <- sum(vars_to_log$infectors == 0)
      s.post$mu[sa] <- pbe0$p$mu
      s.post$mG[sa] <- pbe0$p$gen.mean
      s.post$mS[sa] <- pbe0$p$sample.mean
      s.post$wh.h[sa] <- pbe0$p$wh.history
      s.post$wh.s[sa] <- pbe0$p$wh.slope
      s.post$wh.e[sa] <- pbe0$p$wh.exponent
      s.post$wh.0[sa] <- pbe0$p$wh.level
      s.post$dist.e[sa] <- pbe0$p$dist.exponent
      s.post$dist.s[sa] <- pbe0$p$dist.scale
      s.post$dist.m[sa] <- pbe0$p$dist.mean
      s.post$logLik[sa] <- pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikdist + pbe0$logLikcoal
    }
    
    res <- destroy_pbe(s.post)
    
    
    return(res)
    
}



#' @rdname sample_phybreak
#' @export
sample.phybreak <- function(...) {
  .Deprecated("sample_phybreak")
  sample_phybreak(...)
}

