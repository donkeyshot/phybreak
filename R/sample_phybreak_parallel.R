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
sample_phybreak_parallel <- function(x, nsample, thin = 1, thinswap = 1, classic = 0, keepphylo = 0, withinhost_only = 0, 
                            parameter_frequency = 1, status_interval = 10, histtime = -1e5, history = FALSE,
                            nchains = 1, heats = NULL, all_chains = FALSE, 
                            outfile = "~/phybreak_parallel_output.txt", phybreakdir) {
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
                  hist_dens = c(x$s$hist_dens, rep(NA, nsample)),
                  mG = c(x$s$mG, rep(NA, nsample)), 
                  mS = c(x$s$mS, rep(NA, nsample)), 
                  tG = c(x$s$tG, rep(NA, nsample)),
                  tS = c(x$s$tS, rep(NA, nsample)),
                  wh.h = c(x$s$wh.h, rep(NA, nsample)), 
                  wh.s = c(x$s$wh.s, rep(NA, nsample)), 
                  wh.e = c(x$s$wh.e, rep(NA, nsample)), 
                  wh.0 = c(x$s$wh.0, rep(NA, nsample)), 
                  dist.e = c(x$s$dist.e, rep(NA, nsample)), 
                  dist.s = c(x$s$dist.s, rep(NA, nsample)), 
                  dist.m = c(x$s$dist.m, rep(NA, nsample)), 
                  logLik = with(x, cbind(s$logLik, matrix(NA, nrow = 5, ncol = nsample))),
                  historyinf = c(x$s$historyinf, rep(NA, nsample)),
                  heat = c(x$s$heat, rep(NA, nsample)))
  
  #s.posts <- lapply(1:nchains, function(i) s.post)
  
  ### Check if phybreak object contains the history host, else make one 
  if (!history)
    build_pbe(x, histtime)
  
  ### Set heats for MC3, heats = 1 if there is one chain
  if (is.null(heats))
    heats <- 1/(1+1*(1:nchains-1))
  else if (inherits(heats, "numeric") & length(heats) != nchains)
    stop("length of heats is not the same as number of chains")
  else if (!inherits(heats, "numeric"))
    stop("heats is not a numeric vector")
  
  ### Make swap scheme
  swap <- t(replicate(nsample/thinswap, sample(nchains, 2)))
  
  ### Initialize environments of chains, add heat and chain number 
  ### to the environment
  envirs <- list()
  for (n in 1:nchains){
    heat <- heats[n]
    copy2pbe0("heat", environment())
    chain <- n
    copy2pbe0("chain", environment())
    envirs[[n]] <- as.environment(as.list(pbe0, all.names = TRUE))
  }
  
  ### Set up cluster
  cores <- detectCores()
  cl <- makeCluster(min(nchains, cores-1), outfile = outfile)
  mgrinit(cl)
  makebarr(cl)
  
  ### Make shared variable
  mgrmakevar(cl, "shared_heats", 1, nchains)
  shared_heats[,] <- heats
  
  mgrmakevar(cl, "shared_lik", 1, nchains)
  
  ### Export variables and functions
  #clusterExport(cl, varlist = ls(pos=which(search() == "package:phybreak")))
  clusterExport(cl, varlist = c("x", "envirs", "nsample", "thin", "thinswap", "protocoldistribution",
                                "s.post", "status_interval", "parameter_frequency",
                                "phybreakdir"), 
                envir = environment())
  
  #message(paste0("  sample      logLik  introductions       mu  gen.mean  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
  #print_screen_log(length(x$s$mu))
  print("Intermediate results of MCMC chain are shown in the parallel output file")
  
  ### Process MCMC for each chain in parallel and share heats and likelihoods
  process_mcmc <- function(shared_heats, shared_lik) {
    require(Rdsm)
    devtools::clean_dll()
    devtools::load_all(phybreakdir)

    me <- myinfo$id

    for (i in ls(envir=envirs[[me]])) copy2pbe0(i, envirs[[me]])
    rm(envirs)
    curtime <- Sys.time()

    for(sa in 1:nsample * thin) {
      pbe0$heat <- shared_heats[1,me]

      ### Print intermediate results to outFile
      if( me == 1){
        if(Sys.time() - curtime > status_interval) {
          print_screen_log(floor(sa / thin))
          curtime <- Sys.time()
        }
      }

      ### Propose updates of hosts and parameters
      for(i in  sample(c(rep(-(1:12), parameter_frequency), 1:(x$p$obs+1)))) {
        if(i > 0) {
          which_protocol <- sample(c("edgewise", "classic", "keepphylo", "withinhost"),
                                   1,
                                   prob = protocoldistribution)
          update_host(i, which_protocol)
        }

        if (i == -1)  update_mu()
        if (i == -2 && x$h$est.mG)  update_mG()
        if (i == -3 && x$h$est.mS)  update_mS()
        if (i == -11 && x$h$est.tG) update_tG()
        if (i == -12 && x$h$est.tS) update_tS()
        if (i == -10 && x$h$est.wh.h) update_wh_history()
        if (i == -4 && x$h$est.wh.s)  update_wh_slope()
        if (i == -5 && x$h$est.wh.e)  update_wh_exponent()
        if (i == -6 && x$h$est.wh.0)  update_wh_level()
        if (i == -7 && x$h$est.dist.e)  update_dist_exponent()
        if (i == -8 && x$h$est.dist.s)  update_dist_scale()
        if (i == -9 && x$h$est.dist.m)  update_dist_mean()
      }

      shared_lik[1, me] <- sum(pbe0$logLikcoal, pbe0$logLikgen, pbe0$logLiksam,
                               pbe0$logLikdist, pbe0$logLikseq)
      
      ### Store current state
      if (sa %% thin == 0){
  
        s.post$historyinf[sa] <- pbe0$v$inftimes[1]
        remove_history(keepenv = TRUE)
        vars_to_log <- environment2phybreak(pbe2$v)
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
        s.post$tG[sa] <- pbe0$p$trans.growth
        s.post$tS[sa] <- pbe0$p$trans.sample
        s.post$wh.h[sa] <- pbe0$p$wh.history
        s.post$wh.s[sa] <- pbe0$p$wh.slope
        s.post$wh.e[sa] <- pbe0$p$wh.exponent
        s.post$wh.0[sa] <- pbe0$p$wh.level
        s.post$dist.e[sa] <- pbe0$p$dist.exponent
        s.post$dist.s[sa] <- pbe0$p$dist.scale
        s.post$dist.m[sa] <- pbe0$p$dist.mean
        s.post$logLik[, sa] <- c(pbe0$logLikcoal, pbe0$logLikgen, pbe0$logLiksam, pbe0$logLikdist, pbe0$logLikseq)
        s.post$heat[sa] <- shared_heats[1, me]
      }
  
      if (sa %% thinswap == 0){
        barr()
        if (me == 1){
          shared_heats[,] <- swap_heats(shared_heats[,], shared_lik[,])
        }
        barr()
      }
  
      gc()
    }
    return(s.post)
  }
  
  clusterExport(cl, "process_mcmc", envir = environment())
  #clusterEvalQ(cl, process_mcmc(shared_heats, shared_lik))

  posts <- clusterEvalQ(cl, process_mcmc(shared_heats, shared_lik))
  stopCluster(cl)
  
  ### Sort samples per chain
  heat_index <- matrix(nrow = nsample, ncol = 3)
  for (i in 1:nchains) {
    for (j in 1:length(posts)) {
      heat_index[which(posts[[j]]$heat == heats[i]),j] <- i
    }
  }
  
  s.posts <- lapply(1:nchains, function(i){
    s <- s.post
    for (n in names(s)){
      for (j in 1:length(posts)){
        if (inherits(s[[n]], "matrix"))
          s[[n]][,which(heat_index[,i]==j)] <- posts[[j]][[n]][,which(heat_index[,i]==j)]
        else
          s[[n]][which(heat_index[,i]==j)] <- posts[[j]][[n]][which(heat_index[,i]==j)]
      }
    }
    return(s)
  })
  
  s.posts <- lapply(1:length(heats), function(nheat){
    s <- s.post
    for (chain in 1:nchains){
      smp <- which(posts[[chain]]$heat == heats[nheat])
      for (n in names(s)){
        if(inherits(s[[n]], "matrix"))
          s[[n]][,smp] <- posts[[chain]][[n]][,smp]
        else
          s[[n]][smp] <- posts[[chain]][[n]][smp]
      }
    }
    s$chain <- posts[[nheat]]$heat
    return(s)
  })
  
  
  s.post <- s.posts[[1]]
  
  res <- destroy_pbe(s.post)
  
  if(all_chains)
    return(list(phyb.obj = res, chains = s.posts))
  else
    return(res)
}
  
  
  