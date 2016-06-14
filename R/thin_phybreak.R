### thin samples in phybreak.object ###


thin.phybreak <- function(phybreak.object, thin = 1, nkeep = Inf) {
  if(!inherits(phybreak.object, "phybreak")) stop("object should be of class \"phybreak\"")
  if(nkeep == 0) {
    return(
      within(
        phybreak.object,
        s <- list(
          nodetimes = c(),
          nodehosts = c(),
          nodeparents = c(),
          mu = c(),
          mG = c(),
          mS = c(),
          slope = c(),
          logLik = c()
        )
      )
    )
  }
  tokeep <- seq(thin, length(phybreak.object$s$logLik), thin)
  tokeep <- tail(tokeep, nkeep)
  return(
    within(
      phybreak.object,
      s <- list(
        nodetimes = s$nodetimes[,tokeep],
        nodehosts = s$nodehosts[,tokeep],
        nodeparents = s$nodeparents[,tokeep],
        mu = s$mu[tokeep],
        mG = s$mG[tokeep],
        mS = s$mS[tokeep],
        slope = s$slope[tokeep],
        logLik = s$logLik[tokeep]
      )
    )
  )
}