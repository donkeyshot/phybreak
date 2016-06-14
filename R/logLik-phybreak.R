### calculate the log-likelihood ###

### phybreak functions called ###
# .likseq  ##C++
# .lik.gentimes
# .lik.sampletimes
# .lik.coaltimes


logLik.phybreak <- function(phybreak.object, genetic = TRUE, withinhost = TRUE,
                            sampling = TRUE, generation = TRUE) {
  res <- 0
  if(genetic) {
    res <- res + with(phybreak.object, .likseq(t(d$SNP), d$SNPfr,
                                            v$nodeparents, v$nodetimes, p$mu,p$obs))
  } 
  if(generation) {
    res <- res + with(phybreak.object, .lik.gentimes(p$obs, p$shape.gen, p$mean.gen,
                                                  v$nodetimes, v$nodehosts, v$nodetypes))
  } 
  if(sampling) {
    res <- res + with(phybreak.object, .lik.sampletimes(p$shape.sample, p$mean.sample,
                                                     v$nodetimes, v$nodetypes))
  } 
  if(withinhost) {
    res <- res + with(phybreak.object, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, 
                                                   v$nodetimes, v$nodehosts, v$nodetypes))
  }
  return(res)
}
