iors <- match(simul3$infectors, names(simul3$sample.times))
neworder <- order(simul3$infection.times)
plot(c(0,15),c(0,51))
segments(x0 = simul3$infection.times[neworder], y0 = 1:50, x1 = simul3$sample.times[neworder])
segments(x0 = simul3$infection.times[neworder], y0 = 1:50, y1 = match(iors[neworder],neworder), lty = 3)


plotsimul

inftimes <- simul2$infection.times
samtimes <- simul2$sample.times
hosts <- names(inftimes)
infectors <- simul2$infectors
names(infectors) <- hosts
parmat <- matrix(0,nrow = 50, ncol = 50, dimnames = list(hosts,
                                                         hosts))
for(i in hosts) {
  curhost <- infectors[i]
  while(curhost != "index") {
    parmat[i, curhost] <- parmat[i, curhost] + 1
    curhost <- infectors[curhost]
  }
}
hostweights <- colSums(parmat)
insides <- c()
for(i in names(hostweights)) {
  insides <- c(insides, 
               cumsum(sort(hostweights[which(infectors == i)])) < 0.5*hostweights[i])
}
timedorder <- order(inftimes)
inftimes <- inftimes[timedorder]
samtimes <- samtimes[timedorder]
infectors <- infectors[timedorder]
hosts <- names(inftimes)
plotorder <- 1:2
for(i in 3:50) {
  ior <- which(hosts==infectors[i])
  if(plotorder[ior] == 1) {
    plotorder[i] <- insides[hosts[i]] + 0.5
  } else if (plotorder[ior] == i - 1) {
    plotorder[i] <- i - 0.5 - insides[hosts[i]]
  } else {
    #plotorder[i] <- plotorder[ior] + rank(inftimes[match(c(-1, 1) + plotorder[ior],plotorder)])[2] - 1.5
    #plotorder[i] <- plotorder[ior] + (plotorder[ior]*2 > i)  - 0.5
    plotorder[i] <- plotorder[ior] + 
      rank(inftimes[match(c(-1, 1) + plotorder[ior],plotorder)])[1 + insides[hosts[i]]] - 1.5
  }
  plotorder[1:i] <- rank(plotorder)[1:i]
}
#inftimes <- inftimes[plotorder]
#samtimes <- samtimes[plotorder]
#infectors <- infectors[plotorder]
#hosts <- hosts[plotorder]
infectorpos <- match(infectors, hosts)
infectorpos <- plotorder[infectorpos]
infectorpos[is.na(infectorpos)] <- 0
plot(c(0,17),c(0,51))
segments(x0 = inftimes[1:i], y0 = plotorder, x1 = samtimes[1:i])
segments(x0 = inftimes[1:i], y0 = plotorder, y1 = infectorpos[1:i], lty = 3)



infectors <- simul2$infectors
names(infectors) <- hosts
parmat <- matrix(0,nrow = 50, ncol = 50, dimnames = list(hosts,
                                                         hosts))
for(i in hosts) {
  curhost <- infectors[i]
  while(curhost != "index") {
    parmat[i, curhost] <- parmat[i, curhost] + 1
    curhost <- infectors[curhost]
  }
}
hostweights <- colSums(parmat)
insides <- c()
for(i in names(hostweights)) {
  insides <- c(insides, 
               cumsum(sort(hostweights[which(infectors == i)])) < 0.5*hostweights[i])
}
