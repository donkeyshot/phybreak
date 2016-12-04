maketransplot <- function(x, ...) {
  
  ### sort by time of infection
  timedorder <- order(x$sim.infection.times)
  inftimes <- x$sim.infection.times[timedorder]
  samtimes <- x$sample.times[timedorder]
  infectors <- x$sim.infectors[timedorder]
  hosts <- names(inftimes)

  ### determine rank of each host in the plot (line number)
  plotrank <- rankhostsforplot(hosts, infectors)
  
  infectorpos <- match(infectors, hosts)
  infectorpos <- plotrank[infectorpos]
  infectorpos[is.na(infectorpos)] <- 0
  
  tmin <- min(inftimes)
  tmax <- max(samtimes)
  tstep <- as.numeric((tmax-tmin)/2000)
  tgmean <- as.numeric(mean(inftimes - inftimes[infectors], na.rm = T))
  tgvar <- as.numeric(var(inftimes - inftimes[infectors], na.rm = T))
  tgscale <- tgvar/tgmean
  tgshape <- tgmean/tgscale
  maxwd <- max(dgamma(seq(0, as.numeric(tmax - tmin), tstep), shape = tgshape, scale = tgscale))
  obs <- length(hosts)
  plot.new()
  par(mar = 0.1 + c(5, 0, 0, 0), cex = 1)
  plot.window(xlim = c(tmin, tmax + 0.1*(tmax-tmin)), ylim = c(0,obs + 1))
  axis(1, at = axTicks(1), 
       labels = if(inherits(tmin, "Date")) as.Date(axTicks(1), "1970-01-01") else TRUE)
  for(i in 1:obs) {
    x0s <- seq(inftimes[i], tmax - tstep, tstep)
    widths <- abs(1 - (maxwd - dgamma(x0s - inftimes[i], shape = tgshape, scale = tgscale)) / maxwd)
    polygon(x = c(x0s, rev(x0s)), y = plotrank[i] + 0.3 * c(widths, -rev(widths)), col = "gray", border = NA)
    text(max(samtimes) + 10*tstep, plotrank[i], hosts[i], adj = 0, cex = min(1, 30/obs))
  }
  segments(x0 = inftimes, y0 = plotrank, x1 = max(samtimes), lty = 3)
  arrows(x0 = inftimes[infectors != "index"], 
         y0 = plotrank[infectors != "index"], 
         y1 = infectorpos[infectors != "index"], lwd = 1, length = .05, code = 1)
  points(samtimes, plotrank, pch = 4, lwd = 2, cex = min(1, 30/obs))
}

rankhostsforplot <- function(hosts, infectors) {
  ### extract parameters
  Nhosts <- length(hosts)

  ### calculate branch weights of infection tree per host
  infectormatrix <- matrix(0, nrow = Nhosts, ncol = Nhosts, dimnames = list(hosts,
                                                                            hosts))
  for(i in hosts) {
    curhost <- infectors[i]
    while(curhost != "index") {
      infectormatrix[i, curhost] <- infectormatrix[i, curhost] + 1
      curhost <- infectors[curhost]
    }
  }
  branchweights <- colSums(infectormatrix) + 1
  
  ### determine position of each host relative to infector's infector
  insideYN <- c()
  for(i in names(branchweights)) {
    insideYN <- c(insideYN, 
                  cumsum(sort(branchweights[which(infectors == i)])) <= 0.5*(branchweights[i] - 1))
  }
  
  ### determine plot rank by placing hosts chronologically next to their infector,
  ### either at the side of the infector's infector (insideYN == Y) or at the other side
  plotrank <- 1:2
  for(i in 3:Nhosts) {
    ior <- which(hosts == infectors[i])
    if(plotrank[ior] == 1) {
      plotrank[i] <- insideYN[hosts[i]] + 0.5
    } else if (plotrank[ior] == i - 1) {
      plotrank[i] <- i - 0.5 - insideYN[hosts[i]]
    } else {
      iorior <- which(hosts == infectors[ior])
      aboveYN <- xor((plotrank[iorior] < plotrank[ior]), insideYN[hosts[i]])
      if(infectors[ior] == "index") aboveYN <- insideYN[hosts[i]]
      plotrank[i] <- plotrank[ior] + aboveYN - 0.5    }
    plotrank[1:i] <- rank(plotrank)[1:i]
  }
  
  return(plotrank)
}




