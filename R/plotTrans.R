#' Plotting a phybreak object transmission tree.
#' 
#' Plots a \code{phybreak}-object as transmission tree. The default 
#'   is to plot the current state, but any posterior sample can be chosen, as well as various 
#'   consensus trees; in that case, coloured arrows indicate posterior support.
#' 
#' @param x An object of class \code{phybreak}.
#' @param plot.which Either \code{"sample"} to plot the current state or a selected posterior sample, 
#'   \code{"mpc"} or \code{"mtcc"} to plot a consensus transmission tree (see \code{\link{transtree}}) or \code{"mcc"}
#'   to plot the maximum clade credibility tree (see \code{\link{phylotree}}).
#' @param samplenr If \code{plot.which = "sample"}, this indicates which posterior tree should be plotted: 
#'   \code{samplenr = 0} to plot the current state.
#' @param mar Plot margins.
#' @param label.cex Size of host names, as in \code{par("cex")}. Defaults to a value between 0.5 and 1
#'   depending on outbreak size.
#' @param label.space Scales the space at the right-hand side to place the host names.
#' @param label.adj Left-right adjustment of host names.
#' @param arrow.lwd Arrow width.
#' @param arrow.length Arrow point length, as default automatically scaled with outbreak size.
#' @param arrow.col Arrow colour. Defaults to \code{"black"} if \code{plot.which = sample}, and otherwise to five
#'   colours \code{c("blue", "green", "orange", "red", "purple")} indicating posterior support of infectors
#'   in bins of 0.2 width, from low to high support. Any vector of colours will be divided into equal-sized bins.
#' @param sample.pch Character \code{par("pch")} used for sampling events.
#' @param sample.lwd Line width of sampling event character.
#' @param sample.cex Size of sampling event character.
#' @param polygon.col Color of polygons indicating generation interval distributions.
#' @param polygon.border Border of polygon.
#' @param line.lty Line type of horizontal host lines.
#' @param xlab X-axis title.
#' @param axis.cex Size of tick labels.
#' @param title.cex Size of X-axis title.
#' @param ... Further graphical parameters (see details).
#' @section Details: 
#'   Graphical parameters can be added by using names in the format \code{prefix.parameter} for the
#'   different parts of the plot. The \code{parameter} will then be passed on to the appropriate 
#'   graphics function, e.g. \code{arrow.lty} to change the line type of the arrows. The following 
#'   prefixes can be used: \code{label} for the host labels, \code{arrow} for the arrows, \code{sample}
#'   for the sampling time indicators, \code{polygon} for the generation interval distributions, 
#'   \code{line} for the horizontal host lines, \code{axis} for the X-axis, and \code{title} for the
#'   X-axis title.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' plot(MCMCstate, plot.which = "mpc")
#' @export
plotTrans <- function(x, plot.which = c("sample", "edmonds", "mpc", "mtcc"), samplenr = 0,
                      mar = 0.1 + c(4, 0, 0, 0), label.cex = NULL, 
                      label.space = 0.15, label.adj = 0,
                      arrow.lwd = 1, arrow.length = NULL, arrow.col = NULL, sample.pch = 4,
                      sample.lwd = NULL, sample.cex = label.cex, polygon.col = "gray", 
                      polygon.border = NA, line.lty = 3, xlab = "Time", 
                      axis.cex = 1, title.cex = 1, ...) {
  if(!inherits(x, c("phybreak", "phybreakdata"))) {
    stop("object x must be of class \"phybreak\" or \"phybreakdata\"")
  }
  
  if(inherits(x, "phybreakdata")) {
    if(exists("sim.infection.times", x) && exists("sim.infectors", x)) {
      vars <- x
      names(vars$sample.times) <- vars$sample.hosts
      tg.mean <- NA
      tg.shape <- NA
    } else {
      stop("object x should contain sim.infection.times and sim.infectors")
    }
    
  } else {
    # class "phybreak"
    
    # test for plot.which
    plot.which <- plot.which[which(plot.which %in% c("sample", "edmonds", "mpc", "mtcc"))[1]]
    if(is.na(plot.which)) stop("no valid 'plot.which'")
    if(plot.which == "sample" & samplenr > length(x$s$logLik)) {
      warning("requested 'samplenr' not available; current state used")
      samplenr <- 0
    }
    
    # make vars-list with the tree to plot
    if (plot.which != "sample") {
      tree2plot <- suppressWarnings(transtree(x, plot.which, 
                                              infection.times = "infector", time.quantiles = 0.5))
      vars <- list(sample.times = x$d$sample.times,
                   sample.hosts = x$d$hostnames,
                   sim.infection.times = tree2plot[, 3],
                   sim.infectors = as.character(tree2plot[, 1]),
                   post.support = tree2plot[, 2])
      names(vars$sample.times) <- x$d$hostnames
      names(vars$sim.infection.times) <- x$d$hostnames[1:x$p$obs]
      names(vars$sim.infectors) <- x$d$hostnames[1:x$p$obs]
      names(vars$post.support) <- x$d$hostnames[1:x$p$obs]
      if(is.null(arrow.col)) {
        arrow.col <- c("blue", "green", "orange", "red", "purple")
      }
      tg.mean <- median(x$s$mG)
      tg.shape = x$p$shape.gen
    } else if (samplenr == 0) {
      # plot.which == "sample"
      
      vars <- list(
        nodetimes = x$v$nodetimes,
        nodeparents = x$v$nodeparents,
        nodehosts = x$v$nodehosts,
        nodetypes = x$v$nodetypes
      )
      vars <- phybreak2trans(vars, x$d$hostnames, x$d$reference.date)
      if(is.null(arrow.col)) {
        arrow.col <- "black"
      } else {
        arrow.col <- arrow.col[1]
      }
      tg.mean <- x$p$mean.gen
      tg.shape = x$p$shape.gen
    } else {
      # plot.which == "sample" && samplenr > 0
      
      vars <- list(
        nodetimes = c(x$v$nodetimes[x$v$nodetypes == "s"], x$s$nodetimes[, samplenr]),
        nodeparents = x$s$nodeparents[, samplenr],
        nodehosts = c(x$v$nodehosts[x$v$nodetypes == "s"], x$s$nodehosts[, samplenr]),
        nodetypes = x$v$nodetypes
      )
      vars <- phybreak2trans(vars, unique(x$d$hostnames), x$d$reference.date)
      if(is.null(arrow.col)) {
        arrow.col <- "black"
      } else {
        arrow.col <- arrow.col[1]
      }
      tg.mean <- x$s$mG[samplenr]
      tg.shape = x$p$shape.gen
    }
  }
  maketransplot(vars, tg.mean = tg.mean, tg.shape = tg.shape, mar = mar, label.cex = label.cex, 
                label.space = label.space, label.adj = label.adj,
                arrow.lwd = arrow.lwd, arrow.length = arrow.length, arrow.col = arrow.col, 
                sample.pch = sample.pch, sample.lwd = sample.lwd, sample.cex = sample.cex, 
                polygon.col = polygon.col, polygon.border = polygon.border, line.lty = line.lty,
                xlab = xlab, axis.cex = axis.cex, title.cex = title.cex, ...)
}

maketransplot <- function(x, tg.mean = NA, tg.shape = NA, mar = 0.1 + c(4, 0, 0, 4), label.cex = NULL, 
                          label.space = 0.15, label.adj = 0,
                          arrow.lwd = 1, arrow.length = NULL, arrow.col = par("fg"), sample.pch = 4,
                          sample.lwd = NULL, sample.cex = label.cex, polygon.col = "gray", 
                          polygon.border = NA, line.lty = 3, xlab = "Time", 
                          axis.cex = 1, title.cex = 1, ...) {
  oldmar <- par("mar")
  par(mar = mar)
  on.exit(par(mar = oldmar))

  ### transform support into colours
  if(length(arrow.col) == 1) {
    arrow.colours <- rep(arrow.col, length(x$sample.times))
  } else {
    if(exists("post.support", x)) {
      arrow.colours <- arrow.col[ceiling(length(arrow.col) * x$post.support)]
    } else {
      arrow.colours <- rep(par("fg"), length(x$sample.times))
    }
  }
  
  ### sort by time of infection
  ordertimes <- c(x$sim.infection.times, index = -Inf)
  inconsistenttimes <- which(head(ordertimes, -1) - ordertimes[x$sim.infectors] < 0)
  while(length(inconsistenttimes) > 0) {
    ordertimes[inconsistenttimes] <- ordertimes[x$sim.infectors[inconsistenttimes]] + 0.00001
    inconsistenttimes <- which(head(ordertimes, -1) - ordertimes[x$sim.infectors] < 0)
  }
  
  timedorder <- order(head(ordertimes, -1))
  inftimes <- x$sim.infection.times[timedorder]
  samtimes <- x$sample.times
  infectors <- x$sim.infectors[timedorder]
  arrow.colours <- arrow.colours[timedorder]
  hosts <- names(inftimes)

  ### determine rank of each host in the plot (line number)
  plotrank <- rankhostsforplot(hosts, infectors)
  
  ### determine rank of each host's infector
  infectorpos <- match(infectors, hosts)
  infectorpos <- plotrank[infectorpos]
  infectorpos[is.na(infectorpos)] <- 0
  
  ### calculate parameters needed for plotting
  tmin <- min(inftimes)
  tmax <- max(samtimes)
  tstep <- as.numeric(tmax-tmin)/2000
  tgmean <- if(is.na(tg.mean)) as.numeric(mean(inftimes - inftimes[infectors], na.rm = T)) else tg.mean
  tgvar <- if(is.na(tg.shape)) as.numeric(var(inftimes - inftimes[infectors], na.rm = T)) else NA
  tgscale <- if(is.na(tg.shape)) tgvar/tgmean else tgmean/tg.shape
  tgshape <- if(is.na(tg.shape)) tgmean/tgscale else tg.shape
  maxwd <- dgamma(max(0, tgscale * (tgshape - 1)), shape = tgshape, scale = tgscale)
  obs <- length(hosts)
  
  ### some smart graphical parameters
  if(is.null(label.cex)) label.cex <- max(0.5, min(1, 30/obs))
  if(is.null(sample.cex)) sample.cex <- label.cex
  if(is.null(arrow.length)) arrow.length <- max(0.005, min(0.1, 2.5/obs))
  if(is.null(sample.lwd)) sample.lwd <- min(2, 100/obs)

  ### initialize plot
  plot.new()
  par(cex = 1)
  plot.window(xlim = c(tmin, tmax + label.space * label.cex * as.numeric(tmax - tmin)), ylim = c(0, obs + 1))
  
  ### X-axis
  do.call(axis,
          c(list(side = 1,
                 at = c(((par("xaxp")[3] + 1) * par("xaxp")[1] - par("xaxp")[2]) / par("xaxp")[3], 
                        axTicks(1), 
                        ((par("xaxp")[3] + 1) * par("xaxp")[2] + par("xaxp")[1]) / par("xaxp")[3]),
                 mgp = par("mgp") * par("mar")[1]/4.1,
                 labels = if(inherits(tmin, "Date")) {
                   as.Date(c(((par("xaxp")[3] + 1) * par("xaxp")[1] - par("xaxp")[2]) / par("xaxp")[3],
                             axTicks(1), 
                             ((par("xaxp")[3] + 1) * par("xaxp")[2] + par("xaxp")[1]) / par("xaxp")[3]), "1970-01-01")  
                 } else TRUE, 
                 cex.axis = axis.cex),
            graphicalparameters("axis", 1, ...)
            )
          )

  ### Axis title
  do.call(title,
          c(list(xlab = c(xlab, graphicalparameters("title", ...)),
                 line = par("mar")[1]*5/8),
            graphicalparameters("title", 1, ...)
            )
          )
  
  ### Polygons and labels
  for(i in 1:obs) {
    x0s <- seq(inftimes[i], tmax - tstep, tstep)
    widths <- abs(1 - (maxwd - dgamma(x0s - inftimes[i], shape = tgshape, scale = tgscale)) / maxwd)
    
    do.call(polygon,
            c(list(x = c(x0s, rev(x0s)),
                   y = plotrank[i] + 0.3 * c(widths, -rev(widths)),
                   col = polygon.col, 
                   border = polygon.border),
              graphicalparameters("polygon", 1, ...)
              )
            )
#     do.call(text,
#             c(list(x = max(samtimes) + 10 * tstep,
#                    y = plotrank[i],
#                    labels = hosts[i], 
#                    adj = label.adj, 
#                    cex = label.cex),
#               graphicalparameters("label", timedorder, ...)))
  }
  do.call(text,
          c(list(x = max(samtimes) + 10 * tstep,
                 y = plotrank,
                 labels = hosts, 
                 adj = label.adj, 
                 cex = label.cex),
            graphicalparameters("label", timedorder, ...)))
  
  
  ### Horizontal lines
  do.call(segments,
          c(list(x0 = inftimes,
                 y0 = plotrank,
                 x1 = max(samtimes),
                 lty = line.lty),
            graphicalparameters("line", 1, ...)
            )
          )
  
  ### Arrows
  do.call(arrows,
          c(list(x0 = inftimes[infectors != "index"],
                 y0 = plotrank[infectors != "index"],
                 y1 = infectorpos[infectors != "index"], 
                 lwd = arrow.lwd, 
                 length = arrow.length, 
                 col = arrow.colours[infectors != "index"], 
                 code = 1),
            graphicalparameters("arrow", 1, ...)
            )
          )
  do.call(arrows,
          c(list(x0 = inftimes[infectors == "index"] - (par("xaxp")[2] - par("xaxp")[1])/1000,
                 y0 = plotrank[infectors == "index"],
                 x1 = inftimes[infectors == "index"],
                 lwd = arrow.lwd, 
                 length = arrow.length, 
                 col = arrow.colours[infectors == "index"], 
                 code = 2),
            graphicalparameters("arrow", 1, ...)
            )
          )
  
  ### Samples
  do.call(points,
          c(list(x = samtimes, 
                 y = plotrank[match(names(samtimes), hosts)],
                 pch = sample.pch, 
                 lwd = sample.lwd, 
                 cex = sample.cex),
            graphicalparameters("sample", 1, ...)
            )
          )
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
  plotrank <- 1
  for(i in 2:Nhosts) {
    ior <- which(hosts == infectors[i])
    if(infectors[ior] == "index") {
      plotrank[i] <- plotrank[ior] + insideYN[hosts[i]] - 0.5
    } else if(plotrank[ior] == 1) {
      plotrank[i] <- insideYN[hosts[i]] + 0.5
    } else if (plotrank[ior] == i - 1) {
      plotrank[i] <- i - 0.5 - insideYN[hosts[i]]
    } else {
      iorior <- which(hosts == infectors[ior])
      aboveYN <- xor((plotrank[iorior] < plotrank[ior]), insideYN[hosts[i]])
      plotrank[i] <- plotrank[ior] + aboveYN - 0.5    
    }
    plotrank[1:i] <- rank(plotrank)[1:i]
  }
  
#   plotrank <- 1:2
#   for(i in 3:Nhosts) {
#     ior <- which(hosts == infectors[i])
#     if(plotrank[ior] == 1) {
#       plotrank[i] <- insideYN[hosts[i]] + 0.5
#     } else if (plotrank[ior] == i - 1) {
#       plotrank[i] <- i - 0.5 - insideYN[hosts[i]]
#     } else {
#       iorior <- which(hosts == infectors[ior])
#       aboveYN <- xor((plotrank[iorior] < plotrank[ior]), insideYN[hosts[i]])
#       if(infectors[ior] == "index") aboveYN <- insideYN[hosts[i]]
#       plotrank[i] <- plotrank[ior] + aboveYN - 0.5    }
#     plotrank[1:i] <- rank(plotrank)[1:i]
#   }
  
  return(plotrank)
}


graphicalparameters <- function(which, timedorder, ...) {
  res <- list(...)
  res <- res[grep(paste0(which, "."), names(res))]
  res <- lapply(res, function(x) rep_len(x, length(timedorder))[timedorder])
  names(res) <- substring(names(res), nchar(which) + 2)
  return(res)
}


