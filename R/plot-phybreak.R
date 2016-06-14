### plot functions for phybreak-object ###

### makes use of plot function in 'simmap',
### colouring the branches by host in which they reside,
### but no unique colour per host.
### calls:
# MLtrans
# get.phylo
plot.phybreak <- function(x, plot.which = "sample", samplenr = 0, ...) {
  if(!inherits(x, "phybreak")) {
    stop("object x must be of class \"phybreak\"")
  }

  if(plot.which == "mpc") {
    plotSimmap(transtree(x, "mpc", phylo.class = TRUE),
               colors = setNames(nm = c("black","red","blue")), ...)
  } else if(plot.which == "mtcc") {
    plotSimmap(transtree(x, "mtcc", phylo.class = TRUE),
               colors = setNames(nm = c("black","red","blue")), ...)
  } else if(plot.which == "mcc") {
    plotSimmap(phylotree(x, phylo.class = TRUE),
               colors = setNames(nm = c("black","red","blue")), ...)
  } else {
    plotSimmap(get.phylo(x, samplenr, TRUE),
               colors = setNames(nm = c("black","red","blue")), ...)
  }
}
