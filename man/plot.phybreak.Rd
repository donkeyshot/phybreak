\name{plot.phybreak}
\alias{plot.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plotting a phybreak object
}
\description{
Plots a phybreak object as phylogenetic tree with coloured branches indicating hosts. The default is to plot the current state, but any posterior sample can be chosen, as well as various consensus trees
}
\usage{
plot.phybreak(x, plot.which = "sample", samplenr = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 object of class phybreak
}
  \item{plot.which}{
can be "sample" to plot the current state or a selected posterior sample, "mpc" or "mtcc" to plot a consensus transmission tree (see transtree) or "mcc" to plot a consensus phylogenetic tree (see phylotree).
}
  \item{samplenr}{
if plot.which = "sample", this indicates which posterior tree should be plotted; samplenr = 0 to plot the current state.
}
  \item{\dots}{
additional options for plot.simmap.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
no value returned, but a phylogenetic tree is plotted with coloured branches indicating hosts in which the lineages reside.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Don Klinkenberg
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#First build a phybreak-object.
simulation <- sim.phybreak(obsize = 20)
MCMCstate <- phybreak(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 200)
MCMCstate <- burnin.phybreak(MCMCstate, nsample = 20, thin = 10)

plot(MCMCstate)
plot(MCMCstate, "mpc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
