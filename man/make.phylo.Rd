\name{Transform a phybreak-object.}
\alias{make.phylo.phybreak}
\alias{make.multiPhylo}
\alias{make.phyDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform phybreak-objects for further analysis.
%%  ~~function to do ... ~~
}
\description{
Create an object of class 'phylo' from the current state or any posterior sample, an object of class 'multiPhylo' from the posterior samples, or an object of class 'phyDat' from the data in a phybreak-object. With these objects, functionalities of packages \pkg{ape} and \pkg{phangorn} become available.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.phylo.phybreak(phybreak.object, samplenr = 0)

make.multiPhylo(phybreak.object, samplesize = Inf)

make.phyDat(phybreak.object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object.
%%     ~~Describe \code{phybreak.object} here~~
}
  \item{samplenr}{
  which sample to use (if 0, then use the current state).
%%     ~~Describe \code{samplenr} here~~
}
  \item{samplesize}{
  the number of posterior samples to use. Samples are taken from the tail of the MCMC-chain.
%%     ~~Describe \code{samplenr} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
objects of class 'phylo', 'multiPhylo', or 'phyDat'.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Don Klinkenberg
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#First build a phybreak-object containing samples.
simulation <- sim.phybreak.gentime(obsize = 20)
MCMCstate <- make.phybreak.obkData(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, nburnin = 200)
MCMCstate <- sample.phybreak(MCMCstate, nsample = 100, thin = 10)

tree0 <- make.phylo.phybreak(MCMCstate)
plot(tree0) #in "ape"

treepost <- make.multiPhylo(MCMCstate)
dataphyDat <- make.phyDat(MCMCstate)
pml(tree0, dataphyDat, rate = 0.75*get.parameters(MCMCstate)["mu"]) #in "phangorn"
logLik.phybreak(MCMCstate) #should give the same result as 'pml'
sapply(treepost, parsimony, data=dataphyDat) #in "phangorn"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
