\name{Accessing current state of phybreak-object}
\alias{get.infectors}
\alias{get.parameters}
\alias{get.mcmc}
\alias{get.phylo}
\alias{get.phyDat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Access to the current state of a phybreak-object.
%%  ~~function to do ... ~~
}
\description{
Extracting specific elements of a phybreak-object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.infectors(phybreak.object)

get.parameters(phybreak.object)

get.mcmc(phybreak.object, thin = 1, samplesize = Inf)

get.phylo(phybreak.object, post.tree = 0, simmap = FALSE)

get.phyDat(phybreak.object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  an object of class 'phybreak'.
%%     ~~Describe \code{phybreak.object} here~~
}
  \item{thin}{
  thinning interval for making the mcmc-object.
}
  \item{samplesize}{
  number of samples to keep in the mcmc-object, counting from tail of the chain.
}
  \item{post.tree}{
  the posterior tree sample to choose. If post.tree = 0 (default), the current state is used.
}
  \item{simmap}{
  logical indicating whether to return simmap elements (package \pkg{phytools}), colouring the branches on the tree to indicate host.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{get.infectors }{returns a vector of infectors (position i referring to host i).}
\item{get.parameters }{returns a named vector with all current parameter values.}
\item{get.mcmc }{returns an object of class 'mcmc' (package \pkg{coda}), with sampled parameters, infection times, and infectors.}
\item{get.phylo }{returns an object of class 'phylo' (package \pkg{ape}) and optionally of class 'simmap' (package \pkg{phytools}), with the posterior tree sample as chosen (default = current state).}
\item{get.phyDat }{returns the sequence data in class 'phyDat' (package \pkg{phybreak}).}
%% ...
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

get.infectors(MCMCstate)
get.parameters(MCMCstate)
get.mcmc(MCMCstate)
plot(get.phylo(MCMCstate))
get.phyDat(MCMCstate)

#function from package phangorn:
parsimony(get.phylo(MCMCstate), get.phyDat(MCMCstate))

tree0 <- get.phylo(MCMCstate)
dataphyDat <- get.phyDat(MCMCstate)
pml(tree0, dataphyDat, rate = 0.75*get.parameters(MCMCstate)["mu"]) #in phangorn
logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
sampling = FALSE, generation = FALSE) #should give the same result as 'pml'

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
