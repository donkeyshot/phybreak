\name{Sampling from a phybreak MCMC-chain.}
\alias{sample.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling from a phybreak MCMC-chain.
}
\description{
Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree (plus associated parameters), by use of a phybreak-object.
}
\usage{
sample.phybreak(phybreak.object, nsample, thin, keepphylo = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
an object of class 'phybreak'.
}
  \item{nsample}{
the number of samples to take.
}
  \item{thin}{
the thinning to use (values after every 'thin'th update cycle will be included in the posterior).
}
  \item{keepphylo}{
proportion of transmission tree updates keeping the phylotree intact.
}
}
\details{
For the algorithm used to update the trees, see Klinkenberg et al (in prep)
}
\value{
The same phybreak-object used to call the function, including (additional) samples from the posterior.
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
#First build a burnt-in phybreak-object
simulation <- sim.phybreak.gentime(obsize = 20)
MCMCstate <- make.phybreak.obkData(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, nburnin = 200)

MCMCstate <- sample.phybreak(MCMCstate, nsample = 100, thin = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
