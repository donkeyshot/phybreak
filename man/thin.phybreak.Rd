\name{Thinning a phybreak-object}
\alias{thin.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove posterior samples from a phybreak-object.
}
\description{
The function removes (all or some) posterior samples by thinning and/or removing the first part of the chain, but keeps the state of variables and parameters intact.
}
\usage{
thin.phybreak(phybreak.object, thin = 1, nkeep = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
an object of class 'phybreak'.
}
  \item{thin}{
the thinning interval to use.
}
  \item{nkeep}{
the number of most recent samples to keep. If nkeep = Inf (default), the whole chain is thinned and returned. Otherwise, samples from the tail are kept.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The same phybreak-object as in the argument, but with fewer posterior samples.
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
#First build a phybreak-object with samples
simulation <- sim.phybreak(obsize = 20)
MCMCstate <- phybreak(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 50)
MCMCstate <- sample.phybreak(MCMCstate, nsample = 20, thin = 10)

MCMCstate <- thin.phybreak(MCMCstate, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
