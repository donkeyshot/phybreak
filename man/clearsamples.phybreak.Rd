\name{Clearing a phybreak-object}
\alias{clearsamples.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove posterior samples from a phybreak-object.
}
\description{
The function removes (all or some) posterior samples, but keeps the state of variables and parameters intact.
}
\usage{
clearsamples.phybreak(phybreak.object, clearall = TRUE, nkeep = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
a phybreak-object.
}
  \item{clearall}{
whether to remove all samples.
}
  \item{nkeep}{
the number of most recent samples to keep (if clearall = FALSE).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The same phybreak-object as in the argument, but with posterior samples removed.
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
simulation <- sim.phybreak.gentime(obsize = 20)
MCMCstate <- make.phybreak.obkData(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, nburnin = 50)
MCMCstate <- sample.phybreak(MCMCstate, nsample = 10, thin = 10)

MCMCstate <- clearsamples.phybreak(MCMCstate)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
