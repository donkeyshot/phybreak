\name{phybreak infector-based output analysis}
\alias{MLinfectors}
\alias{MLtree.infector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyse posterior samples of infectors.
%%  ~~function to do ... ~~
}
\description{
MLinfectors returns the most likely infector for each host. Because MLinfectors may not be consistent with a proper tree (multiple roots or cycles), MLtree.infector returns a proper transmission tree based on posterior probabilities of infectors using Edmunds's algorithm to remove inconsistencies.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MLinfectors(phybreak.object, support = FALSE, infectornames = TRUE, samplesize = Inf)

MLtree.infector(phybreak.object, infectornames = TRUE, samplesize = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object.
%%     ~~Describe \code{phybreak.object} here~~
}
  \item{support}{
  whether to return the support (=posterior probability) for each infector.
%%     ~~Describe \code{support} here~~
}
  \item{infectornames}{
  whether to return the infector names, or just the number (=place in the vector of hosts).
%%     ~~Describe \code{infectornames} here~~
}
  \item{samplesize}{
  the number of posterior samples that is used, taken from the tail.
%%     ~~Describe \code{samplesize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe with an infector for each host, and posterior support.
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

MLinfector(MCMCstate, support = TRUE)
MLtree.infector(MCMCstate)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
