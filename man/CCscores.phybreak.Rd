\name{Transmission clade credibility scores}
\alias{CCscores.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate transmission clade credibility scores per posterior (tree) sample.
%%  ~~function to do ... ~~
}
\description{
For each posterior sample in phybreak, the clade credibility score (log-transformed) is calculated. This score is the product of frequencies of each clade (= host + all progeny) among all posterior samples.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CCscores.phybreak(phybreak.object, samplesize = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object containing posterior samples.
%%     ~~Describe \code{phybreak.object} here~~
}
  \item{samplesize}{
  the number of samples to use, taken from the tail of the chain.
%%     ~~Describe \code{samplesize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector of length 'samplesize' with log-transformed clade credibility scores per posterior sample.
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

CCscores.phybreak(MCMCstate)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
