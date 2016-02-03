\name{Making coda output in phybreak}
\alias{make.coda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Making coda output from a phybreak-object.
%%  ~~function to do ... ~~
}
\description{
The function creates coda output from a phybreak-object to use functionality of the package pkg{coda}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.coda(phybreak.object, samplesize = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object.
%%     ~~Describe \code{phybreak.object} here~~
}
  \item{samplesize}{
  the number of samples to include (taken from the tail of the chain).
%%     ~~Describe \code{samplesize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
object of class 'mcmc' [package \pkg{coda}] with samples of all variables in the phybreak-object.
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

codasamples <- make.coda(MCMCstate)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
