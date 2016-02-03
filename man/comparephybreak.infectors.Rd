\name{Compare phybreak-output to (simulated) obkData}
\alias{comparephybreak.infectors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to assess the ability of the method to retrieve the (true simulated) data in an obkData-object.
%%  ~~function to do ... ~~
}
\description{
comparephybreak.infectors counts for posterior sample how many infectors are correct, using the obkData-object as a reference.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comparephybreak.infectors(phybreak.object, sim.object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object with posterior samples.
%%     ~~Describe \code{phybreak.object} here~~
}
  \item{sim.object}{
  an obkData-object containing a vector $infector in the @individuals slot.
%%     ~~Describe \code{sim.object} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector with counts of the numbers of equal infectors in each posterior sample and the obkData-object.
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
Maybe not to include in the final package.
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

comparephybreak.infectors(MCMCstate, simulation)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
