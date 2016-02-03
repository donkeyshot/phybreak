\name{Accessing current state of phybreak-object}
\alias{get.infectors}
\alias{get.parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Access to the current state of a phybreak-object.
%%  ~~function to do ... ~~
}
\description{
Obtaining the current state infectors and parameters in a phybreak-object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.infectors(phybreak.object)

get.parameters(phybreak.object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object.
%%     ~~Describe \code{phybreak.object} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of infectors (position i referring to host i), or a named vector with all parameters.
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
#First build a phybreak-object.
simulation <- sim.phybreak.gentime(obsize = 20)
MCMCstate <- make.phybreak.obkData(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, nburnin = 200)

get.infectors(MCMCstate)
get.parameters(MCMCstate)["mu"]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
