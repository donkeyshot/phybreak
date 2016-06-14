\name{Posterior updating without sampling}
\alias{burnin.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC updating of a phybreak-object without sampling.
}
\description{
This function allows the MCMC chain to burn in. If used after samples have been taken, these samples will be returned unchanged in the output.
}
\usage{
burnin.phybreak(phybreak.object, ncycles, keepphylo = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
an object of class 'phybreak'.
}
  \item{ncycles}{
number of iterations to be carried out.
}
  \item{keepphylo}{
proportion of transmission tree updates keeping the phylotree intact.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
the phybreak-object provided as input, with variables and parameters changed due to the updating.
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
#First create a phybreak-object
simulation <- sim.phybreak(20)
MCMCstate <- phybreak(simulation)

MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
