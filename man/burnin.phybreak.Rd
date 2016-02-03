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
burnin.phybreak(phybreak.object, nburnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
a phybreak-object.
}
  \item{nburnin}{
number of iterations to be carried out.
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
simulation <- sim.phybreak.gentime()
MCMCstate <- make.phybreak.obkData(simulation)

MCMCstate <- burnin.phybreak(MCMCstate, nburnin = 50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
