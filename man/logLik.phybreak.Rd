\name{Log-likelihood calculation}
\alias{logLik.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-likelihood calculation for a phybreak-object.
}
\description{
The likelihood of a phybreak-object is calculated, with the option to include or exclude parts of the likelihood for genetic data, generation times and sampling times.
}
\usage{
logLik.phybreak(phybreak.obj, seq = TRUE, t.gen = FALSE, t.sample = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.obj}{
a phybreak-object.
}
  \item{seq}{
whether to include the likelihood of the sequence data.
}
  \item{t.gen}{
whether to include the likelihood of the generation times
}
  \item{t.sample}{
whether to include the likelihood of the sampling times
}
}
\details{
The sequence likelihood is calculated by Felsenstein's pruning algorithm, assuming a prior probability of 0.25 for each nucleotide. The generation time and sampling time likelihood are log-densities of the gamma distributions for these variables.
}
\value{
The log-likelihood.
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
#First create an MCMC-object
simulation <- sim.phybreak.gentime()
MCMCstate <- make.phybreak.obkData(simulation)

logLik.phybreak(MCMCstate)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
