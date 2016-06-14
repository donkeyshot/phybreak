\name{infectorsets}
\alias{infectorsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summarize sampled infectors for each host in a phybreak-object
}
\description{
The function returns per host a table with posterior infectors, with support per infector.
}
\usage{
infectorsets(phybreak.object, which.hosts = "all", percentile = 0.95, minsupport = 0, samplesize = Inf, infector.name = TRUE, support = "proportion")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
an object of class phybreak, that should contain posterior samples
}
  \item{which.hosts}{
a vector with hosts (positions in the dataset), or "all" for all hosts
}
  \item{percentile}{
return infectors ordered by support, until a cumulative support of 'percentile'
}
  \item{minsupport}{
only return infectors with more support than 'minsupport'
}
  \item{samplesize}{
the number of samples to include (taken from the tail of the chain)
}
  \item{infector.name}{
whether to return the names of the infectors, or their position in the dataset
}
  \item{support}{
give support in "proportion" or "count"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a named list with dataframes, each with a vector of infectors and a vector of supports
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
#First build a phybreak-object containing samples.
simulation <- sim.phybreak(obsize = 20)
MCMCstate <- phybreak(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 200)
MCMCstate <- sample.phybreak(MCMCstate, nsample = 100, thin = 10)

infectorsets(MCMCstate)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
