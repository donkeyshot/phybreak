\name{Outbreak simulation}
\alias{sim.phybreak.gentime}
\alias{sim.phybreak.SIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating outbreaks in package "phybreak".
}
\description{
In sim.phybreak.gentime, gamma distributions are used for generation times and for times between infection and sampling. In sim.phybreak.SIR, a gamma distribution is used for the infectious period duration; infectiousness is constant during the infectious period, and sampling is done at the end of the infectious period.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.phybreak.gentime(obsize = 50, R0 = 1.5, shape.gen = 3, mean.gen = 1, shape.sample = 3, mean.sample = 1, nSUP = 0, multSUP = 10, disp = Inf, wh.model = 3, lambda = 0, rate0 = 1, mutrate = 1, nr.sites = 10000)

sim.phybreak.SIR(obsize = 50, R0 = 1.5, shape.infper = 3, mean.infper = 1, nSUP = 0, multSUP = 10, wh.model = 3, lambda = 0, rate0 = 1, mutrate = 1, nr.sites = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsize}{
the outbreak size (number of cases) to obtain
}
  \item{R0}{
the basic reproduction ratio used for simulation, not including the superspreaders. The offspring distribution is Poisson, unless 'disp' < Inf.
}
  \item{shape.gen}{
the shape parameter of the gamma-distributed generation time.
}
  \item{mean.gen}{
the mean generation time.
}
  \item{shape.sample}{
the shape parameter of the gamma-distributed sampling time, i.e. time from infection to sampling.
}
  \item{mean.sample}{
the mean sampling time.
}
  \item{shape.infper}{
the shape parameter of the gamma-distributed infectious period.
}
  \item{mean.infper}{
the mean infectious period.
}
  \item{nSUP}{
the number of superspreaders in the population (NOTE: not necessarily in the outbreak!)
}
  \item{multSUP}{
the proportional infectivity of the superspreaders.
}
  \item{disp}{
the dispersion parameter of the offspring distribution. The individual expected offspring is drawn from a gamma distribution with mean 'R0' and shape parameter 'disp' (or mean 'R0*multSUP' in case of a superspreader).
}
  \item{wh.model}{
the model for within-host pathogen dynamics (effective pathogen population size = N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Options are:
\enumerate{
\item{effective size = 0, so coalescence occurs 'just before' transmission in the infector}
\item{effective size = Inf, so coalescence occurs 'just after' transmission in the infectee}
\item{effective size at time t after infection = rate0 * exp(lambda * t)}
}
}
  \item{lambda}{
within-host exponential growth rate, used if wh.model = 3
}
  \item{rate0}{
within-host initial effective pathogen population size, used if wh.model = 3
}
  \item{mutrate}{
expected number of mutations per unit of time along each branch
}
  \item{nr.sites}{
number of available sites for mutations, so that the per-site mutation rate mu = mutrate/nr.sites
}
}
\details{
The population size is calculated by solving N from the final size equation << S/N == exp(-Rmean * (1 - S/N)) >> (rounded upwards), in which << Rmean = R0 * (1 + nSUP*(multSUP-1)/N) >>, i.e. the mean basic reproduction ratio.
}
\value{
An object of class 'obkData' [package "OutbreakTools"], containing the following slots:
\itemize{
  \item{@individuals:}{ data frame with two vectors, position i referring to host i
  \itemize{
    \item{$infector:}{ who infected host i?}
    \item{$date:}{ when was host i infected, with outbreak starting on 01-01-2000}
}}
  \item{@records:}{ NULL}
  \item{@dna:}{ object of class 'obkSequences' [package \pkg{OutbreakTools}], with slots
  \itemize{
    \item{@dna:}{ List containing
    \itemize{
      \item{$SNPs:}{ sequence data of class 'DNAbin'}}}
    \item{@meta:}{ data frame with three vectors, position i referring to host i
    \itemize{
      \item{$individualID}
      \item{$date:}{ sampling dates}
      \item{$locus:}{ all 1s}
}}
}}
  \item{@contacts:}{ NULL}
  \item{@context:}{ NULL}
  \item{@trees:}{ object of class 'multiPhylo', containing a single tree of class 'phylo' [package \pkg{ape}]}
}

%% ...
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
simulation <- sim.phybreak.gentime()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
