\name{Outbreak simulation}
\alias{sim.phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating outbreaks of class 'obkData' in package "phybreak".
}
\description{
Simulate outbreaks of class 'obkData' (package OutbreakTools), with the outbreak model of phybreak. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.phybreak(obsize = 50, popsize = NA, R0 = 1.5, shape.gen = 10, mean.gen = 1, shape.sample = 10, mean.sample = 1, wh.model = 3, slope = 1, mutrate = 0.0001, sequence.length = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsize}{
the outbreak size (number of cases) to obtain.
}
  \item{popsize}{
the population size in which to simulate. If it is not defined (default), an optimal population size will be chosen based on R0 and obsize.
}
  \item{R0}{
the basic reproduction ratio used for simulation. The offspring distribution is Poisson.
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
  \item{wh.model}{
the model for within-host pathogen dynamics (effective pathogen population size = N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Options are:
\enumerate{
\item{effective size = 0, so coalescence occurs 'just before' transmission in the infector}
\item{effective size = Inf, so coalescence occurs 'just after' transmission in the infectee}
\item{effective size at time t after infection = slope * t}
}
}
  \item{slope}{
within-host increase of population size, used if wh.model = 3.
}
  \item{mutrate}{
expected number of mutations per nucleotide per unit of time along each lineage.
}
  \item{sequence.length}{
number of available nucleotides for mutations.
}
}
\details{
If popsize = NA, the population size is calculated by solving N from the final size equation << S/N == exp(-R0 * (1 - S/N)) >> (rounded upwards), in which S/N is the proportion of the population not involved in the outbreak.
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
simulation <- sim.phybreak()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
