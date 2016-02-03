\name{Making a phybreak-object}
\alias{make.phybreak}
\alias{make.phybreak.obkData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a phybreak-object from data and model assumptions
}
\description{
make.phybreak takes sampling times and sequences and creates a phybreak-object for further analysis. Mutation rate and mean generation and sampling times will be used for initialization of the chain, the remaining parameter values will be fixed. All variables are initialized by random samples from the prior distribution.

make.phybreak.obkData creates a phybreak-object from an object of class 'obkData' [package \pkg{OutbreakTools}]. This object should contain sequence data in the 'dna' slot, including meta-data with individualID and sampling dates.
}
\usage{
make.phybreak(t.sample, SNP.sample, SNP.freqs = c(), hostnames = NULL, mu = .01, shape.gen = 3, mean.gen = 1, shape.sample = 3, mean.sample = 1, wh.model = 3, wh.lambda= 0, wh.rate0 = 1, est.mean.gen = TRUE, est.mean.sample = TRUE, prior.mean.gen.shape = 0, prior.mean.gen.scale = 0, prior.mean.sample.shape = 0, prior.mean.sample.scale = 0)

make.phybreak.obkData(obk.object, mu = .01, shape.gen = 3, mean.gen = 1, shape.sample = 3, mean.sample = 1, wh.model = 3, wh.lambda= 0, wh.rate0 = 1, est.mean.gen = TRUE, est.mean.sample = TRUE, prior.mean.gen.shape = 0, prior.mean.gen.scale = 0, prior.mean.sample.shape = 0, prior.mean.sample.scale = 0, use.tree = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.sample}{
vector of sampling times. If available, the item names will be used to identify the hosts, unless 'hostnames' is given as well.
}
  \item{SNP.sample}{
matrix of SNPs, each row a host, each vector a SNP.
}
  \item{SNP.freqs}{
vector of SNP-frequencies, indicating per column in 'SNP.sample' how frequent that particular distribution of nucleotides across hosts was observed.
}
  \item{hostnames}{
vector of host names to identify the hosts, used in the same order as the time and SNP data.
}
  \item{obk.object}{
object of class 'obkData'.
}
  \item{mu}{
initial value for mutation rate (defined per site per unit of time).
NOTE: mutation is defined as assignment of a random nucleotide at a particular site; this could be the nucleotide that was there before the mutation event. Therefore, the actual rate of change of nucleotides is 0.75*mu.
}
  \item{shape.gen}{
shape parameter of the generation time distribution (fixed).
}
  \item{mean.gen}{
initial value for the mean generation time.
}
  \item{shape.sample}{
shape parameter of the sampling time distribution (fixed).
}
  \item{mean.sample}{
initial value for mean time between infection and sampling.
}
  \item{wh.model}{
within-host model.
\enumerate{
  \item{indicates population size 0, coalescence 'just before' transmission}
  \item{indicates population size Inf, coalescence 'just after' transmission}
  \item{indicates a within-host population size since infection of rate0 * exp(lambda * t)}
}
}
  \item{wh.lambda}{
parameter lambda of within-host model 3.
}
  \item{wh.rate0}{
parameter rate0 of within-host model 3.
}
  \item{est.mean.gen}{
whether to estimate the mean generation time or keep it fixed.
}
  \item{est.mean.sample}{
whether to estimate the mean sampling time or keep it fixed.
}
  \item{prior.mean.gen.shape}{
shape parameter of the (gamma) prior distribution of mean.gen (only if est.mean.gen = TRUE).
}
  \item{prior.mean.gen.scale}{
scale parameter of the (gamma) prior distribution of mean.gen (only if est.mean.gen = TRUE).
}
  \item{prior.mean.sample.shape}{
shape parameter of the (gamma) prior distribution of mean.sample (only if est.mean.sample = TRUE).
}
  \item{prior.mean.sample.scale}{
scale parameter of the (gamma) prior distribution of mean.sample (only if est.mean.sample = TRUE).
}
  \item{use.tree}{
whether to use the transmission and phylogenetic tree given in 'obk.object', to create a phybreak-object with an exact copy of the outbreak. This requires more data in 'obk.object': the slot @individuals with vectors $infector and $date, and the slot @trees with at least one phylogenetic tree. Such data are created when using sim.phybreak.gentime(...) or sim.phybreak.SIR(...).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of lists is provided:
\item{d }{containing the data}
\item{v }{containing the variables specifying the phylogenetic and transmission tree. All vectors contain information on the sampling nodes (first n items), on the coalescent nodes (items n+1 to 2n-1), and on the transmission nodes (all further items):
\itemize{
  \item{$ nodetimes: times of the nodes.}
  \item{$ nodehosts: host IDs in which the nodes reside.}
  \item{$ nodeparents: parent node of each node.}
  \item{$ nodetypes: either "s"ampling, "c"oalescent, or "t"ransmission.}
}




}
\item{p }{containing the parameter values}
\item{h }{containing 'helper' information for MCMC-updating: vector $si with recent mu samples to be used for efficient sampling of mutation rate, matrix $dist with weights for infector sampling based on sequence distances, logicals $est.mG and $est.mS whether to estimate mean.gen and mean.sample, and shape and scale parameters for priors of mean.gen and mean.sample.}
\item{s }{will contain vectors and matrices for posterior samples; in matrices, the rows are nodes in the phylogenetic tree, the columns are the samples}
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
MCMCstate <- make.phybreak.obkData(simulation)
MCMCstate2 <- make.phybreak.obkData(simulation, use.tree = TRUE)

sampletimedata <- c(0,2,2,4,4)
sampleSNPdata <- matrix(c("a","a","a","a","a",
                          "a","c","c","c","c",
                          "t","t","t","g","g"), nrow = 5)
MCMCstate <- make.phybreak(sampletimedata, sampleSNPdata, wh.model = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
