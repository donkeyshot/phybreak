\name{Constructing a phybreak-object}
\alias{phybreak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a phybreak-object from data and prior distributions.
}
\description{
phybreak takes as data either an obkData-object (package OutbreakTools) or a matrix with sequences (individuals in rows, loci in columns). If only sequences are used, a vector with sampling times is also needed. Parameter values are used as initial values in the MCMC-chain or kept fixed. All variables are initialized by random samples from the prior distribution, unless the tree in the obkData is used (use.tree = TRUE).
}
\usage{
phybreak(data, ti = NULL, freq = NULL, id = NULL,
         mu = .0001, gen.shape = 3, gen.mean = 1,
         sample.shape = 3, sample.mean = 1, 
         wh.model = 3, wh.slope = 1,
         est.mean.gen = TRUE, prior.mean.gen.mean = 1, prior.mean.gen.sd = Inf,
         est.mean.sample = TRUE, prior.mean.sample.mean = 1, prior.mean.sample.sd = Inf,
         est.wh = TRUE, prior.wh.shape = 1, prior.wh.mean = 1,
         use.tree = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
an object of class 'obkData', containing sequences and sampling times as metadata with these sequences; or a matrix with nucleotides, each row a host, each column a locus. All nucleotides that are not "a", "c", "g", or "t", will be turned into "n".
}
  \item{ti}{
vector of sampling times, necessary if data is given as a matrix of sequences. If available, the item names will be used to identify the hosts, unless 'id' is given as well.
}
  \item{freq}{
vector of nucleotide frequencies, optional if data is given as a matrix of sequences. The vector indicates per column in 'data' how frequent that particular distribution of nucleotides across hosts was observed.
}
  \item{id}{
vector of host names to identify the hosts, optional if data is given as a matrix of sequences. The order should match the order of time and sequence data.
}
  \item{mu}{
initial value for mutation rate (defined per site per unit of time).
NOTE: mutation is defined as assignment of a random nucleotide at a particular site; this could be the nucleotide that was there before the mutation event. Therefore, the actual rate of change of nucleotides is 0.75*mu.
}
  \item{gen.shape}{
shape parameter of the generation time distribution (fixed).
}
  \item{gen.mean}{
initial value for the mean generation time.
}
  \item{sample.shape}{
shape parameter of the sampling time distribution (fixed).
}
  \item{sample.mean}{
initial value for mean time between infection and sampling.
}
  \item{wh.model}{
within-host model.
\enumerate{
  \item{indicates population size 0, coalescence 'just before' transmission}
  \item{indicates population size Inf, coalescence 'just after' transmission}
  \item{indicates a within-host population size since infection of slope * t}
}
}
  \item{wh.slope}{
initial value for parameter slope of within-host model 3.
}
  \item{est.mean.gen}{
whether to estimate the mean generation time or keep it fixed.
}
  \item{prior.mean.gen.mean}{
mean of the (gamma) prior distribution of gen.mean (only if est.mean.gen = TRUE).
}
  \item{prior.mean.gen.sd}{
standard deviation of the (gamma) prior distribution of gen.mean (only if est.mean.gen = TRUE).
}
  \item{est.mean.sample}{
whether to estimate the mean sampling time or keep it fixed.
}
  \item{prior.mean.sample.mean}{
mean of the (gamma) prior distribution of sample.mean (only if est.mean.sample = TRUE).
}
  \item{prior.mean.sample.sd}{
standard deviation of the (gamma) prior distribution of sample.mean (only if est.mean.sample = TRUE).
}
  \item{est.wh}{
whether to estimate the slope of the within-host model.
}
  \item{prior.wh.shape}{
shape parameter of the (gamma) prior distribution of slope (only if est.wh = TRUE).
}
  \item{prior.wh.mean}{
mean of the (gamma) prior distribution of slope (only if est.wh = TRUE).
}
  \item{use.tree}{
whether to use the transmission and phylogenetic tree given in data of class 'obkData', to create a phybreak-object with an exact copy of the outbreak. This requires more data in 'obk.object': the slot @individuals with vectors $infector and $date, and the slot @trees with at least one phylogenetic tree. Such data are created when using sim.phybreak(...).
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
\item{h }{containing 'helper' information for MCMC-updating: $si.mu and $si.wh for efficiently proposing mu and slope, matrix $dist with weights for infector sampling based on sequence distances, logicals $est.mG, $est.mS, and $est.wh whether to estimate gen.mean, sample.mean, and slope, and parameters for the priors of gen.mean, sample.mean, and slope.}
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
simulation <- sim.phybreak()
MCMCstate <- phybreak(simulation)
MCMCstate2 <- phybreak(simulation, use.tree = TRUE)

sampletimedata <- c(0,2,2,4,4)
sampleSNPdata <- matrix(c("a","a","a","a","a",
                          "a","c","c","c","c",
                          "t","t","t","g","g"), nrow = 5)
MCMCstate <- phybreak(sampleSNPdata, sampletimedata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
