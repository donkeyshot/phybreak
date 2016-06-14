\name{Summarizing posterior trees}
\alias{transtree}
\alias{phylotree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyse posterior samples of the transmission and phylogenetic trees.%%  ~~function to do ... ~~
}
\description{
Various methods to create summary transmission trees or phylogenetic trees from phybreak mcmc-output.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transtree(phybreak.object, method = c("count", "edmonds", "mpc", "mtcc"), 
                    samplesize = Inf, infector.name = TRUE, support = "proportion",
                    infection.times = "all", time.quantiles = c(.025, 0.5, 0.975),
                    phylo.class = FALSE)

phylotree(phybreak.object,
                    samplesize = Inf, support = "proportion", clade.times = TRUE,
                    time.quantiles = c(.025, 0.5, 0.975), phylo.class = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  an object of class phybreak.
%%     ~~Describe \code{x} here~~
}
  \item{method}{
  the method used to create the tree (see details).
%%     ~~Describe \code{x} here~~
}
  \item{samplesize}{
  the number of posterior samples that is used, taken from the tail.
%%     ~~Describe \code{x} here~~
}
  \item{infector.name}{
  whether to return the infector names, or the position in the vector of hosts.
%%     ~~Describe \code{x} here~~
}
  \item{support}{
  whether to return the support (= posterior probability) for each infector as a "proportion" or as a "count" of posterior trees in which that transmission link, transmission cluster, or clade is present.
%%     ~~Describe \code{x} here~~
}
  \item{infection.times}{
  whether to base the summary infection times on "all" samples, or only on samples in which the "infector" was the same as in the consensus tree. A third option is "infector.sd", in which case only posterior trees with the same infector or transmission cluster as in the consensus tree were used to calculate a mean and standard deviation. In that case, also the infection times in the first sampled tree with consensus tree topology are given.
%%     ~~Describe \code{x} here~~
}
  \item{clade.times}{
  logical whether to include summary infection times, based only on posterior trees with the same clade as in the consensus tree were used to calculate a mean and standard deviation.
%%     ~~Describe \code{x} here~~
}
  \item{time.quantiles}{
  used in transtree only if infection.times is "all" or "infector".
%%     ~~Describe \code{x} here~~
}
  \item{phylo.class}{
  whether to return an object of class 'phylo', in which case a single tree ("mpc" or "mtcc" with transtree, "mcc" with phylotree) from the posterior is returned (not with summary infection times).
%%     ~~Describe \code{x} here~~
}
}
\details{
Five methods are supported for transmission tree reconstruction (who infected whom). 
- "count" gives the most frequent infector from the posterior samples. Note that this may result in an improper tree (multiple roots and/or cycles). Support is measured by the frequency of the infector in the posterior distribution.
- "edmonds" starts from the most frequent infector (method "count"), multiple roots and cycles are removed by selecting one by one the next most frequent option that minimizes the loss in support. Support is measured by the frequency of the infector in the posterior distribution.
- "mpc" gives the maximum parent credibility tree as described in Hall et al (2015). This is the tree in the posterior with maximum support = product of frequencies among all posterior samples. Support is measured by the frequency of the infector in the posterior distribution.
- "mtcc" gives the maximum transmission cluster credibility tree. This is equivalent to the mcc phylogenetic tree, with clades defined as host + all progeny. Support is measured by the frequency of the clade in the posterior distribution.

One method is supported for phylogenetic tree reconstruction.
- "mcc" gives the maximum clade credibility tree.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
transtree: A dataframe with per item (=host) its infector and support per infector/cluster, and summary infection times
phylotree: A dataframe with per item (=node) its parent and support per clade, and optionally summary node times. The first n nodes are the samples, the last n-1 nodes are the internal nodes.
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

transtree(MCMCstate, method = "edmonds")
transtree(MCMCstate, method = "mpc", infection.times = "infector.sd")
phylotree(MCMCstate)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
