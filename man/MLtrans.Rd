\name{Summarizing posterior trees}
\alias{MLtrans}
\alias{MLphylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyse posterior samples of the transmission and phylogenetic trees.%%  ~~function to do ... ~~
}
\description{
Various methods to create summary transmission trees or phylogenetic trees from phybreak mcmc-output.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MLtrans(phybreak.object, method = c("count", "edmunds", "mpc", "mcc", "cc.construct"), support = FALSE, infectornames = FALSE, samplesize = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phybreak.object}{
  a phybreak-object.
%%     ~~Describe \code{x} here~~
}
  \item{method}{
  the method used to create the tree (see details).
%%     ~~Describe \code{x} here~~
}
  \item{support}{
  whether to return the support (=posterior probability) for each infector.
%%     ~~Describe \code{x} here~~
}
  \item{infectornames}{
  whether to return the infector names, or just the number (=place in the vector of hosts).
%%     ~~Describe \code{x} here~~
}
  \item{samplesize}{
  the number of posterior samples that is used, taken from the tail.
%%     ~~Describe \code{x} here~~
}
}
\details{
Five methods are supported for transmission tree reconstruction (who infected whom). 
\item{count }{the most frequent infector from the posterior samples. Note that this may result in an improper tree (multiple roots and/or cycles). Support is measured by the frequency of the infector in the posterior distribution.}
\item{edmunds }{starting from the most frequent infector (method "count"), multiple roots and cycles are removed by selecting one by one the next most frequent option that minimizes the loss in support. Support is measured by the frequency of the infector in the posterior distribution.}
\item{mpc }{maximum parent credibility tree as described in Hall et al (2015). This is the tree in the posterior with maximum support = product of frequencies among all posterior samples. Support is measured by the frequency of the infector in the posterior distribution.}
\item{mtcc }{maximum transmission clade credibility tree as (will be) described in Klinkenberg et al (in prep). Equivalent to the mcc phylogenetic tree, with clades defined as host + all progeny. Support is measured by the frequency of the clade in the posterior distribution.}
\item{cc_construct }{constructed transmission clade credibility tree as (will be) described in Klinkenberg et al (in prep). Posterior transmission clades are combined to create a tree with maximal support, which may not be among the actually sampled trees. Support is measured by the frequency of the clade in the posterior distribution.}

Two methods are supported for phylogenetic tree reconstruction.
\item{mcc }{maximum clade credibility.}
\item{cc_construct }{constructed phylogenetic tree, as described above.}
%%  ~~ If necessary, more details than the description above ~~
}
\value{
MLtrans: A dataframe with per item (=host) its $infector and (optionally) $support per infector/clade
MLphylo: A dataframe with per item (=node) its $parent and (optionally) $support per clade. The first n nodes are the samples, the last n-1 nodes are the internal nodes.
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
simulation <- sim.phybreak.gentime(obsize = 20)
MCMCstate <- make.phybreak.obkData(simulation)
MCMCstate <- burnin.phybreak(MCMCstate, nburnin = 200)
MCMCstate <- sample.phybreak(MCMCstate, nsample = 100, thin = 10)

MLtrans(MCMCstate, support = TRUE)
MLtrans(MCMCstate, method = "mpc")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
