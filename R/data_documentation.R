#' Foot-and-Mouth Disease outbreak (2001)
#' 
#' DNA nucleotide patterns and sampling days of FMD virus from 15 farms in the UK. The nucleotides per sample 
#'   are not ordered as a virus sequence, but do contain the correct numbers of each nucleotide.
#' 
#' @format A matrix with 15 rows (one per farm) and 8196 column (one per nucleotide); the names contain sampling days
#' 
#' @references \href{http://dx.doi.org/10.1098/rspb.2007.1442}{Cottam et al. (2008)}
#'   Integrating genetic and epidemiological data to determine transmission pathways 
#'   of foot-and-mouth disease virus. \emph{Proc R Soc B}, \strong{275}: 887-895. 
#'   
"FMD_2001"

#' Foot-and-Mouth Disease outbreak (2007)
#' 
#' DNA sequences and sampling days of FMD virus from 10 farms in the UK (11 sequences)
#' 
#' @format A list with two elements:
#' \describe{
#'   \item{sequences}{the sequences (class \code{DNAbin})}
#'   \item{dates}{the sampling days (class \code{Date})}
#' }
#' 
#' @references \href{http://dx.doi.org/10.1371/journal.ppat.1000050}{Cottam et al. (2008)} Transmission pathways 
#'   of foot-and-mouth disease virus in the United Kingdom in 2007. \emph{PLoS Pathog}, \strong{4}(4): e1000050. 
"FMD_2007"

#' Mycobacterium tuberculosis outbreak
#' 
#' DNA SNP patterns and sampling days of an M. tuberculosis outbreak with 33 cases in Canada (Didelot et al, 2013). 
#'   The dataset contained only 20 SNP patterns with unspecified nucleotides, so the sequences in this dataset
#'    contain only `a` and `c` on these 20 loci. Additional loci with only `a` are added to a total sequence length
#'    of 440,000, which is 10% of the actual M. tuberculosis genome.
#' 
#' @format A list with two elements:
#' \describe{
#'   \item{sequences_Mtb}{the sequences (class \code{DNAbin})}
#'   \item{dates_Mtb}{the sampling days (class \code{Date})}
#' }
#' 
#' @references \href{http://dx.doi.org/10.1093/molbev/msu121}{Didelot et al. (2013)} Bayesian inference 
#'   of infectious disease transmission from whole-genome sequence data. \emph{Mol Biol Evol}, \strong{31}(7): 1869-1879. 
"M_tuberculosis_2013"

#' Avian influenza (H7N7) epidemic
#' 
#' DNA sequences (base order randomized) of NA, HA, and PB2 genes of samples taken from farms during the Dutch avian
#'   influenza (H7N7) epidemic in 2003. The dataset contains 241 farms, of which 231 were sampled and sequenced 
#'   (for the other 10, identified by the prefix UNK, only the detection day + 2 days is given, as in Klinkenberg 
#'   et al (2017)). The labels of the sequences refer to the labels in GISAID (\url{www.gisaid.org}), where the 
#'   sequences in correct base order can be found.
#'   
#' @format A list with four elements:
#' \describe{
#'   \item{sequences_HA}{the HA gene}
#'   \item{sequences_NA}{the NA gene}
#'   \item{sequences_PB2}{the PB2 gene}
#'   \item{sampledays}{the days at which the samples were taken (or the detection day + 2 days, for the 10 
#'   unsampled farms)}
#' }
#' 
#' @references \href{http://dx.doi.org/10.1371/journal.ppat.1002094}{Bataille et al. (2011)} Evolutionary 
#'   analysis of inter-farm transmission dynamics in a highly pathogenic avian influenza epidemic. 
#'   \emph{PLoS Pathog}, \strong{7}(6): e1002094. 
#'     
#'   \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'     inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'     \emph{PLoS Comput Biol}, \strong{13}(5): e1005495. 
"AvianFluH7N7_2003"