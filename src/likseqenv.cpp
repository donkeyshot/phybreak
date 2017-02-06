#include <Rcpp.h>
using namespace Rcpp;

// Updating the log-likelihood and "likarray" of sequence data
// given an existing "likarray" for a similar tree based on the same
// data (provided in environment "pbenv"). "likarray" contains the
// Felsenstein algorithm results for each node.
// Apart from the environment with current "likarray", the function
// requires the coalescent nodes that require updating,
// and the 'tips' = coalescent and sampling nodes from
// which to start the calculation.
// The function returns the log-likelihood, but more importantly
// changes "likarray" within the environment

// [[Rcpp::export(name=".likseqenv")]]
double likseqenv(Environment pbenv, 
                 IntegerVector nodestochange, IntegerVector tips) {
  List data = pbenv["d"];
  List vars = pbenv["v"];
  List pars = pbenv["p"];
  NumericVector likarray = pbenv["likarray"];
  IntegerVector SNPfreqs = pbenv["likarrayfreq"];
  IntegerVector SNPgenes = pbenv["likarraygene"];
  IntegerVector nodeparents = vars["nodeparents"];
  NumericVector nodetimes = vars["nodetimes"];
  int ngenes = data["ngenes"];
  double mu = pars["mu"];
  int obs = pars["obs"];
  int nSNPs = SNPfreqs.size();
  NumericVector nlens(nodetimes.size());
  NumericVector SNPsums(nSNPs);
  LogicalVector routefree(nodetimes.size());
  int curnode;
  int nextnode;
  int rootnode = nodetimes.size() * ngenes;
  double edgelen;
  double totprob;
  double result;
  
  for(unsigned int i = 0; i < nodestochange.size(); ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      for(int k = 0; k < 4; ++k) {
        likarray[( (nodestochange[i] - 1) * nSNPs + j) * 4 + k] = 1;
      }
    }
  }
  
  
  for(int i = 0; i < 3*obs - 1; ++i) {
    for(int j = 0; j < ngenes; ++j) {
      if(nodeparents[ngenes * i + j] > 0) {
        nlens[ngenes * i + j] = nodetimes[ngenes * i + j] - nodetimes[ngenes * (nodeparents[ngenes * i + j] - 1) + j];
      } else {
        rootnode = i;
      }
    }
  }
  nextnode = 0;
  while(nodeparents[nextnode * ngenes] - 1 != rootnode) {
    ++nextnode;
  }
  rootnode = nextnode;
  for(int i = 0; i < routefree.size(); ++i) {
    routefree[i] = true;
  }
  for(int i = 0; i < nodestochange.size(); ++i) {
    for(int j = 0; j < ngenes; ++j) {
      routefree[ngenes * (nodestochange[i] - 1) + j] = false;
    }
  }
  
  for(int i = 0; i < tips.size(); ++i) {
    for(int j = 0; j < ngenes; ++j) {
      curnode = tips[i] - 1;
      nextnode = nodeparents[ngenes * curnode + j] - 1;
      edgelen = nlens[ngenes * curnode + j];
      while(routefree[ngenes * curnode + j] && nextnode != -1) {
        if(nextnode < 2*obs - 1) {
          for(int k = 0; k < nSNPs; ++k) {
            if(SNPgenes[k] - 1 == j) {
              totprob = likarray[curnode*nSNPs*4 + k*4];
              for(int l = 1; l < 4; ++l) {
                totprob += likarray[curnode*nSNPs*4 + k*4 + l];
              }
              for(int l = 0; l < 4; ++l) {
                likarray[(nextnode * nSNPs + k) * 4 + l] *=
                  0.25*totprob + (likarray[(curnode * nSNPs + k) * 4 + l] -
                  0.25*totprob)*exp(-mu*edgelen);
              }
            }
          }
          curnode = nextnode;
          edgelen = nlens[ngenes * curnode + j];
          nextnode = nodeparents[ngenes * curnode + j] - 1;
        } else {
          edgelen += nlens[ngenes * nextnode + j];
          nextnode = nodeparents[ngenes * nextnode + j] - 1;
        }
      }
      routefree[ngenes * curnode + j] = true;
    }
  }

  for(int j = 0; j < nSNPs; ++j) {
    SNPsums[j] = 0.25*likarray[(rootnode * nSNPs + j) * 4];
    for(int k = 1; k < 4; ++k) {
      SNPsums[j] += 0.25*likarray[(rootnode * nSNPs + j) * 4 + k];
    }
    SNPsums[j] = log(SNPsums[j]) * SNPfreqs[j];
  }
  
  result = SNPsums[0];
  for(int j = 1; j < nSNPs; ++j) {
    result += SNPsums[j];
  }
  
  pbenv["likarray"] = likarray;
  pbenv["logLikseq"] = result;
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
/*** R
#phybreakenv <- new.env()
#phybreakenv.prop <- new.env()
#.build.phybreakenv(curstate)
#.propose.phybreakenv(curstate2)
*/

