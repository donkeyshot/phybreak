#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
//#include "CCtrans_support.h"
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

std::vector<int> ptroottrans(std::vector<int> pars, int ID) {
  std::vector<int> ans;
  ans.push_back(ID);
  int nextID = pars[ID - 1];
  while(nextID > 0) {
    ans.push_back(nextID);
    nextID = pars[nextID - 1];
  }
  return(ans);
}

void makecladearray(std::vector<bool> &ca, const std::vector<int> &pars, 
                    int n, int ss) {
  std::vector<int> parset(n);
  std::vector<int> pos;
  
  for(int i = 0; i < ss; ++i) {
    for(int ii = 0; ii< n; ++ii) {
      parset[ii] = pars[i*n + ii];
    }
    for(int j = 0; j < n; ++j) {
      pos = ptroottrans(parset, j + 1);
      for(unsigned int k = 0; k < pos.size(); ++k) {
        ca[i*n*n + (pos[k]-1)*n + j] = true;
      }
    }
  }
}

void tallyclades(std::vector<int> &pos, std::vector<int> &count,
                 std::vector<int> &refs, const std::vector<bool> &ca,
                 int n, int ss) {
  bool cladeequal = false;
  bool nodeequal = true;
  unsigned int nextclade = 0;
  int nextnode = 0;
  
  pos.push_back(0);
  count.push_back(1);
  
  for(int i = 1; i < n * ss; ++i) {
    while(!cladeequal && nextclade < pos.size()) {
      while(nodeequal && nextnode < n) {
        if(ca[i*n + nextnode] !=
           ca[pos[nextclade]*n + nextnode]) {
          nodeequal = false;
        }
        ++nextnode;
      }
      cladeequal = nodeequal;
      nodeequal = true;
      ++nextclade;
      nextnode = 0;
    }
    if(!cladeequal) {
      pos.push_back(i);
      count.push_back(1);
      refs[i] = pos.size() - 1;
    } else {
      ++count[nextclade - 1];
      refs[i] = nextclade - 1;
    }
    cladeequal = false;
    nodeequal = true;
    nextclade = 0;
    nextnode = 0;
  }
}

void cladetimestats(std::vector<double> &sums, std::vector<double> &ssq,
                    const std::vector<double> &tims, 
                    const std::vector<int> &clIDs, 
                    const std::vector<int> &whclID,
                    const std::vector<bool> &ca, int n, int ss) {
  bool anytime = false;
  double mintime = 0;
  
  for(int i = 0; i < ss; ++i) { //each tree
    for(int j = 0; j < n; ++j) { //each clade
      for(int k = 0; k < n; ++k) { //each host
        if(ca[clIDs[whclID[i*n + j]] * n + k]) { //is host in clade?
          if(!anytime) { //is it first host in clade?
            mintime = tims[i*n + k]; //mintime = time of host
            anytime = true;
          } else if(mintime > tims[i*n + k]) { //is time of host smallest?
            mintime = tims[i*n + k]; //mintime = time of host
          }
        }
      }
      sums[whclID[i*n + j]] += mintime;
      ssq[whclID[i*n + j]] += mintime*mintime;
      anytime = false;
    }
  }
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

