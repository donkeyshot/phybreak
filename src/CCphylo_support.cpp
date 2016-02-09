#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
#include "CCphylo_support.h"
using namespace Rcpp;

std::vector<int> ptrootphylo(std::vector<int> pars, int ID) {
  std::vector<int> ans;
  int nextID = pars[ID - 1];
  while(nextID > 0) {
    ans.push_back(nextID);
    nextID = pars[nextID - 1];
  }
  return(ans);
}

void makephylocladearray(std::vector<bool> &ca, const std::vector<int> &pars, 
                         int n, int ss) {
  std::vector<int> parset(2*n-1);
  std::vector<int> pos;
  
  for(int i = 0; i < ss; ++i) {
    for(int ii = 0; ii< 2*n - 1; ++ii) {
      parset[ii] = pars[i*(2*n-1) + ii];
    }
    for(int j = 0; j < n; ++j) {
      pos = ptrootphylo(parset, j + 1);
      for(unsigned int k = 0; k < pos.size(); ++k) {
        ca[i*n*(n-1) + (pos[k]-1-n)*n + j] = true;
      }
    }
  }
}

void tallyphyloclades(std::vector<int> &pos, std::vector<int> &count,
                      std::vector<int> &refs, const std::vector<bool> &ca,
                      int n, int ss) {
  bool cladeequal = false;
  bool nodeequal = true;
  unsigned int nextclade = 0;
  int nextnode = 0;
  
  pos.push_back(0);
  count.push_back(1);
  
  for(int i = 1; i < (n - 1) * ss; ++i) {
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

void buildphyloclades(std::vector<bool> &uc, std::vector<int> &cf, std::vector<int> &wc, 
                      int &nc, const std::vector<int> &pars, int n, int ss) {
  //make an array with all clades for all trees
  std::vector<bool> cladearray(n * (n - 1) * ss);
  makephylocladearray(cladearray, pars, n, ss);
  
  //find positions of unique clades in the cladearray, plus clade frequency
  std::vector<int> uniquecladepos;
  tallyphyloclades(uniquecladepos, cf, wc, cladearray, n, ss);
  
  //build the array with unique clades
  nc = cf.size();
  uc.resize(nc * n);
  for(int i = 0; i < nc; ++i) {
    for(int j = 0; j < n; ++j) {
      uc[n * i + j] = cladearray[n * uniquecladepos[i] + j];
    }
  }
}

void phylocladetimestats(std::vector<double> &sums, std::vector<double> &ssq,
                         const std::vector<double> &tims, 
                         const std::vector<int> &pc, 
                         const std::vector<int> &wc,
                         const std::vector<bool> &uc, int n, int ss) {
  
  for(int i = 0; i < n - 1; ++i) { //each clade in best tree
    for(int j = 0; j < ss; ++j) { //each post tree
      for(int k = 0; k < n - 1; ++k) { //each clade in post tree
        if(pc[i] == wc[j * (n - 1) + k]) {          //if clades are equal
          sums[i] += tims[j * (n - 1) + k];
          ssq[i] += tims[j * (n - 1) + k] * tims[j * (n - 1) + k];
        }
      }
    }
  }
}


//pc = postclade, will contain n positions of clades in 'uc' that are accepted 
//uc = array with all unique clades, with each n bool to indicate the hosts
//op = ordered pairs of clade indicators (positions in 'uc') and their frequencies (negative) among all posterior trees
//n = outbreak size
void selectphyloclades(std::vector<int> &pc, const std::vector<bool> &uc, 
                       const std::vector<std::pair<int, int> > &op, int n) {
  
  
  //the two most frequent clades are added to pc
  pc[0] = op[0].second;
  pc[1] = op[1].second;
  
  
  
  //'nextclade' by 'nextclade' from 'op' is tested for compatibility with
  //accepted clades in 'pc'. If accepted, it is added to pc.
  int nextclade = 2;  //next clade to test
  int candclade = 0;  //the ID (in 'pc') of the candidate clade
  int cladecount = 2;  //nr of clades accepted thus far
  bool cladeaccept;  //whether to accept the clade undergoing testing
  bool overlap, candmore, accmore;  //do clades overlap, 
  //does the candidate or accepted clade have non-overlapping hosts? 
  
  while(cladecount < n - 1) {
    candclade = op[nextclade].second;
    cladeaccept = true;
    for(int i = 0; i < cladecount && cladeaccept; ++i) { //testing against all accepted clades
      overlap = false;
      candmore = false;
      accmore = false;
      for(int j = 0; j < n; ++j) {
        if(uc[candclade * n + j] && uc[pc[i] * n + j]) {overlap = true;}
        if(uc[candclade * n + j] && !uc[pc[i] * n + j]) {candmore = true;}
        if(!uc[candclade * n + j] && uc[pc[i] * n + j]) {accmore = true;}
      }
      if(overlap && candmore && accmore) {  //candidate clade not nested -> not compatible
        cladeaccept = false;
      } else {cladeaccept = true;}   //candidate clade and accepted clade are not nested
    }
    if(cladeaccept) {
      pc[cladecount] = candclade;  //add candidate to accepted clades
      ++cladecount;
    }
    ++nextclade;
  }
}


//Function to find the node parents given the bestcladetree parents, scores, bestcladetree, besttreescores, cladeparents
void findnodeparents(std::vector<int> &pars, const std::vector<int> &pc, 
                     const std::vector<bool> &uc, int n) {
  std::vector<int> cladesizes(n);
  for(int i = 0; i < n - 1; ++i) {
    for(int j = 0; j < n; ++j) {
      cladesizes[i + 1] += uc[pc[i] * n + j];
    }
  }
  cladesizes[0] = n + 1;
  
  int parclade = 0;
  bool cladeinclade = false;
  for(int i = 0; i < n; ++i) { //all sample nodes
    parclade = 0;
    for(int j = 0; j < n - 1; ++j) {  //for all clades
      //test if sample is in clade and if clade is smallest thus far
      if(uc[pc[j] * n + i] && (cladesizes[j + 1] < cladesizes[parclade]) ) {
        parclade = j + 1;
      }
    }
    pars[i] = parclade + n;
  }
  for(int i = 0; i < n - 1; ++i) {  //all internal nodes i
    parclade = 0;
    for(int j = 0; j < n - 1; ++j) {  //for all clades j
      //test if internal node i is at base of smaller clade than j
      if(cladesizes[i + 1] < cladesizes[j + 1]) {
        cladeinclade = true;
        for(int k = 0; k < n; ++k) {  //for all sample nodes k
          //test if node k is only in smaller clade i
          if(!uc[pc[j] * n + k] && uc[pc[i] * n + k]) {
            cladeinclade = false;
          }
        }
        //test if enclosing clade j is smallest thus far
        if(cladeinclade && cladesizes[j + 1] < cladesizes[parclade]) {
          parclade = j + 1;
        }
      }
    }
    if(parclade > 0) {
      pars[i + n] = parclade + n;
    }
  }
  
}
