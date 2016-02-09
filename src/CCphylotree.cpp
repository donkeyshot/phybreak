#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

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



// [[Rcpp::export(name=".CCphylotree")]]
std::vector<double> CCphylotree(const std::vector<int> &pars, 
                                const std::vector<double> &tims,
                               std::vector<int> dims) {

  
  //The ss*n vector with n parents in each of ss trees is transformed into
  //a n_clade*n uniquecladearray, an ss*(n-1) vector indicating which n-1 clades
  //are in each of the ss trees, and an n_clade long vector with clade frequencies
  std::vector<bool> uniq_clades;
  std::vector<int> which_clades((dims[0] - 1) * dims[1]);
  std::vector<int> clade_freqs;
  int n_clade = 0;
  buildphyloclades(uniq_clades, clade_freqs, which_clades, n_clade,
              pars, dims[0], dims[1]);
  
  
  //For each clade, its log(frequency)
  std::vector<double> cladescores((dims[0] - 1) * dims[1]);
  for(int i = 0; i < (dims[0]-1)*dims[1]; ++i) {
    cladescores[i] = std::log(clade_freqs[which_clades[i]]);
  }

  //Determine the posterior tree with highest sum(log(frequency))
  int posttree = 0;
  std::vector<double> treescores(dims[1]);
  for(int i = 0; i < dims[1]; ++i) {
    for(int j = 0; j < dims[0] - 1; ++j) {
      treescores[i] += cladescores[i * (dims[0] - 1) + j];
    }
    if(treescores[i] > treescores[posttree]) {
      posttree = i;
    }
  }
  
  //Determine the posterior clades
  std::vector<int> postclades(dims[0] - 1);
  for(int i = 0; i < dims[0] - 1; ++i) {
    postclades[i] = which_clades[posttree * (dims[0] - 1) + i];
  }
  
  //For each posterior clade, the sums and sums-of-squares of the infection
  //times of the root hosts are calculated
  std::vector<double> cladetimesums(dims[0] - 1);
  std::vector<double> cladetimesumsqs(dims[0] - 1);
  phylocladetimestats(cladetimesums, cladetimesumsqs, tims,
                 postclades, which_clades, uniq_clades, dims[0], dims[1]);
  

  //Make vector with results: parents, clade support,
  //mean node time, SD(node time)
  //tree infection time,
  std::vector<double> result(5 * (2*dims[0] - 1));
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = pars[(2*dims[0]-1) * posttree + i];
    result[i + 2*dims[0] - 1] = dims[1];
    result[i + 2* (2*dims[0] - 1)] = 0;
    result[i + 3* (2*dims[0] - 1)] = 0;
    result[i + 4* (2*dims[0] - 1)] = 0;
  }
  for(int i = 0; i < dims[0] - 1; ++i) {
    result[i + dims[0]] = pars[(2*dims[0]-1) * posttree + dims[0] + i];
    result[i + dims[0] + 2*dims[0] - 1] = clade_freqs[postclades[i]];
    result[i + dims[0] + 2*(2*dims[0] - 1)] = (
      cladetimesums[i] / result[i + dims[0] + 2*dims[0] - 1]);
    result[i + dims[0] + 3*(2*dims[0] - 1)] = std::sqrt((
      cladetimesumsqs[i] -
        result[i + dims[0] + 2*(2*dims[0] - 1)] * 
        result[i + dims[0] + 2*(2*dims[0] - 1)] * 
        result[i + dims[0] + 2*dims[0] - 1]
    ) / (result[i + dims[0] + 2*dims[0] - 1] - 1) );
    result[i + dims[0] + 4*(2*dims[0] - 1)] = tims[(dims[0]-1) * posttree + i];
    
  }
      
      


  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

