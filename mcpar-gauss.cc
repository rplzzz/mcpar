#include <iostream>
#include "mcpar.hh"
#include "rosenbrock.hh"
#include "mcout.hh"

int main(void)
{
  const int nparam=2;
  const float gmu[2] = {0.0f,0.0f};
  const float gsig2[2] = {1.0f, 2.0f};
  Gaussian L(nparam, gmu, gsig2);
  MCout rslts(nparam);
  //MCPar mcpar(nparam,4);             // 2 parameters, 4 chains
  //MCPar mcpar(2,4,1,0,1.0f);    // 2 param, 4 chain, no MPI, plocal=1.0
  MCPar mcpar(2,4,1,0,1.0e-9f);    // 2 param, 4 chain, no MPI, plocal=0.0

  //float pinit[8] = {0.0f,0.0f, 1.0f,1.0f, 0.0f,1.0f, 1.0f,0.0f};
  float pinit[8] = {0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f};

  mcpar.run(10000000,500, pinit, L, rslts);

  // output
  for(int i=0; i<rslts.size(); ++i) {
    const float *pset = rslts.getpset(i);
    for(int j=0; j<rslts.nparam(); ++j)
      std::cout << pset[j] << "\t";
    std::cout << "\n";
  }

  return 0;
}
