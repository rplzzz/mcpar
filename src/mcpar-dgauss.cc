#include <iostream>
#include "mcpar.hh"
#include "rosenbrock.hh"
#include "mcout.hh"

int main(void)
{
  const int nparam=2;
  const float gmu[2] = {0.0f,0.0f};
  const float gsig2[2] = {1.0f, 2.0f};
  DualGaussian L(5.0f);

  MCout rslts(nparam);
  MCPar mcpar(nparam,4);        // 2 parameters, 4 chains

  float pinit[8] = {0.0f,0.0f, 2.0f,2.0f, 0.0f,1.5f, 0.0f,-2.0f};

  mcpar.run(100000,500, pinit, L, rslts);

    // output
  for(int i=0; i<rslts.size(); ++i) {
    const float *pset = rslts.getpset(i);
    for(int j=0; j<rslts.nparam(); ++j)
      std::cout << pset[j] << "\t";
    std::cout << "\n";
  }

  return 0;
}

