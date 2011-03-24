#include <iostream>
#include "mcpar.hh"
#include "rosenbrock.hh"
#include "mcout.hh"

int main(void)
{
  Gaussian L(2);
  MCout rslts(2);
  MCPar mcpar(2,4);             // 2 parameters, 4 chains
  //MCPar mcpar(2,4,1,0,1.0f);    // 2 param, 4 chain, no MPI, plocal=0

  //float pinit[8] = {0.0f,0.0f, 1.0f,1.0f, 0.0f,1.0f, 1.0f,0.0f};
  float pinit[8] = {0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f};

  mcpar.run(500,50, pinit, L, rslts);

  // output

  return 0;
}
