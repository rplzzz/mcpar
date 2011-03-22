#include <iostream>
#include "mcpar.hh"
#include "rosenbrock.hh"
#include "mcout.hh"

int main(void)
{
  Gaussian L(2);
  MCout rslts(2);
  MCPar mcpar(2,4);             // 2 parameters, 4 chains

  //float pinit[8] = {0.0f,0.0f, 1.0f,1.0f, 0.0f,1.0f, 1.0f,0.0f};
  float pinit[8] = {0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f};

  mcpar.run(50,50, pinit, L, rslts);

  // output

  return 0;
}
