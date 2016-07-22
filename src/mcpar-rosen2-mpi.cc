#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "mpi.h"
#include "mcpar.hh"
#include "rosenbrock.hh"
#include "mcout.hh"

int main(int argc, char *argv[])
{
  const int nparam=2;
  const float gmu[2] = {0.0f,0.0f};
  const float gsig2[2] = {1.0f, 2.0f};
  Rosenbrock1 L(2);
  MCout rslts(nparam);

  // Set up MPI
  int mpistat = MPI_Init(&argc, &argv);
  if(mpistat != MPI_SUCCESS) {
    std::cerr << "Error on MPI_Init.  Exiting.\n";
    return mpistat;
  }
  int size,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  // Set up the Parallel MC
  // 2 parameters, 4 chains per MPI process
  MCPar mcpar(nparam,4,size,rank);        

  float pinit[8] = {0.0f,0.0f, 2.0f,2.0f, 0.0f,1.5f, 0.0f,-2.0f};

  mcpar.run(100000,500, pinit, L, rslts);

  // output
  std::stringstream ofname;
  ofname << "mcpar-dgauss." << std::setfill('0') << std::setw(3) << rank << ".txt";
  //std::string ofn(ofname.str());
  std::ofstream outfile(ofname.str().c_str());
  for(int i=0; i<rslts.size(); ++i) {
    const float *pset = rslts.getpset(i);
    for(int j=0; j<rslts.nparam(); ++j)
      outfile << pset[j] << "\t";
    outfile << "\n";
  }

  MPI_Finalize();
  
  return 0;
}

