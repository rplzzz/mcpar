#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "mcpar.hh"
#include "rosenbrock.hh"
#include "mcout.hh"

int main(int argc, char *argv[])
{
  const int nparam=2;
  const float gmu[2] = {0.0f,0.0f};
  const float gsig2[2] = {1.0f, 2.0f};
  Rosenbrock1 L(2);
  int nsamp = 100000;

  // Set up MPI
  int mpistat = MPI_Init(&argc, &argv);
  if(mpistat != MPI_SUCCESS) {
    std::cerr << "Error on MPI_Init.  Exiting.\n";
    return mpistat;
  }

  // Set up the results object.
  MCout rslts(nparam, &std::cout, MPI_COMM_WORLD);

  int size,rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  if(argc > 1)
    nsamp = atoi(argv[1]);
  if(rank==0)
    std::cout << "nsamp = " << nsamp << "\n";
  
  
  // Set up the Parallel MC
  // 2 parameters, 4 chains per MPI process
  MCPar mcpar(nparam,4,size,rank);        

  float pinit[8] = {0.0f,0.0f, 2.0f,2.0f, 0.0f,1.5f, 0.0f,-2.0f};

  mcpar.run(nsamp ,500, pinit, L, rslts);

  // output
  rslts.output();

  MPI_Finalize();
  
  return 0;
}

