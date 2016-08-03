#ifdef USE_RFUNC
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "mcpar.hh"
#include "mcutil.hh"
#include "rfunc.hh"
#include "mcout.hh"

#define NPARAM 9
#define NPSET 8

int main(int argc, char *argv[])
{

  // Set up MPI
  int mpistat = MPI_Init(&argc, &argv);
  if(mpistat != MPI_SUCCESS) {
    std::cerr << "Error on MPI_Init.  Exiting.\n";
    return 1;
  }
  int size,rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // check arguments
  // usage:  mcpar-rfunc <R-Source> <Input-Data> <nsamp> <Output-file>
  if(argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <R-Source> <Input-Data> <Output-file> [nwarmup] [nsamp]\n";
  }
  std::string Rfile(argv[1]);
  std::string Dfile(argv[2]); 
  std::string Ofile(argv[3]); 
  std::ofstream outfile(Ofile);      // stream for output file
  // set up results object
  MCout rslts(NPARAM, &outfile, MPI_COMM_WORLD);
  // number of samples
  int nsamp = 100;
  int nwarmup = 1; // This default is way too small for a production calculation; it's meant for preliminary testing. 

  if(argc > 4) {
    int nwtmp = atoi(argv[4]);
    if(nwtmp > 0)
      nwarmup = nwtmp;
    else
      std::cerr << "Arg value for nwarmup is bogus: " << argv[4]
                << "\nUsing default value of " << nwarmup << "\n";
  }

  if(argc > 5) {
    int nstmp = atoi(argv[5]);
    if(nstmp > 0)
      nsamp = nstmp;
    else
      std::cerr << "Arg value for nsamp is bogus: " << argv[5]
                << "\nUsing default value of " << nsamp << "\n";
  }

  if(rank==0)
    std::cout << "nwarmup = " << nwarmup <<  "\tnsamp = " << nsamp << "\n";
  
  
  /* Set up vectors for model parameters.  Parameters are:
     As, An, xiss, xins, xisn, xinn, eps1n, lambda_s, k_s
  */
  // Minimum and maximum values of initial guesses for model parameters
  float plo[NPARAM] = {0.001f, 0.001f, -2.0f, -1.0f, -1.0f, -2.0f, 0.05f, 0.0f, 0.001f};
  float phi[NPARAM] = {1.0f,   1.0f,   0.0f,   1.0f,  1.0f,  0.0f, 1.5f,  5.0f, 10.0f};
  float pguess[NPARAM*NPSET];
  mcutil mcu;
  mcu.qriguess(rank, NPSET, NPARAM, plo, phi, pguess);

  // Set up the R likelihood function
  RFunc L(NPARAM, Rfile, Dfile, argc, argv);
  
  // Set up and run the parallel MC
  MCPar mcpar(NPARAM, NPSET, size, rank);
  // Uncomment the next two lines to log progress updates while the code runs.
  // mcpar.logging = true;
  // mcpar.logstep = 10;
  mcpar.run(nsamp, nwarmup, pguess, L, rslts);

  // Write output
  rslts.output();

  // Close down and quit.
  if(rank==0)
    std::cout << "FIN.\n";

  MPI_Finalize();
  return 0;
}

#else  // USE_RFUNC
#error "Building mcpar-rfunc.cc requires the USE_RFUNC to be defined in the environment."
#endif
