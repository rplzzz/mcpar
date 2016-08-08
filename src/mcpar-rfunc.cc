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
  
  
  // Set up the R likelihood function
  RFunc L(Rfile, Dfile, argc, argv);
  const int nparam = L.nparam();
  std::vector<float> plo = L.plo();
  std::vector<float> phi = L.phi();

  // set up results object
  MCout rslts(nparam, &outfile, MPI_COMM_WORLD); 
  float *pguess = new float[nparam*NPSET];
  mcutil mcu;
  mcu.qriguess(rank, NPSET, nparam, &plo[0], &phi[0], pguess);

  // Set up and run the parallel MC
  MCPar mcpar(nparam, NPSET, size, rank);
  // Uncomment the next two lines to log progress updates while the code runs.
  // mcpar.logging = true;
  // mcpar.logstep = 10;
  mcpar.run(nsamp, nwarmup, pguess, L, rslts);

  // Write output
  rslts.output();

  // Close down and quit.
  if(rank==0)
    std::cout << "FIN.\n";

  delete [] pguess;
  
  MPI_Finalize();
  return 0;
}

#else  // USE_RFUNC
#error "Building mcpar-rfunc.cc requires the USE_RFUNC to be defined in the environment."
#endif
