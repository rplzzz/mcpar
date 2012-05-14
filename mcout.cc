#include <mpi.h>
#include "mcout.hh"

void MCout::output(std::ostream &outfile, int mpisize, int mpirank) const
{
  /* collect all the samples from the other processes into one big
     buffer and output them.  This is not a great way to do this, but
     it will serve for now. */
  size_t ntot = 0;
  float *buf = collect(mpisize, mpirank, &ntot); 
  if(mpirank == 0) {

    int npset = ntot / nparam_;
    int indx = 0;
    for(int i=0; i<npset; ++i) {
      for(int j=0; j<nparam_; ++j)
        outfile << buf[indx++] << "  ";
      outfile << "\n";
    }
    delete [] buf;              // buf only needs to be deleted in the rank 0 process
  } 
}

//! Collect the results from the back-end processes, return in a newly-allocated buffer.
//! \remark The buffer will be allocated only in the rank=0 process; other processes will return NULL.
float * MCout::collect(int mpisize, int mpirank, size_t *ntot) const
{
  float *buf = 0;
  if(mpirank == 0) {
    *ntot = mpisize*pvals.size();
    buf = new float[*ntot];
    
    if(mpisize > 1) {
      int mpistat = MPI_Gather((void*) &pvals[0], pvals.size(), MPI_FLOAT,
                               (void*) buf, pvals.size(), MPI_FLOAT,
                               0, MPI_COMM_WORLD);
      if(mpistat != MPI_SUCCESS) {
        std::cerr << "Unable to gather output data.  Aborting.\n";
        MPI_Abort(MPI_COMM_WORLD, mpistat);
      }
    }
    else {
      const float *pv = &pvals[0];
      const int imax = pvals.size();
#pragma ivdep
      for(int i=0; i<imax; ++i)
        buf[i] = pv[i];
    }
  }
  else {                        // mpirank > 0
    int mpistat = MPI_Gather((void*) &pvals[0], pvals.size(), MPI_FLOAT,
                             NULL, pvals.size(), MPI_FLOAT,
                             0, MPI_COMM_WORLD);
    if(mpistat != MPI_SUCCESS) {
      std::cerr << "Unable to gather output data.  Aborting.\n";
      MPI_Abort(MPI_COMM_WORLD, mpistat);
    }
  }
  return buf;
}
