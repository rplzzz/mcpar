#include <mpi.h>
#include <limits>
#include "mcout.hh"

MCout::MCout(int np, std::ostream *aoutstream, MPI_Comm acomm) : nparam_(np), next(0),npset(0), maxsamps(0), nextout(0)
{
  int stat = MPI_Comm_dup(acomm, &mComm);
  if(stat != MPI_SUCCESS) {
    std::cerr << "MCout: unable to duplicate input communicator in constructor." << std::endl;
    MPI_Abort(acomm,1);
  }
  MPI_Comm_rank(mComm, &mrank);
  MPI_Comm_size(mComm, &msize);

  maxlparams.resize(nparam_);
  maxlval = -(std::numeric_limits<float>::infinity());
  
  // for now, do all output through the rank-zero process.  Should really change this to use a parallel output scheme
  if(mrank == 0)
    outstream = aoutstream;
  else
    outstream = 0; 
}


void MCout::output()
{
  /* collect all the samples from the other processes into one big
     buffer and output them.  This is not a great way to do this, but
     it will serve for now. */
  size_t ntot = 0;
  float *buf = collect(&ntot); 
  if(mrank == 0 && ntot>0) {

    int npset = ntot / nparam_;
    int indx = 0;
    for(int i=0; i<npset; ++i) {
      for(int j=0; j<nparam_; ++j)
        (*outstream) << buf[indx++] << "  ";
      (*outstream) << "\n";
    }
    delete [] buf;              // buf only needs to be deleted in the rank 0 process
  } 
}

//! Collect the results from the back-end processes, return in a newly-allocated buffer.
//! \remark The buffer will be allocated only in the rank=0 process; other processes will return NULL.
float * MCout::collect(size_t *ntot)
{
  float *buf = 0;
  int nout = next - nextout;

  if(nout <= 0) {
    *ntot = 0;
    return buf;
  }
  
  if(mrank == 0) {
    *ntot = msize*nout;
    buf = new float[*ntot];
    
    if(msize > 1) {
      int mpistat = MPI_Gather((void*) &pvals[nextout], nout, MPI_FLOAT,
                               (void*) buf, nout, MPI_FLOAT,
                               0, MPI_COMM_WORLD);
      if(mpistat != MPI_SUCCESS) {
        std::cerr << "Unable to gather output data.  Aborting.\n";
        MPI_Abort(MPI_COMM_WORLD, mpistat);
      }
    }
    else {
      const float *pv = &pvals[nextout];
      const int imax = next-nextout;
#pragma ivdep
      for(int i=0; i<imax; ++i)
        buf[i] = pv[i];
    }
  }
  else {                        // mpirank > 0
    int mpistat = MPI_Gather((void*) &pvals[nextout], nout, MPI_FLOAT,
                             NULL, nout, MPI_FLOAT,
                             0, MPI_COMM_WORLD);
    if(mpistat != MPI_SUCCESS) {
      std::cerr << "Unable to gather output data.  Aborting.\n";
      MPI_Abort(MPI_COMM_WORLD, mpistat);
    }
  }
  nextout = next;
  return buf;
}

const std::vector<float> &MCout::maxlike(float *lmax)
{
  struct {
    float val;
    int rank;
  } sendata,rcvdata;
  int root;
  
  sendata.val = maxlval;
  sendata.rank = mrank;

  int stat = MPI_Allreduce(&sendata, &rcvdata, 1, MPI_FLOAT_INT, MPI_MAXLOC, mComm);

  if(stat != MPI_SUCCESS) {
    // MPI failure can't be good.  Best just to bail
    std::cerr << "rank " << mrank << ":  MPI failure in MCout::maxlike().\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  root = rcvdata.rank;
  stat = MPI_Bcast(&maxlparams[0], nparam_, MPI_FLOAT, root, mComm);

  if(stat != MPI_SUCCESS) {
    // MPI failure can't be good.  Best just to bail
    std::cerr << "rank " << mrank << ":  MPI failure in MCout::maxlike().\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  *lmax = maxlval = rcvdata.val;

  return maxlparams;
}

void MCout::add(const float *pv, float lval)
{
  float *strt = &pvals[next];
  for(int i=0; i<nparam_; ++i) {
    assert(strt+i < &pvals[0]+pvals.size());
    strt[i] = pv[i];
  }
  next += nparam_;
  npset++;
  
  if(lval > maxlval) {
    maxlval = lval;
    for(int i=0; i<nparam_; ++i)
      maxlparams[i] = pv[i];
  }
}
