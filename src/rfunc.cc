#ifdef USE_RFUNC
/* Since this likelihood function uses the R libraries, building it is
   optional.  Defining the variable USE_RFUNC in your environment will
   cause it to be built. */
#include "mpi.h"
#include "rfunc.hh"
#include <RInside.h>


RFunc::RFunc(const std::string &asrcfile, const std::string &ainputfile, int argc, char *argv[])
{
  // call the setup function, which will initialize R
  setup(argc, argv, asrcfile, ainputfile);
}


void RFunc::setup(int argc, char *argv[], const std::string &srcfile, const std::string &inputfile)
{
  // start up the embedded R environment
  Rp.reset(new RInside(argc,argv));

  // R command to source the user's code file
  std::string srccmd = "source('" + srcfile + "')";
  Rp->parseEvalQ(srccmd);

  // store the MPI rank, just in case we need it for some reason
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  (*Rp)["input.mpi.rank"] = rank;

  // Call the mc.setup function.  The setup function returns a matrix
  // of recommended limits for the parameter values.
  std::string setupcmd = "mc.setup('" + inputfile + "')";
  Rcpp::NumericMatrix plohi = Rp->parseEval(setupcmd);
  // This matrix gives us the information we need to fill out the parameter info.
  _nparam = plohi.ncol();
  _plo.resize(_nparam);
  _phi.resize(_nparam);

  for(int j=0; j<_nparam; ++j) {
    _plo[j] = plohi(0,j);         // plo is in the first row.
    _phi[j] = plohi(1,j);         // phi is in the second.
  }

  // R should now be ready to go.
}

int RFunc::operator()(int npset, const float *x, float *restrict y)
{
  int ntot = npset*_nparam;
  // TODO:  avoid reallocating these vectors and strings every time we call this function.
  Rcpp::NumericVector xR(x,x+ntot);
  std::string lfcmd = "mc.likelihood(input.params, input.npset)";

  for(int i=0;i<ntot;++i)
    xR[i] = x[i];

  (*Rp)["input.params"] = xR;
  (*Rp)["input.npset"] = npset;
  Rcpp::NumericVector yR = (*Rp).parseEval(lfcmd);

  // copy results back to output vector
  for(int i=0; i<npset; ++i)
    y[i] = yR[i];

  return 0;
}

#endif
