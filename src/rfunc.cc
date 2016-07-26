#ifdef USE_RFUNC
/* Since this likelihood function uses the R libraries, building it is
   optional.  Defining the variable USE_RFUNC in your environment will
   cause it to be built. */
#include "rfunc.hh"
#include <RInside.h> 

RFunc::RFunc(int anparam, const std::string &asrcfile, const std::string &ainputfile, int argc, char *argv[]) :
  nparam(anparam),
  srcfile(asrcfile),
  inputfile(ainputfile)
{
  // call the setup function, which will initialize R
  setup(argc, argv);
}


void RFunc::setup(int argc, char *argv[])
{
  // start up the embedded R environment
  Rp.reset(new RInside(argc,argv));

  // R command to source the user's code file
  std::string srccmd = "source('" + srcfile + "')";
  Rp->parseEvalQ(srccmd);

  // Call the mc.setup function
  std::string setupcmd = "mc.setup('" + inputfile + "')";
  Rp->parseEvalQ(setupcmd);

  // R should now be ready to go.
}

int RFunc::operator()(int npset, const float *x, float *restrict y)
{
  int ntot = npset*nparam;
  // TODO:  avoid reallocating these vectors and strings every time we call this function.
  Rcpp::NumericVector xR(x,x+ntot);
  std::string lfcmd = "mc.likelihood(input.params, input.npset)";

  (*Rp)["input.params"] = xR;
  (*Rp)["input.npset"] = npset;
  Rcpp::NumericVector yR = (*Rp).parseEval(lfcmd);

  // copy results back to output vector
  for(int i=0; i<npset; ++i)
    y[i] = yR[i];

  return 0;
}

#endif
