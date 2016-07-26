#ifndef RFUNC_HH_
#define RFUNC_HH_

#ifdef USE_RFUNC
/* Since this likelihood function uses the R libraries, building it is
   optional.  Defining the variable USE_RFUNC in your environment will
   cause it to be built. */

#include <string>
#include <memory>
#include "vlfunc.hh"
#include <RInside.h>

/*
 * Likelihood function for a model written in R.

 * This class uses the RInside package for R to create an instance of
 * the R interpreter inside the C++ code.  The setup function will
 * cause the embedded R instance to source a specified R source file.
 * The R source should provide two functions:

 * mc.setup(datafile): perform all one-time setup tasks, including
 *         reading the deta file named in its argument.  This function
 *         will be called once when the object is constructed.

 * mc.likelihood(x, npset): evaluate the likelihood function.  This
 *        includes running the model and comparing to whatever
 *        "observed" data you are using.  This function should be able
 *        to accept an arbitrary number of parameter sets concatenated
 *        in the x argument (i.e., as if you had said x <- c(p1, p2,
 *        ... , pn)).  The npset argument should give the number of
 *        parameter sets present.  The return value will be a vector
 *        of likelihood values.
 */

class RFunc : public VLFunc
{
private:
  
  int nparam;                   // need the number of parameters b/c it's not passed to operator()
  std::string srcfile;
  std::string inputfile;

  std::unique_ptr<RInside> Rp;

  void setup(int argc, char *argv[]);

public:
  // constructor: First argument is the name of the file with the R
  // source; second is the file with the input data.  argc and argv
  // are the arguments to main().
  RFunc(int nparam, const std::string &asrcfile, const std::string &ainputfile, int argc, char *argv[]);
  virtual int operator()(int anpset, const float *x, float *restrict y);
};


#endif

#endif
