#ifndef VLFUNC_HH_
#define VLFUNC_HH_

/* Vector likelihood function:
 * npset:   number of full sets of parameters
 *     x:   input values (== npset * [number of function parameters])
 *     y:   output values (== npset)
 */
class VLFunc {
public:
  virtual int operator()(int npset, const float *x, float *restrict y) = 0;
};


#endif
