#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  assert(nrhs==1);
  assert(nlhs==0);

#define SIGNAL_IN   prhs[0]

  // size of array
  assert(mxIsDouble(SIGNAL_IN));
  assert(mxGetNumberOfElements(SIGNAL_IN) == 1);
  int signal = (int)mxGetPr(SIGNAL_IN)[0]; // cast from double to int

  EMPIRE_API_sendConvergenceSignal(signal);

#undef NAME_IN
#undef SIZE_IN
#undef SIGNAL_IN
}
