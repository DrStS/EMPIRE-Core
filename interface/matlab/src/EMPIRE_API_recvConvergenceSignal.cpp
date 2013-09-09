#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  assert(nrhs==0);
  assert(nlhs==1);

  // convergence signal
  plhs[0] = mxCreateDoubleScalar(-1);
  mxGetPr(plhs[0])[0] = EMPIRE_API_recvConvergenceSignal();
}
