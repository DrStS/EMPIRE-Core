#include "mex.h"
#include "matrix.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  assert(nrhs==0);
  assert(nlhs==0);

  EMPIRE_API_Disconnect();
}

