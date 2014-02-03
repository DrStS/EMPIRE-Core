#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==0);
    assert(nlhs==1);

#define SIGNAL_OUT plhs[0]

    // convergence signal
    SIGNAL_OUT = mxCreateDoubleScalar(-1);
    mxGetPr(SIGNAL_OUT)[0] = EMPIRE_API_recvConvergenceSignal();

#undef SIGNAL_OUT
}
