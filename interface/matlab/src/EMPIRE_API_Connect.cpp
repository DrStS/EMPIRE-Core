#include "mex.h"
#include "matrix.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==0);

    assert(mxIsChar(prhs[0]));
    char inputFileName[EMPIRE_API_NAME_STRING_LENGTH];

    mxGetNChars(prhs[0], inputFileName, EMPIRE_API_NAME_STRING_LENGTH);

    EMPIRE_API_Connect(inputFileName);
}
