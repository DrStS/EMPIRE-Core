#include "mex.h"
#include "matrix.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==0);

#define INPUT_FILE_NAME_IN prhs[0]

    assert(INPUT_FILE_NAME_IN);
    char inputFileName[EMPIRE_API_NAME_STRING_LENGTH];

    mxGetNChars(INPUT_FILE_NAME_IN, inputFileName, EMPIRE_API_NAME_STRING_LENGTH);

    EMPIRE_API_Connect(inputFileName);

#undef INPUT_FILE_NAME_IN
}
