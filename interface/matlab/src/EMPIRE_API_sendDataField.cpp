#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==2);
    assert(nlhs==0);

    // size of the datafield array
    assert(mxIsDouble(prhs[0]));
    assert(mxGetNumberOfElements(prhs[0]) == 1);
    int sizeOfArray = (int) mxGetPr(prhs[0])[0]; // cast from double to int

    // array of the datafield
    assert(mxIsDouble(prhs[1]));
    assert(mxGetNumberOfElements(prhs[1]) == sizeOfArray);
    double *dataField = mxGetPr(prhs[1]);

    EMPIRE_API_sendDataField(sizeOfArray, dataField);
}
