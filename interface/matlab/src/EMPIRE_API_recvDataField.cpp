#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==3);
    assert(nlhs==0);

#define NAME_IN       prhs[0]
#define SIZE_IN       prhs[1]
#define DATA_FIELD_IN prhs[2]

    // data field name
    assert(mxIsChar(NAME_IN));
    char name[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(NAME_IN, name, EMPIRE_API_NAME_STRING_LENGTH);

    // size of the datafield array
    assert(mxIsDouble(SIZE_IN));
    assert(mxGetNumberOfElements(SIZE_IN) == 1);
    int sizeOfArray = (int) mxGetPr(SIZE_IN)[0]; // cast from double to int

    // array of the datafield
    assert(mxIsDouble(DATA_FIELD_IN));
    assert(mxGetNumberOfElements(DATA_FIELD_IN) == sizeOfArray);
    double *dataField = mxGetPr(DATA_FIELD_IN);

    EMPIRE_API_recvDataField("", sizeOfArray, dataField);

#undef NAME_IN
#undef SIZE_IN
#undef DATA_FIELD_IN
}
