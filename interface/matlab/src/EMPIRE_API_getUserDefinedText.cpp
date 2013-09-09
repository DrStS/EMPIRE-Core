#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==1);

    // name of the xml element
    assert(mxIsChar(prhs[0]));
    char elementName[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(prhs[0], elementName, EMPIRE_API_NAME_STRING_LENGTH);

    // text
    char *text = EMPIRE_API_getUserDefinedText(elementName);
    plhs[0] = mxCreateString(text);
}

