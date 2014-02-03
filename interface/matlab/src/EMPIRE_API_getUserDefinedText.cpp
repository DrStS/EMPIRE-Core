#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include <string>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==1);

#define ELEMENT_NAME_IN prhs[0]
#define TEXT_OUT plhs[0]

    // name of the xml element
    assert(ELEMENT_NAME_IN);
    char elementName[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(ELEMENT_NAME_IN, elementName, EMPIRE_API_NAME_STRING_LENGTH);

    // text
    char *text = EMPIRE_API_getUserDefinedText(elementName);
    TEXT_OUT = mxCreateString(text);

#undef ELEMENT_NAME_IN
#undef TEXT_OUT
}

