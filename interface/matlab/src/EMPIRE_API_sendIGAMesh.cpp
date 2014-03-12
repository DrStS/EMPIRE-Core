#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

#include <iostream>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==3);
    assert(nlhs==0);

#define NAME_IN                  prhs[0]
#define NUM_PATCHES_IN             prhs[1]
#define NUM_NODES_IN          prhs[2]

    // mesh name
    assert(mxIsChar(NAME_IN));
    char name[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(NAME_IN, name, EMPIRE_API_NAME_STRING_LENGTH);

    // number of patches
    assert(mxIsDouble(NUM_PATCHES_IN));
    assert(mxGetNumberOfElements(NUM_PATCHES_IN) == 1);
    int numPatches = (int) mxGetPr(NUM_PATCHES_IN)[0]; // cast from double to int

    // number of nodes
    assert(mxIsDouble(NUM_NODES_IN));
    assert(mxGetNumberOfElements(NUM_NODES_IN) == 1);
    int numNodes = (int) mxGetPr(NUM_NODES_IN)[0]; // cast from double to int

    EMPIRE_API_sendIGAMesh(name, numPatches, numNodes);


#undef NAME_IN
#undef NUM_PATCHES_IN
#undef NUM_NODES_IN

}
