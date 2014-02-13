#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

#include <iostream>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==5);
    assert(nlhs==0);

#define NAME_IN                  prhs[0]
#define NUM_PATCHES_IN             prhs[1]
#define NUM_CONTROL_PTS_IN          prhs[2]
#define CONTROL_PTS_IN                 prhs[3]
#define CONTROL_PTS_IDS_IN              prhs[4]

    // mesh name
    assert(mxIsChar(NAME_IN));
    char name[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(NAME_IN, name, EMPIRE_API_NAME_STRING_LENGTH);

    // number of patches
    assert(mxIsDouble(NUM_PATCHES_IN));
    assert(mxGetNumberOfElements(NUM_PATCHES_IN) == 1);
    int numPatches = (int) mxGetPr(NUM_PATCHES_IN)[0]; // cast from double to int

    // number of control points
    assert(mxIsDouble(NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(NUM_CONTROL_PTS_IN) == 1);
    int numControlPoints = (int) mxGetPr(NUM_CONTROL_PTS_IN)[0]; // cast from double to int

    // global control points
    assert(mxIsDouble(NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(NUM_CONTROL_PTS_IN) == numControlPoints*4);
    double *globalControlPoints = mxGetPr(NUM_CONTROL_PTS_IN);

    // IDs of control points
    assert(mxIsDouble(CONTROL_PTS_IDS_IN));
    assert(mxGetNumberOfElements(CONTROL_PTS_IDS_IN) == numControlPoints);
    int *controlPointIDs = doubleArrayToIntArray(mxGetPr(CONTROL_PTS_IDS_IN), numControlPoints);

    EMPIRE_API_sendIGAMesh(name, numPatches, numControlPoints, globalControlPoints, controlPointIDs);

    delete[] controlPointIDs;

#undef NAME_IN
#undef NUM_PATCHES_IN
#undef NUM_CONTROL_PTS_IN
#undef CONTROL_PTS_IN
#undef CONTROL_PTS_IDS_IN
}
