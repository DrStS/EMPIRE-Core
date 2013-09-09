#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==6);
    assert(nlhs==0);

    // number of nodes
    assert(mxIsDouble(prhs[0]));
    assert(mxGetNumberOfElements(prhs[0]) == 1);
    int numNodes = (int) mxGetPr(prhs[0])[0]; // cast from double to int

    // number of elements
    assert(mxIsDouble(prhs[1]));
    assert(mxGetNumberOfElements(prhs[1]) == 1);
    int numElems = (int) mxGetPr(prhs[1])[0]; // cast from double to int

    // coordinates of nodes
    assert(mxIsDouble(prhs[2]));
    assert(mxGetNumberOfElements(prhs[2]) == numNodes*3);
    double *nodes = mxGetPr(prhs[2]);

    // IDs of nodes
    assert(mxIsDouble(prhs[3]));
    assert(mxGetNumberOfElements(prhs[3]) == numNodes);
    int *nodeIDs = doubleArrayToIntArray(mxGetPr(prhs[3]), numNodes);

    // number of nodes of each element
    assert(mxIsDouble(prhs[4]));
    assert(mxGetNumberOfElements(prhs[4]) == numElems);
    int *numNodesPerElem = doubleArrayToIntArray(mxGetPr(prhs[4]), numElems);

    // connectivity table of each element
    int elemsArrayLength = 0;
    for (int i = 0; i < numElems; i++)
        elemsArrayLength += numNodesPerElem[i];
    assert(mxIsDouble(prhs[5]));
    assert(mxGetNumberOfElements(prhs[5]) == elemsArrayLength);
    int *elems = doubleArrayToIntArray(mxGetPr(prhs[5]), elemsArrayLength);

    EMPIRE_API_sendMesh(numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);
    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elems;
}
