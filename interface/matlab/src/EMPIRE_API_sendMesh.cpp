#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

#include <iostream>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==7);
    assert(nlhs==0);

#define NAME_IN                  prhs[0]
#define NUM_NODES_IN             prhs[1]
#define NUM_ELEMENTS_IN          prhs[2]
#define NODES_IN                 prhs[3]
#define NODE_IDS_IN              prhs[4]
#define NUM_NODES_PER_ELEMENT_IN prhs[5]
#define ELEMENT_TABLE_IN         prhs[6]

    // mesh name
    assert(mxIsChar(NAME_IN));
    char name[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(NAME_IN, name, EMPIRE_API_NAME_STRING_LENGTH);

    // number of nodes
    assert(mxIsDouble(NUM_NODES_IN));
    assert(mxGetNumberOfElements(NUM_NODES_IN) == 1);
    int numNodes = (int) mxGetPr(NUM_NODES_IN)[0]; // cast from double to int

    // number of elements
    assert(mxIsDouble(NUM_ELEMENTS_IN));
    assert(mxGetNumberOfElements(NUM_ELEMENTS_IN) == 1);
    int numElems = (int) mxGetPr(NUM_ELEMENTS_IN)[0]; // cast from double to int

    // coordinates of nodes
    assert(mxIsDouble(NODES_IN));
    assert(mxGetNumberOfElements(NODES_IN) == numNodes*3);
    double *nodes = mxGetPr(NODES_IN);

    // IDs of nodes
    assert(mxIsDouble(NODE_IDS_IN));
    assert(mxGetNumberOfElements(NODE_IDS_IN) == numNodes);
    int *nodeIDs = doubleArrayToIntArray(mxGetPr(NODE_IDS_IN), numNodes);

    // number of nodes of each element
    assert(mxIsDouble(NUM_NODES_PER_ELEMENT_IN));
    assert(mxGetNumberOfElements(NUM_NODES_PER_ELEMENT_IN) == numElems);
    int *numNodesPerElem = doubleArrayToIntArray(mxGetPr(NUM_NODES_PER_ELEMENT_IN), numElems);

    // connectivity table of each element
    int elemsArrayLength = 0;
    for (int i = 0; i < numElems; i++)
        elemsArrayLength += numNodesPerElem[i];
    assert(mxIsDouble(ELEMENT_TABLE_IN));
    assert(mxGetNumberOfElements(ELEMENT_TABLE_IN) == elemsArrayLength);
    int *elems = doubleArrayToIntArray(mxGetPr(ELEMENT_TABLE_IN), elemsArrayLength);

    EMPIRE_API_sendMesh(name, numNodes, numElems, nodes, nodeIDs, numNodesPerElem, elems);

    delete[] nodeIDs;
    delete[] numNodesPerElem;
    delete[] elems;

#undef NAME_IN
#undef NUM_NODES_IN
#undef NUM_ELEMENTS_IN
#undef NODES_IN
#undef NODE_IDS_IN
#undef NUM_NODES_PER_ELEMENT_IN
#undef ELEMENT_TABLE_IN
}
