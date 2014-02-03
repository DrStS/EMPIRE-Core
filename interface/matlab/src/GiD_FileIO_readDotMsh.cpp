#include "mex.h"
#include "matrix.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "GiDFileIO.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==1);
    assert(nlhs==7);

#define FILE_NAME_IN              prhs[0]
#define NUM_NODES_OUT             plhs[0]
#define NUM_ELEMENTS_OUT          plhs[1]
#define NODES_OUT                 plhs[2]
#define NODE_IDS_OUT              plhs[3]
#define NUM_NODES_PER_ELEMENT_OUT plhs[4]
#define ELEMENT_TABLE_OUT         plhs[5]
#define ELEMENT_IDS_OUT           plhs[6]

    assert(FILE_NAME_IN);
    char meshfile[EMPIRE_API_NAME_STRING_LENGTH];
    mxGetNChars(FILE_NAME_IN, meshfile, EMPIRE_API_NAME_STRING_LENGTH);

    int numNodes = -1;
    int numElems = -1;
    double *nodeCoors;
    int *nodeIDs;
    int *numNodesPerElem;
    int *elemTable;
    int *elemIDs;
    GiDFileIO::readDotMsh(meshfile, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem,
                  elemTable,elemIDs);

    NUM_ELEMENTS_OUT = mxCreateDoubleScalar((double) numElems);
    NUM_NODES_OUT = mxCreateDoubleScalar((double) numNodes);

    // nodeCoors
    NODES_OUT = mxCreateDoubleMatrix(3, numNodes, mxREAL);
    for (int i=0; i<numNodes * 3; i++) {
        mxGetPr(NODES_OUT)[i] = nodeCoors[i];
    }
    delete[] nodeCoors;

    // nodeIDs
    NODE_IDS_OUT = mxCreateDoubleMatrix(1, numNodes, mxREAL);
    for (int i=0; i<numNodes; i++) {
        mxGetPr(NODE_IDS_OUT)[i] = (double) nodeIDs[i];
    }
    delete[] nodeIDs;

    // numNodesPerElem
    NUM_NODES_PER_ELEMENT_OUT = mxCreateDoubleMatrix(1, numElems, mxREAL);
    for (int i=0; i<numElems; i++) {
        mxGetPr(NUM_NODES_PER_ELEMENT_OUT)[i] = (double) numNodesPerElem[i];
    }
    int elemsArrayLength = 0;
    for (int i = 0; i < numElems; i++)
        elemsArrayLength += numNodesPerElem[i];
    delete[] numNodesPerElem;

    // elemTable
    ELEMENT_TABLE_OUT = mxCreateDoubleMatrix(1, elemsArrayLength, mxREAL);
    for (int i=0; i<elemsArrayLength; i++) {
        mxGetPr(ELEMENT_TABLE_OUT)[i] = (double) elemTable[i];
    }
    delete[] elemTable;

    // elemIDs
    ELEMENT_IDS_OUT = mxCreateDoubleMatrix(1, numElems, mxREAL);
    for (int i=0; i<numElems; i++) {
        mxGetPr(ELEMENT_IDS_OUT)[i] = (double) elemIDs[i];
    }
    delete[] elemIDs;

#undef FILE_NAME_IN
#undef NUM_NODES_OUT
#undef NUM_ELEMENTS_OUT
#undef NODES_OUT
#undef NODE_IDS_OU
#undef NUM_NODES_PER_ELEMENT_OUT
#undef ELEMENT_TABLE_OUT
#undef ELEMENT_IDS_OUT
}
