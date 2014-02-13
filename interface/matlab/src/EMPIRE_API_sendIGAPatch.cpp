#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"
#include "HelperFunctions.h"

#include <iostream>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    assert(nrhs==9);
    assert(nlhs==0);

#define P_DEGREE_IN                  prhs[0]
#define U_NUM_KNOTS_IN             prhs[1]
#define U_KNOT_VECTOR_IN          prhs[2]
#define Q_DEGREE_IN                 prhs[3]
#define V_NUM_KNOTS_IN              prhs[4]
#define V_KNOT_VECTOR_IN            prhs[5]
#define U_NUM_CONTROL_PTS_IN         prhs[6]
#define V_NUM_CONTROL_PTS_IN         prhs[7]
#define CONTROL_PTS_IDS_IN         prhs[8]


    // p degree
    assert(mxIsDouble(P_DEGREE_IN));
    assert(mxGetNumberOfElements(P_DEGREE_IN) == 1);
    int pDegree = (int) mxGetPr(P_DEGREE_IN)[0]; // cast from double to int

    // u number of knots
    assert(mxIsDouble(U_NUM_KNOTS_IN));
    assert(mxGetNumberOfElements(U_NUM_KNOTS_IN) == 1);
    int uNoKnots = (int) mxGetPr(U_NUM_KNOTS_IN)[0]; // cast from double to int

    // u knot vector
    assert(mxIsDouble(U_KNOT_VECTOR_IN));
    assert(mxGetNumberOfElements(U_KNOT_VECTOR_IN) == uNoKnots);
    double *uKnotVector = mxGetPr(U_KNOT_VECTOR_IN);

    // q degree
    assert(mxIsDouble(Q_DEGREE_IN));
    assert(mxGetNumberOfElements(Q_DEGREE_IN) == 1);
    int qDegree = (int) mxGetPr(Q_DEGREE_IN)[0]; // cast from double to int

    // v number of knots
    assert(mxIsDouble(V_NUM_KNOTS_IN));
    assert(mxGetNumberOfElements(V_NUM_KNOTS_IN) == 1);
    int vNoKnots = (int) mxGetPr(V_NUM_KNOTS_IN)[0]; // cast from double to int

    // v knot vector
    assert(mxIsDouble(V_KNOT_VECTOR_IN));
    assert(mxGetNumberOfElements(V_KNOT_VECTOR_IN) == vNoKnots);
    double *vKnotVector = mxGetPr(V_KNOT_VECTOR_IN);

    // u number of control points
    assert(mxIsDouble(U_NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(U_NUM_CONTROL_PTS_IN) == 1);
    int uNoControlPoints = (int) mxGetPr(U_NUM_CONTROL_PTS_IN)[0]; // cast from double to int

    // v number of control points
    assert(mxIsDouble(V_NUM_CONTROL_PTS_IN));
    assert(mxGetNumberOfElements(V_NUM_CONTROL_PTS_IN) == 1);
    int vNoControlPoints = (int) mxGetPr(V_NUM_CONTROL_PTS_IN)[0]; // cast from double to int

    // IDs of control points
    assert(mxIsDouble(CONTROL_PTS_IDS_IN));
    assert(mxGetNumberOfElements(CONTROL_PTS_IDS_IN) == uNoControlPoints * vNoControlPoints);
    int *controlPointNetIDs = doubleArrayToIntArray(mxGetPr(CONTROL_PTS_IDS_IN), uNoControlPoints * vNoControlPoints);

    EMPIRE_API_sendIGAPatch(pDegree,  uNoKnots, uKnotVector, qDegree, vNoKnots,
            vKnotVector, uNoControlPoints, vNoControlPoints, controlPointNetIDs);

    delete[] controlPointNetIDs;

#undef P_DEGREE_IN
#undef U_NUM_KNOTS_IN
#undef U_KNOT_VECTOR_IN
#undef Q_DEGREE_IN
#undef V_NUM_KNOTS_IN
#undef V_KNOT_VECTOR_IN
#undef U_NUM_CONTROL_PTS_IN
#undef V_NUM_CONTROL_PTS_IN
#undef CONTROL_PTS_IDS_IN
}
