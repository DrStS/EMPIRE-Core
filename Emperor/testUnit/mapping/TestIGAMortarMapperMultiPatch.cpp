/*
 * TestIGAMortarMapper.cpp
 *
 *  Created on: May 8, 2013
 *      Author: chenshen
 */

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include "IGAMortarMapper.h"
#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "DataField.h"
#include "MathLibrary.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace EMPIRE {
class TestIGAMortarMapperMultiPatch: public CppUnit::TestFixture {
private:
    IGAMortarMapper* mapper;
public:
private:
    IGAMortarMapper* theMapperM;
    IGAMesh* theIGAMeshM;
    FEMesh* theFEMeshM;

    IGAMortarMapper* theMapperS;
    IGAMesh* theIGAMeshS;
    FEMesh* theFEMeshS;
    double Tol;
    int a;
    int b;

public:
    void setUp() {

        Tol = 1e-13;

        // The polynomial degrees
        int pM = 1;
        int qM = 1;

        // Number of knots in both directions
        int uNoKnotsM = 4;
        int vNoKnotsM = 4;

        // The knot vectors in each direction
        double* uKnotVector1M = new double[4];
        uKnotVector1M[0] = 0.0;
        uKnotVector1M[1] = 0.0;
        uKnotVector1M[2] = 1.0;
        uKnotVector1M[3] = 1.0;
        double* vKnotVector1M = new double[4];
        vKnotVector1M[0] = 0.0;
        vKnotVector1M[1] = 0.0;
        vKnotVector1M[2] = 1.0;
        vKnotVector1M[3] = 1.0;

        double* uKnotVector2M = new double[4];
        double* vKnotVector2M = new double[4];
        double* uKnotVector3M = new double[4];
        double* vKnotVector3M = new double[4];
        double* uKnotVector4M = new double[4];
        double* vKnotVector4M = new double[4];

        for (int i = 0; i < 4; i++) {
            uKnotVector2M[i] = uKnotVector1M[i];
            vKnotVector2M[i] = vKnotVector1M[i];
            uKnotVector3M[i] = uKnotVector1M[i];
            vKnotVector3M[i] = vKnotVector1M[i];
            uKnotVector4M[i] = uKnotVector1M[i];
            vKnotVector4M[i] = vKnotVector1M[i];
        }

        // The Control Point net
        int uNoControlPointsM = uNoKnotsM - pM - 1;
        int vNoControlPointsM = vNoKnotsM - qM - 1;
        int numControlPointsM = 9;
        double* globalControlPointsM = new double[numControlPointsM * 4];
        globalControlPointsM[0] = -1.0;
        globalControlPointsM[1] = -1.0;
        globalControlPointsM[2] = 0.0;
        globalControlPointsM[3] = 1.0;

        globalControlPointsM[4] = 0.0;
        globalControlPointsM[5] = -1.0;
        globalControlPointsM[6] = 0.0;
        globalControlPointsM[7] = 1.0;

        globalControlPointsM[8] = 1.0;
        globalControlPointsM[9] = -1.0;
        globalControlPointsM[10] = 0.0;
        globalControlPointsM[11] = 1.0;

        globalControlPointsM[12] = -1.0;
        globalControlPointsM[13] = 0.0;
        globalControlPointsM[14] = 0.0;
        globalControlPointsM[15] = 1.0;

        globalControlPointsM[16] = 0.0;
        globalControlPointsM[17] = 0.0;
        globalControlPointsM[18] = 0.0;
        globalControlPointsM[19] = 1.0;

        globalControlPointsM[20] = 1.0;
        globalControlPointsM[21] = 0.0;
        globalControlPointsM[22] = 0.0;
        globalControlPointsM[23] = 1.0;

        globalControlPointsM[24] = -1.0;
        globalControlPointsM[25] = 1.0;
        globalControlPointsM[26] = 0.0;
        globalControlPointsM[27] = 1.0;

        globalControlPointsM[28] = 0.0;
        globalControlPointsM[29] = 1.0;
        globalControlPointsM[30] = 0.0;
        globalControlPointsM[31] = 1.0;

        globalControlPointsM[32] = 1.0;
        globalControlPointsM[33] = 1.0;
        globalControlPointsM[34] = 0.0;
        globalControlPointsM[35] = 1.0;

        int* controlPointNetIDM = new int[numControlPointsM];
        for (int i = 0; i < numControlPointsM; i++)
            controlPointNetIDM[i] = i;

        int* controlPointNetIDPatch1M = new int[4];
        controlPointNetIDPatch1M[0] = 0;
        controlPointNetIDPatch1M[1] = 1;
        controlPointNetIDPatch1M[2] = 3;
        controlPointNetIDPatch1M[3] = 4;
        int* controlPointNetIDPatch2M = new int[4];
        controlPointNetIDPatch2M[0] = 1;
        controlPointNetIDPatch2M[1] = 2;
        controlPointNetIDPatch2M[2] = 4;
        controlPointNetIDPatch2M[3] = 5;

        int* controlPointNetIDPatch3M = new int[4];
        controlPointNetIDPatch3M[0] = 3;
        controlPointNetIDPatch3M[1] = 4;
        controlPointNetIDPatch3M[2] = 6;
        controlPointNetIDPatch3M[3] = 7;

        int* controlPointNetIDPatch4M = new int[4];
        controlPointNetIDPatch4M[0] = 4;
        controlPointNetIDPatch4M[1] = 5;
        controlPointNetIDPatch4M[2] = 7;
        controlPointNetIDPatch4M[3] = 8;

        // Construct the IGA Surface Patch
        theIGAMeshM = new IGAMesh("The IGA Mesh", numControlPointsM, globalControlPointsM,
                controlPointNetIDM);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector1M, qM, vNoKnotsM, vKnotVector1M,
                uNoControlPointsM, vNoControlPointsM, controlPointNetIDPatch1M);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector2M, qM, vNoKnotsM, vKnotVector2M,
                uNoControlPointsM, vNoControlPointsM, controlPointNetIDPatch2M);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector3M, qM, vNoKnotsM, vKnotVector3M,
                uNoControlPointsM, vNoControlPointsM, controlPointNetIDPatch3M);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector4M, qM, vNoKnotsM, vKnotVector4M,
                uNoControlPointsM, vNoControlPointsM, controlPointNetIDPatch4M);

        int numNodesM = 4;
        int numElemsM = 1;
        theFEMeshM = new FEMesh("Fluid", numNodesM, numElemsM);
        for (int i = 0; i < numElemsM; i++)
            theFEMeshM->numNodesPerElem[i] = 4;
        theFEMeshM->initElems();
        theFEMeshM->elems[0 * 4 + 0] = 1;
        theFEMeshM->elems[0 * 4 + 1] = 2;
        theFEMeshM->elems[0 * 4 + 2] = 4;
        theFEMeshM->elems[0 * 4 + 3] = 3;

        for (int i = 0; i < numNodesM; i++)
            theFEMeshM->nodeIDs[i] = i + 1;
        theFEMeshM->nodes[0 * 3 + 0] = -0.4;
        theFEMeshM->nodes[0 * 3 + 1] = -0.3;
        theFEMeshM->nodes[0 * 3 + 2] = 1.0e-06;

        theFEMeshM->nodes[1 * 3 + 0] = 0.5;
        theFEMeshM->nodes[1 * 3 + 1] = -0.3;
        theFEMeshM->nodes[1 * 3 + 2] = 1.0e-06;

        theFEMeshM->nodes[2 * 3 + 0] = -0.4;
        theFEMeshM->nodes[2 * 3 + 1] = 0.3;
        theFEMeshM->nodes[2 * 3 + 2] = 1.0e-06;

        theFEMeshM->nodes[3 * 3 + 0] = 0.5;
        theFEMeshM->nodes[3 * 3 + 1] = 0.3;
        theFEMeshM->nodes[3 * 3 + 2] = 1.0e-06;

        theMapperM = new IGAMortarMapper("Test IGA Mortar Mapper Multi Patch", theIGAMeshM,
                theFEMeshM, 1e-2, 16, 25);

        // The polynomial degrees
        int pS = 1;
        int qS = 1;

        // Number of knots in both directions
        int uNoKnotsS = 5;
        int vNoKnotsS = 5;

        // The knot vectors in each direction
        double* uKnotVectorS = new double[5];
        uKnotVectorS[0] = 0.000000000000000e+00;
        uKnotVectorS[1] = 0.000000000000000e+00;
        uKnotVectorS[2] = 1.000000000000000e+00;
        uKnotVectorS[3] = 2.000000000000000e+00;
        uKnotVectorS[4] = 2.000000000000000e+00;
        double* vKnotVectorS = new double[5];
        vKnotVectorS[0] = 0.000000000000000e+00;
        vKnotVectorS[1] = 0.000000000000000e+00;
        vKnotVectorS[2] = 1.000000000000000e+00;
        vKnotVectorS[3] = 2.000000000000000e+00;
        vKnotVectorS[4] = 2.000000000000000e+00;

        // The Control Point net
        int uNoControlPoints = uNoKnotsS - pS - 1;
        int vNoControlPoints = vNoKnotsS - qS - 1;
        int numControlPoints = uNoControlPoints * vNoControlPoints;
        double* globalControlPointsS = new double[numControlPoints * 4];
        globalControlPointsS[0] = -1.0;
        globalControlPointsS[1] = -1.0;
        globalControlPointsS[2] = 0.0;
        globalControlPointsS[3] = 1.0;

        globalControlPointsS[4] = 0.0;
        globalControlPointsS[5] = -1.0;
        globalControlPointsS[6] = 0.0;
        globalControlPointsS[7] = 1.0;

        globalControlPointsS[8] = 1.0;
        globalControlPointsS[9] = -1.0;
        globalControlPointsS[10] = 0.0;
        globalControlPointsS[11] = 1.0;

        globalControlPointsS[12] = -1.0;
        globalControlPointsS[13] = 0.0;
        globalControlPointsS[14] = 0.0;
        globalControlPointsS[15] = 1.0;

        globalControlPointsS[16] = 0.0;
        globalControlPointsS[17] = 0.0;
        globalControlPointsS[18] = 0.0;
        globalControlPointsS[19] = 1.0;

        globalControlPointsS[20] = 1.0;
        globalControlPointsS[21] = 0.0;
        globalControlPointsS[22] = 0.0;
        globalControlPointsS[23] = 1.0;

        globalControlPointsS[24] = -1.0;
        globalControlPointsS[25] = 1.0;
        globalControlPointsS[26] = 0.0;
        globalControlPointsS[27] = 1.0;

        globalControlPointsS[28] = 0.0;
        globalControlPointsS[29] = 1.0;
        globalControlPointsS[30] = 0.0;
        globalControlPointsS[31] = 1.0;

        globalControlPointsS[32] = 1.0;
        globalControlPointsS[33] = 1.0;
        globalControlPointsS[34] = 0.0;
        globalControlPointsS[35] = 1.0;

        int* controlPointNetIDS = new int[numControlPoints];
        for (int i = 0; i < numControlPoints; i++)
            controlPointNetIDS[i] = i;

        int* controlPointNetIDPatchS = new int[numControlPoints];
        for (int i = 0; i < numControlPoints; i++)
            controlPointNetIDPatchS[i] = i;

        // Construct the IGA Surface Patch
        theIGAMeshS = new IGAMesh("The IGA Mesh", numControlPoints, globalControlPointsS,
                controlPointNetIDS);
        theIGAMeshS->addPatch(pS, uNoKnotsS, uKnotVectorS, qS, vNoKnotsS, vKnotVectorS,
                uNoControlPoints, vNoControlPoints, controlPointNetIDPatchS);

        int numNodes = 4;
        int numElems = 1;
        theFEMeshS = new FEMesh("Fluid", numNodes, numElems);
        for (int i = 0; i < numElems; i++)
            theFEMeshS->numNodesPerElem[i] = 4;
        theFEMeshS->initElems();
        theFEMeshS->elems[0 * 4 + 0] = 1;
        theFEMeshS->elems[0 * 4 + 1] = 2;
        theFEMeshS->elems[0 * 4 + 2] = 4;
        theFEMeshS->elems[0 * 4 + 3] = 3;

        for (int i = 0; i < numNodes; i++)
            theFEMeshS->nodeIDs[i] = i + 1;
        theFEMeshS->nodes[0 * 3 + 0] = -0.4;
        theFEMeshS->nodes[0 * 3 + 1] = -0.3;
        theFEMeshS->nodes[0 * 3 + 2] = 1.0e-06;

        theFEMeshS->nodes[1 * 3 + 0] = 0.5;
        theFEMeshS->nodes[1 * 3 + 1] = -0.3;
        theFEMeshS->nodes[1 * 3 + 2] = 1.0e-06;

        theFEMeshS->nodes[2 * 3 + 0] = -0.4;
        theFEMeshS->nodes[2 * 3 + 1] = 0.3;
        theFEMeshS->nodes[2 * 3 + 2] = 1.0e-06;

        theFEMeshS->nodes[3 * 3 + 0] = 0.5;
        theFEMeshS->nodes[3 * 3 + 1] = 0.3;
        theFEMeshS->nodes[3 * 3 + 2] = 1.0e-06;

        theMapperS = new IGAMortarMapper("Test IGA Mortar Mapper", theIGAMeshS, theFEMeshS, 1e-2,
                16, 25);

    }

    void tearDown() {

        delete theFEMeshM;
        delete theIGAMeshM;
        delete theMapperM;

        delete theFEMeshS;
        delete theIGAMeshS;
        delete theMapperS;

    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/

    void testMapping() {

        CPPUNIT_ASSERT(fabs((*theMapperS->C_NR)(0, 1) - (*theMapperM->C_NR)(0, 1)) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapperS->C_NR)(1, 6) - (*theMapperM->C_NR)(1, 6)) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapperS->C_NR)(3, 8) - (*theMapperM->C_NR)(3, 8)) < Tol);

    }

    void testLeakage() {
        for (int i = 0; i < 100000000; i++) {
            IGAMortarMapper* theMapper = new IGAMortarMapper("Test IGA Mortar Mapper", theIGAMeshS,
                    theFEMeshS, 1e-2, 16, 25);
            delete theMapper;
        }
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperMultiPatch);
    CPPUNIT_TEST (testMapping);
//    CPPUNIT_TEST (testLeakage);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperMultiPatch);
