/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EMPIRE.
 *
 *  EMPIRE is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EMPIRE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
 */
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
class TestIGAMortarMapperMultiPatchPlanarSurface: public CppUnit::TestFixture {
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
        int numNodes = 9;

        double* controlPoints1 = new double[4 * 4];
        controlPoints1[0] = -1.0;
        controlPoints1[1] = -1.0;
        controlPoints1[2] = 0.0;
        controlPoints1[3] = 1.0;

        controlPoints1[4] = 0.0;
        controlPoints1[5] = -1.0;
        controlPoints1[6] = 0.0;
        controlPoints1[7] = 1.0;

        controlPoints1[8] = -1.0;
        controlPoints1[9] = 0.0;
        controlPoints1[10] = 0.0;
        controlPoints1[11] = 1.0;

        controlPoints1[12] = 0.0;
        controlPoints1[13] = 0.0;
        controlPoints1[14] = 0.0;
        controlPoints1[15] = 1.0;

        int* dofIndex1 = new int[4];
        dofIndex1[0] = 0;
        dofIndex1[1] = 1;
        dofIndex1[2] = 3;
        dofIndex1[3] = 4;

        double* controlPoints2 = new double[4 * 4];
        controlPoints2[0] = 0.0;
        controlPoints2[1] = -1.0;
        controlPoints2[2] = 0.0;
        controlPoints2[3] = 1.0;

        controlPoints2[4] = 1.0;
        controlPoints2[5] = -1.0;
        controlPoints2[6] = 0.0;
        controlPoints2[7] = 1.0;

        controlPoints2[8] = 0.0;
        controlPoints2[9] = 0.0;
        controlPoints2[10] = 0.0;
        controlPoints2[11] = 1.0;

        controlPoints2[12] = 1.0;
        controlPoints2[13] = 0.0;
        controlPoints2[14] = 0.0;
        controlPoints2[15] = 1.0;

        int* dofIndex2 = new int[4];
        dofIndex2[0] = 1;
        dofIndex2[1] = 2;
        dofIndex2[2] = 4;
        dofIndex2[3] = 5;

        double* controlPoints3 = new double[4 * 4];
        controlPoints3[0] = -1.0;
        controlPoints3[1] = 0.0;
        controlPoints3[2] = 0.0;
        controlPoints3[3] = 1.0;

        controlPoints3[4] = 0.0;
        controlPoints3[5] = 0.0;
        controlPoints3[6] = 0.0;
        controlPoints3[7] = 1.0;

        controlPoints3[8] = -1.0;
        controlPoints3[9] = 1.0;
        controlPoints3[10] = 0.0;
        controlPoints3[11] = 1.0;

        controlPoints3[12] = 0.0;
        controlPoints3[13] = 1.0;
        controlPoints3[14] = 0.0;
        controlPoints3[15] = 1.0;

        int* dofIndex3 = new int[4];
        dofIndex3[0] = 3;
        dofIndex3[1] = 4;
        dofIndex3[2] = 6;
        dofIndex3[3] = 7;

        double* controlPoints4 = new double[4 * 4];
        controlPoints4[0] = 0.0;
        controlPoints4[1] = 0.0;
        controlPoints4[2] = 0.0;
        controlPoints4[3] = 1.0;

        controlPoints4[4] = 1.0;
        controlPoints4[5] = 0.0;
        controlPoints4[6] = 0.0;
        controlPoints4[7] = 1.0;

        controlPoints4[8] = 0.0;
        controlPoints4[9] = 1.0;
        controlPoints4[10] = 0.0;
        controlPoints4[11] = 1.0;

        controlPoints4[12] = 1.0;
        controlPoints4[13] = 1.0;
        controlPoints4[14] = 0.0;
        controlPoints4[15] = 1.0;

        int* dofIndex4 = new int[4];
        dofIndex4[0] = 4;
        dofIndex4[1] = 5;
        dofIndex4[2] = 7;
        dofIndex4[3] = 8;

        // Construct the IGA Surface Patch
        theIGAMeshM = new IGAMesh("The IGA Mesh", numNodes);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector1M, qM, vNoKnotsM, vKnotVector1M,
                uNoControlPointsM, vNoControlPointsM, controlPoints1, dofIndex1);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector2M, qM, vNoKnotsM, vKnotVector2M,
                uNoControlPointsM, vNoControlPointsM, controlPoints2, dofIndex2);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector3M, qM, vNoKnotsM, vKnotVector3M,
                uNoControlPointsM, vNoControlPointsM, controlPoints3, dofIndex3);
        theIGAMeshM->addPatch(pM, uNoKnotsM, uKnotVector4M, qM, vNoKnotsM, vKnotVector4M,
                uNoControlPointsM, vNoControlPointsM, controlPoints4, dofIndex4);

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

        int* dofIndexNetS = new int[numControlPoints];
        for (int i = 0; i < numControlPoints; i++)
            dofIndexNetS[i] = i;

        // Construct the IGA Surface Patch
        theIGAMeshS = new IGAMesh("The IGA Mesh", numControlPoints);
        theIGAMeshS->addPatch(pS, uNoKnotsS, uKnotVectorS, qS, vNoKnotsS, vKnotVectorS,
                uNoControlPoints, vNoControlPoints, globalControlPointsS, dofIndexNetS);

        int numNodesS = 4;
        int numElemsS = 1;
        theFEMeshS = new FEMesh("Fluid", numNodesS, numElemsS);
        for (int i = 0; i < numElemsS; i++)
            theFEMeshS->numNodesPerElem[i] = 4;
        theFEMeshS->initElems();
        theFEMeshS->elems[0 * 4 + 0] = 1;
        theFEMeshS->elems[0 * 4 + 1] = 2;
        theFEMeshS->elems[0 * 4 + 2] = 4;
        theFEMeshS->elems[0 * 4 + 3] = 3;

        for (int i = 0; i < numNodesS; i++)
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
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperMultiPatchPlanarSurface);
    CPPUNIT_TEST (testMapping);
//    CPPUNIT_TEST (testLeakage);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperMultiPatchPlanarSurface);
