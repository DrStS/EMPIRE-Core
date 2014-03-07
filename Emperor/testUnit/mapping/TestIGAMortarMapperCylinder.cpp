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
class TestIGAMortarMapperCylinder: public CppUnit::TestFixture {
private:
    IGAMortarMapper* mapper;
public:
private:
    IGAMortarMapper* theMapper;
    IGAMesh* theIGAMesh;
    FEMesh* theFEMesh;

    double Tol;
    int a;
    int b;

public:
    void setUp() {

        Tol = 1e-13;

        // The polynomial degrees
        const int p1 = 2;
        const int q1 = 1;

        // Number of knots in both directions
        const int uNoKnots1 = 8;
        const int vNoKnots1 = 4;

        // The knot vectors in each direction
        double* uKnotVector1 = new double[12];
        uKnotVector1[0] = 0.0;
        uKnotVector1[1] = 0.0;
        uKnotVector1[2] = 0.0;
        uKnotVector1[3] = 0.5;
        uKnotVector1[4] = 0.5;
        uKnotVector1[5] = 1.0;
        uKnotVector1[6] = 1.0;
        uKnotVector1[7] = 1.0;
        double* vKnotVector1 = new double[4];
        vKnotVector1[0] = 0.0;
        vKnotVector1[1] = 0.0;
        vKnotVector1[2] = 1.0;
        vKnotVector1[3] = 1.0;

        // The Control Point net
        const int uNoControlPoints1 = uNoKnots1 - p1 - 1;
        const int vNoControlPoints1 = vNoKnots1 - q1 - 1;
        const int numNodes = 16;

        double* controlPoints1 = new double[uNoControlPoints1 * vNoControlPoints1 * 4];
        controlPoints1[0] = 0.0;
        controlPoints1[1] = 0.0;
        controlPoints1[2] = 0.0;
        controlPoints1[3] = 1.0;

        controlPoints1[4] = 0.0;
        controlPoints1[5] = 0.0;
        controlPoints1[6] = 1.0;
        controlPoints1[7] = 1.0;

        controlPoints1[8] = 0.0;
        controlPoints1[9] = 0.5;
        controlPoints1[10] = 0.0;
        controlPoints1[11] = 0.707106781186548;

        controlPoints1[12] = 0.0;
        controlPoints1[13] = 0.5;
        controlPoints1[14] = 1.0;
        controlPoints1[15] = 0.707106781186548;

        controlPoints1[16] = 0.5;
        controlPoints1[17] = 0.5;
        controlPoints1[18] = 0.0;
        controlPoints1[19] = 1.0;

        controlPoints1[20] = 0.5;
        controlPoints1[21] = 0.5;
        controlPoints1[22] = 1.0;
        controlPoints1[23] = 1.0;

        controlPoints1[24] = 1.0;
        controlPoints1[25] = 0.5;
        controlPoints1[26] = 0.0;
        controlPoints1[27] = 0.707106781186548;

        controlPoints1[28] = 1.0;
        controlPoints1[29] = 0.5;
        controlPoints1[30] = 1.0;
        controlPoints1[31] = 0.707106781186548;

        controlPoints1[32] = 1.0;
        controlPoints1[33] = 0.0;
        controlPoints1[34] = 0.0;
        controlPoints1[35] = 1.0;

        controlPoints1[36] = 1.0;
        controlPoints1[37] = 0.0;
        controlPoints1[38] = 1.0;
        controlPoints1[39] = 1.0;

        const int p2 = 2;
        const int q2 = 1;

        // Number of knots in both directions
        const int uNoKnots2 = 8;
        const int vNoKnots2 = 4;

        // The knot vectors in each direction
        double* uKnotVector2 = new double[12];
        uKnotVector2[0] = 0.0;
        uKnotVector2[1] = 0.0;
        uKnotVector2[2] = 0.0;
        uKnotVector2[3] = 0.5;
        uKnotVector2[4] = 0.5;
        uKnotVector2[5] = 1.0;
        uKnotVector2[6] = 1.0;
        uKnotVector2[7] = 1.0;
        double* vKnotVector2 = new double[4];
        vKnotVector2[0] = 0.0;
        vKnotVector2[1] = 0.0;
        vKnotVector2[2] = 1.0;
        vKnotVector2[3] = 1.0;

        // The Control Point net
        const int uNoControlPoints2 = uNoKnots2 - p2 - 1;
        const int vNoControlPoints2 = vNoKnots2 - q2 - 1;

        double* controlPoints2 = new double[uNoControlPoints2 * vNoControlPoints2 * 4];

        controlPoints2[0] = 1.0;
        controlPoints2[1] = 0.0;
        controlPoints2[2] = 0.0;
        controlPoints2[3] = 1.0;

        controlPoints2[4] = 1.0;
        controlPoints2[5] = 0.0;
        controlPoints2[6] = 1.0;
        controlPoints2[7] = 1.0;

        controlPoints2[8] = 1.0;
        controlPoints2[9] = -0.5;
        controlPoints2[10] = 0.0;
        controlPoints2[11] = 0.707106781186548;

        controlPoints2[12] = 1.0;
        controlPoints2[13] = -0.5;
        controlPoints2[14] = 1.0;
        controlPoints2[15] = 0.707106781186548;

        controlPoints2[16] = 0.5;
        controlPoints2[17] = -0.5;
        controlPoints2[18] = 0.0;
        controlPoints2[19] = 1.0;

        controlPoints2[20] = 0.5;
        controlPoints2[21] = -0.5;
        controlPoints2[22] = 1.0;
        controlPoints2[23] = 1.0;

        controlPoints2[24] = 0.0;
        controlPoints2[25] = -0.5;
        controlPoints2[26] = 0.0;
        controlPoints2[27] = 0.707106781186548;

        controlPoints2[28] = 0.0;
        controlPoints2[29] = -0.5;
        controlPoints2[30] = 1.0;
        controlPoints2[31] = 0.707106781186548;

        controlPoints2[32] = 0.0;
        controlPoints2[33] = 0.0;
        controlPoints2[34] = 0.0;
        controlPoints2[35] = 1.0;

        controlPoints2[36] = 0.0;
        controlPoints2[37] = 0.0;
        controlPoints2[38] = 1.0;
        controlPoints2[39] = 1.0;

        int* dofIndexNet1 = new int(10);
        dofIndexNet1[0] = 0;
        dofIndexNet1[1] = 1;
        dofIndexNet1[2] = 2;
        dofIndexNet1[3] = 3;
        dofIndexNet1[4] = 4;
        dofIndexNet1[5] = 5;
        dofIndexNet1[6] = 6;
        dofIndexNet1[7] = 7;
        dofIndexNet1[8] = 8;
        dofIndexNet1[9] = 9;

        int* dofIndexNet2 = new int(10);
        dofIndexNet2[0] = 8;
        dofIndexNet2[1] = 9;
        dofIndexNet2[2] = 10;
        dofIndexNet2[3] = 11;
        dofIndexNet2[4] = 12;
        dofIndexNet2[5] = 13;
        dofIndexNet2[6] = 14;
        dofIndexNet2[7] = 15;
        dofIndexNet2[8] = 0;
        dofIndexNet2[9] = 1;

        // Construct the IGA Surface Patch
        theIGAMesh = new IGAMesh("The IGA Mesh", numNodes);
        theIGAMesh->addPatch(p1, uNoKnots1, uKnotVector1, q1, vNoKnots1, vKnotVector1,
                uNoControlPoints1, vNoControlPoints1, controlPoints1, dofIndexNet1);
        theIGAMesh->addPatch(p2, uNoKnots2, uKnotVector2, q2, vNoKnots2, vKnotVector2,
                uNoControlPoints2, vNoControlPoints2, controlPoints2, dofIndexNet2);

        int numNodesFE = 36;
        int numElemsFE = 18;
        theFEMesh = new FEMesh("Fluid", numNodesFE, numElemsFE);
        for (int i = 0; i < numElemsFE; i++)
            theFEMesh->numNodesPerElem[i] = 4;
        theFEMesh->initElems();

        int elems[72] = { 15, 16, 1561, 1560, 16, 17, 1562, 1561, 17, 18, 1563, 1562, 18, 19, 1564,
                1563, 19, 20, 1565, 1564, 20, 21, 1566, 1565, 21, 22, 1567, 1566, 22, 23, 1568,
                1567, 23, 24, 1569, 1568, 1560, 2341, 796, 15, 2341, 2342, 797, 796, 2342, 2343,
                798, 797, 2343, 2344, 799, 798, 2344, 2345, 800, 799, 2345, 2346, 801, 800, 2346,
                2347, 802, 801, 2347, 2348, 803, 802, 2348, 1569, 24, 803 };

        double nodes[108] = { 0, 0, 0, 0.0301537, 0.17101, 0, 0.116978, 0.321394, 0, 0.25, 0.433013,
                0, 0.413176, 0.492404, 0, 0.586824, 0.492404, 0, 0.75, 0.433013, 0, 0.883022,
                0.321394, 0, 0.969846, 0.17101, 0, 1, 6.12323e-17, 0, 0.0299778, -0.170939, 0,
                0.117288, -0.321099, 0, 0.249889, -0.433225, 0, 0.413167, -0.492341, 0, 0.586834,
                -0.492341, 0, 0.750111, -0.433225, 0, 0.882712, -0.321098, 0, 0.970022, -0.170939,
                0, 0, 0, 1, 0.0301537, 0.17101, 1, 0.116978, 0.321394, 1, 0.25, 0.433013, 1,
                0.413176, 0.492404, 1, 0.586824, 0.492404, 1, 0.75, 0.433013, 1, 0.883022, 0.321394,
                1, 0.969846, 0.17101, 1, 1, 6.12323e-17, 1, 0.0299778, -0.170939, 1, 0.117288,
                -0.321099, 1, 0.249889, -0.433225, 1, 0.413167, -0.492341, 1, 0.586834, -0.492341,
                1, 0.750111, -0.433225, 1, 0.882712, -0.321098, 1, 0.970022, -0.170939, 1 };

        double nodeIDs[36] = { 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 796, 797, 798, 799, 800, 801,
                802, 803, 1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569, 2341, 2342,
                2343, 2344, 2345, 2346, 2347, 2348 };

        for (int i = 0; i < numNodesFE * 3; i++)
            theFEMesh->nodes[i] = nodes[i];

        for (int i = 0; i < numNodesFE; i++)
            theFEMesh->nodeIDs[i] = nodeIDs[i];

        for (int i = 0; i < numElemsFE * 4; i++)
            theFEMesh->elems[i] = elems[i];


        bool isMappingIGA2FEM = true;
        theMapper = new IGAMortarMapper("Test IGA Mortar Mapper Cylinder", theIGAMesh, theFEMesh,
                 1e-1, 16, 25, isMappingIGA2FEM);

        // The comparison must be extended for all the entries
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(2, 21) - 4.853032548893787e-03) < Tol);


    }

    void tearDown() {

        delete theFEMesh;
        delete theIGAMesh;
        delete theMapper;

    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/

    void testMapping() {

    }

    void testLeakage() {
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperCylinder);
    CPPUNIT_TEST (testMapping);
//    CPPUNIT_TEST (testLeakage);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

//CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperCylinder);
