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
        controlPoints1[5] = 0.5;
        controlPoints1[6] = 0.0;
        controlPoints1[7] = 0.707106781186548;

        controlPoints1[8] = 0.5;
        controlPoints1[9] = 0.5;
        controlPoints1[10] = 0.0;
        controlPoints1[11] = 1.0;

        controlPoints1[12] = 1.0;
        controlPoints1[13] = 0.5;
        controlPoints1[14] = 0.0;
        controlPoints1[15] = 0.707106781186548;

        controlPoints1[16] = 1.0;
        controlPoints1[17] = 0.0;
        controlPoints1[18] = 0.0;
        controlPoints1[19] = 1.0;

        controlPoints1[20] = 0.0;
        controlPoints1[21] = 0.0;
        controlPoints1[22] = 1.0;
        controlPoints1[23] = 1.0;

        controlPoints1[24] = 0.0;
        controlPoints1[25] = 0.5;
        controlPoints1[26] = 1.0;
        controlPoints1[27] = 0.707106781186548;

        controlPoints1[28] = 0.5;
        controlPoints1[29] = 0.5;
        controlPoints1[30] = 1.0;
        controlPoints1[31] = 1.0;

        controlPoints1[32] = 1.0;
        controlPoints1[33] = 0.5;
        controlPoints1[34] = 1.0;
        controlPoints1[35] = 0.707106781186548;

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
        controlPoints2[5] = -0.5;
        controlPoints2[6] = 0.0;
        controlPoints2[7] = 0.707106781186548;

        controlPoints2[8] = 0.5;
        controlPoints2[9] = -0.5;
        controlPoints2[10] = 0.0;
        controlPoints2[11] = 1.0;

        controlPoints2[12] = 0.0;
        controlPoints2[13] = -0.5;
        controlPoints2[14] = 0.0;
        controlPoints2[15] = 0.707106781186548;

        controlPoints2[16] = 0.0;
        controlPoints2[17] = 0.0;
        controlPoints2[18] = 0.0;
        controlPoints2[19] = 1.0;

        controlPoints2[20] = 1.0;
        controlPoints2[21] = 0.0;
        controlPoints2[22] = 1.0;
        controlPoints2[23] = 1.0;

        controlPoints2[24] = 1.0;
        controlPoints2[25] = -0.5;
        controlPoints2[26] = 1.0;
        controlPoints2[27] = 0.707106781186548;

        controlPoints2[28] = 0.5;
        controlPoints2[29] = -0.5;
        controlPoints2[30] = 1.0;
        controlPoints2[31] = 1.0;

        controlPoints2[32] = 0.0;
        controlPoints2[33] = -0.5;
        controlPoints2[34] = 1.0;
        controlPoints2[35] = 0.707106781186548;

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
        dofIndexNet1[5] = 8;
        dofIndexNet1[6] = 9;
        dofIndexNet1[7] = 10;
        dofIndexNet1[8] = 11;
        dofIndexNet1[9] = 12;

        int* dofIndexNet2 = new int(10);
        dofIndexNet2[0] = 4;
        dofIndexNet2[1] = 5;
        dofIndexNet2[2] = 6;
        dofIndexNet2[3] = 7;
        dofIndexNet2[4] = 0;
        dofIndexNet2[5] = 12;
        dofIndexNet2[6] = 13;
        dofIndexNet2[7] = 14;
        dofIndexNet2[8] = 15;
        dofIndexNet2[9] = 8;

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

        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(0, 0) - 0.03768673999069) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(0, 1) - 0.00970316331738) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(0, 10) - 0.00969610014000) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(0, 18) - 0.01884336999535) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(0, 19) - 0.00485158165869) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(0, 28) - 0.00484805007000) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(1, 1) - 0.03908086441244) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(1, 2) - 0.00970563322816) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(1, 18) - 0.00485158165869) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(1, 19) - 0.01954043220622) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(1, 20) - 0.00485281661408) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(2, 2) - 0.03906373003786) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(2, 3) - 0.00970606509779) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(2, 19) - 0.00485281661408) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(2, 20) - 0.01953186501893) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(2, 21) - 0.00485303254889) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(3, 3) - 0.03906927504677) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(3, 4) - 0.00970466300276) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(3, 20) - 0.00485303254889) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(3, 21) - 0.01953463752338) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(3, 22) - 0.00485233150138) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(4, 4) - 0.03847682614662) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(4, 5) - 0.00961956552707) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(4, 21) - 0.00485233150138) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(4, 22) - 0.01923841307331) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(4, 23) - 0.00480978276354) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(5, 5) - 0.03847682614662) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(5, 6) - 0.00970466300276) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(5, 22) - 0.00480978276354) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(5, 23) - 0.01923841307331) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(5, 24) - 0.00485233150138) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(6, 6) - 0.03906927504677) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(6, 7) - 0.00970606509779) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(6, 23) - 0.00485233150138) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(6, 24) - 0.01953463752338) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(6, 25) - 0.00485303254889) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(7, 7) - 0.03906371830076) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(7, 8) - 0.00970562751175) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(7, 24) - 0.00485303254889) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(7, 25) - 0.01953185915038) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(7, 26) - 0.00485281375587) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(8, 8) - 0.03908086530402) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(8, 9) - 0.00970316903091) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(8, 25) - 0.00485281375587) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(8, 26) - 0.01954043265201) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(8, 27) - 0.00485158451546) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(9, 9) - 0.03768675806919) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(9, 17) - 0.00969610394528) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(9, 26) - 0.00485158451546) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(9, 27) - 0.01884337903459) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(9, 35) - 0.00484805197264) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(10, 10) - 0.03907684086726) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(10, 11) - 0.00971120993334) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(10, 18) - 0.00484805007000) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(10, 28) - 0.01953842043363) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(10, 29) - 0.00485560496667) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(11, 11) - 0.03907930557183) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(11, 12) - 0.00970810506566) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(11, 28) - 0.00485560496667) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(11, 29) - 0.01953965278591) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(11, 30) - 0.00485405253283) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(12, 12) - 0.03906988953966) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(12, 13) - 0.00970301035605) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(12, 29) - 0.00485405253283) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(12, 30) - 0.01953494476983) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(12, 31) - 0.00485150517803) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(13, 13) - 0.03847820972050) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(13, 14) - 0.00962178767409) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(13, 30) - 0.00485150517803) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(13, 31) - 0.01923910486025) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(13, 32) - 0.00481089383705) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(14, 14) - 0.03847811543097) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(14, 15) - 0.00970295550313) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(14, 31) - 0.00481089383705) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(14, 32) - 0.01923905771549) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(14, 33) - 0.00485147775157) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(15, 15) - 0.03906986015873) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(15, 16) - 0.00970814779231) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(15, 32) - 0.00485147775157) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(15, 33) - 0.01953493007936) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(15, 34) - 0.00485407389616) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(16, 16) - 0.03907929644105) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(16, 17) - 0.00971116340341) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(16, 33) - 0.00485407389616) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(16, 34) - 0.01953964822052) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(16, 35) - 0.00485558170171) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(17, 17) - 0.03907675770560) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(17, 27) - 0.00484805197264) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(17, 34) - 0.00485558170171) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(17, 35) - 0.01953837885280) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(18, 18) - 0.03768673999069) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(18, 19) - 0.00970316331738) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(18, 28) - 0.00969610014000) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(19, 19) - 0.03908086441244) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(19, 20) - 0.00970563322816) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(20, 20) - 0.03906373003786) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(20, 21) - 0.00970606509779) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(21, 21) - 0.03906927504677) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(21, 22) - 0.00970466300276) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(22, 22) - 0.03847682614662) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(22, 23) - 0.00961956552707) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(23, 23) - 0.03847682614662) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(23, 24) - 0.00970466300276) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(24, 24) - 0.03906927504677) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(24, 25) - 0.00970606509779) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(25, 25) - 0.03906371830076) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(25, 26) - 0.00970562751175) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(26, 26) - 0.03908086530402) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(26, 27) - 0.00970316903091) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(27, 27) - 0.03768675806919) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(27, 35) - 0.00969610394528) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(28, 28) - 0.03907684086726) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(28, 29) - 0.00971120993334) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(29, 29) - 0.03907930557183) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(29, 30) - 0.00970810506566) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(30, 30) - 0.03906988953966) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(30, 31) - 0.00970301035605) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(31, 31) - 0.03847820972050) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(31, 32) - 0.00962178767409) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(32, 32) - 0.03847811543097) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(32, 33) - 0.00970295550313) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(33, 33) - 0.03906986015873) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(33, 34) - 0.00970814779231) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(34, 34) - 0.03907929644105) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(34, 35) - 0.00971116340341) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NN)(35, 35) - 0.03907675770560) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 0) - 0.05054357497080) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 1) - 0.00298874383804) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 2) - 0.00028481210500) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 6) - 0.00028419770032) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 7) - 0.00298467483391) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 8) - 0.02527178748540) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 9) - 0.00149437191902) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 10) - 0.00014240605250) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 14) - 0.00014209885016) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(0, 15) - 0.00149233741695) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(1, 0) - 0.03881230898402) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(1, 1) - 0.01563623482738) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(1, 2) - 0.00404111714657) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(1, 8) - 0.01940615449201) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(1, 9) - 0.00781811741369) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(1, 10) - 0.00202055857328) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(2, 0) - 0.02129006843939) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(2, 1) - 0.02306606730734) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(2, 2) - 0.01411929261708) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(2, 8) - 0.01064503421970) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(2, 9) - 0.01153303365367) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(2, 10) - 0.00705964630854) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(3, 0) - 0.00831052617122) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(3, 1) - 0.02056443790860) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(3, 2) - 0.02960503906748) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(3, 8) - 0.00415526308561) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(3, 9) - 0.01028221895430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(3, 10) - 0.01480251953374) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 0) - 0.00142203648292) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 1) - 0.00888076010430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 2) - 0.04708246019639) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 3) - 0.00039770154847) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 4) - 0.00001809634437) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 8) - 0.00071101824146) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 9) - 0.00444038005215) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 10) - 0.02354123009819) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 11) - 0.00019885077423) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(4, 12) - 0.00000904817219) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 0) - 0.00001809634437) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 1) - 0.00039770154847) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 2) - 0.04708246019639) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 3) - 0.00888076010430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 4) - 0.00142203648292) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 8) - 0.00000904817219) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 9) - 0.00019885077423) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 10) - 0.02354123009819) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 11) - 0.00444038005215) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(5, 12) - 0.00071101824146) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(6, 2) - 0.02960503906748) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(6, 3) - 0.02056443790860) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(6, 4) - 0.00831052617122) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(6, 10) - 0.01480251953374) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(6, 11) - 0.01028221895430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(6, 12) - 0.00415526308561) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(7, 2) - 0.01411929071142) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(7, 3) - 0.02306606141315) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(7, 4) - 0.02129005878573) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(7, 10) - 0.00705964535571) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(7, 11) - 0.01153303070657) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(7, 12) - 0.01064502939287) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(8, 2) - 0.00404111855487) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(8, 3) - 0.01563623742912) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(8, 4) - 0.03881230586269) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(8, 10) - 0.00202055927743) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(8, 11) - 0.00781811871456) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(8, 12) - 0.01940615293135) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 2) - 0.00028481260236) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 3) - 0.00298874713049) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 4) - 0.05054359625599) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 5) - 0.00298467702545) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 6) - 0.00028419803109) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 10) - 0.00014240630118) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 11) - 0.00149437356525) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 12) - 0.02527179812800) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 13) - 0.00149233851272) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(9, 14) - 0.00014209901554) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(10, 0) - 0.03881417007576) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(10, 6) - 0.00403863960302) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(10, 7) - 0.01563134126181) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(10, 8) - 0.01940708503788) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(10, 14) - 0.00201931980151) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(10, 15) - 0.00781567063091) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(11, 0) - 0.02130310926382) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(11, 6) - 0.01412164810958) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(11, 7) - 0.02307386319742) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(11, 8) - 0.01065155463191) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(11, 14) - 0.00706082405479) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(11, 15) - 0.01153693159871) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(12, 0) - 0.00831140168581) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(12, 6) - 0.02960431776919) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(12, 7) - 0.02056528550636) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(12, 8) - 0.00415570084291) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(12, 14) - 0.01480215888460) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(12, 15) - 0.01028264275318) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 0) - 0.00142204135231) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 4) - 0.00001810930205) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 5) - 0.00039788735921) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 6) - 0.04708406330674) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 7) - 0.00888090643033) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 8) - 0.00071102067615) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 12) - 0.00000905465102) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 13) - 0.00019894367961) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 14) - 0.02354203165337) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(13, 15) - 0.00444045321517) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 0) - 0.00001810849456) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 4) - 0.00142204816758) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 5) - 0.00888093402615) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 6) - 0.04708389361561) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 7) - 0.00039787430429) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 8) - 0.00000905424728) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 12) - 0.00071102408379) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 13) - 0.00444046701308) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 14) - 0.02354194680781) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(14, 15) - 0.00019893715214) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(15, 4) - 0.00831142920534) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(15, 5) - 0.02056529709789) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(15, 6) - 0.02960423715095) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(15, 12) - 0.00415571460267) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(15, 13) - 0.01028264854894) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(15, 14) - 0.01480211857547) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(16, 4) - 0.02130312496395) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(16, 5) - 0.02307385510831) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(16, 6) - 0.01412162756451) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(16, 12) - 0.01065156248197) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(16, 13) - 0.01153692755416) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(16, 14) - 0.00706081378226) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(17, 4) - 0.03881411072314) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(17, 5) - 0.01563129491713) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(17, 6) - 0.00403861941403) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(17, 12) - 0.01940705536157) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(17, 13) - 0.00781564745856) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(17, 14) - 0.00201930970701) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 0) - 0.02527178748540) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 1) - 0.00149437191902) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 2) - 0.00014240605250) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 6) - 0.00014209885016) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 7) - 0.00149233741695) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 8) - 0.05054357497080) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 9) - 0.00298874383804) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 10) - 0.00028481210500) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 14) - 0.00028419770032) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(18, 15) - 0.00298467483391) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(19, 0) - 0.01940615449201) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(19, 1) - 0.00781811741369) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(19, 2) - 0.00202055857328) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(19, 8) - 0.03881230898402) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(19, 9) - 0.01563623482738) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(19, 10) - 0.00404111714657) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(20, 0) - 0.01064503421970) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(20, 1) - 0.01153303365367) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(20, 2) - 0.00705964630854) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(20, 8) - 0.02129006843939) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(20, 9) - 0.02306606730734) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(20, 10) - 0.01411929261708) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(21, 0) - 0.00415526308561) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(21, 1) - 0.01028221895430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(21, 2) - 0.01480251953374) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(21, 8) - 0.00831052617122) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(21, 9) - 0.02056443790860) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(21, 10) - 0.02960503906748) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 0) - 0.00071101824146) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 1) - 0.00444038005215) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 2) - 0.02354123009819) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 3) - 0.00019885077423) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 4) - 0.00000904817219) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 8) - 0.00142203648292) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 9) - 0.00888076010430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 10) - 0.04708246019639) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 11) - 0.00039770154847) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(22, 12) - 0.00001809634437) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 0) - 0.00000904817219) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 1) - 0.00019885077423) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 2) - 0.02354123009819) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 3) - 0.00444038005215) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 4) - 0.00071101824146) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 8) - 0.00001809634437) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 9) - 0.00039770154847) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 10) - 0.04708246019639) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 11) - 0.00888076010430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(23, 12) - 0.00142203648292) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(24, 2) - 0.01480251953374) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(24, 3) - 0.01028221895430) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(24, 4) - 0.00415526308561) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(24, 10) - 0.02960503906748) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(24, 11) - 0.02056443790860) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(24, 12) - 0.00831052617122) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(25, 2) - 0.00705964535571) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(25, 3) - 0.01153303070657) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(25, 4) - 0.01064502939287) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(25, 10) - 0.01411929071142) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(25, 11) - 0.02306606141315) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(25, 12) - 0.02129005878573) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(26, 2) - 0.00202055927743) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(26, 3) - 0.00781811871456) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(26, 4) - 0.01940615293135) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(26, 10) - 0.00404111855487) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(26, 11) - 0.01563623742912) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(26, 12) - 0.03881230586269) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 2) - 0.00014240630118) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 3) - 0.00149437356525) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 4) - 0.02527179812800) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 5) - 0.00149233851272) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 6) - 0.00014209901554) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 10) - 0.00028481260236) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 11) - 0.00298874713049) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 12) - 0.05054359625599) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 13) - 0.00298467702545) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(27, 14) - 0.00028419803109) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(28, 0) - 0.01940708503788) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(28, 6) - 0.00201931980151) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(28, 7) - 0.00781567063091) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(28, 8) - 0.03881417007576) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(28, 14) - 0.00403863960302) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(28, 15) - 0.01563134126181) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(29, 0) - 0.01065155463191) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(29, 6) - 0.00706082405479) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(29, 7) - 0.01153693159871) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(29, 8) - 0.02130310926382) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(29, 14) - 0.01412164810958) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(29, 15) - 0.02307386319742) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(30, 0) - 0.00415570084291) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(30, 6) - 0.01480215888460) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(30, 7) - 0.01028264275318) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(30, 8) - 0.00831140168581) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(30, 14) - 0.02960431776919) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(30, 15) - 0.02056528550636) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 0) - 0.00071102067615) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 4) - 0.00000905465102) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 5) - 0.00019894367961) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 6) - 0.02354203165337) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 7) - 0.00444045321517) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 8) - 0.00142204135231) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 12) - 0.00001810930205) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 13) - 0.00039788735921) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 14) - 0.04708406330674) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(31, 15) - 0.00888090643033) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 0) - 0.00000905424728) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 4) - 0.00071102408379) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 5) - 0.00444046701308) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 6) - 0.02354194680781) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 7) - 0.00019893715214) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 8) - 0.00001810849456) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 12) - 0.00142204816758) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 13) - 0.00888093402615) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 14) - 0.04708389361561) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(32, 15) - 0.00039787430429) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(33, 4) - 0.00415571460267) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(33, 5) - 0.01028264854894) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(33, 6) - 0.01480211857547) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(33, 12) - 0.00831142920534) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(33, 13) - 0.02056529709789) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(33, 14) - 0.02960423715095) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(34, 4) - 0.01065156248197) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(34, 5) - 0.01153692755416) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(34, 6) - 0.00706081378226) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(34, 12) - 0.02130312496395) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(34, 13) - 0.02307385510831) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(34, 14) - 0.01412162756451) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(35, 4) - 0.01940705536157) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(35, 5) - 0.00781564745856) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(35, 6) - 0.00201930970701) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(35, 12) - 0.03881411072314) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(35, 13) - 0.01563129491713) < Tol);
        CPPUNIT_ASSERT(fabs((*theMapper->C_NR)(35, 14) - 0.00403861941403) < Tol);
    }

    void testLeakage() {
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperCylinder);
    CPPUNIT_TEST (testMapping);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperCylinder);
