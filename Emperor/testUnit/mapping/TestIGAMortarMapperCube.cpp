/*
 * TestIGAMortarMapper.cpp
 *
 *  Created on: May 8, 2013
 *      Author: Chenshen
 */

#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include "IGAMortarMapper.h"
#include "FEMesh.h"
#include "IGAPatchSurface.h"
#include "IGAMesh.h"
#include "DataField.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

namespace EMPIRE {
/*
 * Test the projection of a point on a multipatch structure. The benchmark consists of a 3-patch structure perpendicular to each other, and two fluid elements matching the two structural patches.
 */
class TestIGAMortarMapperCube: public CppUnit::TestFixture {
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

        // Provide an id for the basis
        int id_basis = 1;

        int numControlPoints = 8;

        int* controlPointID = new int[numControlPoints];
        for (int i = 0; i < numControlPoints; i++)
            controlPointID[i] = i + 1;

        theIGAMesh = new IGAMesh("IGAMesh", numControlPoints);

        // The polynomial degrees
        int p = 1;
        int q = 1;

        // Number of knots in both directions
        int uNoKnots = 4;
        int vNoKnots = 4;

        // The Control Point net
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;

        // The knot vectors in each direction
        double* uKnotVector1 = new double[uNoKnots];
        double* vKnotVector1 = new double[vNoKnots];
        double* uKnotVector2 = new double[uNoKnots];
        double* vKnotVector2 = new double[vNoKnots];
        double* uKnotVector3 = new double[uNoKnots];
        double* vKnotVector3 = new double[vNoKnots];

        uKnotVector1[0] = 0;
        uKnotVector1[1] = 0;
        uKnotVector1[2] = 1;
        uKnotVector1[3] = 1;

        vKnotVector1[0] = 0;
        vKnotVector1[1] = 0;
        vKnotVector1[2] = 1;
        vKnotVector1[3] = 1;
        for (int i = 0; i < 4; i++) {
            uKnotVector2[i] = uKnotVector1[i];
            vKnotVector2[i] = uKnotVector1[i];
            uKnotVector3[i] = uKnotVector1[i];
            vKnotVector3[i] = uKnotVector1[i];
        }

        double* controlPointNet1 = new double[4 * 4];

        controlPointNet1[0] = 0.0;
        controlPointNet1[1] = 0.0;
        controlPointNet1[2] = 0.0;
        controlPointNet1[3] = 1.0;

        controlPointNet1[4] = 0.0;
        controlPointNet1[5] = 0.0;
        controlPointNet1[6] = 1.0;
        controlPointNet1[7] = 1.0;

        controlPointNet1[8] = 0.0;
        controlPointNet1[9] = 1.0;
        controlPointNet1[10] = 0.0;
        controlPointNet1[11] = 1.0;

        controlPointNet1[12] = 0.0;
        controlPointNet1[13] = 1.0;
        controlPointNet1[14] = 1.0;
        controlPointNet1[15] = 1.0;

        double* controlPointNet2 = new double[4 * 4];

        controlPointNet2[0] = 0.0;
        controlPointNet2[1] = 1.0;
        controlPointNet2[2] = 0.0;
        controlPointNet2[3] = 1.0;

        controlPointNet2[4] = 0.0;
        controlPointNet2[5] = 1.0;
        controlPointNet2[6] = 1.0;
        controlPointNet2[7] = 1.0;

        controlPointNet2[8] = 1.0;
        controlPointNet2[9] = 1.0;
        controlPointNet2[10] = 0.0;
        controlPointNet2[11] = 1.0;

        controlPointNet2[12] = 1.0;
        controlPointNet2[13] = 1.0;
        controlPointNet2[14] = 1.0;
        controlPointNet2[15] = 1.0;

        double* controlPointNet3 = new double[4 * 4];

        controlPointNet3[0] = 1.0;
        controlPointNet3[1] = 0.0;
        controlPointNet3[2] = 0.0;
        controlPointNet3[3] = 1.0;

        controlPointNet3[4] = 1.0;
        controlPointNet3[5] = 0.0;
        controlPointNet3[6] = 1.0;
        controlPointNet3[7] = 1.0;

        controlPointNet3[8] = 1.0;
        controlPointNet3[9] = 1.0;
        controlPointNet3[10] = 0.0;
        controlPointNet3[11] = 1.0;

        controlPointNet3[12] = 1.0;
        controlPointNet3[13] = 1.0;
        controlPointNet3[14] = 1.0;
        controlPointNet3[15] = 1.0;

        int* dofIndexNet1 = new int[4];
        int* dofIndexNet2 = new int[4];
        int* dofIndexNet3 = new int[4];
        dofIndexNet1[0] = 0;
        dofIndexNet1[1] = 1;
        dofIndexNet1[2] = 2;
        dofIndexNet1[3] = 3;

        dofIndexNet2[0] = 2;
        dofIndexNet2[1] = 3;
        dofIndexNet2[2] = 6;
        dofIndexNet2[3] = 7;

        dofIndexNet3[0] = 4;
        dofIndexNet3[1] = 5;
        dofIndexNet3[2] = 6;
        dofIndexNet3[3] = 7;

        theIGAMesh->addPatch(p, uNoKnots, uKnotVector1, q, vNoKnots, vKnotVector1, uNoControlPoints,
                vNoControlPoints, controlPointNet1, dofIndexNet1);
        theIGAMesh->addPatch(p, uNoKnots, uKnotVector2, q, vNoKnots, vKnotVector2, uNoControlPoints,
                vNoControlPoints, controlPointNet2, dofIndexNet2);
        theIGAMesh->addPatch(p, uNoKnots, uKnotVector3, q, vNoKnots, vKnotVector3, uNoControlPoints,
                vNoControlPoints, controlPointNet3, dofIndexNet3);

        int numNodes = 6;
        int numElems = 2;
        theFEMesh = new FEMesh("Fluid", numNodes, numElems);
        theFEMesh->numNodesPerElem[0] = 4;
        theFEMesh->numNodesPerElem[1] = 4;
        theFEMesh->initElems();
        theFEMesh->elems[0] = 1;
        theFEMesh->elems[1] = 2;
        theFEMesh->elems[2] = 4;
        theFEMesh->elems[3] = 3;
        theFEMesh->elems[4] = 3;
        theFEMesh->elems[5] = 4;
        theFEMesh->elems[6] = 6;
        theFEMesh->elems[7] = 5;
        for (int i = 0; i < numNodes; i++)
            theFEMesh->nodeIDs[i] = i + 1;

        theFEMesh->nodes[0 * 3 + 0] = 0.05;
        theFEMesh->nodes[0 * 3 + 1] = 0.05;
        theFEMesh->nodes[0 * 3 + 2] = 0.05;

        theFEMesh->nodes[1 * 3 + 0] = 0.05;
        theFEMesh->nodes[1 * 3 + 1] = 0.05;
        theFEMesh->nodes[1 * 3 + 2] = 0.95;

        theFEMesh->nodes[2 * 3 + 0] = 0.05;
        theFEMesh->nodes[2 * 3 + 1] = 0.95;
        theFEMesh->nodes[2 * 3 + 2] = 0.05;

        theFEMesh->nodes[3 * 3 + 0] = 0.05;
        theFEMesh->nodes[3 * 3 + 1] = 0.95;
        theFEMesh->nodes[3 * 3 + 2] = 0.95;

        theFEMesh->nodes[4 * 3 + 0] = 0.95;
        theFEMesh->nodes[4 * 3 + 1] = 0.95;
        theFEMesh->nodes[4 * 3 + 2] = 0.05;

        theFEMesh->nodes[5 * 3 + 0] = 0.95;
        theFEMesh->nodes[5 * 3 + 1] = 0.95;
        theFEMesh->nodes[5 * 3 + 2] = 0.95;

        bool isMappingIGA2FEM = true;
        theMapper = new IGAMortarMapper("Test IGA Mortar Mapper for Cube", theIGAMesh, theFEMesh, 0.5, 16,
                25, isMappingIGA2FEM);

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

        int nS = theIGAMesh->getNumNodes();
        int nF = theFEMesh->numNodes;
        double fieldS[nS];
        double fieldF[nF];

        for (int i = 0; i < nS; i++)
            fieldS[i] = 1.0;
        theMapper->consistentMapping(fieldS, fieldF);
        for (int i = 0; i < nF; i++)
            CPPUNIT_ASSERT(fabs(fieldF[i] - 1.0) < Tol);
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMortarMapperCube);

    CPPUNIT_TEST (testMapping);
//	CPPUNIT_TEST(testMappingPrint);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */
CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMortarMapperCube);
