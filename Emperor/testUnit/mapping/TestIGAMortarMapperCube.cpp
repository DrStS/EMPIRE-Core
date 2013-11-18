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


		double* globalControlPoints = new double[numControlPoints * 4];
		globalControlPoints[0] = 0.0;
		globalControlPoints[1] = 0.0;
		globalControlPoints[2] = 0.0;
		globalControlPoints[3] = 1.0;

		globalControlPoints[4] = 0.0;
		globalControlPoints[5] = 0.0;
		globalControlPoints[6] = 1.0;
		globalControlPoints[7] = 1.0;

		globalControlPoints[8] = 0.0;
		globalControlPoints[9] = 1.0;
		globalControlPoints[10] = 0.0;
		globalControlPoints[11] = 1.0;

		globalControlPoints[12] = 0.0;
		globalControlPoints[13] = 1.0;
		globalControlPoints[14] = 1.0;
		globalControlPoints[15] = 1.0;

		globalControlPoints[16] = 1.0;
		globalControlPoints[17] = 0.0;
		globalControlPoints[18] = 0.0;
		globalControlPoints[19] = 1.0;

		globalControlPoints[20] = 1.0;
		globalControlPoints[21] = 0.0;
		globalControlPoints[22] = 1.0;
		globalControlPoints[23] = 1.0;

		globalControlPoints[24] = 1.0;
		globalControlPoints[25] = 1.0;
		globalControlPoints[26] = 0.0;
		globalControlPoints[27] = 1.0;

		globalControlPoints[28] = 1.0;
		globalControlPoints[29] = 1.0;
		globalControlPoints[30] = 1.0;
		globalControlPoints[31] = 1.0;


		int* controlPointID = new int[numControlPoints];
		for (int i = 0; i < numControlPoints; i++)
			controlPointID[i] = i + 1;


		theIGAMesh = new IGAMesh("IGAMesh", numControlPoints, globalControlPoints, controlPointID);

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
		int* controlPointNet1 = new int[4];
		int* controlPointNet2 = new int[4];
		int* controlPointNet3 = new int[4];
		controlPointNet1[0] = 1;
		controlPointNet1[1] = 2;
		controlPointNet1[2] = 3;
		controlPointNet1[3] = 4;

		controlPointNet2[0] = 3;
		controlPointNet2[1] = 4;
		controlPointNet2[2] = 7;
		controlPointNet2[3] = 8;

		controlPointNet3[0] = 5;
		controlPointNet3[1] = 6;
		controlPointNet3[2] = 7;
		controlPointNet3[3] = 8;

		theIGAMesh->addPatch(p, uNoKnots, uKnotVector1, q, vNoKnots,vKnotVector1, uNoControlPoints, vNoControlPoints, controlPointNet1);
		theIGAMesh->addPatch(p, uNoKnots, uKnotVector2, q, vNoKnots,vKnotVector2, uNoControlPoints, vNoControlPoints, controlPointNet2);
		theIGAMesh->addPatch(p, uNoKnots, uKnotVector3, q, vNoKnots,vKnotVector3, uNoControlPoints, vNoControlPoints, controlPointNet3);

		int numNodes = 6;
		int numElems = 2;
		theFEMesh = new FEMesh("Fluid", numNodes, numElems);
		theFEMesh->initElems();
		theFEMesh->numNodesPerElem[0] = 4;
		theFEMesh->elems[0] = 1;
		theFEMesh->elems[1] = 2;
		theFEMesh->elems[2] = 4;
		theFEMesh->elems[3] = 3;
		theFEMesh->numNodesPerElem[1] = 4;
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

		theMapper = new IGAMortarMapper("Test IGA Mortar Mapper", theIGAMesh, theFEMesh,0.5,16,25);

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

		int nS = theIGAMesh->getNumControlPoints();
		int nF = theFEMesh->numNodes;
		double fieldS[nS];
		double fieldF[nF];

		for (int i = 0; i < nS; i++)
			fieldS[i] = 1.0;
		theMapper->consistentMapping(fieldS, fieldF);
		for (int i = 0; i < nF; i++)
			CPPUNIT_ASSERT(fabs(fieldF[i] - 1.0) < 1e-13);

	}

	void testMappingPrint() {

		std::cout.precision(15);
		int nS = theIGAMesh->getNumControlPoints();
		int nF = theFEMesh->numNodes;
		double fieldS[nS];
		double fieldF[nF];

		cout << "nS : " << nS << endl;
		cout << "1: consistent Mapping" << endl;
		for (int i = 0; i < nS; i++)
			fieldS[i] = 1;
		theMapper->consistentMapping(fieldS, fieldF);
		for (int i = 0; i < nF; i++)
			cout << i << "  :   " << fieldF[i] << endl;
		cout << endl;

		cout << "nS : " << nS << endl;
		cout << "2: consistent Mapping" << endl;
		for (int i = 0; i < nS; i++)
			fieldS[i] = i;
		theMapper->consistentMapping(fieldS, fieldF);
		for (int i = 0; i < nF; i++)
			cout << i << "  :   " << fieldF[i] << endl;
		cout << endl;

		cout << "nS : " << nS << endl;
		cout << "3: Conservative Mapping" << endl;
		for (int i = 0; i < nF; i++)
			fieldF[i] = 2;
		theMapper->conservativeMapping(fieldF, fieldS);

		cout << "nS : " << nS << endl;
		for (int i = 0; i < nS; i++)
			cout << i << "  :   " << fieldS[i] << endl;
		cout << endl;

	}


// Make the tests
CPPUNIT_TEST_SUITE(TestIGAMortarMapperCube);

//	CPPUNIT_TEST(testMapping);
//	CPPUNIT_TEST(testMappingPrint);
	CPPUNIT_TEST_SUITE_END()
	;
}
;

} /* namespace EMPIRE */
CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAMortarMapperCube);
