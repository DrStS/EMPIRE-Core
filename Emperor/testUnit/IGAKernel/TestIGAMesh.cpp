/*  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
 *  Stefan Sicklinger, Tianyang Wang, Andreas Apostolatos, Munich
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
// inclusion of standard libraries   (only if really necessary here in *.h)
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>
#include <string>
#include <math.h>

// Inclusion of user-defined libraries
#include "IGAMesh.h"
#include "IGAPatchSurface.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class IGAPatchSurface
 ***********/

class TestIGAMesh: public CppUnit::TestFixture {

private:
	IGAMesh* theIGAMesh;
	double Tol;
	double relTol;
	double TolDeriv;

public:
	void setUp() {
		// Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
		Tol = 1e-15;

		// Assign a relaxed tolerance value (corresponding to maximum accuracy provided by MATLAB) for the Newton-Rapson iteration error (accumulative error appears here)
		relTol = 1e-14;

		// Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the derivative functional values
		TolDeriv = 1e-13;

		int numControlPoints = 6;
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
		globalControlPoints[17] = 1.0;
		globalControlPoints[18] = 0.0;
		globalControlPoints[19] = 1.0;

		globalControlPoints[20] = 1.0;
		globalControlPoints[21] = 1.0;
		globalControlPoints[22] = 1.0;
		globalControlPoints[23] = 1.0;

		int* controlPointID = new int[6];
		for (int i = 0; i < numControlPoints; i++) {
			controlPointID[i] = i + 1;
		}

		theIGAMesh = new IGAMesh("IGAMesh", numControlPoints,
				globalControlPoints, controlPointID);

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
		double* vKnotVector1 = new double[uNoKnots];
		double* uKnotVector2 = new double[uNoKnots];
		double* vKnotVector2 = new double[uNoKnots];
		uKnotVector1[0] = 0;
		uKnotVector1[1] = 0;
		uKnotVector1[2] = 1;
		uKnotVector1[3] = 1;
		for (int i = 0; i < 4; i++) {
			vKnotVector1[i] = uKnotVector1[i];
			uKnotVector2[i] = uKnotVector1[i];
			vKnotVector2[i] = uKnotVector1[i];
		}
		int* controlPointNet1 = new int[4];
		int* controlPointNet2 = new int[4];
		controlPointNet1[0] = 1;
		controlPointNet1[1] = 2;
		controlPointNet1[2] = 3;
		controlPointNet1[3] = 4;

		controlPointNet2[0] = 3;
		controlPointNet2[1] = 4;
		controlPointNet2[2] = 5;
		controlPointNet2[3] = 6;

		theIGAMesh->addPatch(p, uNoKnots, vKnotVector1, q, vNoKnots,
				vKnotVector1, uNoControlPoints, vNoControlPoints,
				controlPointNet1);
		theIGAMesh->addPatch(p, uNoKnots, vKnotVector2, q, vNoKnots,
				vKnotVector2, uNoControlPoints, vNoControlPoints,
				controlPointNet2);

	}

	void tearDown() {
		delete theIGAMesh;
	}

	/***********************************************************************************************
	 * \brief Test case: Test the constructor
	 ***********/
	void testConstructor() {
		IGAPatchSurface* patch1 = theIGAMesh->getSurfacePatches()[0];
		IGAPatchSurface* patch2 = theIGAMesh->getSurfacePatches()[1];
		CPPUNIT_ASSERT(patch1->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() == 1);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree() == 1);

		CPPUNIT_ASSERT(patch1->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[0] == 0);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[1] == 0);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[2] == 1);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[3] == 1);

		CPPUNIT_ASSERT(patch1->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[0] == 0);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[1] == 0);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[2] == 1);
		CPPUNIT_ASSERT(patch1->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[3] == 1);

		CPPUNIT_ASSERT(patch1->getControlPointNet()[0]->getId() == 	1);
		CPPUNIT_ASSERT(patch1->getControlPointNet()[1]->getId() == 	2);
		CPPUNIT_ASSERT(patch1->getControlPointNet()[2]->getId() == 	3);
		CPPUNIT_ASSERT(patch1->getControlPointNet()[3]->getId() == 	4);

		CPPUNIT_ASSERT(patch2->getControlPointNet()[0]->getId() == 	3);
		CPPUNIT_ASSERT(patch2->getControlPointNet()[1]->getId() == 	4);
		CPPUNIT_ASSERT(patch2->getControlPointNet()[2]->getId() == 	5);
		CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getId() == 	6);

		CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getX() == 1.0);
		CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getY() == 1.0);
		CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getZ() == 1.0);

//		theIGAMesh->getSurfacePatches()[0]->printSelf();
//		theIGAMesh->getSurfacePatches()[1]->printSelf();

	}

// Make the tests
	CPPUNIT_TEST_SUITE (TestIGAMesh);
	CPPUNIT_TEST (testConstructor);CPPUNIT_TEST_SUITE_END()
	;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMesh);

