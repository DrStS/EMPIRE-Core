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

        int numNodes = 6;

        int* controlPointID = new int[6];
        for (int i = 0; i < numNodes; i++) {
            controlPointID[i] = i + 1;
        }

        theIGAMesh = new IGAMesh("IGAMesh", numNodes);

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

        double* controlPoints1 = new double[4 * 4];
        controlPoints1[0] = 0.0;
        controlPoints1[1] = 0.0;
        controlPoints1[2] = 0.0;
        controlPoints1[3] = 1.0;

        controlPoints1[4] = 0.0;
        controlPoints1[5] = 0.0;
        controlPoints1[6] = 1.0;
        controlPoints1[7] = 1.0;

        controlPoints1[8] = 0.0;
        controlPoints1[9] = 1.0;
        controlPoints1[10] = 0.0;
        controlPoints1[11] = 1.0;

        controlPoints1[12] = 0.0;
        controlPoints1[13] = 1.0;
        controlPoints1[14] = 1.0;
        controlPoints1[15] = 1.0;

        int* dofIndexNet1 = new int[4];
        dofIndexNet1[0] = 0;
        dofIndexNet1[1] = 1;
        dofIndexNet1[2] = 2;
        dofIndexNet1[3] = 3;

        double* controlPoints2 = new double[4 * 4];
        controlPoints2[0] = 0.0;
        controlPoints2[1] = 1.0;
        controlPoints2[2] = 0.0;
        controlPoints2[3] = 1.0;

        controlPoints2[4] = 0.0;
        controlPoints2[5] = 1.0;
        controlPoints2[6] = 1.0;
        controlPoints2[7] = 1.0;

        controlPoints2[8] = 1.0;
        controlPoints2[9] = 1.0;
        controlPoints2[10] = 0.0;
        controlPoints2[11] = 1.0;

        controlPoints2[12] = 1.0;
        controlPoints2[13] = 1.0;
        controlPoints2[14] = 1.0;
        controlPoints2[15] = 1.0;

        int* dofIndexNet2 = new int[4];
        dofIndexNet2[0] = 2;
        dofIndexNet2[1] = 3;
        dofIndexNet2[2] = 4;
        dofIndexNet2[3] = 5;

        theIGAMesh->addPatch(p, uNoKnots, uKnotVector1, q, vNoKnots, vKnotVector1, uNoControlPoints,
                vNoControlPoints, controlPoints1, dofIndexNet1);
        theIGAMesh->addPatch(p, uNoKnots, uKnotVector2, q, vNoKnots, vKnotVector2, uNoControlPoints,
                vNoControlPoints, controlPoints2, dofIndexNet2);

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

        CPPUNIT_ASSERT(patch1->getControlPointNet()[0]->getDofIndex() == 0);
        CPPUNIT_ASSERT(patch1->getControlPointNet()[1]->getDofIndex() == 1);
        CPPUNIT_ASSERT(patch1->getControlPointNet()[2]->getDofIndex() == 2);
        CPPUNIT_ASSERT(patch1->getControlPointNet()[3]->getDofIndex() == 3);

        CPPUNIT_ASSERT(patch2->getControlPointNet()[0]->getDofIndex() == 2);
        CPPUNIT_ASSERT(patch2->getControlPointNet()[1]->getDofIndex() == 3);
        CPPUNIT_ASSERT(patch2->getControlPointNet()[2]->getDofIndex() == 4);
        CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getDofIndex() == 5);

        CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getX() == 1.0);
        CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getY() == 1.0);
        CPPUNIT_ASSERT(patch2->getControlPointNet()[3]->getZ() == 1.0);

//		theIGAMesh->getSurfacePatches()[0]->printSelf();
//		theIGAMesh->getSurfacePatches()[1]->printSelf();

    }

    void testLeakage() {
        for (int i = 0; i < 1e10; i++) {

            // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
            Tol = 1e-15;

            // Assign a relaxed tolerance value (corresponding to maximum accuracy provided by MATLAB) for the Newton-Rapson iteration error (accumulative error appears here)
            relTol = 1e-14;

            // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the derivative functional values
            TolDeriv = 1e-13;

            int numNodes = 6;

            int* controlPointID = new int[6];
            for (int i = 0; i < numNodes; i++) {
                controlPointID[i] = i + 1;
            }

            theIGAMesh = new IGAMesh("IGAMesh", numNodes);

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

            double* controlPoints1 = new double[4 * 4];
            controlPoints1[0] = 0.0;
            controlPoints1[1] = 0.0;
            controlPoints1[2] = 0.0;
            controlPoints1[3] = 1.0;

            controlPoints1[4] = 0.0;
            controlPoints1[5] = 0.0;
            controlPoints1[6] = 1.0;
            controlPoints1[7] = 1.0;

            controlPoints1[8] = 0.0;
            controlPoints1[9] = 1.0;
            controlPoints1[10] = 0.0;
            controlPoints1[11] = 1.0;

            controlPoints1[12] = 0.0;
            controlPoints1[13] = 1.0;
            controlPoints1[14] = 1.0;
            controlPoints1[15] = 1.0;

            int* dofIndexNet1 = new int[4];
            dofIndexNet1[0] = 0;
            dofIndexNet1[1] = 1;
            dofIndexNet1[2] = 2;
            dofIndexNet1[3] = 3;

            double* controlPoints2 = new double[4 * 4];
            controlPoints2[0] = 0.0;
            controlPoints2[1] = 1.0;
            controlPoints2[2] = 0.0;
            controlPoints2[3] = 1.0;

            controlPoints2[4] = 0.0;
            controlPoints2[5] = 1.0;
            controlPoints2[6] = 1.0;
            controlPoints2[7] = 1.0;

            controlPoints2[8] = 1.0;
            controlPoints2[9] = 1.0;
            controlPoints2[10] = 0.0;
            controlPoints2[11] = 1.0;

            controlPoints2[12] = 1.0;
            controlPoints2[13] = 1.0;
            controlPoints2[14] = 1.0;
            controlPoints2[15] = 1.0;

            int* dofIndexNet2 = new int[4];
            dofIndexNet2[0] = 2;
            dofIndexNet2[1] = 3;
            dofIndexNet2[2] = 4;
            dofIndexNet2[3] = 5;

            theIGAMesh->addPatch(p, uNoKnots, uKnotVector1, q, vNoKnots, vKnotVector1,
                    uNoControlPoints, vNoControlPoints, controlPoints1, dofIndexNet1);
            theIGAMesh->addPatch(p, uNoKnots, uKnotVector2, q, vNoKnots, vKnotVector2,
                    uNoControlPoints, vNoControlPoints, controlPoints2, dofIndexNet2);

            delete theIGAMesh;

        }
    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestIGAMesh);
    CPPUNIT_TEST (testConstructor);
//    CPPUNIT_TEST (testLeakage);
    CPPUNIT_TEST_SUITE_END();
}
;

} /* namespace EMPIRE */

//CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestIGAMesh);

