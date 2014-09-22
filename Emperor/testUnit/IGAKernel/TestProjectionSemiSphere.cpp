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
#include "IGAPatchSurface.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class IGAPatchSurface
 ***********/

class TestProjectionSemiSphere: public CppUnit::TestFixture {

private:
    IGAPatchSurface* theIGAPatchSurface;
    double Tol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
        Tol = 1e-13;

        // Provide an id for the basis
        int id_basis = 3;

        // Provide an id for the patch itself
        int id_patch = 1;

        // The polynomial degrees
        int p = 2;
        int q = 2;

        // Number of knots in both directions
        int uNoKnots = 6;
        int vNoKnots = 6;

        // The knot vectors in each direction
        double* uKnotVector = new double[uNoKnots];
        for (int i = 0; i < 3; i++)
            uKnotVector[i] = 0;
        for (int i = 3; i < 6; i++)
            uKnotVector[i] = 1;
        double* vKnotVector = new double[vNoKnots];
        for (int i = 0; i < 3; i++)
            vKnotVector[i] = 0;
        for (int i = 3; i < 6; i++)
            vKnotVector[i] = 1;

        // The Control Point net
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;
        IGAControlPoint** controlPointNet =
                new IGAControlPoint*[uNoControlPoints * vNoControlPoints];

// Control Points for a NURBS

        controlPointNet[0] = new IGAControlPoint(0, 0.0, -1.0, 0.0, 1.0);
        controlPointNet[1] = new IGAControlPoint(3, 0.0, -1.0, 1.0, 0.707106781186548);
        controlPointNet[2] = new IGAControlPoint(6, 0.0, 0.0, 1.0, 1.0);
        controlPointNet[3] = new IGAControlPoint(1, 1.0, -1.0, 0.0, 0.707106781186548);
        controlPointNet[4] = new IGAControlPoint(4, 1.0, -1.0, 1.0, 0.5);
        controlPointNet[5] = new IGAControlPoint(7, 0.0, 0.0, 1.0, 0.707106781186548);
        controlPointNet[6] = new IGAControlPoint(2, 1.0, 0.0, 0.0, 1.0);
        controlPointNet[7] = new IGAControlPoint(5, 1.0, 0.0, 1.0, 0.707106781186548);
        controlPointNet[8] = new IGAControlPoint(8, 0.0, 0.0, 1.0, 1.0);

        theIGAPatchSurface = new IGAPatchSurface(id_basis, p, uNoKnots, uKnotVector, q, vNoKnots,
                vKnotVector, uNoControlPoints, vNoControlPoints, controlPointNet);
    }

    void tearDown() {
        delete theIGAPatchSurface;
    }

    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/

    /***********************************************************************************************
     * \brief Test case: Test the projection of an arbitrary point on the 2D IGA patch
     ***********/
    void testProjectionOnIGAPatch() {

        // Initial guesses for the Newton-Rapson iteration
        double u = 0.8;
        double v = 0.8;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { 0.0946463, -0.0608348, 0.992555 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 0;

        // Compute the orthogonal projection of the point on the NURBS patch
        flag = theIGAPatchSurface->computePointProjectionOnPatch(u, v, vertex);

        // Compare the values with the ones from MATLAB
        CPPUNIT_ASSERT(flag);
        CPPUNIT_ASSERT(fabs(u - 9.21928587677731e-01) < Tol);
        CPPUNIT_ASSERT(fabs(v - 6.29734882836685e-01) < Tol);

    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestProjectionSemiSphere);
    CPPUNIT_TEST (testProjectionOnIGAPatch);
    CPPUNIT_TEST_SUITE_END()
    ;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestProjectionSemiSphere);

