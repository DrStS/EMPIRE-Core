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

class TestProjectionSemiSphereFine: public CppUnit::TestFixture {

private:
    IGAPatchSurface* theIGAPatchSurface;
    double Tol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB) for the functional values
        Tol = 1e7;

        // Provide an id for the basis
        int id_basis = 3;

        // Provide an id for the patch itself
        int id_patch = 1;

        // The polynomial degrees
        int p = 2;
        int q = 2;

        // Number of knots in both directions
        int uNoKnots = 10;
        int vNoKnots = 10;

        // The knot vectors in each direction
        double* uKnotVector = new double[uNoKnots];
        for (int i = 0; i < 3; i++)
            uKnotVector[i] = 0;
        uKnotVector[3] = 0.2;
        uKnotVector[4] = 0.4;
        uKnotVector[5] = 0.6;
        uKnotVector[6] = 0.8;
        for (int i = 7; i < 10; i++)
            uKnotVector[i] = 1;

        double* vKnotVector = new double[vNoKnots];
        for (int i = 0; i < 10; i++)
            vKnotVector[i] = uKnotVector[i];

        // The Control Point net
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;
        IGAControlPoint** controlPointNet =
                new IGAControlPoint*[uNoControlPoints * vNoControlPoints];

// Control Points for a NURBS
        controlPointNet[0] = new IGAControlPoint(1, 0.00000000000000, -1.00000000000000,
                0.00000000000000, 1.00000000000000);
        controlPointNet[1] = new IGAControlPoint(2, 0.00000000000000, -1.00000000000000,
                0.15022110482233, 0.94142135623731);
        controlPointNet[2] = new IGAControlPoint(3, 0.00000000000000, -0.90816493864283,
                0.44898963185701, 0.87112698372208);
        controlPointNet[3] = new IGAControlPoint(4, 0.00000000000000, -0.71687947785800,
                0.71687947785799, 0.84769552621700);
        controlPointNet[4] = new IGAControlPoint(5, 0.00000000000000, -0.44898963185701,
                0.90816493864283, 0.87112698372208);
        controlPointNet[5] = new IGAControlPoint(6, 0.00000000000000, -0.15022110482233,
                1.00000000000000, 0.94142135623731);
        controlPointNet[6] = new IGAControlPoint(7, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 1.00000000000000);
        controlPointNet[7] = new IGAControlPoint(8, 0.15022110482233, -1.00000000000000,
                0.00000000000000, 0.94142135623731);
        controlPointNet[8] = new IGAControlPoint(9, 0.15022110482233, -1.00000000000000,
                0.15022110482233, 0.88627416997970);
        controlPointNet[9] = new IGAControlPoint(10, 0.13642554044383, -0.90816493864283,
                0.44898963185701, 0.82009754647056);
        controlPointNet[10] = new IGAControlPoint(11, 0.10769042718829, -0.71687947785800,
                0.71687947785800, 0.79803867196751);
        controlPointNet[11] = new IGAControlPoint(12, 0.06744771855133, -0.44898963185701,
                0.90816493864283, 0.82009754647056);
        controlPointNet[12] = new IGAControlPoint(13, 0.02256638033404, -0.15022110482233,
                1.00000000000000, 0.88627416997970);
        controlPointNet[13] = new IGAControlPoint(14, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 0.94142135623731);
        controlPointNet[14] = new IGAControlPoint(15, 0.44898963185701, -0.90816493864283,
                0.00000000000000, 0.87112698372208);
        controlPointNet[15] = new IGAControlPoint(16, 0.44898963185701, -0.90816493864283,
                0.15022110482233, 0.82009754647056);
        controlPointNet[16] = new IGAControlPoint(17, 0.40775664146669, -0.82476355578014,
                0.44898963185701, 0.75886222176873);
        controlPointNet[17] = new IGAControlPoint(18, 0.32187145284930, -0.65104480702321,
                0.71687947785799, 0.73845044686812);
        controlPointNet[18] = new IGAControlPoint(19, 0.20159168951509, -0.40775664146669,
                0.90816493864283, 0.75886222176873);
        controlPointNet[19] = new IGAControlPoint(20, 0.06744771855133, -0.13642554044383,
                1.00000000000000, 0.82009754647056);
        controlPointNet[20] = new IGAControlPoint(21, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 0.87112698372208);
        controlPointNet[21] = new IGAControlPoint(22, 0.71687947785799, -0.71687947785800,
                0.00000000000000, 0.84769552621700);
        controlPointNet[22] = new IGAControlPoint(23, 0.71687947785799, -0.71687947785800,
                0.15022110482233, 0.79803867196751);
        controlPointNet[23] = new IGAControlPoint(24, 0.65104480702321, -0.65104480702321,
                0.44898963185701, 0.73845044686812);
        controlPointNet[24] = new IGAControlPoint(25, 0.51391618577395, -0.51391618577395,
                0.71687947785800, 0.71858770516832);
        controlPointNet[25] = new IGAControlPoint(26, 0.32187145284930, -0.32187145284930,
                0.90816493864283, 0.73845044686812);
        controlPointNet[26] = new IGAControlPoint(27, 0.10769042718829, -0.10769042718829,
                1.00000000000000, 0.79803867196751);
        controlPointNet[27] = new IGAControlPoint(28, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 0.84769552621700);
        controlPointNet[28] = new IGAControlPoint(29, 0.90816493864283, -0.44898963185701,
                0.00000000000000, 0.87112698372208);
        controlPointNet[29] = new IGAControlPoint(30, 0.90816493864283, -0.44898963185701,
                0.15022110482233, 0.82009754647056);
        controlPointNet[30] = new IGAControlPoint(31, 0.82476355578014, -0.40775664146669,
                0.44898963185701, 0.75886222176873);
        controlPointNet[31] = new IGAControlPoint(32, 0.65104480702321, -0.32187145284930,
                0.71687947785800, 0.73845044686812);
        controlPointNet[32] = new IGAControlPoint(33, 0.40775664146669, -0.20159168951509,
                0.90816493864283, 0.75886222176873);
        controlPointNet[33] = new IGAControlPoint(34, 0.13642554044383, -0.06744771855133,
                1.00000000000000, 0.82009754647056);
        controlPointNet[34] = new IGAControlPoint(35, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 0.87112698372208);
        controlPointNet[35] = new IGAControlPoint(36, 1.00000000000000, -0.15022110482233,
                0.00000000000000, 0.94142135623731);
        controlPointNet[36] = new IGAControlPoint(37, 1.00000000000000, -0.15022110482233,
                0.15022110482233, 0.88627416997970);
        controlPointNet[37] = new IGAControlPoint(38, 0.90816493864283, -0.13642554044383,
                0.44898963185701, 0.82009754647056);
        controlPointNet[38] = new IGAControlPoint(39, 0.71687947785800, -0.10769042718829,
                0.71687947785800, 0.79803867196751);
        controlPointNet[39] = new IGAControlPoint(40, 0.44898963185701, -0.06744771855133,
                0.90816493864283, 0.82009754647056);
        controlPointNet[40] = new IGAControlPoint(41, 0.15022110482233, -0.02256638033404,
                1.00000000000000, 0.88627416997970);
        controlPointNet[41] = new IGAControlPoint(42, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 0.94142135623731);
        controlPointNet[42] = new IGAControlPoint(43, 1.00000000000000, 0.00000000000000,
                0.00000000000000, 1.00000000000000);
        controlPointNet[43] = new IGAControlPoint(44, 1.00000000000000, 0.00000000000000,
                0.15022110482233, 0.94142135623731);
        controlPointNet[44] = new IGAControlPoint(45, 0.90816493864283, 0.00000000000000,
                0.44898963185701, 0.87112698372208);
        controlPointNet[45] = new IGAControlPoint(46, 0.71687947785800, 0.00000000000000,
                0.71687947785799, 0.84769552621700);
        controlPointNet[46] = new IGAControlPoint(47, 0.44898963185701, 0.00000000000000,
                0.90816493864283, 0.87112698372208);
        controlPointNet[47] = new IGAControlPoint(48, 0.15022110482233, 0.00000000000000,
                1.00000000000000, 0.94142135623731);
        controlPointNet[48] = new IGAControlPoint(49, 0.00000000000000, 0.00000000000000,
                1.00000000000000, 1.00000000000000);

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
    void testProjectionOnIGAPatch1() {

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { 0.09464634, -0.0608348, 0.9925547 };
        // Initial guesses for the Newton-Rapson iteration
        double u,v;

        theIGAPatchSurface->findInitialGuess4PointProjection(u, v, vertex);

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch
        flag = theIGAPatchSurface->computePointProjectionOnPatch(u, v, vertex);

        // Compare the values with the ones from MATLAB
        CPPUNIT_ASSERT(flag);
        CPPUNIT_ASSERT(fabs(u - 0.629734882836684) < Tol);
        CPPUNIT_ASSERT(fabs(v - 0.921928587677731) < Tol);

    }

    void testProjectionOnIGAPatch2() {

        // Initial guesses for the Newton-Rapson iteration
        double u = 0.8;
        double v = 0.8;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { 0.0946463, -0.0608348, 0.992555 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch
        flag = theIGAPatchSurface->computePointProjectionOnPatch(u, v, vertex);

        // Compare the values with the ones from MATLAB

        CPPUNIT_ASSERT(flag);
        CPPUNIT_ASSERT(fabs(u - 0.629734882836684) < Tol);
        CPPUNIT_ASSERT(fabs(v - 0.921928587677731) < Tol);

    }

    void testProjectionOnIGAPatch3() {

        // Initial guesses for the Newton-Rapson iteration
        double u = 0.8;
        double v = 0.8;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { 0.095131799857645, -2.0382349911844e-05, 0.99546468558177 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch
        flag = theIGAPatchSurface->computePointProjectionOnPatch(u, v, vertex);

        CPPUNIT_ASSERT(flag);
        CPPUNIT_ASSERT(fabs(u - 9.99848506406883e-01) < Tol);
        CPPUNIT_ASSERT(fabs(v - 9.33884369415096e-01) < Tol);
        // Compare the values with the ones from MATLAB

    }

// Make the tests
    CPPUNIT_TEST_SUITE (TestProjectionSemiSphereFine);

    CPPUNIT_TEST (testProjectionOnIGAPatch1);
    CPPUNIT_TEST (testProjectionOnIGAPatch2);
    CPPUNIT_TEST (testProjectionOnIGAPatch3);CPPUNIT_TEST_SUITE_END()
    ;
}
;

} /* namespace EMPIRE */

//CPPUNIT_TEST_SUITE_REGISTRATION (EMPIRE::TestProjectionSemiSphereFine);

