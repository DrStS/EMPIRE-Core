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
// inclusion of standard libraries   (only if really necessary here in *.h)
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>
#include <string>
#include <math.h>

// Inclusion of user-defined libraries
#include "IGAPatch2D.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class IGAPatch2D
 ***********/

class TestIGAPatch2D: public CppUnit::TestFixture {

private:
    IGAPatch2D* igaPatch2D;
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

        // Provide an id for the basis
        int id_basis = 3;

        // Provide an id for the patch itself
        int id_patch = 1;

        // The polynomial degrees
        int p = 4;
        int q = 3;

        // Number of knots in both directions
        int uNoKnots = 15;
        int vNoKnots = 13;

        // The knot vectors in each direction
        double* uKnotVector = new double[uNoKnots];
        for (int i = 0; i < 5; i++)
            uKnotVector[i] = -4;
        for (int i = 5; i < 7; i++)
            uKnotVector[i] = 5;
        for (int i = 7; i < 10; i++)
            uKnotVector[i] = 11;
        for (int i = 10; i < 15; i++)
            uKnotVector[i] = 21;
        double* vKnotVector = new double[vNoKnots];
        for (int i = 0; i < 4; i++)
            vKnotVector[i] = -11;
        for (int i = 4; i < 7; i++)
            vKnotVector[i] = -0.1;
        for (int i = 7; i < 9; i++)
            vKnotVector[i] = 0.0;
        for (int i = 9; i < 13; i++)
            vKnotVector[i] = 1.0;

        // The Control Point net
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;
        IGAControlPoint* controlPointNet = new IGAControlPoint[uNoControlPoints * vNoControlPoints];
        double* controlPointWeights = new double[uNoControlPoints * vNoControlPoints];

        // Control Points for a NURBS

        // First row
        controlPointNet[0] = IGAControlPoint(0, 0.0, -1.0, 4.500000000000000, 1.0);
        controlPointNet[1] = IGAControlPoint(1, 0.0, -0.500000000000000, 4.500000000000000, 3.0);
        controlPointNet[2] = IGAControlPoint(2, 0.0, -0.250000000000000, 4.500000000000000, 1.0);
        controlPointNet[3] = IGAControlPoint(3, 0.0, -0.125000000000000, 4.500000000000000, 4.0);
        controlPointNet[4] = IGAControlPoint(4, 0.0, 0.0, 4.500000000000000, 5.0);
        controlPointNet[5] = IGAControlPoint(5, 0.0, 0.125000000000000, 4.500000000000000,
                1.200000000000000);
        controlPointNet[6] = IGAControlPoint(6, 0.0, 0.250000000000000, 4.500000000000000, 1.0);
        controlPointNet[7] = IGAControlPoint(7, 0.0, 0.500000000000000, 4.500000000000000,
                3.400000000000000);
        controlPointNet[8] = IGAControlPoint(8, 0.0, 1.0, 4.500000000000000, 5.600000000000000);

        // Second row
        controlPointNet[9] = IGAControlPoint(9, 0.450000000000000, -1.0, 4.500000000000000,
                0.560000000000000);
        controlPointNet[10] = IGAControlPoint(10, 0.450000000000000, -0.500000000000000,
                4.500000000000000, 1.340000000000000);
        controlPointNet[11] = IGAControlPoint(11, 0.450000000000000, -0.250000000000000,
                4.500000000000000, 1.001000000000000);
        controlPointNet[12] = IGAControlPoint(12, 0.450000000000000, -0.125000000000000,
                4.500000000000000, 4.890000000000000);
        controlPointNet[13] = IGAControlPoint(13, 0.450000000000000, 0.0, 4.500000000000000,
                5.210000000000000);
        controlPointNet[14] = IGAControlPoint(14, 0.450000000000000, 0.125000000000000,
                4.500000000000000, 0.200000000000000);
        controlPointNet[15] = IGAControlPoint(15, 0.450000000000000, 0.250000000000000,
                4.500000000000000, 1.0);
        controlPointNet[16] = IGAControlPoint(16, 0.450000000000000, 0.500000000000000,
                4.500000000000000, 3.400000000000000);
        controlPointNet[17] = IGAControlPoint(17, 0.450000000000000, 1.0, 4.500000000000000,
                5.600000000000000);

        // Third row
        controlPointNet[18] = IGAControlPoint(9, 0.500000000000000, -1.0, 4.500000000000000, 2.6);
        controlPointNet[19] = IGAControlPoint(10, 0.500000000000000, -0.500000000000000,
                4.500000000000000, 0.3);
        controlPointNet[20] = IGAControlPoint(11, 0.500000000000000, -0.250000000000000,
                4.500000000000000, 1.2);
        controlPointNet[21] = IGAControlPoint(12, 0.500000000000000, -0.125000000000000,
                4.500000000000000, 5.8);
        controlPointNet[22] = IGAControlPoint(13, 0.500000000000000, 0.0, 4.500000000000000, 3.2);
        controlPointNet[23] = IGAControlPoint(14, 0.500000000000000, 0.125000000000000,
                4.500000000000000, 6.0);
        controlPointNet[24] = IGAControlPoint(15, 0.500000000000000, 0.250000000000000,
                4.500000000000000, 1.0);
        controlPointNet[25] = IGAControlPoint(16, 0.500000000000000, 0.500000000000000,
                4.500000000000000, 1.4);
        controlPointNet[26] = IGAControlPoint(17, 0.500000000000000, 1.0, 4.500000000000000, 5.0);

        // Fourth row
        controlPointNet[27] = IGAControlPoint(27, 0.562500000000000, -1.0, 4.500000000000000, 1.0);
        controlPointNet[28] = IGAControlPoint(28, 0.562500000000000, -0.500000000000000,
                4.500000000000000, 3.0);
        controlPointNet[29] = IGAControlPoint(29, 0.562500000000000, -0.250000000000000,
                4.500000000000000, 1.0);
        controlPointNet[30] = IGAControlPoint(30, 0.562500000000000, -0.125000000000000,
                4.500000000000000, 4.0);
        controlPointNet[31] = IGAControlPoint(31, 0.562500000000000, 0.0, 4.500000000000000, 5.0);
        controlPointNet[32] = IGAControlPoint(32, 0.562500000000000, 0.125000000000000,
                4.500000000000000, 1.2);
        controlPointNet[33] = IGAControlPoint(33, 0.562500000000000, 0.250000000000000,
                4.500000000000000, 1.0);
        controlPointNet[34] = IGAControlPoint(34, 0.562500000000000, 0.500000000000000,
                4.500000000000000, 3.4);
        controlPointNet[35] = IGAControlPoint(35, 0.562500000000000, 1.0, 4.500000000000000, 5.6);

        // Fifth row
        controlPointNet[36] = IGAControlPoint(36, 0.642857142857143, -1.0, 4.500000000000000, 1.5);
        controlPointNet[37] = IGAControlPoint(37, 0.642857142857143, -0.500000000000000,
                4.500000000000000, 2.3);
        controlPointNet[38] = IGAControlPoint(38, 0.642857142857143, -0.250000000000000,
                4.500000000000000, 1.12);
        controlPointNet[39] = IGAControlPoint(39, 0.642857142857143, -0.125000000000000,
                4.500000000000000, 1.89);
        controlPointNet[40] = IGAControlPoint(40, 0.642857142857143, 0.0, 4.500000000000000, 2.21);
        controlPointNet[41] = IGAControlPoint(41, 0.642857142857143, 0.125000000000000,
                4.500000000000000, 3.2);
        controlPointNet[42] = IGAControlPoint(42, 0.642857142857143, 0.250000000000000,
                4.500000000000000, 4.3);
        controlPointNet[43] = IGAControlPoint(43, 0.642857142857143, 0.500000000000000,
                4.500000000000000, 3.4);
        controlPointNet[44] = IGAControlPoint(44, 0.642857142857143, 1.0, 4.500000000000000, 2.6);

        // Sixth row
        controlPointNet[45] = IGAControlPoint(45, 0.750000000000000, -1.0, 4.500000000000000, 3.0);
        controlPointNet[46] = IGAControlPoint(46, 0.750000000000000, -0.500000000000000,
                4.500000000000000, 5.0);
        controlPointNet[47] = IGAControlPoint(47, 0.750000000000000, -0.250000000000000,
                4.500000000000000, 2.0);
        controlPointNet[48] = IGAControlPoint(48, 0.750000000000000, -0.125000000000000,
                4.500000000000000, 5.0);
        controlPointNet[49] = IGAControlPoint(49, 0.750000000000000, 0.0, 4.500000000000000, 6.0);
        controlPointNet[50] = IGAControlPoint(50, 0.750000000000000, 0.125000000000000,
                4.500000000000000, 7.2);
        controlPointNet[51] = IGAControlPoint(51, 0.750000000000000, 0.250000000000000,
                4.500000000000000, 2.0);
        controlPointNet[52] = IGAControlPoint(52, 0.750000000000000, 0.500000000000000,
                4.500000000000000, 1.4);
        controlPointNet[53] = IGAControlPoint(53, 0.750000000000000, 1.0, 4.500000000000000, 3.6);

        // Seventh row
        controlPointNet[54] = IGAControlPoint(54, 0.900000000000000, -1.0, 4.500000000000000,
                1.101);
        controlPointNet[55] = IGAControlPoint(55, 0.900000000000000, -0.500000000000000,
                4.500000000000000, 3.1232);
        controlPointNet[56] = IGAControlPoint(56, 0.900000000000000, -0.250000000000000,
                4.500000000000000, 1.001);
        controlPointNet[57] = IGAControlPoint(57, 0.900000000000000, -0.125000000000000,
                4.500000000000000, 4.12);
        controlPointNet[58] = IGAControlPoint(58, 0.900000000000000, 0.0, 4.500000000000000, 5.8);
        controlPointNet[59] = IGAControlPoint(59, 0.900000000000000, 0.125000000000000,
                4.500000000000000, 1.2);
        controlPointNet[60] = IGAControlPoint(60, 0.900000000000000, 0.250000000000000,
                4.500000000000000, 2.21);
        controlPointNet[61] = IGAControlPoint(61, 0.900000000000000, 0.500000000000000,
                4.500000000000000, 1.12);
        controlPointNet[62] = IGAControlPoint(62, 0.900000000000000, 1.0, 4.500000000000000, 4.8);

        // Eighth row
        controlPointNet[63] = IGAControlPoint(63, 1.125000000000000, -1.0, 4.500000000000000, 9.0);
        controlPointNet[64] = IGAControlPoint(64, 1.125000000000000, -0.500000000000000,
                4.500000000000000, 3.3);
        controlPointNet[65] = IGAControlPoint(65, 1.125000000000000, -0.250000000000000,
                4.500000000000000, 1.1);
        controlPointNet[66] = IGAControlPoint(66, 1.125000000000000, -0.125000000000000,
                4.500000000000000, 4.4);
        controlPointNet[67] = IGAControlPoint(67, 1.125000000000000, 0.0, 4.500000000000000, 5.5);
        controlPointNet[68] = IGAControlPoint(68, 1.125000000000000, 0.125000000000000,
                4.500000000000000, 1.0);
        controlPointNet[69] = IGAControlPoint(69, 1.125000000000000, 0.250000000000000,
                4.500000000000000, 1.1);
        controlPointNet[70] = IGAControlPoint(70, 1.125000000000000, 0.500000000000000,
                4.500000000000000, 3.0);
        controlPointNet[71] = IGAControlPoint(71, 1.125000000000000, 1.0, 4.500000000000000, 5.5);

        // Nineth row
        controlPointNet[72] = IGAControlPoint(72, 1.500000000000000, -1.0, 4.500000000000000, 0.23);
        controlPointNet[73] = IGAControlPoint(73, 1.500000000000000, -0.500000000000000,
                4.500000000000000, 3.032);
        controlPointNet[74] = IGAControlPoint(74, 1.500000000000000, -0.250000000000000,
                4.500000000000000, 1.01);
        controlPointNet[75] = IGAControlPoint(75, 1.500000000000000, -0.125000000000000,
                4.500000000000000, 4.40404);
        controlPointNet[76] = IGAControlPoint(76, 1.500000000000000, 0.0, 4.500000000000000,
                5.0505);
        controlPointNet[77] = IGAControlPoint(77, 1.500000000000000, 0.125000000000000,
                4.500000000000000, 1.01);
        controlPointNet[78] = IGAControlPoint(78, 1.500000000000000, 0.250000000000000,
                4.500000000000000, 1.0);
        controlPointNet[79] = IGAControlPoint(79, 1.500000000000000, 0.500000000000000,
                4.500000000000000, 3.0);
        controlPointNet[80] = IGAControlPoint(80, 1.500000000000000, 1.0, 4.500000000000000, 3.0);

        // 10th row
        controlPointNet[81] = IGAControlPoint(81, 4.500000000000000, -1.0, 4.500000000000000, 1.02);
        controlPointNet[82] = IGAControlPoint(82, 4.500000000000000, -0.500000000000000,
                4.500000000000000, 1.023);
        controlPointNet[83] = IGAControlPoint(83, 4.500000000000000, -0.250000000000000,
                4.500000000000000, 0.902);
        controlPointNet[84] = IGAControlPoint(84, 4.500000000000000, -0.125000000000000,
                4.500000000000000, 1.102);
        controlPointNet[85] = IGAControlPoint(85, 4.500000000000000, 0.0, 4.500000000000000, 3.2);
        controlPointNet[86] = IGAControlPoint(86, 4.500000000000000, 0.125000000000000,
                4.500000000000000, 4.202);
        controlPointNet[87] = IGAControlPoint(87, 4.500000000000000, 0.250000000000000,
                4.500000000000000, 1.01);
        controlPointNet[88] = IGAControlPoint(88, 4.500000000000000, 0.500000000000000,
                4.500000000000000, 3.401);
        controlPointNet[89] = IGAControlPoint(89, 4.500000000000000, 1.0, 4.500000000000000,
                5.9006);

        for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
            controlPointWeights[i] = controlPointNet[i].getW();

        // Test just one object of the class (that works pretty also)
        igaPatch2D = new IGAPatch2D(id_patch, id_basis, p, uNoKnots, uKnotVector, q, vNoKnots,
                vKnotVector, uNoControlPoints, vNoControlPoints, controlPointNet);

        // nurbsBasis1D->printPolynomialDegree();
        // nurbsBasis1D->printNoKnots();
        // nurbsBasis1D->printKnotVector();
        // nurbsBasis1D->printNoBasisFunctions();
        // nurbsBasis1D->printControlPointNet();
        // NurbsBasis1D nurbsBasis1DCopy = *nurbsBasis1D;
        // delete nurbsBasis1D;

        // Make test on memory leakage (that works pretty fine)
        /* for (int j = 1; j < 1e9; j++) {

         // Copy the knot vectors
         double* uKnotVectorCopy = new double[uNoKnots];
         for (int i = 0; i < uNoKnots; i++)
         uKnotVectorCopy[i] = uKnotVector[i];
         double* vKnotVectorCopy = new double[vNoKnots];
         for (int i = 0; i < vNoKnots; i++)
         vKnotVectorCopy[i] = vKnotVector[i];

         // Copy the Control Point net and the Control Point weights
         IGAControlPoint* controlPointNetCopy = new IGAControlPoint[uNoControlPoints * vNoControlPoints];
         for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
         controlPointNetCopy[i] = controlPointNet[i];

         // Create an object of the class IGAPatch2D
         IGAPatch2D* igaPatch2DCopy = new IGAPatch2D(id_patch, id_basis, p, uNoKnots, uKnotVectorCopy, q,
         vNoKnots, vKnotVectorCopy, uNoControlPoints, vNoControlPoints,
         controlPointNetCopy);

         // Free the memory on the heap from the pointer
         delete igaPatch2DCopy;
         }*/
    }

    void tearDown() {
        delete igaPatch2D;
    }

    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        CPPUNIT_ASSERT(igaPatch2D->getId()==1);
        CPPUNIT_ASSERT(igaPatch2D->getIGABasis()->getId()==3);
        CPPUNIT_ASSERT(igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()==4);
        CPPUNIT_ASSERT(igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()==3);
        CPPUNIT_ASSERT(igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getNoKnots()==15);
        CPPUNIT_ASSERT(igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getNoKnots()==13);
        for (int i = 0; i < 5; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == -4);
        for (int i = 5; i < 7; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 5);
        for (int i = 7; i < 10; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 11);
        for (int i = 10; i < 15; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] == 21);
        for (int i = 0; i < 4; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == -11);
        for (int i = 4; i < 7; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == -0.1);
        for (int i = 7; i < 9; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == 0.0);
        for (int i = 9; i < 13; i++)
            CPPUNIT_ASSERT(
                    igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] == 1.0);
        CPPUNIT_ASSERT(igaPatch2D->getUNoControlPoints()==10);
        CPPUNIT_ASSERT(igaPatch2D->getVNoControlPoints()==9);
    }

    /***********************************************************************************************
     * \brief Test case: Test the knot span
     ***********/
    void testIGAPatch2DKnotSpan() {

        // The parametric coordinates on the Spline parameter space
        double u = 13.3300000033001;
        double v = -5.00998989;
        int uCorrectknotSpan = 9;
        int vCorrectknotSpan = 3;

        // Check the find knot span function
        CPPUNIT_ASSERT(
                igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u)==uCorrectknotSpan);
        CPPUNIT_ASSERT(
                igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v)==vCorrectknotSpan);
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constructor for proper functionality
     ***********/
    void testIGAPatch2DCopyConstructor() {

        IGAPatch2D* igaPatch2DCopy = new IGAPatch2D(*igaPatch2D);

        // Test the ID's
        CPPUNIT_ASSERT(igaPatch2DCopy->getId()==igaPatch2D->getId());
        CPPUNIT_ASSERT(igaPatch2DCopy->getIGABasis()->getId()==igaPatch2D->getIGABasis()->getId());

        // Test the polynomial degrees
        CPPUNIT_ASSERT(
                igaPatch2DCopy->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()==igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree());
        CPPUNIT_ASSERT(
                igaPatch2DCopy->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()==igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree());

        // Test the knot vectors
        CPPUNIT_ASSERT(
                igaPatch2DCopy->getIGABasis()->getUBSplineBasis1D()->getNoKnots()==igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < igaPatch2DCopy->getIGABasis()->getUBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i]==igaPatch2DCopy->getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i]);
        CPPUNIT_ASSERT(
                igaPatch2DCopy->getIGABasis()->getVBSplineBasis1D()->getNoKnots()==igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < igaPatch2DCopy->getIGABasis()->getVBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i]==igaPatch2DCopy->getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i]);

        // Test the Control Points
        CPPUNIT_ASSERT(igaPatch2DCopy->getUNoControlPoints()==igaPatch2D->getUNoControlPoints());
        CPPUNIT_ASSERT(igaPatch2DCopy->getVNoControlPoints()==igaPatch2D->getVNoControlPoints());
        int noCPu = igaPatch2D->getUNoControlPoints();
        int noCPv = igaPatch2D->getVNoControlPoints();
        for (int i = 0; i < noCPu * noCPv; i++) {
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getControlPointNet()[i].getId()==igaPatch2D->getControlPointNet()[i].getId());
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getControlPointNet()[i].getX()==igaPatch2D->getControlPointNet()[i].getX());
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getControlPointNet()[i].getY()==igaPatch2D->getControlPointNet()[i].getY());
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getControlPointNet()[i].getZ()==igaPatch2D->getControlPointNet()[i].getZ());
            CPPUNIT_ASSERT(
                    igaPatch2DCopy->getControlPointNet()[i].getW()==igaPatch2D->getControlPointNet()[i].getW());
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constructor for leakage
     ***********/
    void testIGAPatch2DCopyConstructor4Leakage() {

        for (int i = 1; i < 1e9; i++) {
            IGAPatch2D* igaPatch2DCopy = new IGAPatch2D(*igaPatch2D);
            delete igaPatch2DCopy;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local B-Spline basis functions
     ***********/
    void testIGAPatch2DBSplineBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = 17.002333009;
        double v = -0.000001234;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);
        int uNoLocalBasisFunctions =
                igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1;
        int vNoLocalBasisFunctions =
                igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree() + 1;
        int noLocalBasisFunctions = uNoLocalBasisFunctions * vNoLocalBasisFunctions;

        // The polynomial degrees
        int p = 4;
        int q = 3;

        // Number of knots in both directions
        int uNoKnots = 15;
        int vNoKnots = 13;

        // The number of the Control Points at each parametric direction
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;

        /*
         * Make the patch's underlying basis a B-Spline by setting all the weights to 1. First copy everything to a new
         * instance of the class IGAPatch2D in order not to lose the stored Control Point information.
         *
         *              !!!!!!!!!!!!!!!!!!!!!!!!!
         *              !!!!!   ACHTUNG     !!!!!
         *              !!!!!!!!!!!!!!!!!!!!!!!!!
         *
         * By changing the Control Point weights to 1.0 we still call the computeLocalBasisFunctions from the NurbsBasis2D
         * and not this function from BSplineBasis2D cause this is decided in the constructor of the instance. However it should be
         * computed the B-Spline basis functions but in an inefficient way
         *
         */

        IGAPatch2D* igaPatch2DCopy = new IGAPatch2D(*igaPatch2D);
        for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
            igaPatch2DCopy->getControlPointNet()[i].setW(1.0);

        double* localBasisFunctions = new double[noLocalBasisFunctions];
        igaPatch2DCopy->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, uKnotSpan,
                v, vKnotSpan);

        /* cout << endl;
         cout << endl;
         cout << "The non-zero B-Spline basis functions at (u,v) = ( " << u << " , " << v << " ):";
         cout << endl;
         for (int i = 0; i < noLocalBasisFunctions; i++) {
         cout << localBasisFunctions[i] << " ";
         }
         cout << endl;
         cout << endl;*/

        // Value provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.000000000000000, 0.000000000000000,
                0.000000000000001, 0.000000000000001, 0.000000000000000, 0.000000000007292,
                0.000000000074448, 0.000000000157816, 0.000000000157969, 0.000000000059296,
                0.014511603260064, 0.148153855088449, 0.314060685220477, 0.314366157723487,
                0.118001972754625, 0.001451101233533, 0.014814782213839, 0.031404789640423,
                0.031435335646814, 0.011799716761467 };
        for (int i = 0; i < noLocalBasisFunctions; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
        delete igaPatch2DCopy;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local NURBS basis functions
     ***********/
    void testIGAPatch2DNurbsBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = 17.002333009;
        double v = -0.000001234;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);
        int uNoLocalBasisFunctions =
                igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1;
        int vNoLocalBasisFunctions =
                igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree() + 1;
        int noLocalBasisFunctions = uNoLocalBasisFunctions * vNoLocalBasisFunctions;

        // The polynomial degrees
        int p = 4;
        int q = 3;

        // Number of knots in both directions
        int uNoKnots = 15;
        int vNoKnots = 13;

        // The number of the Control Points at each parametric direction
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;

        double* localBasisFunctions = new double[noLocalBasisFunctions];
        igaPatch2D->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, uKnotSpan, v,
                vKnotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero 2D NURBS basis functions at (u,v) = ( " << u << " , " << v << " ):";
         cout << endl;
         for (int i = 0; i < noLocalBasisFunctions; i++) {
         cout << localBasisFunctions[i] << " ";
         }
         cout << endl;
         cout << endl;*/

        // Value provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.000000000000000, 0.000000000000001,
                0.000000000000002, 0.000000000000002, 0.000000000000000, 0.000000000028724,
                0.000000000283475, 0.000000000569836, 0.000000000523774, 0.000000000124570,
                0.068593675533497, 0.116715997109942, 0.206181528917218, 0.208445893008528,
                0.325522863981816, 0.001905302287253, 0.021494320831880, 0.022679044633442,
                0.020637366823345, 0.007824005342697 };
        for (int i = 0; i < noLocalBasisFunctions; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives (Test 1)
     ***********/
    void testIGAPatch2DBasisFunctionsAndDerivativesTest1() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 5.0;
        double v = -0.1;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // The number of local basis functions
        int noBasisFunctions =
                (igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                        * (igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
                                + 1);

        // Return the functional values R, dR/du, dR/dv
        int derivDegree = 1;
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
                * noBasisFunctions / 2];

        igaPatch2D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v, vKnotSpan);

        /* // Printing the results for convenient debugging
         int vDerivOrder = 2;
         cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives in u-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " <<  vDerivOrder << " order in v-direction : " << endl;
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         cout
         << localBasisFunctionsAndDerivatives[vDerivOrder * (derivDegree + 1)
         * noBasisFunctions + i * noBasisFunctions + counter] << " ";
         counter++;
         }
         counter = 0;
         cout << endl;
         }
         cout << endl;

         int uDerivOrder = 1;
         cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives in v-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " <<  uDerivOrder << " order in u-direction : " << endl;
         cout << endl;

         counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         cout
         << localBasisFunctionsAndDerivatives[i * (derivDegree + 1)
         * noBasisFunctions + uDerivOrder * noBasisFunctions + counter] << " ";
         counter++;
         }
         counter = 0;
         cout << endl;
         }
         cout << endl;*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 0.263008729169028, 0.544155991384197,
                0.192835279446775, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -8.290991647209133, -17.153775821811998, -6.078869306854625, 0, 0,
                4.353247931073575, 20.405849676907383, 6.764539167894797, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, -0.135862499631207, 0.021214057279146, 0.114648442352061, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        // The tolerance value is different for the functional and the derivative values due to high gradients which appear
        int indexBasisDeriv = 0;
        int counter = 0;
        for (int i = 0; i <= derivDegree; i++)
            for (int j = 0; j <= derivDegree - i; j++)
                for (int k = 0; k < noBasisFunctions; k++) {
                    // Get the index of the basis function
                    indexBasisDeriv = igaPatch2D->getIGABasis()->indexDerivativeBasisFunction(
                            derivDegree, i, j, k);

                    // Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(localBasisFunctionsAndDerivatives[indexBasisDeriv]-CorrectlocalBasisFunctionsAndDerivatives[counter])<=Tol);
                    counter++;
                }

        /*for (int i = 0; i < noBasisFunctions * (noDerivs - 4); i++)
         CPPUNIT_ASSERT(
         fabs(localBasisFunctionsAndDerivatives[i]-CorrectlocalBasisFunctionsAndDerivatives[i])<=Tol);

         for (int i = noBasisFunctions * (noDerivs - 4); i < noBasisFunctions * noDerivs; i++)
         CPPUNIT_ASSERT(
         fabs(localBasisFunctionsAndDerivatives[i]-CorrectlocalBasisFunctionsAndDerivatives[i])<=TolDeriv);*/

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives (Test 2)
     ***********/
    void testIGAPatch2DBasisFunctionsAndDerivativesTest2() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 10.1100000099;
        double v = -.44449811111;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // The number of local basis functions
        int noBasisFunctions =
                (igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                        * (igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
                                + 1);

        // Return the functional values R, dR/du, dR/dv, d^2R/du^2
        int derivDegree = 2;

        // Compute the local basis functions and their derivatives
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
                * noBasisFunctions / 2];
        igaPatch2D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v, vKnotSpan);

        /* // Printing the results for convenient debugging
         int vDerivOrder = 2;
         cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives in u-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " <<  vDerivOrder << " order in v-direction : " << endl;
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         cout
         << localBasisFunctionsAndDerivatives[vDerivOrder * (derivDegree + 1)
         * noBasisFunctions + i * noBasisFunctions + counter] << " ";
         counter++;
         }
         counter = 0;
         cout << endl;
         }
         cout << endl;

         int uDerivOrder = 1;
         cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives in v-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " <<  uDerivOrder << " order in u-direction : " << endl;
         cout << endl;

         counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         cout
         << localBasisFunctionsAndDerivatives[i * (derivDegree + 1)
         * noBasisFunctions + uDerivOrder * noBasisFunctions + counter] << " ";
         counter++;
         }
         counter = 0;
         cout << endl;
         }
         cout << endl;*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 0.000000001494750, 0.000000034733248,
                0.000001142264686, 0.000015482614677, 0.000001612194963, 0.000000015853674,
                0.000009578112787, 0.000160996570433, 0.002371954492568, 0.000420381247378,
                0.000001943040915, 0.000097825197827, 0.002402142981880, 0.029070893940158,
                0.004128284878567, 0.000095917860106, 0.003996515615600, 0.041401247148402,
                0.742282434457903, 0.173541595299479, -0.000000013254836, -0.000000308000237,
                -0.000010129136093, -0.000137293495221, -0.000014296285614, -0.000000093062335,
                -0.000056224290060, -0.000945062777637, -0.013923563061895, -0.002467671629555,
                -0.000005581519726, -0.000281009662359, -0.006900322240780, -0.083508158143702,
                -0.011858784501478, 0.000011984124746, 0.000499330798590, 0.005172735399946,
                0.092741907308132, 0.021682553430115, 0.000000079700876, 0.000001851994888,
                0.000060906148819, 0.000825541090056, 0.000085963076331, 0.000000278411293,
                0.000168204218230, 0.002827310856599, 0.041654630717609, 0.007382445858471,
                -0.000000442540616, -0.000022280345704, -0.000547104194573, -0.006621091306643,
                -0.000940244602629, -0.000004477619159, -0.000186564575884, -0.001932685083262,
                -0.034651086318766, -0.008101235485938, -0.000000006878969, -0.000000116297524,
                -0.000002345944302, -0.000000314936181, 0.000001088353830, -0.000000072959955,
                -0.000032070447573, -0.000330649272242, -0.000048248587454, 0.000283789213634,
                -0.000008942039681, -0.000327548646394, -0.004933439430687, -0.000591339156391,
                0.002786905283366, -0.000441422156723, -0.013381555154059, -0.085028471120272,
                -0.015098973891198, 0.117153734068774, 0.000000060984946, 0.000001030931773,
                0.000020791459376, 0.000002637998707, -0.000009667188256, 0.000000428122352,
                0.000188160379250, 0.001939328810523, 0.000259518791082, -0.001670066547248,
                0.000025667210968, 0.000939928570973, 0.014147640792019, 0.001408140796047,
                -0.008046834916190, -0.000056110523722, -0.001711851497747, -0.011037333897326,
                -0.009304538180043, 0.012903067902514, 0.000000024293786, 0.000000261548917,
                0.000002025164269, -0.000003727924282, 0.000000686727844, 0.000000257665580,
                0.000072125274329, 0.000285436909738, -0.000571122315798, 0.000179064886343,
                0.000031579732070, 0.000736645035069, 0.004258850158427, -0.006999727997080,
                0.001758477256494, 0.001558927709635, 0.030094632581332, 0.073401837154189,
                -0.178727738779250, 0.073921484918389 };

        // The tolerance value is different for the functional and the derivative values due to high gradients which appear
        int indexBasisDeriv = 0;
        int counter = 0;
        for (int i = 0; i <= derivDegree; i++)
            for (int j = 0; j <= derivDegree - i; j++)
                for (int k = 0; k < noBasisFunctions; k++) {
                    // Get the index of the basis function
                    indexBasisDeriv = igaPatch2D->getIGABasis()->indexDerivativeBasisFunction(
                            derivDegree, i, j, k);

                    // Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(localBasisFunctionsAndDerivatives[indexBasisDeriv]-CorrectlocalBasisFunctionsAndDerivatives[counter])<=Tol);
                    counter++;
                }

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of a point on the IGA 2D patch
     ***********/
    void testIGAPatch2DPointOnSurface() {

        // The parameters on the IGA surface and their knot span indices
        double u = 8.000000000021;
        double v = -2.11111111231;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // Initialize the Cartesian Coordinates
        double cartesianCoordinatesPointOnPatch[3];

        // Compute the Cartesian Coordinates of the point with the above surface parameters
        igaPatch2D->computeCartesianCoordinates(cartesianCoordinatesPointOnPatch, u, uKnotSpan, v,
                vKnotSpan);

        /*cout << endl;
         cout << "Point on Patch 2D: ";
         for (int i = 0; i < 3; i++) {
         cout << cartesianCoordinatesPointOnPatch[i] << " ";
         }
         cout << endl;*/

        // Values provided by MATLAB
        double CorrectcartesianCoordinatesPointOnPatch[] = { 0.684365187483793, -0.194078057047477,
                4.500000000000000 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(
                    fabs(cartesianCoordinatesPointOnPatch[i]-CorrectcartesianCoordinatesPointOnPatch[i])<=Tol);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of a point on the IGA 2D patch by pre-computing the basis functions
     ***********/
    void testIGAPatch2DPointOnSurfaceMethod2() {

        // The parameters on the IGA surface and their knot span indices
        double u = 8.000000000021;
        double v = -2.11111111231;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // Initialize the Cartesian Coordinates
        double cartesianCoordinatesPointOnPatch[3];

        // Compute the local basis functions
        int pDegree = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        int qDegree = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
        double* localBasisFunctions = new double[noLocalBasisFunctions];
        igaPatch2D->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, uKnotSpan, v,
                vKnotSpan);

        // Compute the Cartesian Coordinates of the point with the above surface parameters
        igaPatch2D->computeCartesianCoordinates(cartesianCoordinatesPointOnPatch,
                localBasisFunctions, uKnotSpan, vKnotSpan);

        /*cout << endl;
         cout << "Point on Patch 2D: ";
         for (int i = 0; i < 3; i++) {
         cout << cartesianCoordinatesPointOnPatch[i] << " ";
         }
         cout << endl;*/

        // Values provided by MATLAB
        double CorrectcartesianCoordinatesPointOnPatch[] = { 0.684365187483793, -0.194078057047477,
                4.500000000000000 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(
                    fabs(cartesianCoordinatesPointOnPatch[i]-CorrectcartesianCoordinatesPointOnPatch[i])<=Tol);

        // Free the memory from the heap
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vectors for the surface patch
     ***********/
    void testIGAPatch2DBaseVectors() {

        // The parameters on the IGA surface and their knot span indices
        double u = 2.222122;
        double v = -3.3339333;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // Modify the Control Point net for getting some curvature into the structure
        igaPatch2D->getControlPointNet()[igaPatch2D->getVNoControlPoints() + 4].setZ(
                102.78054 * 4.500000000000000);
        igaPatch2D->getControlPointNet()[4 * igaPatch2D->getVNoControlPoints()].setZ(
                102.78054 * 4.500000000000000 / 1.5672);

        // Get the polynomial degree of the basis in each direction
        int pDegree = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        int qDegree = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();
        int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

        // Derivative order of the basis functions needed for the computation of the base vectors
        int derivDegree = 1;

        // Compute the local basis functions and their derivatives
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
                * noLocalBasisFunctions / 2];
        igaPatch2D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v, vKnotSpan);

        // Compute the base vectors at the given surface parameters
        double baseVectors[6];
        igaPatch2D->computeBaseVectors(baseVectors, localBasisFunctionsAndDerivatives, uKnotSpan,
                vKnotSpan);

        /* cout << endl;
         cout << "The base vectors at (u,v) = ( " << u << " , " << v << " ) are:" << endl;
         int counter = 0;
         for (int i = 0; i < 2; i++) {
         for (int j = 0; j < 3; j++) {
         cout << baseVectors[counter] << " ";
         counter++;
         }
         cout << endl;
         }
         cout << endl;*/

        // Values provided by MATLAB
        double correctBaseVectors[] = { 0.021609817639182, -0.008848771788876, 0.252711306398794,
                -0.003022870519747, 0.052433901953327, -0.419630625272796 };

        for (int i = 0; i < 6; i++)
            CPPUNIT_ASSERT( fabs(baseVectors[i]-correctBaseVectors[i])<=Tol);

        // Free the memory from the heap
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vectors and their derivatives for the surface patch
     ***********/
    void testIGAPatch2DBaseVectorsAndDerivatives() {

        // Testing if the implemented index works correctly
        /* int derivDegree = 2;
         int noBaseVcts = 2;
         int noCoord = 3;

         // Initialize index
         int indexBaseVec = 0;>

         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= derivDegree - i; j++) {
         for (int k = 0; k < 3; k++) {
         for (int l = 0; l < 2; l++) {
         indexBaseVec = igaPatch2D->indexDerivativeBaseVector(derivDegree, i, j, k,
         l);
         cout << "Indexing quadruple (i,j,k,l) = ( " << i << " , " << j << " , " << k
         << " , " << l << ") and indexing number indexNum = " << indexBaseVec
         << endl;
         }
         }
         }
         }*/

        //
        // The parameters on the IGA surface and their knot span indices
        double u = 9.0012100121;
        double v = -1.00010001;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // Modify the Control Point net for getting some curvature into the structure
        igaPatch2D->getControlPointNet()[igaPatch2D->getVNoControlPoints() + 4].setZ(
                102.78054 * 4.500000000000000);
        igaPatch2D->getControlPointNet()[4 * igaPatch2D->getVNoControlPoints()].setZ(
                102.78054 * 4.500000000000000 / 1.5672);

        // Get the polynomial degrees of the patch
        int pDegree = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree();
        int qDegree = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree();

        // Decide on absolute degree of the partial derivatives for the base vectors
        int derivDegreeBaseVec = 1;

        // Initialize the index of the base vector
        int indexBaseVct = 0;

        // Initialize the indices related to the derivatives of each component
        int uDeriv = 1;
        int vDeriv = 1;

        // Indices related to the base vectors
        int indexUBaseVec = 0;
        int indexVBaseVec = 1;

        // Compute the local basis functions and their derivatives
        int noBasisFcts = (pDegree + 1) * (qDegree + 1);
        int derivDegree = derivDegreeBaseVec + 1;
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
                * noBasisFcts / 2];
        igaPatch2D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v, vKnotSpan);

        // Compute the base vectors and their derivatives
        int noBaseVec = 2;
        int noCoord = 3;
        double* baseVctsAndDerivs = new double[(derivDegreeBaseVec + 1) * (derivDegreeBaseVec + 2)
                * noCoord * noBaseVec / 2];
        igaPatch2D->computeBaseVectorsAndDerivatives(baseVctsAndDerivs,
                localBasisFunctionsAndDerivatives, derivDegreeBaseVec, uKnotSpan, vKnotSpan);

        // Print the 1st base vector and its derivatives
        /*cout << endl;
         cout << "All partial derivatives of the base vector g1 = dX/du up to order "
         << derivDegreeBaseVec << ":" << endl;
         cout << "_________________________________________________________________ " << endl;
         cout << endl;
         for (int i = 0; i <= derivDegreeBaseVec; i++) {
         for (int j = 0; j <= derivDegreeBaseVec - i; j++) {
         cout << "The " << i << "-th w.r.t u partial derivative and " << j
         << "-th w.r.t v partial derivative: \t";
         for (int k = 0; k < noCoord; k++) {
         // Find the index of the derivatives to the base vector
         indexBaseVct = igaPatch2D->indexDerivativeBaseVector(derivDegreeBaseVec, i, j,
         k, indexUBaseVec);

         // print the respective functional value
         cout << baseVctsAndDerivs[indexBaseVct] << " ";
         }
         cout << endl;
         cout << endl;
         }
         }

         // Print the 2nd base vector and its derivatives
         cout << endl;
         cout << "All partial derivatives of the base vector g2 = dX/dv up to order "
         << derivDegreeBaseVec << ":" << endl;
         cout << "_________________________________________________________________ " << endl;
         cout << endl;
         for (int i = 0; i <= derivDegreeBaseVec; i++) {
         for (int j = 0; j <= derivDegreeBaseVec - i; j++) {
         cout << "The " << i << "-th w.r.t u partial derivative and " << j
         << "-th w.r.t v partial derivative: \t";
         for (int k = 0; k < noCoord; k++) {
         // Find the index of the derivatives to the base vector
         indexBaseVct = igaPatch2D->indexDerivativeBaseVector(derivDegreeBaseVec, i, j,
         k, indexVBaseVec);

         // print the respective functional value
         cout << baseVctsAndDerivs[indexBaseVct] << " ";
         }
         cout << endl;
         cout << endl;
         }
         }*/

        // Tests on the pointer baseVctsAndDerivs
        /* cout << endl;
         cout << "Debugging the pointer baseVctsAndDerivs" << endl;
         cout << "_______________________________________" << endl;
         cout << endl;
         for (int i = 0;
         i < (derivDegreeBaseVec + 1) * (derivDegreeBaseVec + 2) * noCoord * noBaseVec / 2 ;
         i++) {
         cout << baseVctsAndDerivs[i] << endl;
         }
         cout << endl;
         cout << endl;*/

        // Analytical values provided by MATLAB
        double correctBaseVctGuAndDerivs[] = { 0.040872970530277, 0.000856207886370,
                -0.021988101482595, 0.000713760263798, -0.000916996592231, 0.075879725575611,
                -0.014073318985134, 0.000629989330548, 0.003707201239999 };
        double correctBaseVctGvAndDerivs[] = { 0.001165581187192, 0.032039673000424,
                -0.092289772055987, -0.000045712613992, -0.020522438760949, 0.219508056298072,
                0.000713760263798, -0.000916996592231, 0.075879725575611 };

        // Testing the base vector Gu and its derivatives
        int counterBaseVec = 0;
        for (int i = 0; i <= derivDegreeBaseVec; i++) {
            for (int j = 0; j <= derivDegreeBaseVec - i; j++) {
                for (int k = 0; k < noCoord; k++) {
                    // Find the index of the derivatives to the base vector
                    indexBaseVct = igaPatch2D->indexDerivativeBaseVector(derivDegreeBaseVec, i, j,
                            k, indexUBaseVec);

                    // Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(baseVctsAndDerivs[indexBaseVct]-correctBaseVctGuAndDerivs[counterBaseVec])<=Tol);

                    // Update counter
                    counterBaseVec++;
                }
            }
        }

        // Testing the base vector Gv and its derivatives
        counterBaseVec = 0;
        for (int i = 0; i <= derivDegreeBaseVec; i++) {
            for (int j = 0; j <= derivDegreeBaseVec - i; j++) {
                for (int k = 0; k < noCoord; k++) {
                    // Find the index of the derivatives to the base vector
                    indexBaseVct = igaPatch2D->indexDerivativeBaseVector(derivDegreeBaseVec, i, j,
                            k, indexVBaseVec);

                    // Assert a message if failure occurs
                    CPPUNIT_ASSERT(
                            fabs(baseVctsAndDerivs[indexBaseVct]-correctBaseVctGvAndDerivs[counterBaseVec])<=Tol);

                    // Update counter
                    counterBaseVec++;
                }
            }
        }

        // Free the memory from the heap
        delete[] localBasisFunctionsAndDerivatives;
        delete[] baseVctsAndDerivs;
    }

    /***********************************************************************************************
     * \brief Test case: Test the projection of an arbitrary point on the 2D IGA patch
     ***********/
    void testProjectionOnIGAPatch() {
        // Modify the Control Point net for getting some curvature into the structure
        igaPatch2D->getControlPointNet()[igaPatch2D->getVNoControlPoints() + 4].setZ(
                102.78054 * 4.500000000000000);
        igaPatch2D->getControlPointNet()[4 * igaPatch2D->getVNoControlPoints()].setZ(
                102.78054 * 4.500000000000000 / 1.5672);

        // Initial guesses for the Newton-Rapson iteration
        double u = 10.0002;
        double v = -6.2365;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { .5, -1.72, 3.8 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch
        flag = igaPatch2D->computePointProjectionOnPatch(u, v, vertex);

        // Compare the values with the ones from MATLAB

        // On the return flag
        CPPUNIT_ASSERT(flag==1);

        // On the Cartesian components of the orthogonal projection
        double orthogonalProjection[] = { 0.735728398615148, -0.142426205779032, 4.516956877569002 };
        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(fabs(orthogonalProjection[i] - vertex[i])<=relTol);

        // On the surface parametric values of the orthogonal projection
        CPPUNIT_ASSERT(fabs(u - 9.089448501939843)<=TolDeriv);
        CPPUNIT_ASSERT(fabs(v + 0.897104716006744)<=TolDeriv);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions for memory leakage
     ***********/
    void testIGAPatch2DBasisFunctions4Leakage() {

        // Initialize variables
        double u = 17.002333009;
        double v = -0.000001234;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);
        int uNoLocalBasisFunctions =
                igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1;
        int vNoLocalBasisFunctions =
                igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree() + 1;
        int noLocalBasisFunctions = uNoLocalBasisFunctions * vNoLocalBasisFunctions;

        // The polynomial degrees
        int p = 4;
        int q = 3;

        // Number of knots in both directions
        int uNoKnots = 15;
        int vNoKnots = 13;

        // The number of the Control Points at each parametric direction
        int uNoControlPoints = uNoKnots - p - 1;
        int vNoControlPoints = vNoKnots - q - 1;

        for (int i = 1; i < 1e9; i++) {
            // Compute the non-zero basis functions at another parametric location
            double* localBasisFunctions = new double[noLocalBasisFunctions];
            igaPatch2D->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, uKnotSpan,
                    v, vKnotSpan);

            // Free the memory from the heap
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (Test 1)
     ***********/
    void testIGAPatch2DBasisFunctionsAndDerivativesTest14Leakage() {

        // Initialize variables
        double u = 5.0;
        double v = -0.1;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // The number of local basis functions
        int noBasisFunctions =
                (igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                        * (igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
                                + 1);

        // Return the functional values R, dR/du, dR/dv
        int derivDegree = 1;

        // Number of iterations
        int noIterations = 1e9;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions and their derivatives at another parametric location
            double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                    * (derivDegree + 1) * noBasisFunctions];

            igaPatch2D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v, vKnotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctionsAndDerivatives;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (Test 2)
     ***********/
    void testIGAPatch2DBasisFunctionsAndDerivativesTest24Leakage() {

        // Initialize variables
        double u = 10.1100000099;
        double v = -.44449811111;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // The number of local basis functions
        int noBasisFunctions =
                (igaPatch2D->getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                        * (igaPatch2D->getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
                                + 1);

        // Return the functional values R, dR/du, dR/dv, d^2R/du^2
        int derivDegree = 2;

        // Number of iterations
        int noIterations = 1e9;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions and their derivatives at another parametric location
            double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                    * (derivDegree + 1) * noBasisFunctions];

            igaPatch2D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, u, uKnotSpan, v, vKnotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctionsAndDerivatives;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of a point on the IGA patch for leakage
     ***********/
    void testIGAPatch2DPointOnSurface4Leakage() {

        // The parameters on the IGA surface and their knot span indices
        double u = 8.000000000021;
        double v = -2.11111111231;
        int uKnotSpan = igaPatch2D->getIGABasis()->getUBSplineBasis1D()->findKnotSpan(u);
        int vKnotSpan = igaPatch2D->getIGABasis()->getVBSplineBasis1D()->findKnotSpan(v);

        // Initialize the Cartesian Coordinates
        double cartesianCoordinatesPointOnPatch[3];

        for (int i = 0; i < 1e9; i++) {
            // Compute the Cartesian Coordinates of the point with the above surface parameters
            igaPatch2D->computeCartesianCoordinates(cartesianCoordinatesPointOnPatch, u, uKnotSpan,
                    v, vKnotSpan);
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the projection of a point on the IGA patch for leakage
     ***********/
    void testIGAPatch2DPointProjection4Leakage() {
        // Modify the Control Point net for getting some curvature into the structure
        igaPatch2D->getControlPointNet()[igaPatch2D->getVNoControlPoints() + 4].setZ(
                102.78054 * 4.500000000000000);
        igaPatch2D->getControlPointNet()[4 * igaPatch2D->getVNoControlPoints()].setZ(
                102.78054 * 4.500000000000000 / 1.5672);

        // Initial guesses for the Newton-Rapson iteration
        double u = 10.0002;
        double v = -6.2365;

        // The vertex to be projected onto the NURBS patch
        double vertex[] = { .5, -1.72, 3.8 };

        // Flag on the convergence of the Newton-Rapson iterations
        bool flag = 1;

        // Compute the orthogonal projection of the point on the NURBS patch iteratively
        int noIterations = 1e9;
        for (int i = 0; i < noIterations; i++) {
            flag = igaPatch2D->computePointProjectionOnPatch(u, v, vertex);
        }
    }

// Make the tests
CPPUNIT_TEST_SUITE(TestIGAPatch2D);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testIGAPatch2DKnotSpan);
    CPPUNIT_TEST(testIGAPatch2DCopyConstructor);
    CPPUNIT_TEST(testIGAPatch2DBSplineBasisFunctions);
    CPPUNIT_TEST(testIGAPatch2DNurbsBasisFunctions);
    CPPUNIT_TEST(testIGAPatch2DBasisFunctionsAndDerivativesTest1);
    CPPUNIT_TEST(testIGAPatch2DBasisFunctionsAndDerivativesTest2);
    CPPUNIT_TEST(testIGAPatch2DPointOnSurface);
    CPPUNIT_TEST(testIGAPatch2DPointOnSurfaceMethod2);
    CPPUNIT_TEST(testIGAPatch2DBaseVectors);
    CPPUNIT_TEST(testIGAPatch2DBaseVectorsAndDerivatives);
    CPPUNIT_TEST(testProjectionOnIGAPatch);

// Make the tests for leakage
    // CPPUNIT_TEST(testIGAPatch2DCopyConstructor4Leakage);
    // CPPUNIT_TEST(testIGAPatch2DBasisFunctions4Leakage);
    // CPPUNIT_TEST(testIGAPatch2DBasisFunctionsAndDerivativesTest14Leakage);
    // CPPUNIT_TEST(testIGAPatch2DBasisFunctionsAndDerivativesTest24Leakage);
    // CPPUNIT_TEST(testIGAPatch2DPointOnSurface4Leakage);
    // CPPUNIT_TEST(testIGAPatch2DPointProjection4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatch2D);
