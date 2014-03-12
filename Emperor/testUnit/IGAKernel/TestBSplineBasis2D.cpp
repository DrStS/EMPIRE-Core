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
#include "BSplineBasis2D.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class BSplineBasis2D
 ***********/

class TestBSplineBasis2D: public CppUnit::TestFixture {

private:
    BSplineBasis2D* bSplineBasis2D;
    double Tol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;

        // Provide an id
        int id = 1;

        // The polynomial degrees
        int p = 3;
        int q = 2;

        // Number of knots
        int noKnotsU = 14;
        int noKnotsV = 8;

        // The knot vectors
        double* knotVectorU = new double[noKnotsU];
        for (int i = 0; i < 4; i++)
            knotVectorU[i] = -1;
        for (int i = 4; i < 5; i++)
            knotVectorU[i] = 2;
        for (int i = 5; i < 7; i++)
            knotVectorU[i] = 3;
        for (int i = 7; i < 10; i++)
            knotVectorU[i] = 5;
        for (int i = 10; i < 14; i++)
            knotVectorU[i] = 10;

        double* knotVectorV = new double[noKnotsV];
        for (int i = 0; i < 3; i++)
            knotVectorV[i] = -1;
        for (int i = 3; i < 5; i++)
            knotVectorV[i] = -.7;
        for (int i = 5; i < 8; i++)
            knotVectorV[i] = 1;

        // Test just one object of the class (that works pretty also)
        bSplineBasis2D = new BSplineBasis2D(id, p, noKnotsU, knotVectorU, q, noKnotsV, knotVectorV);

        // bSplineBasis2D->printPolynomialDegrees();
        // bSplineBasis2D->printNoKnots();
        // bSplineBasis2D->printKnotVectors();
        // bSplineBasis2D->printNoBasisFunctions();

        // Make test on memory leakage (that works pretty fine)
        /* for (int i = 1; i < 1000000000000; i++) {

         // First define the knot vectors
         double* knotVectorU = new double[noKnotsU];
         for (int i = 0; i < 4; i++)
         knotVectorU[i] = -1;
         for (int i = 4; i < 5; i++)
         knotVectorU[i] = 2;
         for (int i = 5; i < 7; i++)
         knotVectorU[i] = 3;
         for (int i = 7; i < 10; i++)
         knotVectorU[i] = 5;
         for (int i = 10; i < 14; i++)
         knotVectorU[i] = 10;

         double* knotVectorV = new double[noKnotsV];
         for (int i = 0; i < 3; i++)
         knotVectorV[i] = -1;
         for (int i = 3; i < 5; i++)
         knotVectorV[i] = -.7;
         for (int i = 5; i < 8; i++)
         knotVectorV[i] = 1;

         // Create an object of the class BSplineBasis2D
         bSplineBasis2D = new BSplineBasis2D(id, p, noKnotsU, knotVectorU, q, noKnotsV,
         knotVectorV);

         // Free the memory
         delete bSplineBasis2D;
         }*/

    }
    void tearDown() {
        delete bSplineBasis2D;
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {

        // Test the ID's
        CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getId()==1);
        CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->getId()==1);
        // Test the polynomial degrees
        CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree()==3);
        CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree()==2);

        // Test the number of knots
        CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getNoKnots()==14);
        CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->getNoKnots()==8);

        // Test the knot vectors
        for (int i = 0; i < 4; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==-1);
        for (int i = 4; i < 5; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==2);
        for (int i = 5; i < 7; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==3);
        for (int i = 7; i < 10; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==5);
        for (int i = 10; i < 14; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==10);
        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->getKnotVector()[i]==-1);
        for (int i = 3; i < 5; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->getKnotVector()[i]==-.7);
        for (int i = 5; i < 8; i++)
            CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->getKnotVector()[i]==1);
    }

    /***********************************************************************************************
     * \brief Test case: Test the knot span
     ***********/
    void testBSpline2DKnotSpans() {

        // The parametric coordinate on the Spline parameter space (from MATLAB this must be 6)
        double u = 28.2 / 4.5;
        int uCorrectknotSpan = 9;
        double v = 2.2 / 3.3;
        int vCorrectknotSpan = 4;

        // Check the find knot span function
        CPPUNIT_ASSERT(bSplineBasis2D->getUBSplineBasis1D()->findKnotSpan(u)==uCorrectknotSpan);
        CPPUNIT_ASSERT(bSplineBasis2D->getVBSplineBasis1D()->findKnotSpan(v)==vCorrectknotSpan);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions
     ***********/
    void testBSpline2DBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = 22.2 / 5.5;
        int uknotSpan = bSplineBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = -1.2 / 8.45;
        int vknotSpan = bSplineBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);
        double* localBasisFunctions = new double[noBasisFunctions];
        bSplineBasis2D->computeLocalBasisFunctions(localBasisFunctions, u, uknotSpan, v, vknotSpan);

        /* cout << endl;
         cout << endl;
         cout << "The non-zero basis functions at ( u , v ) = ( " << u << " , " << v << " )" << endl;
         cout << endl;
         for (int i=0; i<noBasisFunctions; i++) {
         cout << localBasisFunctions[i] << " " ;
         }*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.033651285585847, 0.179685166430089,
                0.175150808439971, 0.062789912459612, 0.032884106018087, 0.175588717039973,
                0.171157733014396, 0.061358432590067, 0.008033604138616, 0.042896414551475,
                0.041813922858180, 0.014989896873687 };
        for (int i = 0; i < noBasisFunctions; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives (Super-efficient algorithm)
     ***********/
    void testBSplineBasisFunctionsAndDerivatives2D() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 30.78 / 3.46;
        int uknotSpan = bSplineBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 0.56 / 4.32;
        int vknotSpan = bSplineBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFcts = (bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);
        int derivDegree = 2;

        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
                * noBasisFcts / 2];
        bSplineBasis2D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                derivDegree, u, uknotSpan, v, vknotSpan);

        // On the derivative order to be displayed
        int vDeriv = 0;
        int uDeriv = 0;

        /* cout << endl;
         cout << endl;

         cout << "The non-zero " << vDeriv
         << "-th w.r.t. v derivatives of the basis functions and their derivatives dR/du at ( u , v ) = ( "
         << u << " , " << v << " ): " << endl;
         cout << endl;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFcts; k++) {
         cout
         << localBasisFunctionsAndDerivatives[(derivDegree - vDeriv)
         * (derivDegree - vDeriv + 1) * noBasisFcts / 2 + i * noBasisFcts + k]
         << " ";
         }
         cout << endl;
         }
         cout << endl;

         cout << endl;
         cout << endl;

         cout << "The non-zero " << uDeriv
         << "-th w.r.t. u derivatives of the basis functions and their derivatives dR/dv at ( u , v ) = ( "
         << u << " , " << v << " ): " << endl;
         cout << endl;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFcts; k++) {
         cout
         << localBasisFunctionsAndDerivatives[(derivDegree - i)
         * (derivDegree - i + 1) / 2 * noBasisFcts + uDeriv * noBasisFcts + k]
         << " ";
         }
         cout << endl;
         }
         cout << endl;

         cout << "The mixed derivative functional ddR/dudv at ( u , v ) = ( " << u << " , " << v
         << " ): " << endl;
         uDeriv = 1;
         vDeriv = 1;
         for (int k = 0; k < noBasisFcts; k++) {
         cout
         << localBasisFunctionsAndDerivatives[(derivDegree - vDeriv)
         * (derivDegree - vDeriv + 1) / 2 * noBasisFcts + uDeriv * noBasisFcts
         + k] << " ";
         }
         cout << endl; */

        // ---------------------------------- The indexing Scheme -------------------------------------------//
        /*cout << endl;
         cout << "Test for the 3D case" << endl;
         cout << "____________________" << endl;
         cout << endl;
         cout << endl;

         int index = 0;
         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= derivDegree - i; j++) {
         for (int k = 0; k < noBasisFunctions; k++) {

         index = (derivDegree - j) * (derivDegree - j + 1) * noBasisFunctions / 2
         + i * noBasisFunctions + k;

         cout << "index triple (i,j,k) = ( " << i << " , " << j << " , " << k
         << " ) and index id = " << index << endl;
         counter++;

         }
         }
         }*/

        // ---------------------------------- The indexing Scheme -------------------------------------------//
        // Values provided by MATLAB
        // For the functional values R, dR/du, d^2R/du^2
        double CorrectlocaldBasisFunctionsdu[] = { 0.00282203329150813, 0.0298751377771175,
                0.105423261056425, 0.124005720684172, 0.00537987623232188, 0.0569534541453134,
                0.200977110439483, 0.236402395176635, 0.00256402611923426, 0.0271437738905324,
                0.0957848356137112, 0.112668375573545, -0.00766824753032839, -0.0464511120030888,
                -0.0413687250883714, 0.0954880846217885, -0.0146186165684558, -0.0885536092654629,
                -0.0788646333599590, 0.182036859193878, -0.00696717044964702, -0.0422042733520504,
                -0.0375865486651720, 0.0867579924668694, 0.0138911709188148, 0.0212367639177693,
                -0.0841470405919829, 0.0490191057553988, 0.0264818917941661, 0.0404854052560027,
                -0.160416485894504, 0.0934491888443348, 0.0126211569401983, 0.0192951718666906,
                -0.0764538145539762, 0.0445374857470872 };

        // For the functional values R, dR/dv, d^2R/dv^2
        double CorrectlocaldBasisFunctionsdv[] = { 0.00282203329150813, 0.0298751377771175,
                0.105423261056425, 0.124005720684172, 0.00537987623232188, 0.0569534541453134,
                0.200977110439483, 0.236402395176635, 0.00256402611923426, 0.0271437738905324,
                0.0957848356137112, 0.112668375573545, -0.00648467224431656, -0.0686492527644403,
                -0.242249195619020, -0.284949315614694, 0.000303537849733967, 0.00321336927833551,
                0.0113393240502520, 0.0133380530713261, 0.00618113439458259, 0.0654358834861048,
                0.230909871568768, 0.271611262543368, 0.00745047449347009, 0.0788736095591442,
                0.278328863051640, 0.327388575387095, -0.0149009489869402, -0.157747219118288,
                -0.556657726103280, -0.654777150774191, 0.00745047449347009, 0.0788736095591442,
                0.278328863051640, 0.327388575387095 };

        // For the functional values d^2R/dudv
        double CorrectlocalddBasisFunctionsdudv[] = { 0.0176206538994780, 0.106738725453906,
                0.0950600491392363, -0.219419428492620, -0.000824796565507483, -0.00499628076592753,
                -0.00444961932141108, 0.0102706966528461, -0.0167958573339705, -0.101742444687979,
                -0.0906104298178253, 0.209148731839774 };

        // Compare the values and assert message if they are not the same up to MATLAB precision

        // For the functional values R, dR/du, d^2R/du^2
        int index = 0;
        int counter = 0;
        vDeriv = 0;
        for (int i = 0; i <= derivDegree; i++) {
            for (int k = 0; k < noBasisFcts; k++) {
                // Find the correct index
                index = bSplineBasis2D->indexDerivativeBasisFunction(derivDegree, i, vDeriv, k);

                // Check the result
                CPPUNIT_ASSERT(
                        fabs( localBasisFunctionsAndDerivatives[index] - CorrectlocaldBasisFunctionsdu[counter]) <= Tol);

                // Update counter for the pre-computed data
                counter++;
            }
        }

        // For the functional values R, dR/dv, d^2R/dv^2
        uDeriv = 0;
        counter = 0;
        for (int i = 0; i <= derivDegree; i++) {
            for (int k = 0; k < noBasisFcts; k++) {
                // Find the correct index
                index = bSplineBasis2D->indexDerivativeBasisFunction(derivDegree, uDeriv, i, k);

                // Check the result
                CPPUNIT_ASSERT(
                        fabs( localBasisFunctionsAndDerivatives[index] - CorrectlocaldBasisFunctionsdv[counter]) <= Tol);

                // Update counter for the pre-computed data
                counter++;
            }
        }

        // For the functional values d^2R/dudv
        uDeriv = 1;
        vDeriv = 1;
        counter = 0;
        for (int k = 0; k < noBasisFcts; k++) {
            // Find the correct index
            index = bSplineBasis2D->indexDerivativeBasisFunction(derivDegree, uDeriv, vDeriv, k);

            // Check the result
            CPPUNIT_ASSERT(
                    fabs( localBasisFunctionsAndDerivatives[index] - CorrectlocalddBasisFunctionsdudv[counter]) <= Tol);

            // Update counter for the pre-computed data
            counter++;
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constuctor of the class BSplineBasis2D
     ***********/
    void testBSplineBasis2DCopyConstructor() {

        BSplineBasis2D* bSplineBasis2DCopy = new BSplineBasis2D(*bSplineBasis2D);

        // Test the ID's
        CPPUNIT_ASSERT(bSplineBasis2DCopy->getId()==bSplineBasis2D->getId());
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getUBSplineBasis1D()->getId()==bSplineBasis2D->getUBSplineBasis1D()->getId());
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getVBSplineBasis1D()->getId()==bSplineBasis2D->getVBSplineBasis1D()->getId());

        // Test the polynomial degrees
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getUBSplineBasis1D()->getPolynomialDegree()==bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree());
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getVBSplineBasis1D()->getPolynomialDegree()==bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree());

        // Test the knot vectors
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getUBSplineBasis1D()->getNoKnots()==bSplineBasis2D->getUBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < bSplineBasis2DCopy->getUBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    bSplineBasis2DCopy->getUBSplineBasis1D()->getKnotVector()[i]==bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]);
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getVBSplineBasis1D()->getNoKnots()==bSplineBasis2D->getVBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < bSplineBasis2DCopy->getVBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    bSplineBasis2DCopy->getVBSplineBasis1D()->getKnotVector()[i]==bSplineBasis2D->getVBSplineBasis1D()->getKnotVector()[i]);

        // Free the memory on the heap from the pointers
        delete bSplineBasis2DCopy;

        /* // Test the copy constructor for memory leakage
         int noIterations = 1e9;

         for (int i=0; i<noIterations; i++) {
         BSplineBasis2D* bSplineBasis2DCopy2 = new BSplineBasis2D(*bSplineBasis2D);
         delete bSplineBasis2DCopy2;
         } */
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy assignement of the class BSplineBasis2D
     ***********/
    void testBSplineBasis2DCopyAssignment() {

        BSplineBasis2D* bSplineBasis2DCopy = bSplineBasis2D;

        // Test the ID's
        CPPUNIT_ASSERT(bSplineBasis2DCopy->getId()==bSplineBasis2D->getId());
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getUBSplineBasis1D()->getId()==bSplineBasis2D->getUBSplineBasis1D()->getId());
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getVBSplineBasis1D()->getId()==bSplineBasis2D->getVBSplineBasis1D()->getId());

        // Test the polynomial degrees
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getUBSplineBasis1D()->getPolynomialDegree()==bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree());
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getVBSplineBasis1D()->getPolynomialDegree()==bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree());

        // Test the knot vectors
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getUBSplineBasis1D()->getNoKnots()==bSplineBasis2D->getUBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < bSplineBasis2DCopy->getUBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    bSplineBasis2DCopy->getUBSplineBasis1D()->getKnotVector()[i]==bSplineBasis2D->getUBSplineBasis1D()->getKnotVector()[i]);
        CPPUNIT_ASSERT(
                bSplineBasis2DCopy->getVBSplineBasis1D()->getNoKnots()==bSplineBasis2D->getVBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < bSplineBasis2DCopy->getVBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    bSplineBasis2DCopy->getVBSplineBasis1D()->getKnotVector()[i]==bSplineBasis2D->getVBSplineBasis1D()->getKnotVector()[i]);

        /* // Test the copy assignement for memory leakage
         int noIterations = 1e9;

         for (int i=0; i<noIterations; i++) {
         BSplineBasis2D* bSplineBasis2DCopy2 = bSplineBasis2D;
         }*/
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the 2D B-Spline basis functions for memory leakage
     ***********/
    void testBSpline2DBasisFunctions4Leakage() {
        // initialize variables
        int noIterations = 1e8;
        double u = 22.2 / 5.5;
        int uknotSpan = bSplineBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = -1.2 / 8.45;
        int vknotSpan = bSplineBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions at another parametric location
            double* localBasisFunctions = new double[noBasisFunctions];
            bSplineBasis2D->computeLocalBasisFunctions(localBasisFunctions, u, uknotSpan, v,
                    vknotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the 2D B-Spline basis functions and their derivatives for memory leakage
     ***********/
    void testBSpline2DBasisFunctionsAndDerivatives4Leakage() {
        // initialize variables
        int noIterations = 1e8;
        double u = 30.78 / 3.46;
        int uknotSpan = bSplineBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 0.56 / 4.32;
        int vknotSpan = bSplineBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (bSplineBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (bSplineBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);
        int derivDegree = 8;

        for (int i = 0; i < noIterations; i++) {
            // Initialize the array holding the B-Spline basis functions and their derivatives
            double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                    * (derivDegree + 1) * noBasisFunctions];

            // Fill up the array with appropriate values
            bSplineBasis2D->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, u, uknotSpan, v, vknotSpan);

            // Clean up the array from the heap
            delete[] localBasisFunctionsAndDerivatives;
        }
    }

    // Make the tests
CPPUNIT_TEST_SUITE(TestBSplineBasis2D);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testBSpline2DKnotSpans);
    CPPUNIT_TEST(testBSpline2DBasisFunctions);
    CPPUNIT_TEST(testBSplineBasisFunctionsAndDerivatives2D);
    CPPUNIT_TEST(testBSplineBasis2DCopyConstructor);
    CPPUNIT_TEST(testBSplineBasis2DCopyAssignment);

    // Make the tests for leakage
    // CPPUNIT_TEST(testBSpline2DBasisFunctions4Leakage);
    // CPPUNIT_TEST(testBSpline2DBasisFunctionsAndDerivatives4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestBSplineBasis2D);
