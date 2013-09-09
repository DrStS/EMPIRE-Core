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
#include "BSplineBasis1D.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class BSplineBasis1D
 ***********/

class TestBSplineBasis1D: public CppUnit::TestFixture {

private:
    BSplineBasis1D* bSplineBasis1D;
    double Tol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;

        // Provide an id
        int id = 1;

        // The polynomial degree
        int p = 3;

        // Number of knots
        int noKnots = 14;

        // The knot vector
        double* knotVector = new double[noKnots];

        for (int i = 0; i < 4; i++) {
            knotVector[i] = -1;
        }
        for (int i = 4; i < 5; i++) {
            knotVector[i] = 2;
        }
        for (int i = 5; i < 7; i++) {
            knotVector[i] = 3;
        }
        for (int i = 7; i < 10; i++) {
            knotVector[i] = 5;
        }
        for (int i = 10; i < 14; i++) {
            knotVector[i] = 10;
        }

        // Test just one object of the class (that works pretty also)
        bSplineBasis1D = new BSplineBasis1D(id, p, noKnots, knotVector);
        // bSplineBasis1D->printPolynomialDegree();
        // bSplineBasis1D->printNoKnots();
        // bSplineBasis1D->printKnotVector();
        // bSplineBasis1D->printNoBasisFunctions();

        // Make test on memory leakage (that works pretty fine)
        /*for (int i = 1; i < 1e9; i++) {
         // Create the knot vector
         double* knotVector = new double[noKnots];

         for (int i = 0; i < 4; i++) {
         knotVector[i] = -1;
         }
         for (int i = 4; i < 5; i++) {
         knotVector[i] = 2;
         }
         for (int i = 5; i < 7; i++) {
         knotVector[i] = 3;
         }
         for (int i = 7; i < 10; i++) {
         knotVector[i] = 5;
         }
         for (int i = 10; i < 14; i++) {
         knotVector[i] = 10;
         }

         bSplineBasis1D = new BSplineBasis1D(id, p, noKnots, knotVector);

         // free the memory
         delete bSplineBasis1D;
         }*/

    }
    void tearDown() {
        delete bSplineBasis1D;
    }
    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        CPPUNIT_ASSERT(bSplineBasis1D->getId()==1);
        CPPUNIT_ASSERT(bSplineBasis1D->getPolynomialDegree()==3);
        CPPUNIT_ASSERT(bSplineBasis1D->getNoKnots()==14);
        for (int i = 0; i < 4; i++) {
            CPPUNIT_ASSERT(bSplineBasis1D->getKnotVector()[i]==-1);
        }
        for (int i = 4; i < 5; i++) {
            CPPUNIT_ASSERT(bSplineBasis1D->getKnotVector()[i]==2);
        }
        for (int i = 5; i < 7; i++) {
            CPPUNIT_ASSERT(bSplineBasis1D->getKnotVector()[i]==3);
        }
        for (int i = 7; i < 10; i++) {
            CPPUNIT_ASSERT(bSplineBasis1D->getKnotVector()[i]==5);
        }
        for (int i = 10; i < 14; i++) {
            CPPUNIT_ASSERT(bSplineBasis1D->getKnotVector()[i]==10);
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the knot span
     ***********/
    void testBSpline1DKnotSpan() {

        // The parametric coordinate on the Spline parameter space (from MATLAB this must be 6)
        double u = 7.0 / 2.0;
        int CorrectknotSpan = 6;

        // Check the find knot span function
        CPPUNIT_ASSERT(bSplineBasis1D->findKnotSpan(u)==CorrectknotSpan);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions
     ***********/
    void testBSpline1DBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = 17.0 / 3.0;
        int knotSpan = bSplineBasis1D->findKnotSpan(u);
        double* localBasisFunctions = new double[bSplineBasis1D->getPolynomialDegree() + 1];
        bSplineBasis1D->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero basis functions at u = " << u << " :";
         cout << endl;
         for (int i=0; i<bSplineBasis1D->getPolynomialDegree()+1; i++) {
         cout << localBasisFunctions[i] << " " ;
         }*/

        // Value provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.650962962962963, 0.300444444444445,
                0.046222222222222, 0.002370370370370 };
        for (int i = 0; i < bSplineBasis1D->getPolynomialDegree() + 1; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives
     ***********/
    void testBSpline1DBasisFunctionsAndDerivativesInefficient() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 39.0 / 5.0;
        int knotSpan = bSplineBasis1D->findKnotSpan(u);
        int derivDegree = 2;
        double** localBasisFunctionsAndDerivatives = new double*[derivDegree + 1];
        for (int i = 0; i <= derivDegree; i++)
            localBasisFunctionsAndDerivatives[i] = new double[bSplineBasis1D->getPolynomialDegree()
                    + 1];

        bSplineBasis1D->computeLocalBasisFunctionsAndDerivativesInefficient(
                localBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);
        /*cout << endl;
         cout << endl;
         cout << "The non-zero basis functions and their derivatives at u = " << u << " :";
         cout << endl;
         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= bSplineBasis1D->getPolynomialDegree(); j++)
         cout << localBasisFunctionsAndDerivatives[i][j] << " ";
         cout << endl;
         }*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[][4] = { { 0.085184000000000,
                0.325248000000000, 0.413952000000000, 0.175616000000000 }, { -0.116160000000000,
                -0.179520000000000, 0.107520000000000, 0.188160000000000 }, { 0.105600000000000,
                -0.076800000000000, -0.163200000000000, 0.134400000000000 } };

        for (int i = 0; i <= derivDegree; i++)
            for (int j = 0; j <= bSplineBasis1D->getPolynomialDegree(); j++)
                CPPUNIT_ASSERT(
                        fabs(localBasisFunctionsAndDerivatives[i][j]-CorrectlocalBasisFunctionsAndDerivatives[i][j])<=Tol);

        // Clear the heap from the pointer
        for (int i = 0; i <= derivDegree; i++)
            delete[] localBasisFunctionsAndDerivatives[i];
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives (Efficient algorithm)
     ***********/
    void testBSpline1DBasisFunctionsAndDerivatives() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 39.0 / 5.0;
        int knotSpan = bSplineBasis1D->findKnotSpan(u);
        int derivDegree = 2;
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                * (bSplineBasis1D->getPolynomialDegree() + 1)];

        bSplineBasis1D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                derivDegree, u, knotSpan);
        /*cout << endl;
         cout << endl;
         cout << "The non-zero basis functions and their derivatives at u = " << u << " :";
         cout << endl;

         // Initialize counter
         int counter = 0;

         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= bSplineBasis1D->getPolynomialDegree(); j++) {
         cout << localBasisFunctionsAndDerivatives[counter] << " ";
         counter += 1;
         }
         cout << endl;
         }

         cout << endl;
         cout << endl;
         cout << "The non-zero basis functions and their derivatives at u = " << u << " :";
         cout << endl;
         for (int i = 0; i < (derivDegree + 1) * (bSplineBasis1D->getPolynomialDegree() + 1); i++) {
         cout << localBasisFunctionsAndDerivatives[i] << " ";
         counter += 1;
         }
         cout << endl;*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 0.085184000000000, 0.325248000000000,
                0.413952000000000, 0.175616000000000, -0.116160000000000, -0.179520000000000,
                0.107520000000000, 0.188160000000000, 0.105600000000000, -0.076800000000000,
                -0.163200000000000, 0.134400000000000 };

        for (int i = 0; i < (derivDegree + 1) * (bSplineBasis1D->getPolynomialDegree() + 1); i++)
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[i]-CorrectlocalBasisFunctionsAndDerivatives[i])<=Tol);

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the B-Spline basis functions for memory leakage
     ***********/
    void testBSpline1DBasisFunctions4Leakage() {
        // initialize variables
        int noIterations = 1000000000;
        int knotSpan = 0;
        double u = 17.0 / 3.0;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions at another parametric location
            knotSpan = bSplineBasis1D->findKnotSpan(u);
            double* localBasisFunctions = new double[bSplineBasis1D->getPolynomialDegree() + 1];
            bSplineBasis1D->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the B-Spline basis functions and their derivatives for memory leakage (Inefficient algorithm)
     ***********/
    void testBSpline1DBasisFunctionsAndDerivativesInefficient4Leakage() {
        // initialize variables
        int noIterations = 1000000000;
        int knotSpan = 0;
        double u = 17.0 / 3.0;
        int derivDegree = 2;

        for (int i = 0; i < noIterations; i++) {
            // Initialize the array holding the B-Spline basis functions and their derivatives
            double** localBasisFunctionsAndDerivatives = new double*[derivDegree + 1];
            for (int i = 0; i <= derivDegree; i++)
                localBasisFunctionsAndDerivatives[i] =
                        new double[bSplineBasis1D->getPolynomialDegree() + 1];

            // Fill up the array with appropriate values
            bSplineBasis1D->computeLocalBasisFunctionsAndDerivativesInefficient(
                    localBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);

            // Clean up the array from the heap
            for (int i = 0; i <= derivDegree; i++)
                delete[] localBasisFunctionsAndDerivatives[i];
            delete[] localBasisFunctionsAndDerivatives;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the B-Spline basis functions and their derivatives for memory leakage (Inefficient algorithm)
     ***********/
    void testBSpline1DBasisFunctionsAndDerivatives4Leakage() {
        // initialize variables
        int noIterations = 1000000000;
        int knotSpan = 0;
        double u = 17.0 / 3.0;
        int derivDegree = 2;

        for (int i = 0; i < noIterations; i++) {
            // Initialize the array holding the B-Spline basis functions and their derivatives
            double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                    * (bSplineBasis1D->getPolynomialDegree() + 1)];

            // Fill up the array with appropriate values
            bSplineBasis1D->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);

            // Clear the memory
            delete[] localBasisFunctionsAndDerivatives;
        }
    }

    // Make the tests
CPPUNIT_TEST_SUITE(TestBSplineBasis1D);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testBSpline1DKnotSpan);
    CPPUNIT_TEST(testBSpline1DBasisFunctions);
    CPPUNIT_TEST(testBSpline1DBasisFunctionsAndDerivativesInefficient);
    CPPUNIT_TEST(testBSpline1DBasisFunctionsAndDerivatives);

    // Make the tests for leakage
    // CPPUNIT_TEST(testBSpline1DBasisFunctions4Leakage);
    // CPPUNIT_TEST(testBSpline1DBasisFunctionsAndDerivativesInefficient4Leakage);
    // CPPUNIT_TEST(testBSpline1DBasisFunctionsAndDerivatives4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestBSplineBasis1D);
