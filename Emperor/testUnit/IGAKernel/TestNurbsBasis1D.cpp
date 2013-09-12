// inclusion of standard libraries   (only if really necessary here in *.h)
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"
#include <iostream>
#include <string>
#include <math.h>

// Inclusion of user-defined libraries
#include "NurbsBasis1D.h"
#include "IGAControlPoint.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class NurbsBasis1D
 ***********/

class TestNurbsBasis1D: public CppUnit::TestFixture {

private:
    NurbsBasis1D* nurbsBasis1D;
    double Tol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;

        // Provide an id
        int id = 1;

        // The polynomial degree
        int p = 5;

        // Number of knots
        int noKnots = 16;

        // The knot vector
        double* knotVector = new double[noKnots];

        for (int i = 0; i < 6; i++) {
            knotVector[i] = -4;
        }
        for (int i = 6; i < 7; i++) {
            knotVector[i] = 0;
        }
        for (int i = 7; i < 10; i++) {
            knotVector[i] = 1;
        }
        for (int i = 10; i < 16; i++) {
            knotVector[i] = 21;
        }

        // The Control Point net
        int noControlPoints = noKnots - p - 1;
        IGAControlPoint** controlPointNet = new IGAControlPoint*[noControlPoints];
        double* controlPointWeights = new double[noControlPoints];
        for (int i = 0; i < noControlPoints; i++) {
            controlPointNet[i] = new IGAControlPoint(i, 3 * i - i, 1 - 2 * i, i * i, 1 + i);
            controlPointWeights[i] = controlPointNet[i]->getW();
        }

        // Test just one object of the class (that works pretty also)
        nurbsBasis1D = new NurbsBasis1D(id, p, noKnots, knotVector, noControlPoints,
                controlPointWeights);
        // nurbsBasis1D->printPolynomialDegree();
        // nurbsBasis1D->printNoKnots();
        // nurbsBasis1D->printKnotVector();
        // nurbsBasis1D->printNoBasisFunctions();
        // nurbsBasis1D->printControlPointNet();
        // NurbsBasis1D nurbsBasis1DCopy = *nurbsBasis1D;
        // delete nurbsBasis1D;

        // Make test on memory leakage (that works pretty fine)
        /* for (int j = 1; j < 1e9; j++) {

         // The knot vector
         double* knotVectorCopy = new double[noKnots];

         for (int i = 0; i < 6; i++) {
         knotVectorCopy[i] = -4;
         }
         for (int i = 6; i < 7; i++) {
         knotVectorCopy[i] = 0;
         }
         for (int i = 7; i < 10; i++) {
         knotVectorCopy[i] = 1;
         }
         for (int i = 10; i < 16; i++) {
         knotVectorCopy[i] = 21;
         }

         // Fill up the array of the Control Point weights
         IGAControlPoint* controlPointNetCopy = new IGAControlPoint[noControlPoints];
         double* controlPointWeightsCopy = new double[noControlPoints];

         for (int i = 0; i < noControlPoints; i++) {
         controlPointNetCopy[i] = IGAControlPoint(i, 3 * i - i, 1 - 2 * i, i * i, 1 + i);
         controlPointWeightsCopy[i] = controlPointNetCopy[i].getW();
         }

         // Create an object of the class NurbsBasis1D
         NurbsBasis1D* nurbsBasis1Dcopy = new NurbsBasis1D(id, p, noKnots, knotVectorCopy,
         noControlPoints, controlPointWeightsCopy);

         // Free the memory on the heap from the pointer
         delete nurbsBasis1Dcopy;
         delete[] controlPointNetCopy;
         }*/
    }

    void tearDown() {
        delete nurbsBasis1D;
    }

    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        CPPUNIT_ASSERT(nurbsBasis1D->getId()==1);
        CPPUNIT_ASSERT(nurbsBasis1D->getPolynomialDegree()==5);
        CPPUNIT_ASSERT(nurbsBasis1D->getNoKnots()==16);
        for (int i = 0; i < 6; i++) {
            CPPUNIT_ASSERT(nurbsBasis1D->getKnotVector()[i]==-4);
        }
        for (int i = 6; i < 7; i++) {
            CPPUNIT_ASSERT(nurbsBasis1D->getKnotVector()[i]==0);
        }
        for (int i = 7; i < 10; i++) {
            CPPUNIT_ASSERT(nurbsBasis1D->getKnotVector()[i]==1);
        }
        for (int i = 10; i < 16; i++) {
            CPPUNIT_ASSERT(nurbsBasis1D->getKnotVector()[i]==21);
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the knot span
     ***********/
    void testNurbs1DKnotSpan() {

        // The parametric coordinate on the Spline parameter space (from MATLAB this must be 6)
        double u = 7.0 / 2.0;
        int CorrectknotSpan = 9;

        // Check the find knot span function
        CPPUNIT_ASSERT(nurbsBasis1D->findKnotSpan(u)==CorrectknotSpan);
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constructor for proper functionality
     ***********/
    void testNurbs1DCopyConstructor() {

        NurbsBasis1D nurbsBasis1DCopy(*nurbsBasis1D);

        CPPUNIT_ASSERT(nurbsBasis1DCopy.getId()==nurbsBasis1D->getId());
        CPPUNIT_ASSERT(nurbsBasis1DCopy.getPolynomialDegree()==nurbsBasis1D->getPolynomialDegree());
        CPPUNIT_ASSERT(nurbsBasis1DCopy.getNoKnots()==nurbsBasis1D->getNoKnots());
        for (int i = 1; i < nurbsBasis1DCopy.getNoKnots(); i++)
            CPPUNIT_ASSERT(nurbsBasis1DCopy.getKnotVector()[i]==nurbsBasis1D->getKnotVector()[i]);
        CPPUNIT_ASSERT(nurbsBasis1DCopy.getNoControlPoints()==nurbsBasis1D->getNoControlPoints());
        for (int i = 1; i < nurbsBasis1DCopy.getNoControlPoints(); i++)
            CPPUNIT_ASSERT(
                    nurbsBasis1DCopy.getControlPointWeights()[i]==nurbsBasis1D->getControlPointWeights()[i]);
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constructor for leakage
     ***********/
    void testNurbs1DCopyConstructor4Leakage() {

        for (int i = 1; i < 1000000000; i++) {
            NurbsBasis1D nurbsBasis1DCopy(*nurbsBasis1D);
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local NURBS basis functions
     ***********/
    void testNurbs1DBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = 23.0 / 3.0;
        int knotSpan = nurbsBasis1D->findKnotSpan(u);
        int noLocalBasisFunctions = nurbsBasis1D->getPolynomialDegree() + 1;

        double* localBasisFunctions = new double[noLocalBasisFunctions];
        nurbsBasis1D->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero basis functions at u = " << u << " :";
         cout << endl;
         for (int i=0; i<noLocalBasisFunctions; i++) {
         cout << localBasisFunctions[i] << " " ;
         }
         cout << endl;
         cout << endl;*/

        // Value provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.074719385344104, 0.307950609596771,
                0.359898372740767, 0.196138386528273, 0.055163921211077, 0.006129324579009 };
        for (int i = 0; i < nurbsBasis1D->getPolynomialDegree() + 1; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives
     ***********/
    void testNurbs1DBasisFunctionsAndDerivativesInefficient() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 51.0 / 4.5;
        int knotSpan = nurbsBasis1D->findKnotSpan(u);
        int derivDegree = 2;
        double** localBasisFunctionsAndDerivatives = new double*[derivDegree + 1];
        for (int i = 0; i <= derivDegree; i++)
            localBasisFunctionsAndDerivatives[i] = new double[nurbsBasis1D->getPolynomialDegree()
                    + 1];

        nurbsBasis1D->computeLocalBasisFunctionsAndDerivativesInefficient(
                localBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);

        /* cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives at u = " << u << " :";
         cout << endl;
         for (int i=0; i<=derivDegree; i++) {
         for (int j=0; j<=nurbsBasis1D->getPolynomialDegree(); j++)
         cout << localBasisFunctionsAndDerivatives[i][j] << " " ;
         cout << endl;
         }*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[][6] = { { 0.013228096315987,
                0.110966134071401, 0.283990266789861, 0.339318650389149, 0.204029964350372,
                0.048466888083230 }, { -0.007268043699170, -0.039735058036026, -0.043573844884038,
                0.017382454544698, 0.051303333334689, 0.021891158739847 }, { 0.003296333811469,
                0.008711340524660, -0.007276560245087, -0.015628735809817, 0.003239989772855,
                0.007657631945920 } };

        for (int i = 0; i <= derivDegree; i++)
            for (int j = 0; j <= nurbsBasis1D->getPolynomialDegree(); j++)
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
    void testNurbs1DBasisFunctionsAndDerivatives() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 51.0 / 4.5;
        int knotSpan = nurbsBasis1D->findKnotSpan(u);
        int derivDegree = 2;
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                * (nurbsBasis1D->getPolynomialDegree() + 1)];

        nurbsBasis1D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                derivDegree, u, knotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives at u = " << u << " :";
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= nurbsBasis1D->getPolynomialDegree(); j++) {
         cout << localBasisFunctionsAndDerivatives[counter] << " ";
         counter += 1;
         }
         cout << endl;
         }*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 0.013228096315987, 0.110966134071401,
                0.283990266789861, 0.339318650389149, 0.204029964350372, 0.048466888083230,
                -0.007268043699170, -0.039735058036026, -0.043573844884038, 0.017382454544698,
                0.051303333334689, 0.021891158739847, 0.003296333811469, 0.008711340524660,
                -0.007276560245087, -0.015628735809817, 0.003239989772855, 0.007657631945920 };

        for (int i = 0; i < (derivDegree + 1) * (nurbsBasis1D->getPolynomialDegree() + 1); i++)
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[i]-CorrectlocalBasisFunctionsAndDerivatives[i])<=Tol);

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives (super-efficient algorithm)
     ***********/
    void testNurbs1DBasisFunctionsAndDerivativesMethod2() {
        // The NURBS parameter
        double u = -4 / 1.410098;

        // Find the correct knot span
        int knotSpan = nurbsBasis1D->findKnotSpan(u);

        // The derivative order
        int derivDegree = 8;

        // Compute the non-zero B-Spline basis functions and their derivatives at another parametric location
        double* localBSplineBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                * (nurbsBasis1D->getPolynomialDegree() + 1)];
        nurbsBasis1D->BSplineBasis1D::computeLocalBasisFunctionsAndDerivatives(
                localBSplineBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);

        // Compute the denominator function and its derivatives
        double* denominatorFctAndDerivatives = new double[derivDegree + 1];
        nurbsBasis1D->computeDenominatorFunctionAndDerivatives(denominatorFctAndDerivatives,
                localBSplineBasisFunctionsAndDerivatives, derivDegree, knotSpan);

        // Compute the non-zero NURBS basis functions and their derivatives at another parametric location
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                * (nurbsBasis1D->getPolynomialDegree() + 1)];
        nurbsBasis1D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                localBSplineBasisFunctionsAndDerivatives, denominatorFctAndDerivatives, derivDegree,
                u, knotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives at u = " << u << " :";
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= nurbsBasis1D->getPolynomialDegree(); j++) {
         cout << localBasisFunctionsAndDerivatives[counter] << " ";
         counter += 1;
         }
         cout << endl;
         }*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 0.077090634624221, 0.372438555973010,
                0.365414676500290, 0.153410040390891, 0.031206570053090, 0.000439522458499,
                -0.171011798859359, -0.300920645476598, 0.144989207625777, 0.239178643404263,
                0.086075792823857, 0.001688800482060, 0.350754290149719, -0.072477526516815,
                -0.481579118557142, 0.046907015651030, 0.151420234970259, 0.004975104302949,
                -0.706298121006933, 0.861346531588323, 0.157121067475015, -0.414792072834333,
                0.092469989762257, 0.010152605015672, 1.543834785719358, -2.330978351186587,
                0.852654740695628, 0.079755404332027, -0.156654737487208, 0.011388157926784,
                -4.012566946121711, 5.968908813513032, -2.491292772116841, 0.619194104575780,
                -0.086326551867316, 0.002083352017056, 12.618331912156258, -18.380785733445034,
                6.925417334828507, -1.355560243364866, 0.193338553002237, -0.000741823177098,
                -46.565171137529695, 67.804747080717675, -25.001262224398399, 4.109516334234632,
                -0.358007335321444, 0.010177282297207, 196.3327196083553, -286.4383376009845,
                106.4190989610515, -17.5012152855986, 1.1874085785556, 0.0003257386208 };

        int counter = 0;
        for (int i = 0; i <= derivDegree; i++) {
            for (int j = 0; j < nurbsBasis1D->getPolynomialDegree(); j++) {
                if (i < derivDegree-1) {
                    CPPUNIT_ASSERT(
                            fabs(localBasisFunctionsAndDerivatives[counter]-CorrectlocalBasisFunctionsAndDerivatives[counter])<=Tol);
                } else {
                    CPPUNIT_ASSERT(
                            fabs(localBasisFunctionsAndDerivatives[counter]-CorrectlocalBasisFunctionsAndDerivatives[counter])<=Tol*1e2);
                }
                counter++;
            }
        }

        // Clear the heap from the pointer
        delete[] localBSplineBasisFunctionsAndDerivatives;
        delete[] denominatorFctAndDerivatives;
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions for memory leakage
     ***********/
    void testNurbs1DBasisFunctions4Leakage() {

        // Compute the non-zero basis functions at another parametric location
        double u = 23.0 / 3.0;
        int knotSpan = nurbsBasis1D->findKnotSpan(u);
        int noLocalBasisFunctions = nurbsBasis1D->getPolynomialDegree() + 1;

        for (int i = 1; i < 1000000000; i++) {
            double* localBasisFunctions = new double[noLocalBasisFunctions];
            nurbsBasis1D->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

            delete[] localBasisFunctions;
        }

    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (inefficient algorithm)
     ***********/
    void testNurbs1DBasisFunctionsAndDerivativesInefficient4Leakage() {
        // initialize variables
        int noIterations = 1000000000;
        int knotSpan = 0;
        int derivDegree = 2;
        double u = 17.0 / 3.0;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions at another parametric location
            knotSpan = nurbsBasis1D->findKnotSpan(u);
            double** localBasisFunctions = new double*[derivDegree + 1];
            for (int i = 0; i <= derivDegree; i++)
                localBasisFunctions[i] = new double[nurbsBasis1D->getPolynomialDegree() + 1];

            nurbsBasis1D->computeLocalBasisFunctionsAndDerivativesInefficient(localBasisFunctions,
                    derivDegree, u, knotSpan);

            // Clear the heap from the pointer
            for (int i = 0; i <= derivDegree; i++)
                delete[] localBasisFunctions[i];
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (efficient algorithm)
     ***********/
    void testNurbs1DBasisFunctionsAndDerivatives4Leakage() {
        // initialize variables
        int noIterations = 1000000000;
        int knotSpan = 0;
        int derivDegree = 2;
        double u = 17.0 / 3.0;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions at another parametric location
            knotSpan = nurbsBasis1D->findKnotSpan(u);
            double* localBasisFunctions = new double[(derivDegree + 1)
                    * (nurbsBasis1D->getPolynomialDegree() + 1)];

            nurbsBasis1D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctions, derivDegree,
                    u, knotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctions;
        }
    }

    // Make the tests
CPPUNIT_TEST_SUITE(TestNurbsBasis1D);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testNurbs1DKnotSpan);
    CPPUNIT_TEST(testNurbs1DCopyConstructor);
    CPPUNIT_TEST(testNurbs1DBasisFunctions);
    CPPUNIT_TEST(testNurbs1DBasisFunctionsAndDerivativesInefficient);
    CPPUNIT_TEST(testNurbs1DBasisFunctionsAndDerivatives);
    CPPUNIT_TEST(testNurbs1DBasisFunctionsAndDerivativesMethod2);

    // Make the tests for leakage
    // CPPUNIT_TEST(testNurbs1DCopyConstructor4Leakage);
    // CPPUNIT_TEST(testNurbs1DBasisFunctions4Leakage);
    // CPPUNIT_TEST(testNurbs1DBasisFunctionsAndDerivativesInefficient4Leakage);
    // CPPUNIT_TEST(testNurbs1DBasisFunctionsAndDerivatives4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestNurbsBasis1D);
