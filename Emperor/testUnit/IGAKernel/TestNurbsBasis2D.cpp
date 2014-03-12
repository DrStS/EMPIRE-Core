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
#include "NurbsBasis2D.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class BSplineBasis2D
 ***********/

class TestNurbsBasis2D: public CppUnit::TestFixture {

private:
    NurbsBasis2D* nurbsBasis2D;
    double Tol;
    double rlxTol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;
        rlxTol = 2 * 1e-15;

        // Provide an id
        int id = 1;

        // The polynomial degrees
        int p = 5;
        int q = 2;

        // Number of knots
        int noKnotsU = 16;
        int noKnotsV = 6;

        // The knot vectors
        double* knotVectorU = new double[noKnotsU];
        for (int i = 0; i < 6; i++)
            knotVectorU[i] = -4;
        for (int i = 6; i < 7; i++)
            knotVectorU[i] = 0;
        for (int i = 7; i < 10; i++)
            knotVectorU[i] = 1;
        for (int i = 10; i < 16; i++)
            knotVectorU[i] = 21;

        double* knotVectorV = new double[noKnotsV];
        for (int i = 0; i < 3; i++)
            knotVectorV[i] = 0;
        for (int i = 3; i < 6; i++)
            knotVectorV[i] = 1;

        // The Control Points
        double mid_Radius_NURBS = 4.5;
        double width_NURBS = 2.0;
        double weight_NURBS = 1;

        int noUControlPoints = noKnotsU - p - 1;
        int noVControlPoints = noKnotsV - q - 1;

        double* controlPointWeights = new double[noUControlPoints * noVControlPoints];

        controlPointWeights[0] = 1;
        controlPointWeights[1] = 2;
        controlPointWeights[2] = 3.0;
        controlPointWeights[3] = 4.0;
        controlPointWeights[4] = 9.0;
        controlPointWeights[5] = 10.0;
        controlPointWeights[6] = 0.5;
        controlPointWeights[7] = 2.0;
        controlPointWeights[8] = 4.0;
        controlPointWeights[9] = 21.0;
        controlPointWeights[10] = 3;
        controlPointWeights[11] = 1.5;
        controlPointWeights[12] = 3.7;
        controlPointWeights[13] = 1.9;
        controlPointWeights[14] = 4.4;
        controlPointWeights[15] = 1.6;
        controlPointWeights[16] = 0.34;
        controlPointWeights[17] = 2.6;
        controlPointWeights[18] = 1.0;
        controlPointWeights[19] = 1.4;
        controlPointWeights[20] = 1;
        controlPointWeights[21] = 4;
        controlPointWeights[22] = 2.1;
        controlPointWeights[23] = 1.0;
        controlPointWeights[24] = 3.2;
        controlPointWeights[25] = 1.02;
        controlPointWeights[26] = 3.63;
        controlPointWeights[27] = 2.0;
        controlPointWeights[28] = 8.0;
        controlPointWeights[29] = 10.54;

        // Test just one object of the class (that works pretty also)
        nurbsBasis2D = new NurbsBasis2D(id, p, noKnotsU, knotVectorU, q, noKnotsV, knotVectorV,
                noUControlPoints, noVControlPoints, controlPointWeights);

        // nurbsBasis2D->printPolynomialDegrees();
        // nurbsBasis2D->printNoKnots();
        // nurbsBasis2D->printKnotVectors();
        // nurbsBasis2D->printNoBasisFunctions();
        // nurbsBasis2D->printControlPointNet();
        // nurbsBasis2D->getUBSplineBasis1D()->printNoBasisFunctions();
        // nurbsBasis2D->getVBSplineBasis1D()->printNoBasisFunctions();

        // Make test on memory leakage (that works pretty fine)
        /* for (int i = 1; i < 1e9; i++) {

         // On the knot vectors
         double* knotVectorUCopy = new double[noKnotsU];
         for (int i = 0; i < 6; i++)
         knotVectorUCopy[i] = -4;
         for (int i = 6; i < 7; i++)
         knotVectorUCopy[i] = 0;
         for (int i = 7; i < 10; i++)
         knotVectorUCopy[i] = 1;
         for (int i = 10; i < 16; i++)
         knotVectorUCopy[i] = 21;

         double* knotVectorVCopy = new double[noKnotsV];
         for (int i = 0; i < 3; i++)
         knotVectorVCopy[i] = 0;
         for (int i = 3; i < 6; i++)
         knotVectorVCopy[i] = 1;

         // On the control point weights
         double* controlPointWeightsCopy = new double[noUControlPoints * noVControlPoints];
         int counter = 0;
         for (int i = 0; i < noUControlPoints; i++) {
         for (int j = 0; j < noVControlPoints; j++) {
         controlPointWeightsCopy[counter] = ControlPointNet[i][j].getW();
         counter++;
         }
         }

         // Create an instance of the class NurbsBasis2D
         NurbsBasis2D* nurbsBasis2DCopy = new NurbsBasis2D(id, p, noKnotsU, knotVectorUCopy, q,
         noKnotsV, knotVectorVCopy, noUControlPoints, noVControlPoints, controlPointWeightsCopy);

         // Destroy the object
         delete nurbsBasis2DCopy;
         }*/

    }

    void tearDown() {
        delete nurbsBasis2D;
    }

    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        // Test the ID's
        CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getId()==1);
        CPPUNIT_ASSERT(nurbsBasis2D->getVBSplineBasis1D()->getId()==1);
        // Test the polynomial degrees
        CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree()==5);
        CPPUNIT_ASSERT(nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree()==2);

        // Test the number of knots
        CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getNoKnots()==16);
        CPPUNIT_ASSERT(nurbsBasis2D->getVBSplineBasis1D()->getNoKnots()==6);

        // Test the knot vectors
        for (int i = 0; i < 6; i++)
            CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==-4);
        for (int i = 6; i < 7; i++)
            CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==0);
        for (int i = 7; i < 10; i++)
            CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==1);
        for (int i = 10; i < 16; i++)
            CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->getKnotVector()[i]==21);
        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(nurbsBasis2D->getVBSplineBasis1D()->getKnotVector()[i]==0);
        for (int i = 3; i < 6; i++)
            CPPUNIT_ASSERT(nurbsBasis2D->getVBSplineBasis1D()->getKnotVector()[i]==1);
    }

    /***********************************************************************************************
     * \brief Test case: Test the knot span
     ***********/
    void testNurbs2DKnotSpans() {

        // The parametric coordinate on the Spline parameter space (from MATLAB this must be 6)
        double u = 43.999 / 7.0001;
        int uCorrectknotSpan = 9;
        double v = 0.12 / 10.31;
        int vCorrectknotSpan = 2;

        // Check the find knot span function
        CPPUNIT_ASSERT(nurbsBasis2D->getUBSplineBasis1D()->findKnotSpan(u)==uCorrectknotSpan);
        CPPUNIT_ASSERT(nurbsBasis2D->getVBSplineBasis1D()->findKnotSpan(v)==vCorrectknotSpan);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions
     ***********/
    void testNurbs2DBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = 49.999 / 4.0001;
        int uknotSpan = nurbsBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 2.12 / 10.86;
        int vknotSpan = nurbsBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);
        double* localBasisFunctions = new double[noBasisFunctions];
        nurbsBasis2D->computeLocalBasisFunctions(localBasisFunctions, u, uknotSpan, v, vknotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero basis functions at ( u , v ) = ( " << u << " , " << v << " )" << endl;
         cout << endl;
         for (int i=0; i<noBasisFunctions; i++) {
         cout << localBasisFunctions[i] << " " ;
         }*/

        // Values provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.018786001445944, 0.182942946982802,
                0.025512273410215, 0.135646299494851, 0.183500344551929, 0.260648432024145,
                0.004455525590841, 0.014200056662830, 0.008416131200953, 0.085547185676844,
                0.022255190529181, 0.008429819616952, 0.000392998304663, 0.001097904838662,
                0.010897694721708, 0.007980990742603, 0.021593136806345, 0.007697067398532 };

        // Assert messages if any error exists
        for (int i = 0; i < noBasisFunctions; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives with the efficient algorithm (Test 1)
     ***********/
    void testNurbsBasisFunctionsAndDerivatives2DTest1() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 129.9999 / 144.019;
        int uknotSpan = nurbsBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 0.012 / 10.86;
        int vknotSpan = nurbsBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);

        // On the basis function index
        int indexBasis = 0;

        // Return the functional values R, dR/du, dR/dv, d^2R/dv^2
        int derivDegree = 2;

        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
                * noBasisFunctions / 2];
        nurbsBasis2D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                derivDegree, u, uknotSpan, v, vknotSpan);

        /* // Printing the results for convenient debugging
         int vDerivOrder = 0;
         cout << endl;
         cout << endl;
         cout
         << "The non-zero NURBS basis functions and their derivatives in u-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " << vDerivOrder << " order in v-direction : "
         << endl;
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         // Compute the index of the basis function
         indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, i, vDerivOrder,
         k);

         // Print the functional value at (u,v)
         cout << localBasisFunctionsAndDerivatives[indexBasis] << " ";
         counter++;
         }
         counter = 0;
         cout << endl;
         }
         cout << endl;

         int uDerivOrder = 1;
         cout << endl;
         cout << endl;
         cout
         << "The non-zero NURBS basis functions and their derivatives in v-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " << uDerivOrder << " order in u-direction : "
         << endl;
         cout << endl;

         counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         // Compute the index of the basis function
         indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, uDerivOrder, i,
         k);

         // Print the functional value at (u,v)
         cout << localBasisFunctionsAndDerivatives[indexBasis] << " ";
         }
         cout << endl;
         }
         cout << endl;*/

        // Values provided by MATLAB
        double correctlocalBasisFunctions[] = { 0.000000003035577, 0.000001128444274,
                0.000147986265981, 0.762141830179219, 0.236726535712056, 0.000073745775615,
                0.000000000005037, 0.000000003079088, 0.000000155516541, 0.000824342195867,
                0.000083797003792, 0.000000110944972, 0.000000000000007, 0.000000000000967,
                0.000000000045271, 0.000000331593804, 0.000000029546786, 0.000000000655144 };
        double correctLocaldBasisFunctionsdu[] = { -0.000000155979579, -0.000045970222359,
                -0.004441151717863, -0.191153704793974, 0.195375689810006, 0.000407124528988,
                -0.000000000258816, -0.000000125434973, -0.000004667139527, -0.000206754252973,
                0.000069159536216, 0.000000612488229, -0.000000000000382, -0.000000000039377,
                -0.000000001358622, -0.000000083167439, 0.000000024385622, 0.000000003616821 };
        double correctLocalddBasisFunctionsddu[] = { 0.000006413456436, 0.001396361951662,
                0.086406995271943, -0.094747824711368, 0.005135704815084, 0.001805667526682,
                0.000000010641797, 0.000003810132169, 0.000090803811403, -0.000102480439707,
                0.000001817948607, 0.000002716490969, 0.000000000015696, 0.000000001196077,
                0.000000026433341, -0.000000041223025, 0.000000000641008, 0.000000016041226 };
        double correctLocaldBasisFunctionsdv[] = { -0.000000002500329, -0.000000929471323,
                -0.000121892585646, -0.627757161748128, -0.194985726126753, -0.000060742550740,
                0.000000004559296, 0.000002787121292, 0.000140770062710, 0.746175951291699,
                0.075851156635268, 0.000100424884909, 0.000000000013455, 0.000000001750663,
                0.000000081994708, 0.000600575580302, 0.000053514505483, 0.000001186583133 };
        double correctLocalddBasisFunctionsddv[] = { -0.000000002690489, -0.000001000161142,
                -0.000131162979069, -0.675500474701346, -0.209815130092382, -0.000065362252093,
                0.000000001614931, 0.000000987215533, 0.000049861623513, 0.264300119228149,
                0.026866947008397, 0.000035571112965, 0.000000012200896, 0.000001587444625,
                0.000074350149840, 0.544582516857665, 0.048525223203035, 0.001075955216972 };
        double correctLocalddBasisFunctionsdudv[] = { 0.000000128866703, 0.000038009584762,
                0.003677089342771, 0.255422509759903, -0.130494669921288, -0.000325858208384,
                -0.000000234273458, -0.000113540510385, -0.004224569375917, -0.187043313702530,
                0.062612423512034, 0.000554424784251, -0.000000000691388, -0.000000071317858,
                -0.000002460701837, -0.000150588453583, 0.000044170513214, 0.000006550792988 };

        // Check if the basis functions R themselves are correct
        int counterCheck = 0;

        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 0, 0, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctlocalBasisFunctions[i])<=Tol);
        }

        // Check if the derivatives dR/du are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 1, 0, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocaldBasisFunctionsdu[i])<=Tol);
        }

        // Check if the derivatives ddR/ddu are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 2, 0, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocalddBasisFunctionsddu[i])<=Tol);
        }

        // Check if the derivatives dR/dv are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 0, 1, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocaldBasisFunctionsdv[i])<=Tol);
        }

        // Check if the derivatives ddR/ddv are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 0, 2, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocalddBasisFunctionsddv[i])<=Tol);
        }

        // Check if the derivatives ddR/dudv are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 1, 1, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocalddBasisFunctionsdudv[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives with the efficient algorithm (Test 2)
     ***********/
    void testNurbsBasisFunctionsAndDerivatives2DTest2() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 169.9999 / 144.019;
        int uknotSpan = nurbsBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 5.67 / 10.86;
        int vknotSpan = nurbsBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);

        // Initialize the basis index
        int indexBasis = 0;

        // Return all derivatives up to second order
        int derivDegree = 2;

        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 1)
                * noBasisFunctions];
        nurbsBasis2D->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
                derivDegree, u, uknotSpan, v, vknotSpan);

        // On the derivative orders to be printed and compared
        int vDerivOrder = 0;
        int uDerivOrder = 1;

        /* // Print the results
         cout << endl;
         cout << endl;
         cout
         << "The non-zero NURBS basis functions and their derivatives in u-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " << vDerivOrder << " order in v-direction : "
         << endl;
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         // Compute the index of the basis function
         indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, i, vDerivOrder,
         k);

         // Print the functional value at (u,v)
         cout << localBasisFunctionsAndDerivatives[indexBasis] << " ";
         counter++;
         }
         counter = 0;
         cout << endl;
         }
         cout << endl;

         cout << endl;
         cout << endl;
         cout
         << "The non-zero NURBS basis functions and their derivatives in v-direction at ( u , v ) = ( "
         << u << " , " << v << " ) for " << uDerivOrder << " order in u-direction : "
         << endl;
         cout << endl;

         counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int k = 0; k < noBasisFunctions; k++) {
         // Compute the index of the basis function
         indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, uDerivOrder, i,
         k);

         // Print the functional value at (u,v)
         cout << localBasisFunctionsAndDerivatives[indexBasis] << " ";
         }
         cout << endl;
         }
         cout << endl;*/

        // Values provided by MATLAB
        double correctlocalBasisFunctions[] = { 0.323010060508582, 0.131505168880388,
                0.000123963281455, 0.000000710430068, 0.000000006466374, 0.000000000061800,
                0.345041960589514, 0.045973598924196, 0.000184182207313, 0.000002017949916,
                0.000000003532210, 0.000000000009002, 0.137073947559471, 0.016009375723639,
                0.001074140487682, 0.000000847915816, 0.000000015435554, 0.000000000037021 };
        double correctLocaldBasisFunctionsdu[] = { -0.056111911871781, 0.091362758492980,
                0.000482283200559, 0.000011798426899, 0.000000143561024, 0.000000001717728,
                -0.059939198346272, 0.031939997882401, 0.000716566900993, 0.000033512988314,
                0.000000078419172, 0.000000000250212, -0.023811893825424, 0.011122458077618,
                0.004178978695715, 0.000014081713631, 0.000000342687249, 0.000000001028981 };
        double correctLocalddBasisFunctionsddu[] = { 0.009214604046781, -0.024403020221538,
                0.001100440254561, 0.000130456137708, 0.000002392487507, 0.000000038249830,
                0.009843114611822, -0.008531183023113, 0.001635012503077, 0.000370555757628,
                0.000001306879014, 0.000000005571652, 0.003910349262500, -0.002970811891611,
                0.009535302856053, 0.000155702619368, 0.000005710985748, 0.000000022913014 };
        double correctLocaldBasisFunctionsdv[] = { -0.905575024629653, -0.368681385218402,
                -0.000347537322768, -0.000001991726590, -0.000000018128807, -0.000000000173260,
                0.415526766469325, 0.055365036969099, 0.000221806753344, 0.000002430174585,
                0.000000004253766, 0.000000000010841, 0.714444274075796, 0.083442601755665,
                0.005598536663172, 0.000004419429151, 0.000000080451782, 0.000000000192955 };
        double correctLocalddBasisFunctionsddv[] = { -0.311828690165540, -0.126952963933704,
                -0.000119672147746, -0.000000685837702, -0.000000006242533, -0.000000000059661,
                -2.299717831785261, -0.306415790869841, -0.001227581438879, -0.000013449713180,
                -0.000000023542315, -0.000000000059999, 2.708677369209219, 0.316356495817515,
                0.021225769609308, 0.000016755411389, 0.000000305017380, 0.000000000731551 };
        double correctLocalddBasisFunctionsdudv[] = { 0.214585753745614, -0.232822822188890,
                -0.001330125290797, -0.000032951518997, -0.000000401334027, -0.000000004804780,
                -0.011003846526137, 0.046616282608509, 0.000895603923191, 0.000040716790358,
                0.000000095064857, 0.000000000302922, -0.099805479974777, 0.060810092515316,
                0.021971746695334, 0.000073545762598, 0.000001788859983, 0.000000005369722 };

        // Check if the basis functions R themselves are correct
        int counterCheck = 0;

        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 0, 0, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctlocalBasisFunctions[i])<=Tol);
        }

        // Check if the derivatives dR/du are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 1, 0, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocaldBasisFunctionsdu[i])<=Tol);
        }

        // Check if the derivatives ddR/ddu are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 2, 0, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocalddBasisFunctionsddu[i])<=Tol);
        }

        // Check if the derivatives dR/dv are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 0, 1, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocaldBasisFunctionsdv[i])<=Tol);
        }

        // Check if the derivatives ddR/ddv are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 0, 2, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocalddBasisFunctionsddv[i])<=Tol);
        }

        // Check if the derivatives ddR/dudv are correct
        for (int i = 0; i < noBasisFunctions; i++) {
            // Compute the index
            indexBasis = nurbsBasis2D->indexDerivativeBasisFunction(derivDegree, 1, 1, i);

            // Assert message if failure occurs
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[indexBasis]-correctLocalddBasisFunctionsdudv[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constuctor of the class NurbsBasis2D
     ***********/
    void testNurbsBasis2DCopyConstructor() {

        NurbsBasis2D* nurbsBasis2DCopy = new NurbsBasis2D(*nurbsBasis2D);

        // Test the ID's
        CPPUNIT_ASSERT(nurbsBasis2DCopy->getId()==nurbsBasis2D->getId());
        CPPUNIT_ASSERT(
                nurbsBasis2DCopy->getUBSplineBasis1D()->getId()==nurbsBasis2D->getUBSplineBasis1D()->getId());
        CPPUNIT_ASSERT(
                nurbsBasis2DCopy->getVBSplineBasis1D()->getId()==nurbsBasis2D->getVBSplineBasis1D()->getId());

        // Test the polynomial degrees
        CPPUNIT_ASSERT(
                nurbsBasis2DCopy->getUBSplineBasis1D()->getPolynomialDegree()==nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree());
        CPPUNIT_ASSERT(
                nurbsBasis2DCopy->getVBSplineBasis1D()->getPolynomialDegree()==nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree());

        // Test the knot vectors
        CPPUNIT_ASSERT(
                nurbsBasis2DCopy->getUBSplineBasis1D()->getNoKnots()==nurbsBasis2D->getUBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < nurbsBasis2DCopy->getUBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    nurbsBasis2DCopy->getUBSplineBasis1D()->getKnotVector()[i]==nurbsBasis2D->getUBSplineBasis1D()->getKnotVector()[i]);
        CPPUNIT_ASSERT(
                nurbsBasis2DCopy->getVBSplineBasis1D()->getNoKnots()==nurbsBasis2D->getVBSplineBasis1D()->getNoKnots());
        for (int i = 0; i < nurbsBasis2DCopy->getVBSplineBasis1D()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    nurbsBasis2DCopy->getVBSplineBasis1D()->getKnotVector()[i]==nurbsBasis2D->getVBSplineBasis1D()->getKnotVector()[i]);

        // Test the Control Points Weights
        int counter = 0;
        for (int i = 0; i < nurbsBasis2DCopy->getUNoBasisFnc(); i++) {
            for (int j = 0; j < nurbsBasis2DCopy->getVNoBasisFnc(); j++) {
                CPPUNIT_ASSERT(
                        nurbsBasis2DCopy->getIGAControlPointWeights()[counter]==nurbsBasis2D->getIGAControlPointWeights()[counter]);
                counter++;
            }
        }

        // Free the memory on the heap from the pointers
        delete nurbsBasis2DCopy;

        // Test the copy constructor for memory leakage
        /* int noIterations = 1e9;

         for (int i = 0; i < noIterations; i++) {
         NurbsBasis2D* nurbsBasis2DCopy2 = new NurbsBasis2D(*nurbsBasis2D);
         delete nurbsBasis2DCopy2;
         }*/
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the 2D B-Spline basis functions for memory leakage
     ***********/
    void testNurbs2DBasisFunctions4Leakage() {
        // Auxiliary variables
        double u = 49.999 / 4.0001;
        int uknotSpan = nurbsBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 2.12 / 10.86;
        int vknotSpan = nurbsBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);
        int noIterations = 1e9;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions at another parametric location
            double* localBasisFunctions = new double[noBasisFunctions];
            nurbsBasis2D->computeLocalBasisFunctions(localBasisFunctions, u, uknotSpan, v,
                    vknotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the 2D NURBS basis functions and their derivatives for memory leakage
     ***********/
    void testNurbs2DBasisFunctionsAndDerivatives4Leakage() {
        // initialize variables
        double u = 169.9999 / 144.019;
        int uknotSpan = nurbsBasis2D->getUBSplineBasis1D()->findKnotSpan(u);
        double v = 5.67 / 10.86;
        int vknotSpan = nurbsBasis2D->getVBSplineBasis1D()->findKnotSpan(v);
        int noBasisFunctions = (nurbsBasis2D->getUBSplineBasis1D()->getPolynomialDegree() + 1)
                * (nurbsBasis2D->getVBSplineBasis1D()->getPolynomialDegree() + 1);
        int derivDegree = 2;
        int noIterations = 1e9;

        for (int i = 0; i < noIterations; i++) {
            // Initialize the array holding the NURBS basis functions and their derivatives
            double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                    * (derivDegree + 2) * noBasisFunctions / 2];

            // Compute the non-zero basis functions and their derivatives at another parametric location
            nurbsBasis2D->computeLocalBasisFunctionsAndDerivatives(
                    localBasisFunctionsAndDerivatives, derivDegree, u, uknotSpan, v, vknotSpan);

            // Clean up the array from the heap
            delete[] localBasisFunctionsAndDerivatives;
        }
    }

// Make the tests
CPPUNIT_TEST_SUITE(TestNurbsBasis2D);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testNurbs2DKnotSpans);
    CPPUNIT_TEST(testNurbs2DBasisFunctions);
    CPPUNIT_TEST(testNurbsBasisFunctionsAndDerivatives2DTest1);
    CPPUNIT_TEST(testNurbsBasisFunctionsAndDerivatives2DTest2);
    CPPUNIT_TEST(testNurbsBasis2DCopyConstructor);

// Make the tests for leakage
    // CPPUNIT_TEST(testNurbs2DBasisFunctions4Leakage);
    // CPPUNIT_TEST(testNurbs2DBasisFunctionsAndDerivatives4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;
}
;

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestNurbsBasis2D);
