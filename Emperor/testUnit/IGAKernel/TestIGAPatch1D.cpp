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
#include "IGAPatch1D.h"

using namespace std;

namespace EMPIRE {

/********//**
 * \brief Test the class IGAPatch1D
 ***********/

class TestIGAPatch1D: public CppUnit::TestFixture {

private:
    IGAPatch1D* igaPatch1D;
    double Tol;

public:
    void setUp() {
        // Assign a tolerance value (corresponding to maximum accuracy provided by MATLAB)
        Tol = 1e-15;

        // Provide an id for the basis
        int id_basis = 3;

        // Provide an id for the patch itself
        int id_patch = 1;

        // The polynomial degree
        int p = 8;

        // Number of knots
        int noKnots = 26;

        // The knot vector
        double* knotVector = new double[noKnots];

        for (int i = 0; i < 9; i++)
            knotVector[i] = -13.4;
        for (int i = 9; i < 12; i++)
            knotVector[i] = -3.34;
        for (int i = 12; i < 14; i++)
            knotVector[i] = 1.01099;
        for (int i = 14; i < 15; i++)
            knotVector[i] = 2.000001;
        for (int i = 15; i < 17; i++)
            knotVector[i] = 9;
        for (int i = 17; i < 26; i++)
            knotVector[i] = 2.34;

        // The Control Point net
        int noControlPoints = noKnots - p - 1;
        IGAControlPoint* controlPointNet = new IGAControlPoint[noControlPoints];
        double* controlPointWeights = new double[noControlPoints];

        controlPointNet[0] = IGAControlPoint(0, 0, 0, 1.698550724637681, 15.731999999999999);
        controlPointNet[1] = IGAControlPoint(1, 1.000000000000000, 1.000000000000000,
                0.751282051282051, 17.327999999999999);
        controlPointNet[2] = IGAControlPoint(2, 1.111111111111111, 1.111111111111111,
                0.265158371040724, 9.575999999999999);
        controlPointNet[3] = IGAControlPoint(3, 1.250000000000000, 1.250000000000000,
                4.842975206611571, 5.517599999999999);
        controlPointNet[4] = IGAControlPoint(4, 1.333333333333333, 1.333333333333333,
                1.347126436781609, 19.835999999999995);
        controlPointNet[5] = IGAControlPoint(5, 1.408450704225352, 3.549295774647887,
                0.266254713980644, 9.161039999999998);
        controlPointNet[6] = IGAControlPoint(6, 1.451378809869376, 3.657474600870827,
                6.584269662921349, 4.058400000000000);
        controlPointNet[7] = IGAControlPoint(7, 1.694915254237288, 7.544910179640719,
                0.710303030303030, 19.379999999999999);
        controlPointNet[8] = IGAControlPoint(8, 1.724137931034483, 11.351351351351351,
                0.180307692307692, 11.399999999999999);
        controlPointNet[9] = IGAControlPoint(9, 1.754385964912281, 4.421052631578947,
                5.854145854145855, 4.564559999999999);
        controlPointNet[10] = IGAControlPoint(10, 1.999960000799984, 5.039899202015960,
                1.171976560468791, 11.400456000000000);
        controlPointNet[11] = IGAControlPoint(11, 2.000040000800016, 4.344827586206897,
                1.010344827586207, 11.764799999999999);
        controlPointNet[12] = IGAControlPoint(12, 2.222222222222222, 3.485477178423236,
                0.810511756569848, 10.168799999999999);
        controlPointNet[13] = IGAControlPoint(13, 2.369668246445498, 2.802802802802803,
                0.651762873985096, 18.198960000000000);
        controlPointNet[14] = IGAControlPoint(14, 2.500000000000000, 2.566191446028513,
                0.651111111111111, 20.520000000000000);
        controlPointNet[15] = IGAControlPoint(15, 2.898550724637681, 2.468168462291870,
                5.274527452745275, 5.066159999999999);
        controlPointNet[16] = IGAControlPoint(16, 9.999999909000000, 0, 0, 4.560000000000000);

        for (int i = 0; i < noControlPoints; i++)
            controlPointWeights[1] = controlPointNet[i].getW();

        // Test just one object of the class (that works pretty also)
        igaPatch1D = new IGAPatch1D(id_patch, id_basis, p, noKnots, knotVector, noControlPoints,
                controlPointNet);

        // nurbsBasis1D->printPolynomialDegree();
        // nurbsBasis1D->printNoKnots();
        // nurbsBasis1D->printKnotVector();
        // nurbsBasis1D->printNoBasisFunctions();
        // nurbsBasis1D->printControlPointNet();
        // NurbsBasis1D nurbsBasis1DCopy = *nurbsBasis1D;
        // delete nurbsBasis1D;

        // Make test on memory leakage (that works pretty fine)
        /* for (int j = 1; j < 1e9; j++) {

         // The knot vector to feed the instance igaPatch1DCopy
         double* knotVectorCopy = new double[noKnots];

         for (int i = 0; i < 9; i++)
         knotVectorCopy[i] = -13.4;
         for (int i = 9; i < 12; i++)
         knotVectorCopy[i] = -3.34;
         for (int i = 12; i < 14; i++)
         knotVectorCopy[i] = 1.01099;
         for (int i = 14; i < 15; i++)
         knotVectorCopy[i] = 2.000001;
         for (int i = 15; i < 17; i++)
         knotVectorCopy[i] = 9;
         for (int i = 17; i < 26; i++)
         knotVectorCopy[i] = 2.34;

         // Fill up the array of the Control Point weights and feed them to the instance igaPatch1DCopy
         IGAControlPoint* controlPointNetCopy = new IGAControlPoint[noControlPoints];

         for (int i = 0; i < noControlPoints; i++) {
         controlPointNetCopy[i] = controlPointNet[i];
         }

         // Create an object of the class NurbsBasis1D
         IGAPatch1D* igaPatch1DCopy = new IGAPatch1D(id_patch, id_basis, p, noKnots,
         knotVectorCopy, noControlPoints, controlPointNetCopy);

         // Free the memory on the heap from the pointer
         delete igaPatch1DCopy;
         }*/
    }

    void tearDown() {
        delete igaPatch1D;
    }

    /***********************************************************************************************
     * \brief Test case: Test the constructor
     ***********/
    void testConstructor() {
        CPPUNIT_ASSERT(igaPatch1D->getId()==1);
        CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getId()==3);
        CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getPolynomialDegree()==8);
        CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getNoKnots()==26);
        for (int i = 0; i < 9; i++)
            CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getKnotVector()[i]==-13.4);
        for (int i = 9; i < 12; i++)
            CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getKnotVector()[i]==-3.34);
        for (int i = 12; i < 14; i++)
            CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getKnotVector()[i]==1.01099);
        for (int i = 14; i < 15; i++)
            CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getKnotVector()[i]==2.000001);
        for (int i = 15; i < 17; i++)
            CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getKnotVector()[i]==9);
        for (int i = 17; i < 26; i++)
            CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->getKnotVector()[i]==2.34);
        CPPUNIT_ASSERT(igaPatch1D->getNoControlPoints()==17);
    }

    /***********************************************************************************************
     * \brief Test case: Test the knot span
     ***********/
    void testIGAPatch1DKnotSpan() {

        // The parametric coordinate on the Spline parameter space (from MATLAB this must be 6)
        double u = -1.45601090807;
        int CorrectknotSpan = 11;

        // Check the find knot span function
        CPPUNIT_ASSERT(igaPatch1D->getIGABasis()->findKnotSpan(u)==CorrectknotSpan);
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constructor for proper functionality
     ***********/
    void testIGAPatch1DCopyConstructor() {

        IGAPatch1D igaPatch1DCopy(*igaPatch1D);

        CPPUNIT_ASSERT(igaPatch1DCopy.getId()==igaPatch1D->getId());
        CPPUNIT_ASSERT(igaPatch1DCopy.getIGABasis()->getId()==igaPatch1D->getIGABasis()->getId());
        CPPUNIT_ASSERT(
                igaPatch1DCopy.getIGABasis()->getPolynomialDegree()==igaPatch1D->getIGABasis()->getPolynomialDegree());
        CPPUNIT_ASSERT(
                igaPatch1DCopy.getIGABasis()->getNoKnots()==igaPatch1D->getIGABasis()->getNoKnots());
        for (int i = 1; i < igaPatch1DCopy.getIGABasis()->getNoKnots(); i++)
            CPPUNIT_ASSERT(
                    igaPatch1DCopy.getIGABasis()->getKnotVector()[i]==igaPatch1D->getIGABasis()->getKnotVector()[i]);
        CPPUNIT_ASSERT(igaPatch1DCopy.getNoControlPoints()==igaPatch1D->getNoControlPoints());
        for (int i = 1; i < igaPatch1DCopy.getNoControlPoints(); i++) {
            CPPUNIT_ASSERT(
                    igaPatch1DCopy.getControlPointNet()[i].getId()==igaPatch1D->getControlPointNet()[i].getId());
            CPPUNIT_ASSERT(
                    igaPatch1DCopy.getControlPointNet()[i].getX()==igaPatch1D->getControlPointNet()[i].getX());
            CPPUNIT_ASSERT(
                    igaPatch1DCopy.getControlPointNet()[i].getY()==igaPatch1D->getControlPointNet()[i].getY());
            CPPUNIT_ASSERT(
                    igaPatch1DCopy.getControlPointNet()[i].getZ()==igaPatch1D->getControlPointNet()[i].getZ());
            CPPUNIT_ASSERT(
                    igaPatch1DCopy.getControlPointNet()[i].getW()==igaPatch1D->getControlPointNet()[i].getW());
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the copy constructor for leakage
     ***********/
    void testIGAPatch1DCopyConstructor4Leakage() {

        for (int i = 1; i < 1e9; i++) {
            IGAPatch1D* igaPatch1DCopy = new IGAPatch1D(*igaPatch1D);
            delete igaPatch1DCopy;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local NURBS basis functions
     ***********/
    void testIGAPatch1DBasisFunctions() {
        // Compute the non-zero basis functions at another parametric location
        double u = -4.38900001;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        int noLocalBasisFunctions = igaPatch1D->getIGABasis()->getPolynomialDegree() + 1;

        double* localBasisFunctions = new double[noLocalBasisFunctions];
        igaPatch1D->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

        /* cout << endl;
         cout << endl;
         cout << "The non-zero basis functions at u = " << u << " :";
         cout << endl;
         for (int i = 0; i < noLocalBasisFunctions; i++) {
         cout << localBasisFunctions[i] << " ";
         }
         cout << endl;
         cout << endl;*/

        // Value provided by MATLAB
        double CorrectlocalBasisFunctions[] = { 0.000000022469071, 0.000001700736804,
                0.000028257800147, 0.007902201809734, 0.184134927644108, 0.251399766092729,
                0.175494512813427, 0.350045048202155, 0.030993562431824 };
        for (int i = 0; i < igaPatch1D->getIGABasis()->getPolynomialDegree() + 1; i++) {
            CPPUNIT_ASSERT(fabs(localBasisFunctions[i]-CorrectlocalBasisFunctions[i])<=Tol);
        }

        // Clear the heap from the pointer
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the local basis functions and their derivatives
     ***********/
    void testIGAPatch1DBasisFunctionsAndDerivatives() {
        // Compute the non-zero basis functions and their derivatives at another parametric location
        double u = 1.900000000009;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        int derivDegree = 2;
        double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1)
                * (igaPatch1D->getIGABasis()->getPolynomialDegree() + 1)];

        igaPatch1D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);

        /*cout << endl;
         cout << endl;
         cout << "The non-zero NURBS basis functions and their derivatives at u = " << u << " :";
         cout << endl;

         int counter = 0;
         for (int i = 0; i <= derivDegree; i++) {
         for (int j = 0; j <= igaPatch1D->getIGABasis()->getPolynomialDegree(); j++) {
         cout << localBasisFunctionsAndDerivatives[counter] << " ";
         counter += 1;
         }
         cout << endl;
         }*/

        // igaPatch1D->printControlPointNet();
        // Values provided by MATLAB
        double CorrectlocalBasisFunctionsAndDerivatives[] = { 0.000000000000015, 0.005566997924493,
                0.238884366508277, 0.244756827621970, 0.076507240832072, 0.226075930036191,
                0.198614581661947, 0.007166348317694, 0.002427707097340, -0.000000000001205,
                -0.006176882341997, -0.170552935072810, -0.140197747972416, -0.021169016995988,
                0.036497288219350, 0.233347855685965, 0.046363277599723, 0.021888160879378,
                0.000000000084348, 0.005332585745941, 0.045735656197950, -0.019998956676438,
                -0.048694089492731, -0.271481856770363, -0.104561170920791, 0.221176444249632,
                0.172491387582452 };

        for (int i = 0;
                i < (derivDegree + 1) * (igaPatch1D->getIGABasis()->getPolynomialDegree() + 1); i++)
            CPPUNIT_ASSERT(
                    fabs(localBasisFunctionsAndDerivatives[i]-CorrectlocalBasisFunctionsAndDerivatives[i])<=Tol);

        // Clear the heap from the pointer
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of a point on the IGA patch
     ***********/
    void testIGAPatch1DPointOnCurve() {
        // Compute the Cartesian coordinates of a point on the 1D IGA patch
        double u = .0000009876;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        double cartesianCoordinatesPointOnPatch[3];
        igaPatch1D->computeCartesianCoordinates(cartesianCoordinatesPointOnPatch, u, knotSpan);

        /* cout << endl;
         cout << "Point on Patch 1D: ";
         for (int i = 0; i < 3; i++) {
         cout << cartesianCoordinatesPointOnPatch[i] << " ";
         }
         cout << endl; */

        // Values provided by MATLAB
        double correctCartesianCoordinatesPointOnPatch[] = { 1.709407917665832, 8.541357215107709,
                0.889957939293058 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(
                    fabs(cartesianCoordinatesPointOnPatch[i]-correctCartesianCoordinatesPointOnPatch[i])<=Tol);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of a point on the IGA patch with precomputing the basis functions
     ***********/
    void testIGAPatch1DPointOnCurveMethod2() {

        // Find the knot span of the given parametric coordinate
        double u = .0000009876;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);

        // Compute the local basis functions at u
        int noLocalBasisFunctions = igaPatch1D->getIGABasis()->getPolynomialDegree() + 1;
        double* localBasisFunctions = new double[noLocalBasisFunctions];
        igaPatch1D->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

        // Compute the Cartesian Coordinates of the point with parametric coordinate u
        double cartesianCoordinatesPointOnPatch[3];
        igaPatch1D->computeCartesianCoordinates(cartesianCoordinatesPointOnPatch, u, knotSpan);

        /* cout << endl;
         cout << "Point on Patch 1D: ";
         for (int i = 0; i < 3; i++) {
         cout << cartesianCoordinatesPointOnPatch[i] << " ";
         }
         cout << endl; */

        // Values provided by MATLAB
        double correctCartesianCoordinatesPointOnPatch[] = { 1.709407917665832, 8.541357215107709,
                0.889957939293058 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT(
                    fabs(cartesianCoordinatesPointOnPatch[i]-correctCartesianCoordinatesPointOnPatch[i])<=Tol);

        // Free the memory from the heap
        delete[] localBasisFunctions;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vector on the IGA 1D patch
     ***********/
    void testIGAPatch1DBaseVector() {
        // Compute the Cartesian coordinates of a point on the 1D IGA patch
        double u = -7.0123404321;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        double baseVector[3];
        igaPatch1D->computeBaseVector(baseVector, u, knotSpan);

        /*cout << endl;
         cout << "Base vector on Patch 1D at u = " << u << ": ";
         for (int i = 0; i < 3; i++) {
         cout << baseVector[i] << " ";
         }
         cout << endl;*/

        // Values provided by MATLAB
        double correctBaseVector[] = { 0.040891913044100, 0.529867297482033, -0.054352652178514 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT( fabs(baseVector[i]-correctBaseVector[i])<=Tol);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vector on the IGA 1D patch with precomputing the basis functions and their derivatives
     ***********/
    void testIGAPatch1DBaseVectorMethod2() {

        // Find the correct knot span of the given parametric coordinate
        double u = -7.0123404321;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);

        // Compute the basis functions and their derivatives at u
        int noLocalBasisFunctions = igaPatch1D->getIGABasis()->getPolynomialDegree() + 1;
        int derivDegree = 1;
        double* localBasisFunctionsAndDerivatives = new double[noLocalBasisFunctions
                * (derivDegree + 1)];
        igaPatch1D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree, u, knotSpan);

        // Compute the Cartesian coordinates of a point on the 1D IGA patch
        double baseVector[3];
        igaPatch1D->computeBaseVector(baseVector, localBasisFunctionsAndDerivatives, u, knotSpan);

        /*cout << endl;
         cout << "Base vector on Patch 1D at u = " << u << ": ";
         for (int i = 0; i < 3; i++) {
         cout << baseVector[i] << " ";
         }
         cout << endl;*/

        // Values provided by MATLAB
        double correctBaseVector[] = { 0.040891913044100, 0.529867297482033, -0.054352652178514 };

        for (int i = 0; i < 3; i++)
            CPPUNIT_ASSERT( fabs(baseVector[i]-correctBaseVector[i])<=Tol);

        // Free the memory from the heap
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of the base vector and its derivatives on the 1D IGA patch
     ***********/
    void testIGAPatch1DBaseVectorAndDerivatives() {
        // Compute the Cartesian coordinates of a point on the 1D IGA patch
        double u = -4.2000009;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        int derivDegree = 1;
        double baseVectorAndDerivatives[3 * (derivDegree + 1)];
        igaPatch1D->computeBaseVectorAndDerivatives(baseVectorAndDerivatives, derivDegree, u,
                knotSpan);

        /*cout << endl;
         cout << "Base vector and its derivatives on Patch 1D at u = " << u << ": ";
         for (int i = 0; i < 3 * (derivDegree + 1); i++) {
         cout << baseVectorAndDerivatives[i] << " ";
         }
         cout << endl;*/

        // Values provided by MATLAB
        double correctBaseVectorAndDerivatives[] = { 0.081305824604765, 1.404342898611633,
                0.028677665858097, -0.000047908589758, 0.012782932929461, -0.196076895166805 };

        for (int i = 0; i < 3 * (derivDegree + 1); i++)
            CPPUNIT_ASSERT(
                    fabs(baseVectorAndDerivatives[i]-correctBaseVectorAndDerivatives[i])<=Tol);
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of the base vector and its derivatives on the 1D IGA patch with precomputing the basis functions and their derivatives
     ***********/
    void testIGAPatch1DBaseVectorAndDerivativesMethod2() {

        // Auxiliary variables and knot span index
        double u = -4.2000009;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        int derivDegree = 9;

        // Compute the basis functions and their derivatives (one order higher than the derivatives for the base vector)
        int noLocalBasisFunctions = igaPatch1D->getIGABasis()->getPolynomialDegree() + 1;
        double* localBasisFunctionsAndDerivatives = new double[noLocalBasisFunctions
                * (derivDegree + 2)];
        igaPatch1D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(
                localBasisFunctionsAndDerivatives, derivDegree + 1, u, knotSpan);

        // Compute the Cartesian coordinates of a point on the 1D IGA patch
        double baseVectorAndDerivatives[3 * (derivDegree + 1)];
        igaPatch1D->computeBaseVectorAndDerivatives(baseVectorAndDerivatives, derivDegree, u,
                knotSpan);

        /* cout << endl;
        cout << "Base vector and its derivatives on Patch 1D at u = " << u << ": " << endl;
        int counter = 0;
        for (int i = 0; i < derivDegree + 1; i++) {
            for (int j = 0; j < 3; j++) {
                cout << baseVectorAndDerivatives[counter] << " ";
                counter++;
            }
            cout << endl;
        }
        cout << endl; */

        // Values provided by MATLAB
        double correctBaseVectorAndDerivatives[] = { 0.081305824604765, 1.404342898611633,
                0.028677665858097, -0.000047908589758, 0.012782932929461, -0.196076895166805,
                -0.034075576858149, -0.527243963435928, -0.210608175225937, -0.007400924424824,
                -0.029021307006914, 0.208657290269712, 0.069539364954965, 1.107063596675528,
                0.531181458629499, 0.041927191455165, 0.503453072815761, -0.373441518908986,
                -0.281334035129307, -4.585170551232989, -2.302870899623261, -0.416157910904029,
                -5.634626204065550, 1.002653371302840, 1.758843800367119, 29.321563626545071,
                17.713143261104175, 4.898682692325160, 69.935671579710402, 3.429763873177842 };

        for (int i = 0; i < 3 * (derivDegree + 1); i++)
            CPPUNIT_ASSERT(
                    fabs(baseVectorAndDerivatives[i]-correctBaseVectorAndDerivatives[i])<=Tol*1e2);

        // Free the memory on the heap
        delete[] localBasisFunctionsAndDerivatives;
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions for memory leakage
     ***********/
    void testIGAPatch1DBasisFunctions4Leakage() {

        // Compute the non-zero basis functions at another parametric location
        double u = 1.900000000009;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        int noLocalBasisFunctions = igaPatch1D->getIGABasis()->getPolynomialDegree() + 1;

        for (int i = 1; i < 1e9; i++) {
            double* localBasisFunctions = new double[noLocalBasisFunctions];
            igaPatch1D->getIGABasis()->computeLocalBasisFunctions(localBasisFunctions, u, knotSpan);

            // Free the memory from the heap
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the NURBS basis functions and their derivatives for memory leakage (efficient algorithm)
     ***********/
    void testIGAPatch1DBasisFunctionsAndDerivatives4Leakage() {
        // initialize variables
        int noIterations = 1000000000;
        int knotSpan = 0;
        int derivDegree = 2;
        double u = 0.0;

        for (int i = 0; i < noIterations; i++) {
            // Compute the non-zero basis functions at another parametric location
            knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
            double* localBasisFunctions = new double[(derivDegree + 1)
                    * (igaPatch1D->getIGABasis()->getPolynomialDegree() + 1)];

            igaPatch1D->getIGABasis()->computeLocalBasisFunctionsAndDerivatives(localBasisFunctions,
                    derivDegree, u, knotSpan);

            // Clear the heap from the pointer
            delete[] localBasisFunctions;
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the Coordinates of a point on the IGA patch for leakage
     ***********/
    void testIGAPatch1DPointOnCurve4Leakage() {
        // Compute the Cartesian coordinates of a point on the 1D IGA patch
        double u = .0000009876;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);

        for (int i = 0; i < 1e9; i++) {
            double cartesianCoordinatesPointOnPatch[3];
            igaPatch1D->computeCartesianCoordinates(cartesianCoordinatesPointOnPatch, u, knotSpan);
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vector for memory leakage
     ***********/
    void testIGAPatch1DBaseVector4Leakage() {
        // Compute the Cartesian coordinates of a point on the 1D IGA patch iteratively
        double u = -7.0123404321;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);

        for (int i = 0; i < 1e9; i++) {
            double baseVector[3];
            igaPatch1D->computeBaseVector(baseVector, u, knotSpan);
        }
    }

    /***********************************************************************************************
     * \brief Test case: Test the computation of the base vector and its derivatives for memory leakage
     ***********/
    void testIGAPatch1DBaseVectorAndDerivatives4Leakage() {

        // Initialize variables
        double u = -4.2000009;
        int knotSpan = igaPatch1D->getIGABasis()->findKnotSpan(u);
        int derivDegree = 9;

        for (int i = 0; i < 1e9; i++) {
            // Compute the Cartesian coordinates of the base vector and its derivatives at the parametric location u
            double baseVectorAndDerivatives[3 * (derivDegree + 1)];
            igaPatch1D->computeBaseVectorAndDerivatives(baseVectorAndDerivatives, derivDegree, u,
                    knotSpan);
        }
    }

// Make the tests
CPPUNIT_TEST_SUITE(TestIGAPatch1D);
    CPPUNIT_TEST(testConstructor);
    CPPUNIT_TEST(testIGAPatch1DKnotSpan);
    CPPUNIT_TEST(testIGAPatch1DCopyConstructor);
    CPPUNIT_TEST(testIGAPatch1DBasisFunctions);
    CPPUNIT_TEST(testIGAPatch1DBasisFunctionsAndDerivatives);
    CPPUNIT_TEST(testIGAPatch1DPointOnCurve);
    CPPUNIT_TEST(testIGAPatch1DPointOnCurveMethod2);
    CPPUNIT_TEST(testIGAPatch1DBaseVector);
    CPPUNIT_TEST(testIGAPatch1DBaseVectorMethod2);
    CPPUNIT_TEST(testIGAPatch1DBaseVectorAndDerivatives);
    CPPUNIT_TEST(testIGAPatch1DBaseVectorAndDerivativesMethod2);

// Make the tests for leakage
    // CPPUNIT_TEST(testIGAPatch1DCopyConstructor4Leakage);
    // CPPUNIT_TEST(testIGAPatch1DBasisFunctions4Leakage);
    // CPPUNIT_TEST(testIGAPatch1DBasisFunctionsAndDerivatives4Leakage);
    // CPPUNIT_TEST(testIGAPatch1DPointOnCurve4Leakage);
    // CPPUNIT_TEST(testIGAPatch1DBaseVector4Leakage);
    // CPPUNIT_TEST(testIGAPatch1DBaseVectorAndDerivatives4Leakage);

    CPPUNIT_TEST_SUITE_END()
    ;
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION(EMPIRE::TestIGAPatch1D);
