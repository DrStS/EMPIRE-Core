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
#include <string>
#include "cppunit/TestFixture.h"
#include "cppunit/TestAssert.h"
#include "cppunit/extensions/HelperMacros.h"

#include <math.h>
#include <iostream>
#include <vector>

#include "MortarMath.h"

namespace EMPIRE {
using namespace std;
using namespace MortarMath;

class IntegrandFunctionConstant1: public IntegrandFunction {
public:
    IntegrandFunctionConstant1() {
    }
    virtual ~IntegrandFunctionConstant1() {
    }
    double operator()(double *gaussPoint) {
        return 1.0;
    }
};
class IntegrandFunctionXY: public IntegrandFunction {
public:
    IntegrandFunctionXY() {
    }
    virtual ~IntegrandFunctionXY() {
    }
    double operator()(double *gaussPoint) {
        return gaussPoint[0] * gaussPoint[1];
    }
};

/********//**
 * \brief Test functions in MortarMath.h
 ***********/
class TestMortarMath: public CppUnit::TestFixture {
public:
    void setUp() {
    }

    void tearDown() {
    }

    void testGaussQuadratureOnTriangle() {
        double triangle[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
        vector<int> numGaussPoints;
        numGaussPoints.push_back(3);
        numGaussPoints.push_back(6);
        numGaussPoints.push_back(7);
        numGaussPoints.push_back(12);
        const double EPS = 1E-10;
        for (int i = 0; i < numGaussPoints.size(); i++) {
            int numGaussPointsI = numGaussPoints[i];
            GaussQuadratureOnTriangle gaussQuadratureOnTriangle(triangle, numGaussPointsI);
            {
                IntegrandFunctionConstant1 integrandFunc;
                gaussQuadratureOnTriangle.setIntegrandFunc(&integrandFunc);
                // intigrate(f(x,y,z)=1) on the triangle should equal to 0.5
                CPPUNIT_ASSERT(fabs(gaussQuadratureOnTriangle.computeIntegral() - 0.5) < EPS);
            }
            {
                if (i == 0)
                    continue; // in this case, with 3 Gauss points is not correct
                IntegrandFunctionXY integrandFunc;
                gaussQuadratureOnTriangle.setIntegrandFunc(&integrandFunc);
                // intigrate(f(x,y,z)=xy) on the triangle should equal to 1.0/24.0
                CPPUNIT_ASSERT(fabs(gaussQuadratureOnTriangle.computeIntegral() - 1.0/24.0) < EPS);
            }
        }
    }

    void testGaussQuadratureOnQuad() {
        double quad[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 1, 0 };
        vector<int> numGaussPoints;
        numGaussPoints.push_back(1);
        numGaussPoints.push_back(4);
        numGaussPoints.push_back(9);
        const double EPS = 1E-10;
        for (int i = 0; i < numGaussPoints.size(); i++) {
            int numGaussPointsI = numGaussPoints[i];
            GaussQuadratureOnQuad GaussQuadratureOnQuad(quad, numGaussPointsI);
            {
                IntegrandFunctionConstant1 integrandFunc;
                GaussQuadratureOnQuad.setIntegrandFunc(&integrandFunc);
                // intigrate(f(x,y,z)=1) on the quad should equal to 3.0
                CPPUNIT_ASSERT(fabs(GaussQuadratureOnQuad.computeIntegral() - 3.0) < EPS);
            }
            {
                if (i == 0)
                    continue; // in this case, with 1 Gauss point is not correct
                IntegrandFunctionXY integrandFunc;
                GaussQuadratureOnQuad.setIntegrandFunc(&integrandFunc);
                // intigrate(f(x,y,z)=xy) on the quad should equal to 17.0/6.0
                CPPUNIT_ASSERT(fabs(GaussQuadratureOnQuad.computeIntegral() - 17.0/6.0) < EPS);
            }
        }
    }

    void testMassMatrixOfTriangle() {
        { // 1. verify that the computation of the mass matrix is correct
            double triangle[] = { 0, 0, 0, 2, 0, 0, 1, 2, 0 };

            vector<int> numGaussPoints;
            //numGaussPoints.push_back(3); // not correct in this case
            numGaussPoints.push_back(6);
            numGaussPoints.push_back(7);
            numGaussPoints.push_back(12);
            double massMatrix[3 * 3];
            const double EPS = 1E-10;
            // 1. test on normal mass matrix
            for (int i = 0; i < numGaussPoints.size(); i++) {
                computeMassMatrixOfTrianlge(triangle, numGaussPoints[i], false, massMatrix);
                //printGeneralMatrix(massMatrix, 3, 3);

                // testing with the mass matrix of the triangle computed by hand
                CPPUNIT_ASSERT(massMatrix[0] - 1.0/3.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[1] - 1.0/6.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[2] - 1.0/6.0 < EPS);

                CPPUNIT_ASSERT(massMatrix[3] - 1.0/6.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[4] - 1.0/3.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[5] - 1.0/6.0 < EPS);

                CPPUNIT_ASSERT(massMatrix[6] - 1.0/6.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[7] - 1.0/6.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[8] - 1.0/3.0 < EPS);
            }
            // 2. test on dual mass matrix
            for (int i = 0; i < numGaussPoints.size(); i++) {
                computeMassMatrixOfTrianlge(triangle, numGaussPoints[i], true, massMatrix);
                //printGeneralMatrix(massMatrix, 3, 3);

                // testing with the mass matrix of the triangle computed by hand
                CPPUNIT_ASSERT(massMatrix[0] - 2.0/3.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[1] - 0.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[2] - 0.0 < EPS);

                CPPUNIT_ASSERT(massMatrix[3] - 0.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[4] - 2.0/3.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[5] - 0.0 < EPS);

                CPPUNIT_ASSERT(massMatrix[6] - 0.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[7] - 0.0 < EPS);
                CPPUNIT_ASSERT(massMatrix[8] - 2.0/3.0 < EPS);
            }
        }
        { // 2. check how many gauss points are enough --- The result shows at least 6 gauss points are required
          // The triangle has a singular mass matrix if using only 3 gauss points.
          // However, for other triangles, we may not get a singular mass matrix with 3 gauss points, and it is
          // still consistent because C_BB*1==C_BA*1 even both C_BB and C_BA are wrong!
            double triangle[] = { -4.76182, 0.0112458, -0.156187, -4.75928, 0.0175115, -0.134099,
                    -4.78426, 0.0134955, -0.147207 };

            const double EPS = 1E-10;
            { // 1. test on normal mass matrix
                vector<int> numGaussPoints;
                numGaussPoints.push_back(7);
                numGaussPoints.push_back(6);
                //numGaussPoints.push_back(3); // not correct in this case
                double massMatrixCorrect[3 * 3]; // the one computed by 12 gauss points is regarded as the correct one
                computeMassMatrixOfTrianlge(triangle, 12, false, massMatrixCorrect);
                double massMatrix[3 * 3];

                for (int i = 0; i < numGaussPoints.size(); i++) {
                    computeMassMatrixOfTrianlge(triangle, numGaussPoints[i], false, massMatrix);
                    //printGeneralMatrix(massMatrix, 3, 3);

                    CPPUNIT_ASSERT(massMatrix[0] - massMatrixCorrect[0] < EPS);
                    CPPUNIT_ASSERT(massMatrix[1] - massMatrixCorrect[1] < EPS);
                    CPPUNIT_ASSERT(massMatrix[2] - massMatrixCorrect[2] < EPS);

                    CPPUNIT_ASSERT(massMatrix[3] - massMatrixCorrect[3] < EPS);
                    CPPUNIT_ASSERT(massMatrix[4] - massMatrixCorrect[4] < EPS);
                    CPPUNIT_ASSERT(massMatrix[5] - massMatrixCorrect[5] < EPS);

                    CPPUNIT_ASSERT(massMatrix[6] - massMatrixCorrect[6] < EPS);
                    CPPUNIT_ASSERT(massMatrix[7] - massMatrixCorrect[7] < EPS);
                    CPPUNIT_ASSERT(massMatrix[8] - massMatrixCorrect[8] < EPS);
                }
            }
            { // 2. test on dual mass matrix
                vector<int> numGaussPoints;
                numGaussPoints.push_back(7);
                numGaussPoints.push_back(6);
                //numGaussPoints.push_back(3); // not correct in this case
                double massMatrixCorrect[3 * 3]; // the one computed by 12 gauss points is regarded as the correct one
                computeMassMatrixOfTrianlge(triangle, 12, true, massMatrixCorrect);
                double massMatrix[3 * 3];

                for (int i = 0; i < numGaussPoints.size(); i++) {
                    computeMassMatrixOfTrianlge(triangle, numGaussPoints[i], true, massMatrix);
                    //printGeneralMatrix(massMatrix, 3, 3);

                    CPPUNIT_ASSERT(massMatrix[0] - massMatrixCorrect[0] < EPS);
                    CPPUNIT_ASSERT(massMatrix[1] - massMatrixCorrect[1] < EPS);
                    CPPUNIT_ASSERT(massMatrix[2] - massMatrixCorrect[2] < EPS);

                    CPPUNIT_ASSERT(massMatrix[3] - massMatrixCorrect[3] < EPS);
                    CPPUNIT_ASSERT(massMatrix[4] - massMatrixCorrect[4] < EPS);
                    CPPUNIT_ASSERT(massMatrix[5] - massMatrixCorrect[5] < EPS);

                    CPPUNIT_ASSERT(massMatrix[6] - massMatrixCorrect[6] < EPS);
                    CPPUNIT_ASSERT(massMatrix[7] - massMatrixCorrect[7] < EPS);
                    CPPUNIT_ASSERT(massMatrix[8] - massMatrixCorrect[8] < EPS);
                }
            }
        }
    }

    void testMassMatrixOfQuad() {
        double quad[] = { 0, 0, 0, 1, 0, 0, 2, 2, 0, 0, 1, 0 };

        vector<int> numGaussPoints;
        //numGaussPoints.push_back(1); // not correct in this case
        numGaussPoints.push_back(4);
        numGaussPoints.push_back(9);
        double massMatrix[4 * 4];
        const double EPS = 1E-10;
        // 1. test on normal mass matrix
        for (int i = 0; i < numGaussPoints.size(); i++) {
            computeMassMatrixOfQuad(quad, numGaussPoints[i], false, massMatrix);
            //printGeneralMatrix(massMatrix, 4, 4);

            // testing with the mass matrix of the triangle computed by hand.
            // not by hand actually, by computing the numbers and use
            // http://www.mindspring.com/~alanh/fracs.html
            CPPUNIT_ASSERT(massMatrix[0] - 1.0/6.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[1] - 7.0/72.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[2] - 1.0/18.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[3] - 7.0/72.0 < EPS);

            CPPUNIT_ASSERT(massMatrix[4] - 7.0/72.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[5] - 2.0/9.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[6] - 1.0/8.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[7] - 1.0/18.0 < EPS);

            CPPUNIT_ASSERT(massMatrix[8] - 1.0/18.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[9] - 1.0/8.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[10] - 5.0/18.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[11] - 1.0/8.0 < EPS);

            CPPUNIT_ASSERT(massMatrix[12] - 7.0/72.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[13] - 1.0/18.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[14] - 1.0/8.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[15] - 2.0/9.0 < EPS);
        }
        // 2. test on dual mass matrix
        for (int i = 0; i < numGaussPoints.size(); i++) {
            computeMassMatrixOfQuad(quad, numGaussPoints[i], true, massMatrix);
            //printGeneralMatrix(massMatrix, 4, 4);

            // testing with the mass matrix of the triangle computed by hand.
            // not by hand actually, by computing the numbers and use
            // http://www.mindspring.com/~alanh/fracs.html
            CPPUNIT_ASSERT(massMatrix[0] - 5.0/12.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[1] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[2] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[3] - 0.0 < EPS);

            CPPUNIT_ASSERT(massMatrix[4] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[5] - 0.5 < EPS);
            CPPUNIT_ASSERT(massMatrix[6] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[7] - 0.0 < EPS);

            CPPUNIT_ASSERT(massMatrix[8] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[9] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[10] - 7.0/12.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[11] - 0.0 < EPS);

            CPPUNIT_ASSERT(massMatrix[12] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[13] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[14] - 0.0 < EPS);
            CPPUNIT_ASSERT(massMatrix[15] - 0.5 < EPS);
        }
    }

    void testIntersect() {
        {
            double P1[] = { 0, 0, 0 };
            double P2[] = { 1, 0, 0 };
            double P3[] = { 1, 1, 0 };
            double P4[] = { 0, 1, 0 };

            int planeToProject = 2;
            {
                double intersection[3];
                bool hasIntersection = PolygonClipper::intersect(P1, P2, P3, P4, planeToProject,
                        intersection);
                CPPUNIT_ASSERT(hasIntersection == false);
            }
            {
                double intersection[3];
                bool hasIntersection = PolygonClipper::intersect(P1, P3, P2, P4, planeToProject,
                        intersection);
                CPPUNIT_ASSERT(hasIntersection == true);
                const double EPS = 1E-10;
                CPPUNIT_ASSERT(fabs(intersection[0]-0.5)<EPS);
                CPPUNIT_ASSERT(fabs(intersection[1]-0.5)<EPS);
                CPPUNIT_ASSERT(fabs(intersection[2]-0.0)<EPS);
            }
        }
        { // used to cause problem (solve2x2LinearSystem fails because det0==det1==detA==0)
            double P1[] = { 2.9999, 0.999703, -0.334495 };
            double P2[] = { 3.0001, 1.37398, -0.00033278 };
            double P3[] = { 3.0001, 1.37398, -0.00033278 };
            double P4[] = { 2.9999, 0.999703, -0.334495 };

            int planeToProject = 2;
            {
                double intersection[3];
                bool hasIntersection = PolygonClipper::intersect(P1, P2, P3, P4, planeToProject,
                        intersection);
                //MortarMath::printPoint(intersection);
                CPPUNIT_ASSERT(hasIntersection == false);
            }
        }
    }

    void testComputePolygonCenter() {
        double quad[] = { 0, 0, 1000000, 1, 0, 1000000, 1, 1, 1000000, 0, 1, 1000000 };
        double center[3];
        computePolygonCenter(quad, 4, center);
        const double EPS = 1E-10;
        CPPUNIT_ASSERT(fabs(center[0]-0.5)<EPS);
        CPPUNIT_ASSERT(fabs(center[1]-0.5)<EPS);
        CPPUNIT_ASSERT(fabs(center[2]-1000000.0)<EPS);
    }

    void testComputeElemNormal() {
        {
            double triangle[] = { 0, 0, 1000000, 1, 0, 1000000, 100, 1, 1000000 };
            double normal[3];
            computeNormalOfTriangle(triangle, true, normal);
            const double EPS = 1E-10;
            CPPUNIT_ASSERT(fabs(normal[0]-0.0)<EPS);
            CPPUNIT_ASSERT(fabs(normal[1]-0.0)<EPS);
            CPPUNIT_ASSERT(fabs(normal[2]-1.0)<EPS);
        }
        {
            double quad[] = { 0, 0, 1000000, 1, 0, 1000000, 100, 100, 1000000, 0, 1, 1000000 };
            double normal[3];
            computeNormalOfQuad(quad, true, normal);
            const double EPS = 1E-10;
            CPPUNIT_ASSERT(fabs(normal[0]-0.0)<EPS);
            CPPUNIT_ASSERT(fabs(normal[1]-0.0)<EPS);
            CPPUNIT_ASSERT(fabs(normal[2]-1.0)<EPS);
        }
    }

    void testLinearAlgebra() {
        { // test solve2x2LinearSystem
            double A[] = { 1, 2, 2, 1 };
            double b[] = { 3, 3 };
            solve2x2LinearSystem(A, b, 1e-15);
            const double EPS = 1E-10;
            CPPUNIT_ASSERT(fabs(b[0]-1.0)<EPS);
            CPPUNIT_ASSERT(fabs(b[1]-1.0)<EPS);
        }
        { // test solve3x3LinearSystem
            double A[] = { 1, 1, 1, 1, 0, 0, 0, 1, 0 };
            double A_T[9];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    A_T[i * 3 + j] = A[j * 3 + i];
                }
            }
            double b[] = { 1.0, 0.3, 0.3 };
            solve3x3LinearSystem(A_T, 0, b);
            const double EPS = 1E-10;
            CPPUNIT_ASSERT(fabs(b[0]-0.3)<EPS);
            CPPUNIT_ASSERT(fabs(b[1]-0.3)<EPS);
            CPPUNIT_ASSERT(fabs(b[2]-0.4)<EPS);
        }
        { // test computeMatrixProduct
            double A[] = { 1, 1, 1, 2, 2, 2, 3, 3, 3 };
            double B[] = { 4, 5, 6, 4, 5, 6, 4, 5, 6 };
            computeMatrixProduct(3, 3, A, B);
            const double EPS = 1E-10;
            CPPUNIT_ASSERT(fabs(B[0]-12)<EPS);
            CPPUNIT_ASSERT(fabs(B[1]-15)<EPS);
            CPPUNIT_ASSERT(fabs(B[2]-18)<EPS);
            CPPUNIT_ASSERT(fabs(B[3]-24)<EPS);
            CPPUNIT_ASSERT(fabs(B[4]-30)<EPS);
            CPPUNIT_ASSERT(fabs(B[5]-36)<EPS);
            CPPUNIT_ASSERT(fabs(B[6]-36)<EPS);
            CPPUNIT_ASSERT(fabs(B[7]-45)<EPS);
            CPPUNIT_ASSERT(fabs(B[8]-54)<EPS);
        }
        { // test dcsrmv
          // get the matrix from PARDISO manual
            int numRows = 8;
            int numCols = 8;
            double A[20] = { 7, 1, 2, 7, -4, 8, 2, 1, 5, 7, 9, -4, 7, 3, 5, 17, 11, -3, 2, 5 };
            int JA[20] = { 1, 3, 6, 7, 2, 3, 5, 3, 8, 4, 7, 2, 3, 6, 8, 2, 7, 3, 7, 8 };
            int IA[9] = { 1, 5, 8, 10, 12, 13, 16, 18, 21 };
            double x[] = { 1, 1, 1, 1, 1, 1, 1, 1 };
            double y[numRows];
            {
                char trans = 'N';
                dcsrmv(trans, numRows, numCols, A, JA, IA, x, y);
                const double EPS = 1E-10;
                CPPUNIT_ASSERT(fabs(y[0]-17)<EPS);
                CPPUNIT_ASSERT(fabs(y[1]-6)<EPS);
                CPPUNIT_ASSERT(fabs(y[2]-6)<EPS);
                CPPUNIT_ASSERT(fabs(y[3]-16)<EPS);
                CPPUNIT_ASSERT(fabs(y[4]-(-4))<EPS);
                CPPUNIT_ASSERT(fabs(y[5]-15)<EPS);
                CPPUNIT_ASSERT(fabs(y[6]-28)<EPS);
                CPPUNIT_ASSERT(fabs(y[7]-4)<EPS);
            }
            {
                char trans = 'T';
                dcsrmv(trans, numRows, numCols, A, JA, IA, x, y);
                const double EPS = 1E-10;
                CPPUNIT_ASSERT(fabs(y[0]-7)<EPS);
                CPPUNIT_ASSERT(fabs(y[1]-9)<EPS);
                CPPUNIT_ASSERT(fabs(y[2]-14)<EPS);
                CPPUNIT_ASSERT(fabs(y[3]-7)<EPS);
                CPPUNIT_ASSERT(fabs(y[4]-2)<EPS);
                CPPUNIT_ASSERT(fabs(y[5]-5)<EPS);
                CPPUNIT_ASSERT(fabs(y[6]-29)<EPS);
                CPPUNIT_ASSERT(fabs(y[7]-15)<EPS);
            }
        }
        { // test dcsrsymv
          // get the matrix from PARDISO manual
            int n = 8;
            double A[18] = { 7, 1, 2, 7, -4, 8, 2, 1, 5, 7, 9, 5, -1, 5, 0, 5, 11, 5 };
            int JA[20] = { 1, 3, 6, 7, 2, 3, 5, 3, 8, 4, 7, 5, 6, 7, 6, 8, 7, 8 };
            int IA[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
            double x[] = { 1, 1, 1, 1, 1, 1, 1, 1 };
            double y[n];

            dcsrsymv(n, A, IA, JA, x, y);
            const double EPS = 1E-10;
            CPPUNIT_ASSERT(fabs(y[0]-17)<EPS);
            CPPUNIT_ASSERT(fabs(y[1]-6)<EPS);
            CPPUNIT_ASSERT(fabs(y[2]-15)<EPS);
            CPPUNIT_ASSERT(fabs(y[3]-16)<EPS);
            CPPUNIT_ASSERT(fabs(y[4]-11)<EPS);
            CPPUNIT_ASSERT(fabs(y[5]-6)<EPS);
            CPPUNIT_ASSERT(fabs(y[6]-32)<EPS);
            CPPUNIT_ASSERT(fabs(y[7]-15)<EPS);
        }
    }

    void testClipping() {
        const double EPS = 1e-10;
        { // overlap
            double window[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0 };
            double triangle[] = { 1, 1, 0, 3, 1, 0, 1, 3, 0 };
            PolygonClipper clipper(window, 4, 2);
            vector<double*> polygonResult;
            clipper.clip(triangle, 3, &polygonResult);
            double tmp[polygonResult.size() * 3];
            for (int i = 0; i < polygonResult.size(); i++)
                for (int j = 0; j < 3; j++)
                    tmp[i * 3 + j] = polygonResult[i][j];
            for (int i = 0; i < polygonResult.size(); i++)
                delete[] polygonResult[i];
            //MortarMath::printElem(tmp, polygonResult.size());
            CPPUNIT_ASSERT(fabs(MortarMath::computePolygonArea(tmp, polygonResult.size())-1.0)<EPS);
        }
        { // window is inside triangle
            double window[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0 };
            double triangle[] = { -10, -10, 0, 10, -10, 0, 0, 50, 0 };
            PolygonClipper clipper(window, 4, 2);
            vector<double*> polygonResult;
            clipper.clip(triangle, 3, &polygonResult);
            double tmp[polygonResult.size() * 3];
            for (int i = 0; i < polygonResult.size(); i++)
                for (int j = 0; j < 3; j++)
                    tmp[i * 3 + j] = polygonResult[i][j];
            for (int i = 0; i < polygonResult.size(); i++)
                delete[] polygonResult[i];
            //MortarMath::printElem(tmp, polygonResult.size());
            CPPUNIT_ASSERT(fabs(MortarMath::computePolygonArea(tmp, polygonResult.size())-4.0)<EPS);
        }
        { // window and triangle have overlapped points
            double window[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0 };
            double triangle[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0 };
            PolygonClipper clipper(window, 4, 2);
            vector<double*> polygonResult;
            clipper.clip(triangle, 3, &polygonResult);
            double tmp[polygonResult.size() * 3];
            for (int i = 0; i < polygonResult.size(); i++)
                for (int j = 0; j < 3; j++)
                    tmp[i * 3 + j] = polygonResult[i][j];
            for (int i = 0; i < polygonResult.size(); i++)
                delete[] polygonResult[i];
            //MortarMath::printElem(tmp, polygonResult.size());
            CPPUNIT_ASSERT(fabs(MortarMath::computePolygonArea(tmp, polygonResult.size())-2.0)<EPS);
        }
        { // window and triangle have overlapped points, but no clipped area
            double window[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0 };
            double triangle[] = { 0, 0, 0, 2, 0, 0, 2, -2, 0 };
            PolygonClipper clipper(window, 4, 2);
            vector<double*> polygonResult;
            clipper.clip(triangle, 3, &polygonResult);
            double tmp[polygonResult.size() * 3];
            for (int i = 0; i < polygonResult.size(); i++)
                for (int j = 0; j < 3; j++)
                    tmp[i * 3 + j] = polygonResult[i][j];
            for (int i = 0; i < polygonResult.size(); i++)
                delete[] polygonResult[i];
            //MortarMath::printElem(tmp, polygonResult.size());
            CPPUNIT_ASSERT(fabs(MortarMath::computePolygonArea(tmp, polygonResult.size())-0.0)<EPS);
        }
        { // window and triangle have overlapped points, with numerical error
            double EPS2 = 1e-5;
            double disturb = 1e-5;
            streamsize pcs = cout.precision();
            cout.precision(20);
            for (int i = 0; i < 10; i++) {
                double window[] = { 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0 };
                double triangle[] = { 0, 0, 0, 2.0 - disturb, 0, 0, 2.0 + disturb, 2, 0 };
                PolygonClipper clipper(window, 4, 2);
                vector<double*> polygonResult;
                clipper.clip(triangle, 3, &polygonResult);
                double tmp[polygonResult.size() * 3];
                for (int i = 0; i < polygonResult.size(); i++)
                    for (int j = 0; j < 3; j++)
                        tmp[i * 3 + j] = polygonResult[i][j];
                for (int i = 0; i < polygonResult.size(); i++)
                    delete[] polygonResult[i];
                //MortarMath::printElem(tmp, polygonResult.size());
                CPPUNIT_ASSERT(
                        fabs(MortarMath::computePolygonArea(tmp, polygonResult.size())-(2.0-disturb))<EPS2);
                disturb /= 10.0;
                EPS2 /= 10.0;
            }
            cout.precision(pcs);
        }
        { // test case which used to cause problem
            double window[] = { 2.0 / 30.0, 0, 0, 0.1, 0, 0, 0.1, 0, 0.1, 2.0 / 30.0, 0, 0.1 };
            double triangle[] = { 0.1, 0, 0.1, 0.11, 0, 0, 0.11, 0, 0.1 };
            PolygonClipper clipper(window, 4, 1);
            vector<double*> polygonResult;
            clipper.clip(triangle, 3, &polygonResult);
            double tmp[polygonResult.size() * 3];
            for (int i = 0; i < polygonResult.size(); i++)
                for (int j = 0; j < 3; j++)
                    tmp[i * 3 + j] = polygonResult[i][j];
            for (int i = 0; i < polygonResult.size(); i++)
                delete[] polygonResult[i];
            //MortarMath::printElem(tmp, polygonResult.size());
            CPPUNIT_ASSERT(fabs(MortarMath::computePolygonArea(tmp, polygonResult.size())-0.0)<EPS);
        }
        /*{ // test memory leak, comment it except when checking memory leak
         for (int i = 0; i < 10000000; i++) {
         double window[] = { 2.0 / 30.0, 0, 0, 0.1, 0, 0, 0.1, 0, 0.1, 2.0 / 30.0, 0, 0.1 };
         double triangle[] = { 0.1, 0, 0.1, 0.11, 0, 0, 0.11, 0, 0.1 };
         PolygonClipper clipper(window, 4, 1);
         vector<double*> polygonResult;
         clipper.clip(triangle, 3, &polygonResult);
         double tmp[polygonResult.size() * 3];
         for (int i = 0; i < polygonResult.size(); i++)
         for (int j = 0; j < 3; j++)
         tmp[i * 3 + j] = polygonResult[i][j];
         for (int i = 0; i < polygonResult.size(); i++)
         delete[] polygonResult[i];
         }
         }*/
    }
    void testComputeLocalCoorInQuad() { // data come from some big test
        double quad[12] = { 0.587856, 0.249597, 0.768143, 0.587862, 0.65312, 0.474549, 0.866025,
                0.404234, 0.293715, 0.866025, 0.154482, 0.475424 };
        double point[3] = { 1.24731, -0.366352, 0.358272 };
        double localCoor[2];
        double normal[3];
        computeNormalOfQuad(quad, false, normal);
        int planeToProject = MortarMath::computePlaneToProject(normal);
        MortarMath::computeLocalCoorInQuad(quad, planeToProject, point, localCoor);
        printElem(quad, 4);
        printPoint(point);
        cout << "plane to project:  " << planeToProject << endl;
        cout << "xi: " << localCoor[0] << "  ita:  " << localCoor[1] << endl;

        MortarMath::computeGlobalCoorInQuad(quad, localCoor, point);
        MortarMath::printPoint(point);
    }

CPPUNIT_TEST_SUITE( TestMortarMath );
        CPPUNIT_TEST(testGaussQuadratureOnTriangle);
        CPPUNIT_TEST(testGaussQuadratureOnQuad);
        CPPUNIT_TEST(testMassMatrixOfTriangle);
        CPPUNIT_TEST(testMassMatrixOfQuad);
        CPPUNIT_TEST(testIntersect);
        CPPUNIT_TEST(testComputePolygonCenter);
        CPPUNIT_TEST(testComputeElemNormal);
        CPPUNIT_TEST(testLinearAlgebra);
        CPPUNIT_TEST(testClipping);
        //CPPUNIT_TEST(testComputeLocalCoorInQuad);
    CPPUNIT_TEST_SUITE_END();
};

} /* namespace EMPIRE */

CPPUNIT_TEST_SUITE_REGISTRATION( EMPIRE::TestMortarMath);
