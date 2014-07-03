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
#include "MortarMath.h"
#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <list>

using namespace std;

namespace EMPIRE {
namespace MortarMath {

const double triGaussPoints3[9] = { 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0
        / 6.0, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0 };
const double triWeights3[3] = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };

const double triGaussPoints6[18] = { 0.816847572980459, 0.091576213509771, 0.091576213509771,
        0.091576213509771, 0.816847572980459, 0.09157621350977, 0.09157621350977, 0.09157621350977,
        0.816847572980459, 0.108103018168070, 0.445948490915965, 0.445948490915965,
        0.445948490915965, 0.108103018168070, 0.445948490915965, 0.445948490915965,
        0.445948490915965, 0.108103018168070 };
const double triWeights6[6] = { 0.109951743655322, 0.109951743655322, 0.109951743655322,
        0.223381589678011, 0.223381589678011, 0.223381589678011 };

const double triGaussPoints7[21] = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.797426985353087,
        0.101286507323456, 0.101286507323456, 0.101286507323456, 0.797426985353087,
        0.101286507323456, 0.101286507323456, 0.101286507323456, 0.797426985353087,
        0.059715871789770, 0.470142064105115, 0.470142064105115, 0.470142064105115,
        0.059715871789770, 0.470142064105115, 0.470142064105115, 0.470142064105115,
        0.059715871789770 };
const double triWeights7[7] = { 0.225000000000000, 0.125939180544827, 0.125939180544827,
        0.125939180544827, 0.132394152788506, 0.132394152788506, 0.132394152788506 };

const double triGaussPoints12[36] = { 0.873821971016996, 0.063089014491502, 0.063089014491502,
        0.063089014491502, 0.873821971016996, 0.063089014491502, 0.063089014491502,
        0.063089014491502, 0.873821971016996, 0.501426509658179, 0.249286745170910,
        0.249286745170910, 0.249286745170910, 0.501426509658179, 0.249286745170910,
        0.249286745170910, 0.249286745170910, 0.501426509658179, 0.636502499121399,
        0.310352451033785, 0.053145049844816, 0.310352451033785, 0.053145049844816,
        0.636502499121399, 0.053145049844816, 0.636502499121399, 0.310352451033785,
        0.636502499121399, 0.053145049844816, 0.310352451033785, 0.310352451033785,
        0.636502499121399, 0.053145049844816, 0.053145049844816, 0.310352451033785,
        0.636502499121399 };
const double triWeights12[12] = { 0.050844906370207, 0.050844906370207, 0.050844906370207,
        0.116786275726379, 0.116786275726379, 0.116786275726379, 0.082851075618374,
        0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374,
        0.082851075618374 };

const double quadGaussPoints1[2] = { 0.0, 0.0 };
const double quadWeights1[1] = { 4.0 };

static const double tmpSqrt13 = sqrt(1.0 / 3.0);
const double quadGaussPoints4[8] = { tmpSqrt13, tmpSqrt13, -tmpSqrt13, tmpSqrt13, -tmpSqrt13,
        -tmpSqrt13, tmpSqrt13, -tmpSqrt13 };
const double quadWeights4[4] = { 1.0, 1.0, 1.0, 1.0 };

static const double tmpSqrt35 = sqrt(3.0 / 5.0);
const double quadGaussPoints9[18] = { 0.0, 0.0, -tmpSqrt35, 0, 0, -tmpSqrt35, tmpSqrt35, 0, 0,
        tmpSqrt35, -tmpSqrt35, -tmpSqrt35, tmpSqrt35, -tmpSqrt35, -tmpSqrt35, tmpSqrt35, tmpSqrt35,
        tmpSqrt35 };
const double quadWeights9[9] = { 64.0 / 81.0, 40.0 / 81.0, 40.0 / 81.0, 40.0 / 81.0, 40.0 / 81.0,
        25.0 / 81.0, 25.0 / 81.0, 25.0 / 81.0, 25.0 / 81.0 };

GaussQuadratureOnTriangle::GaussQuadratureOnTriangle(double *_triangle, int _numGaussPoints) :
        triangle(_triangle), numGaussPoints(_numGaussPoints) {
    if (numGaussPoints == 3) {
        gaussPointsLocal = triGaussPoints3;
        weights = triWeights3;
    } else if (numGaussPoints == 6) {
        gaussPointsLocal = triGaussPoints6;
        weights = triWeights6;
    } else if (numGaussPoints == 7) {
        gaussPointsLocal = triGaussPoints7;
        weights = triWeights7;
    } else if (numGaussPoints == 12) {
        gaussPointsLocal = triGaussPoints12;
        weights = triWeights12;
    } else {
        assert(false);
    }
    gaussPointsGlobal = new double[numGaussPoints * 3];
    for (int i = 0; i < numGaussPoints; i++) {
        computeGlobalCoorInTriangle(triangle, &gaussPointsLocal[i * 3], &gaussPointsGlobal[i * 3]);
    }
    area = computeAreaOfTriangle(triangle);
}

GaussQuadratureOnTriangle::~GaussQuadratureOnTriangle() {
    delete[] gaussPointsGlobal;
}

void GaussQuadratureOnTriangle::setIntegrandFunc(IntegrandFunction *_integrandFunc) {
    integrandFunc = _integrandFunc;
}

double GaussQuadratureOnTriangle::computeIntegral() {
    double toReturn = 0;
    for (int i = 0; i < numGaussPoints; i++)
        toReturn += weights[i] * (*integrandFunc)(&gaussPointsGlobal[i * 3]);
    toReturn *= area;
    return toReturn;
}

GaussQuadratureOnQuad::GaussQuadratureOnQuad(double *_quad, int _numGaussPoints) :
        quad(_quad), numGaussPoints(_numGaussPoints) {
    if (numGaussPoints == 1) {
        gaussPointsLocal = quadGaussPoints1;
        weights = quadWeights1;
    } else if (numGaussPoints == 4) {
        gaussPointsLocal = quadGaussPoints4;
        weights = quadWeights4;
    } else if (numGaussPoints == 9) {
        gaussPointsLocal = quadGaussPoints9;
        weights = quadWeights9;
    } else {
        assert(false);
    }
    gaussPointsGlobal = new double[numGaussPoints * 3];
    detJ = new double[numGaussPoints];
    for (int i = 0; i < numGaussPoints; i++) {
        computeGlobalCoorInQuad(quad, &gaussPointsLocal[i * 2], &gaussPointsGlobal[i * 3]);
    }
    for (int i = 0; i < numGaussPoints; i++) {
        detJ[i] = computeDetJOfQuad(quad, &gaussPointsLocal[i * 2]);
    }
}

GaussQuadratureOnQuad::~GaussQuadratureOnQuad() {
    delete[] gaussPointsGlobal;
    delete[] detJ;
}

void GaussQuadratureOnQuad::setIntegrandFunc(IntegrandFunction *_integrandFunc) {
    integrandFunc = _integrandFunc;
}

double GaussQuadratureOnQuad::computeIntegral() {
    double toReturn = 0;
    for (int i = 0; i < numGaussPoints; i++)
        toReturn += detJ[i] * weights[i] * (*integrandFunc)(&gaussPointsGlobal[i * 3]);
    return toReturn;
}

void computeMassMatrixOfTrianlge(const double *triangle, int numGaussPoints, bool dual,
        double *massMatrix) {
    const double *gaussPointsLocal;
    const double *weights;
    if (numGaussPoints == 3) {
        gaussPointsLocal = triGaussPoints3;
        weights = triWeights3;
    } else if (numGaussPoints == 6) {
        gaussPointsLocal = triGaussPoints6;
        weights = triWeights6;
    } else if (numGaussPoints == 7) {
        gaussPointsLocal = triGaussPoints7;
        weights = triWeights7;
    } else if (numGaussPoints == 12) {
        gaussPointsLocal = triGaussPoints12;
        weights = triWeights12;
    } else {
        assert(false);
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            massMatrix[i * 3 + j] = 0.0;
        }
    }
    double area = computeAreaOfTriangle(triangle);
    if (!dual) {
        // since the mass matrix is symmetric, only calculate the upper part
        for (int i = 0; i < 3; i++) {
            for (int j = i; j < 3; j++) {
                for (int k = 0; k < numGaussPoints; k++) {
                    massMatrix[i * 3 + j] += weights[k] * gaussPointsLocal[k * 3 + i]
                            * gaussPointsLocal[k * 3 + j];
                }
                massMatrix[i * 3 + j] *= area;
            }
        }
        // set the lower part
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < i; j++) {
                massMatrix[i * 3 + j] = massMatrix[j * 3 + i];
            }
        }
    } else {
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < numGaussPoints; k++) {
                massMatrix[i * 3 + i] += weights[k] * gaussPointsLocal[k * 3 + i];
            }
            massMatrix[i * 3 + i] *= area;
        }
    }
}

void computeMassMatrixOfQuad(const double *quad, int numGaussPoints, bool dual,
        double *massMatrix) {
    const double *gaussPointsLocal;
    const double *weights;
    if (numGaussPoints == 1) {
        gaussPointsLocal = quadGaussPoints1;
        weights = quadWeights1;
    } else if (numGaussPoints == 4) {
        gaussPointsLocal = quadGaussPoints4;
        weights = quadWeights4;
    } else if (numGaussPoints == 9) {
        gaussPointsLocal = quadGaussPoints9;
        weights = quadWeights9;
    } else {
        assert(false);
    }
    double GPShapeFunc[numGaussPoints * 4];
    for (int i = 0; i < numGaussPoints; i++) {
        for (int j = 0; j < 4; j++) {
            computeShapeFuncOfQuad(&gaussPointsLocal[i * 2], &GPShapeFunc[i * 4]);
        }
    }
    double detJ[numGaussPoints];
    for (int i = 0; i < numGaussPoints; i++) {
        detJ[i] = computeDetJOfQuad(quad, &gaussPointsLocal[i * 2]);
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            massMatrix[i * 4 + j] = 0.0;
        }
    }
    if (!dual) {
        // since the mass matrix is symmetric, only calculate the upper part
        for (int i = 0; i < 4; i++) {
            for (int j = i; j < 4; j++) {
                for (int k = 0; k < numGaussPoints; k++) {
                    massMatrix[i * 4 + j] += weights[k] * detJ[k] * GPShapeFunc[k * 4 + i]
                            * GPShapeFunc[k * 4 + j];
                }
            }
        }
        // set the lower part
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < i; j++) {
                massMatrix[i * 4 + j] = massMatrix[j * 4 + i];
            }
        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < numGaussPoints; k++) {
                massMatrix[i * 4 + i] += weights[k] * detJ[k] * GPShapeFunc[k * 4 + i];
            }
        }
    }
}

void computeShapeFuncOfQuad(const double *xi_eta, double *shapeFuncValues) {
    shapeFuncValues[0] = 0.25 * (1.0 - xi_eta[0]) * (1.0 - xi_eta[1]);
    shapeFuncValues[1] = 0.25 * (1.0 + xi_eta[0]) * (1.0 - xi_eta[1]);
    shapeFuncValues[2] = 0.25 * (1.0 + xi_eta[0]) * (1.0 + xi_eta[1]);
    shapeFuncValues[3] = 0.25 * (1.0 - xi_eta[0]) * (1.0 + xi_eta[1]);
}

double computeDetJOfQuad(const double *quad, const double *xi_eta) {
    // d_N_d_xi[4] contains the partial derivative w.r.t. xi of the shape functions
    double d_N_d_xi[4];
    // d_N_d_eta[4] contains the partial derivative w.r.t. eta of the shape functions
    double d_N_d_eta[4];

    d_N_d_xi[0] = -0.25 * (1 - xi_eta[1]);
    d_N_d_xi[1] = -d_N_d_xi[0];
    d_N_d_xi[2] = 0.25 * (1 + xi_eta[1]);
    d_N_d_xi[3] = -d_N_d_xi[2];

    d_N_d_eta[0] = -0.25 * (1 - xi_eta[0]);
    d_N_d_eta[1] = -0.25 * (1 + xi_eta[0]);
    d_N_d_eta[2] = -d_N_d_eta[1];
    d_N_d_eta[3] = -d_N_d_eta[0];

    // g1 and g2 are the local basis vectors, and det(J)=||g1 x g2||
    double g1[3] = { 0, 0, 0 };
    double g2[3] = { 0, 0, 0 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            g1[i] += quad[j * 3 + i] * d_N_d_xi[j];
            g2[i] += quad[j * 3 + i] * d_N_d_eta[j];
        }
    }

    double crossProduct[3];
    computeVectorCrossProduct(g1, g2, crossProduct);
    return computeVectorLength(crossProduct);
}

void computeGlobalCoorInTriangle(const double *triangle, const double *localCoor,
        double *globalCoor) {
    for (int i = 0; i < 3; i++) {
        globalCoor[i] = 0;
        for (int j = 0; j < 3; j++)
            globalCoor[i] += localCoor[j] * triangle[i + j * 3];
    }
}

void computeGlobalCoorInQuad(const double *quad, const double *localCoor, double *globalCoor) {
    double shapeFuncValues[4];
    computeShapeFuncOfQuad(localCoor, shapeFuncValues);
    for (int i = 0; i < 3; i++)
        globalCoor[i] = 0.0;

    for (int i = 0; i < 3; i++) { // x, y, z of globalCoor
        for (int j = 0; j < 4; j++) { // 4 end nodes of quad
            globalCoor[i] += shapeFuncValues[j] * quad[j * 3 + i];
        }
    }
}

bool computeLocalCoorInTriangle(const double *triangle, int planeToProject, const double *point,
        double *localCoor) {
    /* Normally, the local coordinates can be solved by:
     *  |x_1 x_2 x_3|   |xi_1|    |x_4|
     *  |y_1 y_2 y_3| * |xi_2| =  |y_4|
     *  |z_1 z_2 z_3|   |xi_3|    |z_4|
     *
     * But if x_1 = x_2 = x_3 = 0 than the system cannot be solved.
     * So we remove one equation by xi_1 + xi_2 + xi_3 = 1.
     * This indicates projection.
     * Choose among planes (xy or yz or zx) the one which has the smallest angle
     * with the triangle normal.
     * For example, if the angle between of normals of yz-plane and the triangle
     * is 0, then the system is replaced by
     *  |1.0 1.0 1.0|   |xi_1|    |1.0|
     *  |y_1 y_2 y_3| * |xi_2| =  |y_4|
     *  |z_1 z_2 z_3|   |xi_3|    |z_4|
     */
    double A[9]; // Attention! A is the transpose of the above matrix
    for (int i = 0; i < 9; i++)
        A[i] = triangle[i];

    for (int i = 0; i < 3; i++)
        localCoor[i] = point[i];

    A[planeToProject] = A[planeToProject + 3] = A[planeToProject + 6] = localCoor[planeToProject] =
            1.0;

    solve3x3LinearSystem(A, planeToProject, localCoor);

    // make sure the sum is 1.0
    assert(fabs(localCoor[0] + localCoor[1] + localCoor[2] -1.0) < 1E-12);
    localCoor[0] = 1.0 - localCoor[1] - localCoor[2];

    for (int i = 0; i < 3; i++) {
        if (localCoor[i] > 1.0)
        	return false;
        if (localCoor[i] < 0.0)
            return false;
    }
    return true;
}

bool computeLocalCoorInQuad(const double *quad, int planeToProject, const double *point,
        double *localCoor) {
    /*
     * So we use two coordinates among x, y, z.
     * This indicates projection.
     * Choose among planes (xy or yz or zx) the one which has the smallest angle
     * with the quad normal.
     */
    double x[4];
    double y[4];
    int x_direc = (planeToProject + 1) % 3;
    int y_direc = (planeToProject + 2) % 3;
    for (int i = 0; i < 4; i++) {
        x[i] = quad[i * 3 + x_direc];
        y[i] = quad[i * 3 + y_direc];
    }
    double x0 = point[x_direc];
    double y0 = point[y_direc];

    double a1 = x[0] + x[1] + x[2] + x[3] - 4.0 * x0;
    double b1 = -x[0] + x[1] + x[2] - x[3];
    double c1 = -x[0] - x[1] + x[2] + x[3];
    double d1 = x[0] - x[1] + x[2] - x[3];

    double a2 = y[0] + y[1] + y[2] + y[3] - 4.0 * y0;
    double b2 = -y[0] + y[1] + y[2] - y[3];
    double c2 = -y[0] - y[1] + y[2] + y[3];
    double d2 = y[0] - y[1] + y[2] - y[3];

    double delta[2];
    double J_T[4]; // transpose of Jacobian --- to use column major in lapack
    double F[2]; // -F
    localCoor[0] = 0;
    localCoor[1] = 0;
    //int dummy[2];
    const double EPS = 1E-13; // should be smaller than 1E-15 from experience

    const int MAX_ITER_NUM = 100;
    for (int i = 0; i < MAX_ITER_NUM; i++) {
        J_T[0] = b1 + d1 * localCoor[1];
        J_T[2] = c1 + d1 * localCoor[0];
        J_T[1] = b2 + d2 * localCoor[1];
        J_T[3] = c2 + d2 * localCoor[0];
        F[0] = a1 + b1 * localCoor[0] + c1 * localCoor[1] + d1 * localCoor[0] * localCoor[1];
        F[1] = a2 + b2 * localCoor[0] + c2 * localCoor[1] + d2 * localCoor[0] * localCoor[1];
        delta[0] = -F[0];
        delta[1] = -F[1];
        /*int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 2, 1, J_T, 2, dummy,
         delta, 2);
         if (info != 0) {
         cerr << "ERROR in MortarMath::computeLocalCoorInQuad()!" << endl;
         exit(EXIT_FAILURE);
         }*/
        solve2x2LinearSystem(J_T, delta, EPS);
        if (i >= 10) { // normally should find a solution within 10 iterations
            cout<<"WARNING(MortarMath::computeLocalCoorInQuad): More than 10 iterations are necessary for computing local coordinates in quad"<<endl;
            cout << "iteration #: " << i << endl;
            printElem(quad, 4);
            printPoint(point);
            cout << "plane to project:  " << planeToProject << endl;
            cout << "xi:  " << localCoor[0] << " ita: " << localCoor[1] << endl;
            cout << "delta-xi:  " << delta[0] << " delta-ita: " << delta[1] << endl;
            // if point is far out of quad, return false
            for (int i = 0; i < 2; i++) { // do not care accuracy if point is far outside the quad
                if (localCoor[i] > 2.0)
                    return false;
                if (localCoor[i] < -2.0)
                    return false;
            }
            //assert(false);
        }
        if (fabs(delta[0]) < EPS && fabs(delta[1]) < EPS) {
            break;
        }

        localCoor[0] += delta[0];
        localCoor[1] += delta[1];
    }
    for (int i = 0; i < 2; i++) {
        if (localCoor[i] > 1.0)
            return false;
        if (localCoor[i] < -1.0)
            return false;
    }
    return true;
}

const double PolygonClipper::EPS1 = 1e-10;
const double PolygonClipper::EPS2 = 1e-10;

PolygonClipper::PolygonClipper(const double *_polygonWindow, int _sizePolygonWindow,
        int _planeToProject) :
        planeToProject(_planeToProject) {
    assert(_sizePolygonWindow>2);
    // Michael Andre begin: correcting collapsed nodes so clipping algorithm works properly
    const double eps = 1e-10;
    double lastNoSkip[3];
    int iLastNoSkip;
    bool *skip = new bool[_sizePolygonWindow];
    int numSkip;
    int i, ii;

    sizePolygonWindow = _sizePolygonWindow;
    skip[0] = false;
    iLastNoSkip = 0;
    copyPoint(&_polygonWindow[0], lastNoSkip);

    for (i = 1; i < _sizePolygonWindow; i++) {
        skip[i] = (distanceSquare(lastNoSkip, &_polygonWindow[i * 3]) < eps);
        if (!skip[i]) {
            copyPoint(&_polygonWindow[i * 3], lastNoSkip);
            iLastNoSkip = i;
        }
    }

    if (distanceSquare(&_polygonWindow[0], &_polygonWindow[iLastNoSkip * 3]) < eps)
        skip[iLastNoSkip] = true;

    numSkip = 0;
    for (i = 0; i < _sizePolygonWindow; i++)
        if (skip[i])
            numSkip++;

    sizePolygonWindow = _sizePolygonWindow - numSkip;
    assert(sizePolygonWindow > 2);

    polygonWindow = new double[3 * sizePolygonWindow];

    ii = 0;
    for (i = 0; i < _sizePolygonWindow; i++) {
        if (!skip[i]) {
            copyPoint(&_polygonWindow[i * 3], &polygonWindow[ii * 3]);
            ii++;
        }
    }
    delete[] skip;
    // Michael Andre end

    // numbering the edge and point in the following way:
    // The edge i is the edge between point i and i+1
    insideFlag = new bool[sizePolygonWindow];
    edgeSlopes = new double[sizePolygonWindow];
    reverseXY = new bool[sizePolygonWindow];

    const int XX = (planeToProject + 1) % 3;
    const int YY = (planeToProject + 2) % 3;

    for (int i = 0; i < sizePolygonWindow; i++) {
        // i is the ID of the edge
        double x1, y1, x2, y2, x3, y3;
        int p1_pos, p2_pos, p3_pos; // index of points of the edge in the triangle
        p1_pos = i;
        p2_pos = (i + 1) % sizePolygonWindow;
        p3_pos = (i + 2) % sizePolygonWindow;
        x1 = polygonWindow[p1_pos * 3 + XX];
        y1 = polygonWindow[p1_pos * 3 + YY];
        x2 = polygonWindow[p2_pos * 3 + XX];
        y2 = polygonWindow[p2_pos * 3 + YY];
        x3 = polygonWindow[p3_pos * 3 + XX];
        y3 = polygonWindow[p3_pos * 3 + YY];
        if (fabs(x2 - x1) < fabs(y2 - y1)) {
            reverseXY[i] = true;
            edgeSlopes[i] = (x2 - x1) / (y2 - y1);
            insideFlag[i] = (x3 - x1) > edgeSlopes[i] * (y3 - y1);
        } else {
            reverseXY[i] = false;
            edgeSlopes[i] = (y2 - y1) / (x2 - x1);
            insideFlag[i] = (y3 - y1) > edgeSlopes[i] * (x3 - x1);
        }
    }
}

PolygonClipper::~PolygonClipper() {
    delete[] insideFlag;
    delete[] edgeSlopes;
    delete[] reverseXY;
    delete[] polygonWindow;
}

bool PolygonClipper::clip(const double *polygonToBeClipped, int sizePolygonToBeClipped,
        std::vector<double*> *polygonResult) {
    // numbering the edge and point in the following way:
    // The edge i is the edge between point i and i+1.
    // The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice, 2nd edition, Addison-Wesley" P237 - P240.
    list<double*> *listInput = new list<double*>;
    list<double*> *listOutput = new list<double*>;

    // initialize listOutput to be the slave nodes
    for (int i = 0; i < sizePolygonToBeClipped; i++) {
        double *point = new double[3];
        copyPoint(&polygonToBeClipped[i * 3], point);
        listOutput->push_back(point);
    }
    //printElem(polygonWindow, sizePolygonWindow);
    //printElem(polygonToBeClipped, sizePolygonToBeClipped);
    const double TOL_SQR = EPS1 * EPS1 * longestEdgeLengthSquare(polygonWindow, sizePolygonWindow);

    for (int i = 0; i < sizePolygonWindow; i++) {
        // i is the edge ID
        // 1. make the output the new input
        for (list<double*>::iterator it = listInput->begin(); it != listInput->end(); it++)
            delete[] *it;
        listInput->clear();
        for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) {
            double *point = new double[3];
            copyPoint(*it, point);
            listInput->push_back(point);
        }
        for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++)
            delete[] *it;
        listOutput->clear();

        const double *edgeP1 = &polygonWindow[i * 3];
        const double *edgeP2 = &polygonWindow[(i + 1) % sizePolygonWindow * 3];
        int size = listInput->size();
        // 2. clip the input list by edge, create new output list
        for (list<double*>::iterator it = listInput->begin(); it != listInput->end(); it++) {
            double *p1 = *it;
            list<double*>::iterator it2 = it;
            advance(it2, 1);
            if (it2 == listInput->end())
                it2 = listInput->begin();
            double *p2 = *it2;
            bool inside1 = inside(i, p1);
            bool inside2 = inside(i, p2);
            if (!inside1 && inside2) { // there should be intersection, unless two lines are overlapped
                double *point1 = new double[3];
                if (intersect(edgeP1, edgeP2, p1, p2, planeToProject, point1))
                    listOutput->push_back(point1);
                double *point2 = new double[3];
                copyPoint(p2, point2);
                listOutput->push_back(point2);
            } else if (inside1 && inside2) {
                double *point = new double[3];
                copyPoint(p2, point);
                listOutput->push_back(point);
            } else if (inside1 && !inside2) { // there should be intersection, unless two lines are overlapped
                double *point = new double[3];
                if (intersect(edgeP1, edgeP2, p1, p2, planeToProject, point))
                    listOutput->push_back(point);
            } else {
                // do nothing
            }
        }

        /*cout << "===================================" << i << endl;
         printPoint(edgeP1);
         printPoint(edgeP2);
         double tmp[listOutput->size() * 3];
         int count = 0;
         for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) {
         for (int jj = 0; jj < 3; jj++)
         tmp[count * 3 + jj] = (*it)[jj];
         count++;
         }
         printElem(tmp, listOutput->size());
         cout << "===================================" << endl;*/

        for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) { // remove overlapped points in listOutput
            if (listOutput->size() == 1)
                break; // otherwise, segmentation fault
            double *p1 = *it;
            list<double*>::iterator it2 = it;
            advance(it2, 1);
            if (it2 == listOutput->end())
                it2 = listOutput->begin();
            double *p2 = *it2;
            if (distanceSquare(p1, p2) < TOL_SQR) {
                delete[] *it2; // otherwise, memory leak
                listOutput->erase(it2); // listOutput->erase(it) will cause segmentation fault
            }
        }
    }

    for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++) {
        double *point = new double[3];
        copyPoint(*it, point);
        polygonResult->push_back(point);
    }

    // release the storage
    for (list<double*>::iterator it = listInput->begin(); it != listInput->end(); it++)
        delete[] *it;
    delete listInput;

    for (list<double*>::iterator it = listOutput->begin(); it != listOutput->end(); it++)
        delete[] *it;
    delete listOutput;

    // judge whether there is really overlapped area or not
    int size = polygonResult->size();

    if (size < 3)
        return false;

    double *tmp = new double[size * 3];
    for (int i = 0; i < size; i++)
        for (int j = 0; j < 3; j++)
            tmp[i * 3 + j] = polygonResult->at(i)[j];
    if (computePolygonArea(tmp, size) < TOL_SQR) {
        delete[] tmp;
        return false; // it could happen the points are almost on the same line
    }
    delete[] tmp;

    return true;
    // since there is error when computing the intersection, now we have to check whether the center
    // of the points is inside both elements
    /*double polygon[polygonResult->size() * 3];
     for (int i = 0; i < polygonResult->size(); i++)
     for (int j = 0; j < 3; j++)
     polygon[i * 3 + j] = polygonResult->at(i)[j];
     double polygonCenter[3];
     computePolygonCenter(polygon, polygonResult->size(), polygonCenter);
     if (sizePolygonWindow == 3) {
     double localCoor[3];
     if (!computeLocalCoorInTriangle(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else if (sizePolygonWindow == 4) {
     double localCoor[2];
     if (!computeLocalCoorInQuad(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else {
     assert(false);
     }
     if (sizePolygonWindow == 3) {
     double localCoor[3];
     if (!computeLocalCoorInTriangle(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else if (sizePolygonWindow == 4) {
     double localCoor[2];
     if (!computeLocalCoorInQuad(polygonWindow, planeToProject, polygonCenter, localCoor))
     return false;
     } else {
     assert(false);
     }*/
}

bool PolygonClipper::inside(int edgeID, double *point) {
    const int XX = (planeToProject + 1) % 3;
    const int YY = (planeToProject + 2) % 3;
    int p1_pos = edgeID;
    double x1 = polygonWindow[p1_pos * 3 + XX];
    double y1 = polygonWindow[p1_pos * 3 + YY];
    double x4 = point[XX];
    double y4 = point[YY];
    if (reverseXY[edgeID])
        return ((x4 - x1) > edgeSlopes[edgeID] * (y4 - y1)) == insideFlag[edgeID];
    else
        return ((y4 - y1) > edgeSlopes[edgeID] * (x4 - x1)) == insideFlag[edgeID];
}

bool PolygonClipper::intersect(const double *la0, const double *la1, const double *lb0,
        const double *lb1, int _planeToProject, double *intersection) {
    // equation (in vector): alpha a0a1 - beta b0b1 = a0b0
    double unknown[2];
    double a0a1[3];
    double b0b1[3];
    double a0b0[3];

    for (int i = 0; i < 3; i++) {
        a0a1[i] = la1[i] - la0[i];
        b0b1[i] = lb1[i] - lb0[i];
        a0b0[i] = lb0[i] - la0[i];
    }

    double A[4];
    double a0a1_2D[2];
    double b0b1_2D[2];

    int XX = (_planeToProject + 1) % 3;
    int YY = (_planeToProject + 2) % 3;

    A[0] = a0a1_2D[0] = a0a1[XX];
    A[2] = b0b1_2D[0] = -b0b1[XX];
    unknown[0] = a0b0[XX];
    A[1] = a0a1_2D[1] = a0a1[YY];
    A[3] = b0b1_2D[1] = -b0b1[YY];
    unknown[1] = a0b0[YY];

    // if the system is not solvable, it means la and lb are parallel
    if (!solve2x2LinearSystem(A, unknown, EPS2))
        return false;

    /* intersection point should be on lb (0=<unknown[1]<=1), except:
     * 1. if la and lb are 'very' parallel, the numerical error could be large, therefore
     *    the intersection point may be outside of lb.
     *   (if la and lb are 'very' parallel, we can return any point on lb, or do not return
     *    any point, which won't cause problems)
     * 2. if end point of lb locates exactly on la, due to numerical error, the intersection
     *    point may locate 'a little bit' out of lb.
     * For both cases, it is enough to say, if unknown[1] is out of range, pick up one end
     * point of lb as the intersection point.
     */

    if (unknown[1] < 0.0)
        unknown[1] = 0.0;
    else if (unknown[1] > 1.0)
        unknown[1] = 1.0;

    for (int i = 0; i < 3; i++)
        intersection[i] = lb0[i] + unknown[1] * b0b1[i];
    /*for (int i = 0; i < 3; i++) {
     double tmp = intersection[i];
     assert(tmp==intersection[i]);//check whether intersection[i] is NAN
     }*/

    return true;
}

int computePlaneToProject(const double *planeNormal) {
    // find out the biggest entry in planeNormal, which indicates the plane to project
    double max = fabs(planeNormal[0]);
    int imax = 0;
    for (int i = 1; i < 3; i++) {
        if (fabs(planeNormal[i]) > max) {
            max = fabs(planeNormal[i]);
            imax = i;
        }
    }
    return imax;
}

void computePolygonCenter(const double *polygon, int num, double *center) {
    for (int i = 0; i < 3; i++)
        center[i] = 0;
    for (int i = 0; i < num; i++)
        for (int j = 0; j < 3; j++)
            center[j] += polygon[i * 3 + j];
    for (int i = 0; i < 3; i++)
        center[i] /= (double) num;
}

double computePolygonArea(const double *polygon, int num) {
    double center[3];
    computePolygonCenter(polygon, num, center);

    double area = 0.0;
    for (int i = 0; i < num; i++) {
        double clipTriangle[9];
        buildTrianagle(center, &polygon[i * 3], &polygon[(i + 1) % num * 3], clipTriangle);
        area += computeAreaOfTriangle(clipTriangle);
    }
    return area;
}

void projectToPlane(const double *pointOnPlane, const double *unitNormal, const double *points,
        int num, double *projections) {
    double ridge[3];
    double distance;
    double path[3];
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < 3; j++)
            ridge[j] = pointOnPlane[j] - points[i * 3 + j];
        distance = computeVectorDotProduct(unitNormal, ridge);
        for (int j = 0; j < 3; j++)
            path[j] = distance * unitNormal[j];

        for (int j = 0; j < 3; j++)
            projections[i * 3 + j] = points[i * 3 + j] + path[j];
    }
}

void computeNormalOfTriangle(const double *triangle, bool unitLength, double *normal) {
    const int NUM = 3; // three points in a triangle
    double length = 0; // length of the normal vector
    normal[0] = (triangle[2 * NUM + 1] - triangle[0 * NUM + 1])
            * (triangle[2 * NUM + 2] - triangle[1 * NUM + 2])
            - (triangle[2 * NUM + 1] - triangle[1 * NUM + 1])
                    * (triangle[2 * NUM + 2] - triangle[0 * NUM + 2]);
    normal[1] = -((triangle[2 * NUM + 0] - triangle[0 * NUM + 0])
            * (triangle[2 * NUM + 2] - triangle[1 * NUM + 2])
            - (triangle[2 * NUM + 0] - triangle[1 * NUM + 0])
                    * (triangle[2 * NUM + 2] - triangle[0 * NUM + 2]));
    normal[2] = (triangle[2 * NUM + 0] - triangle[0 * NUM + 0])
            * (triangle[2 * NUM + 1] - triangle[1 * NUM + 1])
            - (triangle[2 * NUM + 0] - triangle[1 * NUM + 0])
                    * (triangle[2 * NUM + 1] - triangle[0 * NUM + 1]);

    if (!unitLength)
        return;

    length = computeVectorLength(normal);

    for (int i = 0; i < 3; i++)
        normal[i] /= length;
}

void computeNormalOfQuad(const double *quad, bool unitLength, double *normal) {
    // d_N_d_xi[4] contains the partial derivative w.r.t. xi of the shape functions at (0, 0)
    double d_N_d_xi[4];
    // d_N_d_eta[4] contains the partial derivative w.r.t. eta of the shape functions at (0, 0)
    double d_N_d_eta[4];

    d_N_d_xi[0] = -0.25;
    d_N_d_xi[1] = 0.25;
    d_N_d_xi[2] = 0.25;
    d_N_d_xi[3] = -0.25;

    d_N_d_eta[0] = -0.25;
    d_N_d_eta[1] = -0.25;
    d_N_d_eta[2] = 0.25;
    d_N_d_eta[3] = 0.25;
    // g1 and g2 are the local basis vectors, and det(J)=||g1 x g2||
    double g1[3] = { 0, 0, 0 };
    double g2[3] = { 0, 0, 0 };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            g1[i] += quad[j * 3 + i] * d_N_d_xi[j];
            g2[i] += quad[j * 3 + i] * d_N_d_eta[j];
        }
    }
    double length;
    normal[0] = g1[1] * g2[2] - g2[1] * g1[2];
    normal[1] = -(g1[0] * g2[2] - g2[0] * g1[2]);
    normal[2] = g1[0] * g2[1] - g2[0] * g1[1];

    if (!unitLength)
        return;

    length = computeVectorLength(normal);

    for (int i = 0; i < 3; i++)
        normal[i] /= length;
}

double computeAreaOfTriangle(const double *triangle) {
    double area = 0.0;
    double normal[3];
    computeNormalOfTriangle(triangle, false, normal); // 2 times the length of the normal vector equals the area
    area = computeVectorLength(normal);
    area = 0.5 * area;
    return area;
}

double computeVectorLength(const double *vec) {
    double length = 0;
    for (int i = 0; i < 3; i++)
        length += vec[i] * vec[i];
    length = sqrt(length);
    return length;
}

double computeVectorLength2D(const double *vec) {
    double length = 0;
    for (int i = 0; i < 2; i++)
        length += vec[i] * vec[i];
    length = sqrt(length);
    return length;
}

void computeVectorCrossProduct(const double *vec1, const double *vec2, double *crossProduct) {
    crossProduct[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
    crossProduct[1] = -(vec1[0] * vec2[2] - vec2[0] * vec1[2]);
    crossProduct[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];
}

double computeVectorDotProduct(const double *vec1, const double *vec2) {
    double product = 0;
    for (int i = 0; i < 3; i++)
        product += vec1[i] * vec2[i];
    return product;
}

double distanceSquare(const double *p1, const double *p2) {
    double d;
    d = (p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
            + (p1[2] - p2[2]) * (p1[2] - p2[2]);
    return d;
}

void copyPoint(const double *origin, double *copy) {
    for (int i = 0; i < 3; i++)
        copy[i] = origin[i];
}

void copyElem(const double *origin, int size, double *copy) {
    for (int i = 0; i < 3 * size; i++)
        copy[i] = origin[i];
}

void buildTrianagle(const double *p0, const double *p1, const double *p2, double *triangle) {
    for (int i = 0; i < 3; i++) {
        triangle[i] = p0[i];
        triangle[3 + i] = p1[i];
        triangle[6 + i] = p2[i];
    }
}

double longestEdgeLengthSquare(const double *elem, int size) {
    double l[size];
    for (int i = 0; i < size; i++) {
        l[i] = distanceSquare(&elem[i * 3], &elem[(i + 1) % size * 3]);
    }
    double max = l[0];
    for (int i = 1; i < size; i++) {
        if (l[i] > max)
            max = l[i];
    }
    return max;
}

bool solve2x2LinearSystem(const double *A, double *b, double EPS) {
    /*
     * The function is called by PolygonClipper::intersect() and computeLocalCoorInQuad()
     * understand the system in this way:
     * A = [v1, v2], Ax = b means solving alpha and beta in:
     * alpha * v1 + beta * v2 = b
     * detA is the area of the quad formed by v1 and v2
     */
    // A is column major
    double detA = A[0] * A[3] - A[2] * A[1];
    double det0 = A[3] * b[0] - A[2] * b[1];
    double det1 = A[0] * b[1] - A[1] * b[0];

    double v1LengthSquare = A[0] * A[0] + A[1] * A[1];
    double v2LengthSquare = A[2] * A[2] + A[3] * A[3];
    double max = (v1LengthSquare > v2LengthSquare) ? v1LengthSquare : v2LengthSquare;
    assert(max > 1E-20); // assume the length of the vector cannot be too small

    if (fabs(detA) < EPS * fabs(det0))
        return false;
    if (fabs(detA) < EPS * fabs(det1))
        return false;
    if (fabs(detA) < 1E-15 * max) // Check detA relative to max. 1E-15 is chosen to pass TestMortarMath::testClipping()
        return false; // det0 and det1 could be 0!!! Without this, it could fail!!!!!

    b[0] = det0 / detA;
    b[1] = det1 / detA;
    return true;
}

void solve3x3LinearSystem(const double *A, int planeToProject, double *b) {
    int XX = (planeToProject + 1) % 3;
    int YY = (planeToProject + 2) % 3;
    double detA1A2 = A[XX + 3] * A[YY + 6] - A[YY + 3] * A[XX + 6];
    double detA0A2 = A[XX] * A[YY + 6] - A[YY] * A[XX + 6];
    double detA0A1 = A[XX] * A[YY + 3] - A[YY] * A[XX + 3];
    double detA0b = A[XX] * b[YY] - A[YY] * b[XX];
    double detbA2 = b[XX] * A[YY + 6] - b[YY] * A[XX + 6];
    double detA1b = A[XX + 3] * b[YY] - A[YY + 3] * b[XX];
    double detA = detA1A2 - detA0A2 + detA0A1;
    double det0 = detA1A2 - detbA2 - detA1b;
    double det1 = detbA2 - detA0A2 + detA0b;
    double det2 = detA1b - detA0b + detA0A1;
    b[0] = det0 / detA;
    b[1] = det1 / detA;
    b[2] = det2 / detA;
}

void computeMatrixProduct(int n, int m, const double *A, double *B) {
    double *C = new double[n * m];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double tmp = 0.0;
            for (int k = 0; k < n; k++) {
                tmp += A[i * n + k] * B[j + k * m];
            }
            C[i * m + j] = tmp;
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            B[i * m + j] = C[i * m + j];
    delete[] C;
}

void dcsrmv(char trans, int numRows, int numCols, const double *A, const int *JA, const int *IA,
        const double *x, double *y) {
    if (trans == 'N') {
        for (int i = 0; i < numRows; i++)
            y[i] = 0.0;
        for (int i = 0; i < numRows; i++) {
            int JA_row_begin = IA[i] - 1;
            int JA_row_end = IA[i + 1] - 1;
            for (int j = JA_row_begin; j < JA_row_end; j++) {
                int col = JA[j] - 1;
                y[i] += A[j] * x[col];
            }
        }
    } else if (trans == 'T') {
        for (int i = 0; i < numCols; i++)
            y[i] = 0.0;
        for (int i = 0; i < numRows; i++) {
            int JA_row_begin = IA[i] - 1;
            int JA_row_end = IA[i + 1] - 1;
            for (int j = JA_row_begin; j < JA_row_end; j++) {
                int col = JA[j] - 1;
                y[col] += A[j] * x[i];
            }
        }
    } else {
        assert(false);
    }
}

void dcsrsymv(int n, const double *A, const int *IA, const int *JA, const double *x, double *y) {
    for (int i = 0; i < n; i++)
        y[i] = 0.0;
    // use the upper part with diagonal
    for (int i = 0; i < n; i++) {
        int JA_row_begin = IA[i] - 1;
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            y[i] += A[j] * x[col];
        }
    }
    // use the lower part without diagonal
    for (int i = 0; i < n; i++) {
        int JA_row_begin = IA[i];
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            y[col] += A[j] * x[i];
        }
    }
}

void printElem(const double *elem, int size) {
    cout << endl << "Element: " << endl;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++)
            cout << elem[i * 3 + j] << "  ";
        cout << endl;
    }
}

void printPoint(const double *p) {
    cout << endl << "Point: ";
    for (int i = 0; i < 3; i++) {
        cout << p[i] << "  ";
    }
    cout << endl;
}

void printDiagonalMatrix(const double *A, int numRows) {
    cout << "======================================================" << endl;
    cout << "diagonal matrix:" << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numRows; j++) {
            if (i==j)
                cout << setw(15) << A[i];
            else
                cout << setw(15) << 0.0;
        }
        cout << endl;
    }
    cout << "======================================================" << endl;
}

void printGeneralMatrix(const double *A, int numRows, int numCols) {
    cout << "======================================================" << endl;
    cout << "general matrix:" << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            cout << setw(15) << A[i * numCols + j];
        }
        cout << endl;
    }
    cout << "======================================================" << endl;
}

void printCSRMatrixUnsymmetric(const double *A, const int *IA, const int *JA, int numRows,
        int numCols) {
    double *matrix = new double[numRows * numCols];
    for (int i = 0; i < numRows * numCols; i++) {
        matrix[i] = 0.0;
    }
    for (int i = 0; i < numRows; i++) {
        int JA_row_begin = IA[i] - 1;
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            matrix[i * numCols + col] = A[j];
        }
    }
    cout << "======================================================" << endl;
    cout << "unsymmetric CSR matrix:" << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            cout << setw(15) << matrix[i * numCols + j];
        }
        cout << endl;
    }
    cout << "======================================================" << endl;
    delete[] matrix;
}

void printCSRMatrixSymmetric(const double *A, const int *IA, const int *JA, int n) {
    double *matrix = new double[n * n];
    for (int i = 0; i < n * n; i++) {
        matrix[i] = 0.0;
    }
    // 1. set up upper part of the matrix (same as unsymmetric matrix)
    for (int i = 0; i < n; i++) {
        int JA_row_begin = IA[i] - 1;
        int JA_row_end = IA[i + 1] - 1;
        for (int j = JA_row_begin; j < JA_row_end; j++) {
            int col = JA[j] - 1;
            matrix[i * n + col] = A[j];
        }
    }
    // 2. copy the upper part to the lower part
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            matrix[i * n + j] = matrix[j * n + i];
        }
    }
    cout << "======================================================" << endl;
    cout << "symmetric CSR matrix:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(15) << matrix[i * n + j];
        }
        cout << endl;
    }
    cout << "======================================================" << endl;
    delete[] matrix;
}

} /* namespace MortarMath */
} /* namespace EMPIRE */
