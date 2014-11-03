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
#include <assert.h>
#include <math.h>
#include <iostream>
#include "stdlib.h"

#include "IGAMortarMath.h"
#include "MortarMath.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {
namespace IGAMortarMath {

const double EPS_IVERTIBILITYOFSQUAREMATRICES = 1e-14;

IGAPolygonClipper::IGAPolygonClipper(const double _u1, const double _u2, const double _v1,
        const double _v2) {

    polygonWindow[0] = _u1;
    polygonWindow[1] = _v1;
    polygonWindow[2] = 0;
    polygonWindow[3] = _u2;
    polygonWindow[4] = _v1;
    polygonWindow[5] = 0;
    polygonWindow[6] = _u2;
    polygonWindow[7] = _v2;
    polygonWindow[8] = 0;
    polygonWindow[9] = _u1;
    polygonWindow[10] = _v2;
    polygonWindow[11] = 0;
    clipper = new MortarMath::PolygonClipper(polygonWindow, 4, 2);
}

IGAPolygonClipper::~IGAPolygonClipper() {
    delete clipper;
}

bool IGAPolygonClipper::clip(const double *_polygonToBeClipped, int _sizePolygonToBeClipped,
        vector<double*> *_polygonResult) {
    double polygonToBeClipped[_sizePolygonToBeClipped * 3];
    for (int i = 0; i < _sizePolygonToBeClipped; i++) {
        polygonToBeClipped[i * 3] = _polygonToBeClipped[i * 2];
        polygonToBeClipped[i * 3 + 1] = _polygonToBeClipped[i * 2 + 1];
        polygonToBeClipped[i * 3 + 2] = 0;
    }
    return clipper->clip(polygonToBeClipped, _sizePolygonToBeClipped, _polygonResult);
}

void computeLinearCombinationValueFromVerticesValues(int _nNodes, int _nValue,
        const double *_values, const double* _coords, double *_returnValue) {
    double shapeFuncs[4];
    computeLowOrderShapeFunc(_nNodes, _coords, shapeFuncs);
    computeLinearCombination(_nNodes, _nValue, _values, shapeFuncs, _returnValue);
}

void computeLinearCombination(int _nNodes, int _nValue, const double *_values,
        const double *_shapeFuncs, double *_returnValue) {

    for (int i = 0; i < _nValue; i++) {
        _returnValue[i] = 0;
        for (int j = 0; j < _nNodes; j++)
            _returnValue[i] += _values[j * _nValue + i] * _shapeFuncs[j];
    }

}

void computeLowOrderShapeFunc(int _nNodes, const double *_coords, double *_shapeFuncs) {
    assert(_coords!=NULL);
    assert(_shapeFuncs!=NULL);
    if (_nNodes == 3) {
        _shapeFuncs[0] = 1 - _coords[0] - _coords[1];
        _shapeFuncs[1] = _coords[0];
        _shapeFuncs[2] = _coords[1];
    } else {
        _shapeFuncs[0] = (1 - _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[1] = (1 + _coords[0]) / 2 * (1 - _coords[1]) / 2;
        _shapeFuncs[2] = (1 + _coords[0]) / 2 * (1 + _coords[1]) / 2;
        _shapeFuncs[3] = (1 - _coords[0]) / 2 * (1 + _coords[1]) / 2;
    }

}

const double triGaussPoints1[2] = { 0.33333333333333, 0.33333333333333 };
const double triWeights1[1] = { 1.0 };

const double triGaussPoints3[6] = { 0.16666666666667, 0.16666666666667, 0.16666666666667,
        0.66666666666667, 0.66666666666667, 0.16666666666667 };
const double triWeights3[3] = { 0.33333333333333, 0.33333333333333, 0.33333333333333 };

const double triGaussPoints4[8] = { 0.33333333333333, 0.33333333333333, 0.20000000000000,
        0.20000000000000, 0.20000000000000, 0.60000000000000, 0.60000000000000, 0.20000000000000 };
const double triWeights4[4] = { -0.56250000000000, 0.52083333333333, 0.52083333333333,
        0.52083333333333, };

const double triGaussPoints6[12] = { 0.44594849091597, 0.44594849091597, 0.44594849091597,
        0.10810301816807, 0.10810301816807, 0.44594849091597, 0.09157621350977, 0.09157621350977,
        0.09157621350977, 0.81684757298046, 0.81684757298046, 0.09157621350977, };
const double triWeights6[6] = { 0.22338158967801, 0.22338158967801, 0.22338158967801,
        0.10995174365532, 0.10995174365532, 0.10995174365532 };

const double triGaussPoints7[14] = { 0.33333333333333, 0.33333333333333, 0.47014206410511,
        0.47014206410511, 0.47014206410511, 0.05971587178977, 0.05971587178977, 0.47014206410511,
        0.10128650732346, 0.10128650732346, 0.10128650732346, 0.79742698535309, 0.79742698535309,
        0.10128650732346, };
const double triWeights7[7] = { 0.22500000000000, 0.13239415278851, 0.13239415278851,
        0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483 };

const double triGaussPoints12[24] = { 0.24928674517091, 0.24928674517091, 0.24928674517091,
        0.50142650965818, 0.50142650965818, 0.24928674517091, 0.06308901449150, 0.06308901449150,
        0.06308901449150, 0.87382197101700, 0.87382197101700, 0.06308901449150, 0.31035245103378,
        0.63650249912140, 0.63650249912140, 0.05314504984482, 0.05314504984482, 0.31035245103378,
        0.63650249912140, 0.31035245103378, 0.31035245103378, 0.05314504984482, 0.05314504984482,
        0.63650249912140 };
const double triWeights12[12] = { 0.11678627572638, 0.11678627572638, 0.11678627572638,
        0.05084490637021, 0.05084490637021, 0.05084490637021, 0.08285107561837, 0.08285107561837,
        0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837 };

const double triGaussPoints13[26] = { 0.33333333333333, 0.33333333333333,
        0.26034596607904, 0.26034596607904,
        0.26034596607904, 0.47930806784192,
        0.47930806784192, 0.26034596607904,

        0.06513010290222, 0.06513010290222,
        0.06513010290222, 0.86973979419557,
        0.86973979419557, 0.06513010290222,

        0.31286549600487, 0.63844418856981,
        0.63844418856981, 0.04869031542532,
        0.04869031542532, 0.31286549600487,
        0.63844418856981, 0.31286549600487,
        0.31286549600487, 0.04869031542532,
        0.04869031542532, 0.63844418856981 };
const double triWeights13[13] = { -0.14957004446768,
        0.17561525743321, 0.17561525743321, 0.17561525743321,
        0.05334723560884, 0.05334723560884, 0.05334723560884,
        0.07711376089026, 0.07711376089026, 0.07711376089026, 0.07711376089026, 0.07711376089026, 0.07711376089026 };

const double triGaussPoints16[32] = { 0.33333333333333, 0.33333333333333, 0.45929258829272,
        0.45929258829272, 0.45929258829272, 0.08141482341455, 0.08141482341455, 0.45929258829272,
        0.17056930775176, 0.17056930775176, 0.17056930775176, 0.65886138449648, 0.65886138449648,
        0.17056930775176, 0.05054722831703, 0.05054722831703, 0.05054722831703, 0.89890554336594,
        0.89890554336594, 0.05054722831703, 0.26311282963464, 0.72849239295540, 0.72849239295540,
        0.00839477740996, 0.00839477740996, 0.26311282963464, 0.72849239295540, 0.26311282963464,
        0.26311282963464, 0.00839477740996, 0.00839477740996, 0.72849239295540 };
const double triWeights16[16] = { 0.14431560767779, 0.09509163426728, 0.09509163426728,
        0.09509163426728, 0.10321737053472, 0.10321737053472, 0.10321737053472, 0.03245849762320,
        0.03245849762320, 0.03245849762320, 0.02723031417443, 0.02723031417443, 0.02723031417443,
        0.02723031417443, 0.02723031417443, 0.02723031417443 };

GaussQuadratureOnTriangle::GaussQuadratureOnTriangle(int _numGaussPoints) :
        GaussQuadrature(_numGaussPoints) {
    switch (_numGaussPoints) {
    case 1:
        gaussPoints = triGaussPoints1;
        weights = triWeights1;
        break;
    case 3:
        gaussPoints = triGaussPoints3;
        weights = triWeights3;
        break;
    case 4:
        gaussPoints = triGaussPoints4;
        weights = triWeights4;
        break;
    case 6:
        gaussPoints = triGaussPoints6;
        weights = triWeights6;
        break;
    case 7:
        gaussPoints = triGaussPoints7;
        weights = triWeights7;
        break;
    case 12:
        gaussPoints = triGaussPoints12;
        weights = triWeights12;
        break;
    case 13:
        gaussPoints = triGaussPoints13;
        weights = triWeights13;
        break;
    case 16:
        gaussPoints = triGaussPoints16;
        weights = triWeights16;
        break;
    default:
        ERROR_OUT()<<"Number of Gauss Points for Triangle = " << numGaussPoints << "doesn't exist! Please choose from 1,3,4,6,7,12,13,16." << endl;
        exit(EXIT_FAILURE);

    }
}

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

static const double tmpG41 = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
static const double tmpW41 = 0.5 + sqrt(30.0) / 36.0;
static const double tmpG42 = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
static const double tmpW42 = 0.5 - sqrt(30.0) / 36.0;
const double GP4[4] = { -tmpG42, -tmpG41, tmpG41, tmpG42 };
const double W4[4] = { tmpW42, tmpW41, tmpW41, tmpW42 };
const double quadGaussPoints16[32] = { GP4[0], GP4[0], GP4[0], GP4[1], GP4[0], GP4[2], GP4[0],
        GP4[3], GP4[1], GP4[0], GP4[1], GP4[1], GP4[1], GP4[2], GP4[1], GP4[3], GP4[2], GP4[0],
        GP4[2], GP4[1], GP4[2], GP4[2], GP4[2], GP4[3], GP4[3], GP4[0], GP4[3], GP4[1], GP4[3],
        GP4[2], GP4[3], GP4[3] };
const double quadWeights16[16] = { W4[0] * W4[0], W4[0] * W4[1], W4[0] * W4[2], W4[0] * W4[3], W4[1]
        * W4[0], W4[1] * W4[1], W4[1] * W4[2], W4[1] * W4[3], W4[2] * W4[0], W4[2] * W4[1], W4[2]
        * W4[2], W4[2] * W4[3], W4[3] * W4[0], W4[3] * W4[1], W4[3] * W4[2], W4[3] * W4[3] };

static const double tmpG51 = sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3;
static const double tmpG52 = sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3;
static const double tmpW51 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
static const double tmpW52 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
const double GP5[5] = { -tmpG52, -tmpG51, 0.0, tmpG51, tmpG52 };
const double W5[5] = { tmpW52, tmpW51, 128.0 / 225.0, tmpW51, tmpW52 };
const double quadGaussPoints25[50] = { GP5[0], GP5[0], GP5[0], GP5[1], GP5[0], GP5[2], GP5[0],
        GP5[3], GP5[0], GP5[4], GP5[1], GP5[0], GP5[1], GP5[1], GP5[1], GP5[2], GP5[1], GP5[3],
        GP5[1], GP5[4], GP5[2], GP5[0], GP5[2], GP5[1], GP5[2], GP5[2], GP5[2], GP5[3], GP5[2],
        GP5[4], GP5[3], GP5[0], GP5[3], GP5[1], GP5[3], GP5[2], GP5[3], GP5[3], GP5[3], GP5[4],
        GP5[4], GP5[0], GP5[4], GP5[1], GP5[4], GP5[2], GP5[4], GP5[3], GP5[4], GP5[4] };
const double quadWeights25[25] = { W5[0] * W5[0], W5[0] * W5[1], W5[0] * W5[2], W5[0] * W5[3], W5[0]
        * W5[4], W5[1] * W5[0], W5[1] * W5[1], W5[1] * W5[2], W5[1] * W5[3], W5[1] * W5[4], W5[2]
        * W5[0], W5[2] * W5[1], W5[2] * W5[2], W5[2] * W5[3], W5[2] * W5[4], W5[3] * W5[0], W5[3]
        * W5[1], W5[3] * W5[2], W5[3] * W5[3], W5[3] * W5[4], W5[4] * W5[0], W5[4] * W5[1], W5[4]
        * W5[2], W5[4] * W5[3], W5[4] * W5[4] };

GaussQuadratureOnQuad::GaussQuadratureOnQuad(int _numGaussPoints) :
        GaussQuadrature(_numGaussPoints) {
    switch (_numGaussPoints) {
    case 1:
        gaussPoints = quadGaussPoints1;
        weights = quadWeights1;
        break;
    case 4:
        gaussPoints = quadGaussPoints4;
        weights = quadWeights4;
        break;
    case 9:
        gaussPoints = quadGaussPoints9;
        weights = quadWeights9;
        break;
    case 16:
        gaussPoints = quadGaussPoints16;
        weights = quadWeights16;
        break;
    case 25:
        gaussPoints = quadGaussPoints25;
        weights = quadWeights25;
        break;
    default:
        ERROR_OUT()<<"Number of Gauss Points for Quadrilateral = " << numGaussPoints << "doesn't exist! Please choose from 1,4,9,16,25." << endl;
        exit(EXIT_FAILURE);
    }
}

bool computeLocalCoordsInTriangle(const double *_coordsTri, const double *_coordsNode,
        double* _localCoords) {
    assert(_coordsTri!=NULL);
    assert(_coordsNode!=NULL);

    double area = computeAreaTriangle(_coordsTri[2] - _coordsTri[0], _coordsTri[3] - _coordsTri[1],
            0, _coordsTri[4] - _coordsTri[0], _coordsTri[5] - _coordsTri[1], 0);
    double area1 = computeAreaTriangle(_coordsTri[2] - _coordsNode[0],
            _coordsTri[3] - _coordsNode[1], 0, _coordsTri[4] - _coordsNode[0],
            _coordsTri[5] - _coordsNode[1], 0);
    double area2 = computeAreaTriangle(_coordsTri[0] - _coordsNode[0],
            _coordsTri[1] - _coordsNode[1], 0, _coordsTri[4] - _coordsNode[0],
            _coordsTri[5] - _coordsNode[1], 0);
    _localCoords[0] = area1 / area;
    _localCoords[1] = area2 / area;
    if (_localCoords[0] < 0 || _localCoords[0] > 1 || _localCoords[1] < 0 || _localCoords[1] > 1
            || _localCoords[0] + _localCoords[1] > 1)
        return false;
    return true;
}

bool computeLocalCoordsInQuad(const double *_coordsQuad, const double *_coordsNode,
        double* _localCoords) {
    /*
     * So we use two coordinates among x, y, z.
     * This indicates projection.
     * Choose among planes (xy or yz or zx) the one which has the smallest angle
     * with the quad normal.
     */
    assert(_coordsQuad!=NULL);
    assert(_coordsNode!=NULL);
    double x[4];
    double y[4];
    for (int i = 0; i < 4; i++) {
        x[i] = _coordsQuad[i * 2];
        y[i] = _coordsQuad[i * 2 + 1];
        if (x[i] == _coordsNode[0] & y[i] == _coordsNode[1]) {
            _localCoords[0] = (i == 0 || i == 3) ? (-1) : (1);
            _localCoords[1] = (i == 0 || i == 1) ? (-1) : (1);
            return true;
        }
    }
    double x0 = _coordsNode[0];
    double y0 = _coordsNode[1];

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
    _localCoords[0] = 0;
    _localCoords[1] = 0;
//int dummy[2];
    const double EPS = 1E-15;

    const int MAX_ITER_NUM = 100;
    for (int i = 0; i < MAX_ITER_NUM; i++) {
        J_T[0] = b1 + d1 * _localCoords[1];
        J_T[2] = c1 + d1 * _localCoords[0];
        J_T[1] = b2 + d2 * _localCoords[1];
        J_T[3] = c2 + d2 * _localCoords[0];
        F[0] = a1 + b1 * _localCoords[0] + c1 * _localCoords[1]
                + d1 * _localCoords[0] * _localCoords[1];
        F[1] = a2 + b2 * _localCoords[0] + c2 * _localCoords[1]
                + d2 * _localCoords[0] * _localCoords[1];
        delta[0] = -F[0];
        delta[1] = -F[1];

        solve2x2LinearSystem(J_T, delta, EPS);
        if (fabs(delta[0]) < EPS && fabs(delta[1]) < EPS) {
            assert(i < 100);
            break;
        }
        _localCoords[0] += delta[0];
        _localCoords[1] += delta[1];
    }
    for (int i = 0; i < 2; i++) {
        if (_localCoords[i] > 1.0)
            return false;
        if (_localCoords[i] < -1.0)
            return false;
    }
    return true;
}

bool computeIntersectionBetweenLineAndTriangle(const double *_X, const double* _X0,
        const double* _n, double* _localCoords) {
    double A[9];
    //  A(1:3,1) = X1-X3;
    for (int i = 0; i < 3; i++)
        A[i * 3] = _X[i] - _X[i + 6];

    //  A(1:3,2) = X2-X3;
    for (int i = 0; i < 3; i++)
        A[i * 3 + 1] = _X[i + 3] - _X[i + 6];

    //  A(1:3,3) = -n;
    for (int i = 0; i < 3; i++)
        A[i * 3 + 2] = -_n[i];

    double b[3];
    // b = X0 - X3
    b[0] = _X0[0] - _X[6];
    b[1] = _X0[1] - _X[7];
    b[2] = _X0[2] - _X[8];

    solve3x3LinearSystem(A, b, EPS_IVERTIBILITYOFSQUAREMATRICES);

    _localCoords[0] = b[0];
    _localCoords[1] = b[1];
    if (fabs(_localCoords[0] - 0.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[0] = 0.0;
    if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[0] = 1.0;
    if (fabs(_localCoords[0] - 0.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[1] = 0.0;
    if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
        _localCoords[1] = 1.0;

    if (_localCoords[0] >= 0.0 && _localCoords[0] <= 1.0 && _localCoords[1] >= 0.0
            && _localCoords[1] <= 1.0)
        return true;
    else
        return false;

}

bool computeIntersectionBetweenLineAndQuad(const double *_X, const double* _X0, const double* _n,
        double* _localCoords) {
    double x[3] = { 0.0, 0.0, 0.0 };
    double f[3];
    double df[9];
    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 3; j++) {
            // f = 0.25*((-X1+X2+X3-X4)*x(1)+(-X1-X2+X3+X4)*x(2)+(X1-X2+X3-X4)*x(1)*x(2)+X1+X2+X3+X4)-x(3)*n-X0;
            f[j] = 0.25
                    * ((-_X[j] + _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[0]
                            + (-_X[j] - _X[j + 3] + _X[j + 6] + _X[j + 9]) * x[1]
                            + (_X[j] - _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[0] * x[1] + _X[j]
                            + _X[j + 3] + _X[j + 6] + _X[j + 9]) - x[2] * _n[j] - _X0[j];
            //df = [.25*((-X1+X2+X3-X4)+(X1-X2+X3-X4)*x(2)) 0.25*((-X1-X2+X3+X4)+(X1-X2+X3-X4)*x(1)) -n];
            df[j * 3] = 0.25
                    * ((-_X[j] + _X[j + 3] + _X[j + 6] - _X[j + 9])
                            + (_X[j] - _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[1]);
            df[j * 3 + 1] = 0.25
                    * ((-_X[j] - _X[j + 3] + _X[j + 6] + _X[j + 9])
                            + (_X[j] - _X[j + 3] + _X[j + 6] - _X[j + 9]) * x[0]);
            df[j * 3 + 2] = -_n[j];
        }
        if (sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]) > 1e-12) {
            solve3x3LinearSystem(df, f, EPS_IVERTIBILITYOFSQUAREMATRICES);
            for (int j = 0; j < 3; j++)
                x[j] -= f[j];
        } else {
            _localCoords[0] = x[0];
            _localCoords[1] = x[1];
            if (fabs(_localCoords[0] - 0.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[0] = 0.0;
            if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[0] = 1.0;
            if (fabs(_localCoords[0] - 0.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[1] = 0.0;
            if (fabs(_localCoords[0] - 1.0) < EPS_IVERTIBILITYOFSQUAREMATRICES)
                _localCoords[1] = 1.0;

            if (_localCoords[0] >= 0.0 && _localCoords[0] <= 1.0 && _localCoords[1] >= 0.0
                    && _localCoords[1] <= 1.0)
                return true;
            else
                return false;
        }
    }
    return false;

}

bool solve2x2LinearSystem(const double* _A, double* _b, double _EPS) {
// A is column major
    double detA = _A[0] * _A[3] - _A[2] * _A[1];
    double det0 = _A[3] * _b[0] - _A[2] * _b[1];
    double det1 = _A[0] * _b[1] - _A[1] * _b[0];
    if (fabs(detA) < _EPS * fabs(det0))
        return false;
    if (fabs(detA) < _EPS * fabs(det1))
        return false;
    if (detA == 0) // det0 and det1 could be 0!!!
        return false; // without this, it could fail!!!!!
    _b[0] = det0 / detA;
    _b[1] = det1 / detA;
    return true;
}

bool solve3x3LinearSystem(const double* _A, double* _b, double _EPS) {
// A is column major
    double A[9];
    double b[3];
    double detA = det3x3(_A);
    if (fabs(detA) < _EPS)
        return false;
    for (int i = 0; i < 3; i++)
        b[i] = _b[i];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 9; j++)
            A[j] = _A[j];
        for (int j = 0; j < 3; j++)
            A[j * 3 + i] = b[j];
        _b[i] = det3x3(A) / detA;
    }
    return true;

}

double det3x3(const double* _A) {
    return _A[0] * _A[4] * _A[8] + _A[1] * _A[5] * _A[6] + _A[2] * _A[3] * _A[7]
            - _A[0] * _A[5] * _A[7] - _A[1] * _A[3] * _A[8] - _A[2] * _A[4] * _A[6];
}

double computeAreaTriangle(double _x1, double _y1, double _z1, double _x2, double _y2, double _z2) {
    double x = _y1 * _z2 - _y2 * _z1;
    double y = _z1 * _x2 - _z2 * _x1;
    double z = _x1 * _y2 - _x2 * _y1;
    return sqrt(x * x + y * y + z * z) / 2;

}

double computeCrossProduct2D(double _x1, double _y1, double _x2, double _y2) {
    return _x1 * _y2 - _x2 * _y1;
}

double computePointDistance(double* _x1, double* _x2) {
    return sqrt(
            (_x1[0] - _x2[0]) * (_x1[0] - _x2[0]) + (_x1[1] - _x2[1]) * (_x1[1] - _x2[1])
                    + (_x1[2] - _x2[2]) * (_x1[2] - _x2[2]));
}

void cleanPolygon(std::vector<double>& polygon) {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polygon.size()/2;
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polygon[p1*2]==polygon[p2*2] && polygon[p1*2+1]==polygon[p2*2+1];

		double v1x=polygon[p2*2]-polygon[p1*2];
		double v1y=polygon[p2*2+1]-polygon[p1*2+1];
		double v2x=polygon[p3*2]-polygon[p1*2];
		double v2y=polygon[p3*2+1]-polygon[p1*2+1];
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		double v=computeCrossProduct2D(v1x,v1y,v2x,v2y);
		// Result of cross product only in Z direction
		bool isColinear=fabs(v)<1e-9?true:false;
		//double n1=v1x*v1x+v1y*v1y;
		//double n2=v2x*v2x+v2y*v2y;
		//bool isWrong=n1>3*fmin(n1,n2);
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polygon.erase(polygon.begin()+p2*2+1);
			polygon.erase(polygon.begin()+p2*2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polygon.size()/2;
			k--;
		}
	}
}
void cleanPolygon(std::vector<pair<double,double> >& polygon) {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polygon.size();
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polygon[p1]==polygon[p2];

		double v1x=polygon[p2].first-polygon[p1].first;
		double v1y=polygon[p2].second-polygon[p1].second;
		double v2x=polygon[p3].first-polygon[p1].first;
		double v2y=polygon[p3].second-polygon[p1].second;
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		double v=computeCrossProduct2D(v1x,v1y,v2x,v2y);
		// Result of cross product only in Z direction
		bool isColinear=fabs(v)<1e-9?true:false;
		//double n1=v1x*v1x+v1y*v1y;
		//double n2=v2x*v2x+v2y*v2y;
		//bool isWrong=n1>3*fmin(n1,n2);
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polygon.erase(polygon.begin()+p2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polygon.size();
			k--;
		}
	}
}

void cleanPolygon(std::vector<pair<double,double> >& polygon,std::vector<pair<double,double> >& polygonSlave) {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polygon.size();
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polygon[p1]==polygon[p2];

		double v1x=polygon[p2].first-polygon[p1].first;
		double v1y=polygon[p2].second-polygon[p1].second;
		double v2x=polygon[p3].first-polygon[p1].first;
		double v2y=polygon[p3].second-polygon[p1].second;
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		double v=computeCrossProduct2D(v1x,v1y,v2x,v2y);
		// Result of cross product only in Z direction
		bool isColinear=fabs(v)<1e-9?true:false;
		//double n1=v1x*v1x+v1y*v1y;
		//double n2=v2x*v2x+v2y*v2y;
		//bool isWrong=n1>3*fmin(n1,n2);
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polygon.erase(polygon.begin()+p2);
			polygonSlave.erase(polygonSlave.begin()+p2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polygon.size();
			k--;
		}
	}
}

} // end of IGAMortarMath
} // end of EMPIRE

