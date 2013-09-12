#include "IGAMortarMath.h"
#include "MortarMath.h"
#include <math.h>
#include <assert.h>
#include <iostream>

using namespace std;

namespace EMPIRE {
namespace IGAMortarMath {

IGAPolygonClipper::IGAPolygonClipper(const double _u1, const double _u2,
		const double _v1, const double _v2) {

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

bool IGAPolygonClipper::clip(const double *_polygonToBeClipped,
		int _sizePolygonToBeClipped, vector<double*> *_polygonResult) {
	double polygonToBeClipped[_sizePolygonToBeClipped * 3];
	for (int i = 0; i < _sizePolygonToBeClipped; i++) {
		polygonToBeClipped[i * 3] = _polygonToBeClipped[i * 2];
		polygonToBeClipped[i * 3 + 1] = _polygonToBeClipped[i * 2 + 1];
		polygonToBeClipped[i * 3 + 2] = 0;
	}
	return clipper->clip(polygonToBeClipped, _sizePolygonToBeClipped,
			_polygonResult);
}

void calcValueFromShapeFuncs(const double *_values, const double *_shapeFuncs,
		int _nNodes, int _nValue, double *_returnValue) {

	for (int i = 0; i < _nValue; i++) {
		_returnValue[i] = 0;
		for (int j = 0; j < _nNodes; j++)
			_returnValue[i] += _values[j * _nValue + i] * _shapeFuncs[j];
	}

}

void calcCoords(const double *_coordsQuadrature, const double *_gaussPoint,
		int _nNodes, double *_coords) {
	assert(_gaussPoint!=NULL);
	assert(_coords!=NULL);
	assert(_coordsQuadrature!=NULL);

	double *shapeFuncs;
	if (_nNodes == 3) {
		shapeFuncs = new double(3);
		calcLinearShapeFunc(_gaussPoint, 3, shapeFuncs);
	} else {
		shapeFuncs = new double(4);
		calcLinearShapeFunc(_gaussPoint, 4, shapeFuncs);
	}
	calcValueFromShapeFuncs(_coordsQuadrature, shapeFuncs, _nNodes, 2, _coords);
	delete shapeFuncs;
}

void calcLinearShapeFunc(const double *_coords, int _nNodes,
		double *_shapeFuncs) {
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

const double triGaussPoints3[6] = { 0.16666666666667, 0.16666666666667,
		0.16666666666667, 0.66666666666667, 0.66666666666667, 0.16666666666667 };
const double triWeights3[3] = { 0.33333333333333, 0.33333333333333,
		0.33333333333333 };

const double triGaussPoints4[8] = { 0.33333333333333, 0.33333333333333,
		0.20000000000000, 0.20000000000000, 0.20000000000000, 0.60000000000000,
		0.60000000000000, 0.20000000000000 };
const double triWeights4[4] = { 0.56250000000000, 0.52083333333333,
		0.52083333333333, 0.52083333333333, };

const double triGaussPoints6[12] = { 0.44594849091597, 0.44594849091597,
		0.44594849091597, 0.10810301816807, 0.10810301816807, 0.44594849091597,
		0.09157621350977, 0.09157621350977, 0.09157621350977, 0.81684757298046,
		0.81684757298046, 0.09157621350977, };
const double triWeights6[6] = { 0.22338158967801, 0.22338158967801,
		0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532 };

const double triGaussPoints7[14] =
		{ 0.33333333333333, 0.33333333333333, 0.47014206410511,
				0.47014206410511, 0.47014206410511, 0.05971587178977,
				0.05971587178977, 0.47014206410511, 0.10128650732346,
				0.10128650732346, 0.10128650732346, 0.79742698535309,
				0.79742698535309, 0.10128650732346, };
const double triWeights7[7] = { 0.22500000000000, 0.13239415278851,
		0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483,
		0.12593918054483 };

const double triGaussPoints12[24] = { 0.24928674517091, 0.24928674517091,
		0.24928674517091, 0.50142650965818, 0.50142650965818, 0.24928674517091,
		0.06308901449150, 0.06308901449150, 0.06308901449150, 0.87382197101700,
		0.87382197101700, 0.06308901449150, 0.31035245103378, 0.63650249912140,
		0.63650249912140, 0.05314504984482, 0.05314504984482, 0.31035245103378,
		0.63650249912140, 0.31035245103378, 0.31035245103378, 0.05314504984482,
		0.05314504984482, 0.63650249912140 };
const double triWeights12[12] = { 0.11678627572638, 0.11678627572638,
		0.11678627572638, 0.05084490637021, 0.05084490637021, 0.05084490637021,
		0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837,
		0.08285107561837, 0.08285107561837 };

const double triGaussPoints13[26] = { 0.33333333333333, 0.33333333333333,
		0.26034596607904, 0.26034596607904, 0.26034596607904, 0.47930806784192,
		0.47930806784192, 0.26034596607904, 0.06513010290222, 0.06513010290222,
		0.06513010290222, 0.86973979419557, 0.86973979419557, 0.06513010290222,
		0.31286549600487, 0.63844418856981, 0.63844418856981, 0.04869031542532,
		0.04869031542532, 0.31286549600487, 0.63844418856981, 0.31286549600487,
		0.31286549600487, 0.04869031542532, 0.04869031542532, 0.63844418856981 };
const double triWeights13[13] = { -0.14957004446768, 0.17561525743321,
		0.17561525743321, 0.17561525743321, 0.05334723560884, 0.05334723560884,
		0.05334723560884, 0.07711376089026, 0.07711376089026, 0.07711376089026,
		0.07711376089026, 0.07711376089026, 0.07711376089026 };

const double triGaussPoints16[32] = { 0.33333333333333, 0.33333333333333,
		0.45929258829272, 0.45929258829272, 0.45929258829272, 0.08141482341455,
		0.08141482341455, 0.45929258829272, 0.17056930775176, 0.17056930775176,
		0.17056930775176, 0.65886138449648, 0.65886138449648, 0.17056930775176,
		0.05054722831703, 0.05054722831703, 0.05054722831703, 0.89890554336594,
		0.89890554336594, 0.05054722831703, 0.26311282963464, 0.72849239295540,
		0.72849239295540, 0.00839477740996, 0.00839477740996, 0.26311282963464,
		0.72849239295540, 0.26311282963464, 0.26311282963464, 0.00839477740996,
		0.00839477740996, 0.72849239295540 };
const double triWeights16[16] = { 0.14431560767779, 0.09509163426728,
		0.09509163426728, 0.09509163426728, 0.10321737053472, 0.10321737053472,
		0.10321737053472, 0.03245849762320, 0.03245849762320, 0.03245849762320,
		0.02723031417443, 0.02723031417443, 0.02723031417443, 0.02723031417443,
		0.02723031417443, 0.02723031417443 };

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
		assert(0);
		break;
	}
}

const double quadGaussPoints1[2] = { 0.0, 0.0 };
const double quadWeights1[1] = { 4.0 };

static const double tmpSqrt13 = sqrt(1.0 / 3.0);
const double quadGaussPoints4[8] = { tmpSqrt13, tmpSqrt13, -tmpSqrt13,
		tmpSqrt13, -tmpSqrt13, -tmpSqrt13, tmpSqrt13, -tmpSqrt13 };
const double quadWeights4[4] = { 1.0, 1.0, 1.0, 1.0 };

static const double tmpSqrt35 = sqrt(3.0 / 5.0);
const double quadGaussPoints9[18] = { 0.0, 0.0, -tmpSqrt35, 0, 0, -tmpSqrt35,
		tmpSqrt35, 0, 0, tmpSqrt35, -tmpSqrt35, -tmpSqrt35, tmpSqrt35,
		-tmpSqrt35, -tmpSqrt35, tmpSqrt35, tmpSqrt35, tmpSqrt35 };
const double quadWeights9[9] =
		{ 64.0 / 81.0, 40.0 / 81.0, 40.0 / 81.0, 40.0 / 81.0, 40.0 / 81.0, 25.0
				/ 81.0, 25.0 / 81.0, 25.0 / 81.0, 25.0 / 81.0 };

static const double tmpG41 = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
static const double tmpW41 = 0.5 + sqrt(30.0) / 36.0;
static const double tmpG42 = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
static const double tmpW42 = 0.5 - sqrt(30.0) / 36.0;
const double GP4[4] = { -tmpG42, -tmpG41, tmpG41, tmpG42 };
const double W4[4] = { tmpW42, tmpW41, tmpW41, tmpW42 };
const double quadGaussPoints16[32] = { GP4[0], GP4[0], GP4[0], GP4[1], GP4[0],
		GP4[2], GP4[0], GP4[3], GP4[1], GP4[0], GP4[1], GP4[1], GP4[1], GP4[2],
		GP4[1], GP4[3], GP4[2], GP4[0], GP4[2], GP4[1], GP4[2], GP4[2], GP4[2],
		GP4[3], GP4[3], GP4[0], GP4[3], GP4[1], GP4[3], GP4[2], GP4[3], GP4[3] };
const double quadWeights16[16] = { W4[0] * W4[0], W4[0] * W4[1], W4[0] * W4[2],
		W4[0] * W4[3], W4[1] * W4[0], W4[1] * W4[1], W4[1] * W4[2], W4[1]
				* W4[3], W4[2] * W4[0], W4[2] * W4[1], W4[2] * W4[2], W4[2]
				* W4[3], W4[3] * W4[0], W4[3] * W4[1], W4[3] * W4[2], W4[3]
				* W4[3] };

static const double tmpG51 = sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3;
static const double tmpG52 = sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3;
static const double tmpW51 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
static const double tmpW52 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
const double GP5[5] = { -tmpG52, -tmpG51, 0.0, tmpG51, tmpG52 };
const double W5[5] = { tmpW52, tmpW51, 128.0 / 225.0, tmpW51, tmpW52 };
const double quadGaussPoints25[50] = { GP5[0], GP5[0], GP5[0], GP5[1], GP5[0],
		GP5[2], GP5[0], GP5[3], GP5[0], GP5[4], GP5[1], GP5[0], GP5[1], GP5[1],
		GP5[1], GP5[2], GP5[1], GP5[3], GP5[1], GP5[4], GP5[2], GP5[0], GP5[2],
		GP5[1], GP5[2], GP5[2], GP5[2], GP5[3], GP5[2], GP5[4], GP5[3], GP5[0],
		GP5[3], GP5[1], GP5[3], GP5[2], GP5[3], GP5[3], GP5[3], GP5[4], GP5[4],
		GP5[0], GP5[4], GP5[1], GP5[4], GP5[2], GP5[4], GP5[3], GP5[4], GP5[4] };
const double quadWeights25[25] = { W5[0] * W5[0], W5[0] * W5[1], W5[0] * W5[2],
		W5[0] * W5[3], W5[0] * W5[4], W5[1] * W5[0], W5[1] * W5[1], W5[1]
				* W5[2], W5[1] * W5[3], W5[1] * W5[4], W5[2] * W5[0], W5[2]
				* W5[1], W5[2] * W5[2], W5[2] * W5[3], W5[2] * W5[4], W5[3]
				* W5[0], W5[3] * W5[1], W5[3] * W5[2], W5[3] * W5[3], W5[3]
				* W5[4], W5[4] * W5[0], W5[4] * W5[1], W5[4] * W5[2], W5[4]
				* W5[3], W5[4] * W5[4] };

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
		assert(0);
		break;
	}
}

bool computeLocalCoordsInTriangle(const double *_coordsTri,
		const double *_coordsNode, double* _localCoords) {
	assert(_coordsTri!=NULL);
	assert(_coordsNode!=NULL);

	double area = areaTriangle(_coordsTri[2] - _coordsTri[0],
			_coordsTri[3] - _coordsTri[1], 0, _coordsTri[4] - _coordsTri[0],
			_coordsTri[5] - _coordsTri[1], 0);
	double area1 = areaTriangle(_coordsTri[2] - _coordsNode[0],
			_coordsTri[3] - _coordsNode[1], 0, _coordsTri[4] - _coordsNode[0],
			_coordsTri[5] - _coordsNode[1], 0);
	double area2 = areaTriangle(_coordsTri[0] - _coordsNode[0],
			_coordsTri[1] - _coordsNode[1], 0, _coordsTri[4] - _coordsNode[0],
			_coordsTri[5] - _coordsNode[1], 0);
	_localCoords[0] = area1 / area;
	_localCoords[1] = area2 / area;
	if (_localCoords[0] < 0 || _localCoords[0] > 1 || _localCoords[1] < 0
			|| _localCoords[1] > 1)
		return false;
	return true;
}

bool computeLocalCoordsInQuad(const double *_coordsQuad,
		const double *_coordsNode, double* _localCoords) {
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
			/*if (i>=10) {
			 printElem(quad,4);
			 printPoint(point);
			 cout << "plane to project:  "<<planeToProject<<endl;
			 }*/
			assert(i < 10);
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

bool solve2x2LinearSystem(const double *A, double *b, double EPS) {
// A is column major
	double detA = A[0] * A[3] - A[2] * A[1];
	double det0 = A[3] * b[0] - A[2] * b[1];
	double det1 = A[0] * b[1] - A[1] * b[0];
	if (fabs(detA) < EPS * fabs(det0))
		return false;
	if (fabs(detA) < EPS * fabs(det1))
		return false;
	if (detA == 0) // det0 and det1 could be 0!!!
		return false; // without this, it could fail!!!!!
	b[0] = det0 / detA;
	b[1] = det1 / detA;
	return true;
}

double areaTriangle(double x1, double y1, double z1, double x2, double y2,
		double z2) {
	double x = y1 * z2 - y2 * z1;
	double y = z1 * x2 - z2 * x1;
	double z = x1 * y2 - x2 * y1;
	return sqrt(x * x + y * y + z * z) / 2;

}

double distance(double* x1, double* x2){
	return sqrt((x1[0]-x2[0]) * (x1[0]-x2[0]) + (x1[1]-x2[1]) * (x1[1]-x2[1]) + (x1[2]-x2[2]) * (x1[2]-x2[2]));
}


} // end of IGAMortarMath
} // end of EMPIRE

