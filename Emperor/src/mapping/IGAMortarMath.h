/******************************************************************************//**
 * \file IGAMortarMath.h
 * The header file of class IGAMortarMath.
 * \date 3/6/2013
 *********************************************************************************/

#ifndef IGAMORTARMATH_H_
#define IGAMORTARMATH_H_

#include <vector>

namespace EMPIRE {

namespace MortarMath{
	class PolygonClipper;
}

namespace IGAMortarMath {

class GaussQuadrature {
public:
	GaussQuadrature(int _numGaussPoints) :
		numGaussPoints(_numGaussPoints) {
	}

	virtual ~GaussQuadrature(){}

	int numGaussPoints;
	const double *gaussPoints;
	const double *weights;
	virtual const double* getGaussPoint(int i) = 0;
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on triangle
 * \author Chenshen Wu
 ***********/
class GaussQuadratureOnTriangle: public GaussQuadrature {
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * param[in] _numGaussPoints, number of Gauss points
	 * \author Chenshen Wu
	 ***********/
	GaussQuadratureOnTriangle(int _numGaussPoints);
	virtual ~GaussQuadratureOnTriangle() {
	}
	;
	const double* getGaussPoint(int i) {
		return &gaussPoints[i * 3];
	}

};

/********//**
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Chenshen Wu
 ***********/
class GaussQuadratureOnQuad: public GaussQuadrature {
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * param[in] _numGaussPoints, number of Gauss points
	 * \author Chenshen Wu
	 ***********/
	GaussQuadratureOnQuad(int _numGaussPoints);

	virtual ~GaussQuadratureOnQuad() {
	}
	;
	const double* getGaussPoint(int i) {
		return &gaussPoints[i * 2];
	}

};

/***********************************************************************************************
 * \brief This class is used to clip fluid element on the IGA patch. The fluid nodes are already
 *  				projected to the IGA in the previous step.
 *        The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice,
 *        2nd edition, Addison-Wesley" P237 - P240.
 * ***********/
class IGAPolygonClipper {
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _u1, _u2, _v1, _v2 define the IGA element that clips other the fluid element.
	 * \author Chenshen Wu
	 ***********/
	IGAPolygonClipper(const double _u1, const double _u2, const double _v1,
			const double _v2);
	/***********************************************************************************************
	 * \brief Destructor
	 * \author Chenshen Wu
	 ***********/
	virtual ~IGAPolygonClipper();
	/***********************************************************************************************
	 * \brief Clip.
	 * \param[in] polygonToBeClipped polygon to be clipped by the window polygon
	 * \param[in] sizePolygonToBeClipped size of the polygon
	 * \param[out] the polygonResult results from the clipping
	 * \return if there is a clip, return true, otherwise, return false
	 * \author Tianyang Wang
	 ***********/
	bool clip(const double *_polygonToBeClipped, int _sizePolygonToBeClipped,
			std::vector<double*> *_polygonResult);

private:
	/// use the polygon clipper in MortarMath temporarily
	MortarMath::PolygonClipper *clipper;

	double polygonWindow[12];

	/// unit test class
	friend class TestIGAMortarMath;

};

/***********************************************************************************************
 * \brief Solve a 2x2 linear system by close form formula
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[out] b the solution is written to b
 * \return whether the determinant is zero or not
 * \author Tianyang Wang
 ***********/
bool solve2x2LinearSystem(const double *A, double *b, double EPS);

void calcValueFromShapeFuncs(const double *values, const double *_shapeFuncs,
		int _nNodes, int _nValue, double *returnValue);

void calcCoords(const double *_coordsQuadrature, const double *_gaussPoint,
		int _nNodes, int _nCoords, double *_coords);

void calcLinearShapeFunc(const double *_coords, int _nCoords,
		double *_shapeFuncs);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle in a 2D space
 * \param[in] _coordsTriangle, coordinates of the triangle. double[6].
 * \param[in] _coordsNode, coordinates of the point. double[2].
 * \param[out] _localCoords, local coordinates of the point. double[3]
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInTriangle(const double *_coordsTriangle,
		const double *_coordsNode, double* _localCoords);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadriliteral in a 2D space
 * \param[in] _coordsQuad, coordinates of the quadriliteral. double[8].
 * \param[in] _coordsNode, coordinates of the point. double[2].
 * \param[out] _localCoords, local coordinates of the point. double[2]
 * \return a boolean saying whether the point is inside the quadriliteral or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInQuad(const double *_coordsQuad,
		const double *_coordsNode, double* _localCoords);

/***********************************************************************************************
 * \brief Compute the area of a given triangle in a 2D space
 * \author Chenshen Wu
 ***********/
double areaTriangle(double x1, double y1, double z1, double x2, double y2,
		double z2);

double distance(double* _x1, double* _x2);

}
}

#endif /* IGAMORTARMATH_H_ */
