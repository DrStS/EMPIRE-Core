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
/******************************************************************************//**
 * \file IGAMortarMath.h
 * The header file of class IGAMortarMath.
 * \date 3/6/2013
 *********************************************************************************/

#ifndef IGAMORTARMATH_H_
#define IGAMORTARMATH_H_

#include <vector>

namespace EMPIRE {

namespace MortarMath {
class PolygonClipper;
}

namespace IGAMortarMath {

/// Tolerance up to which square matrices are assumed regular
extern const double EPS_IVERTIBILITYOFSQUAREMATRICES;

/********//**
 * \brief Class GaussQuadrature base class for general quadrature
 * \author Chenshen Wu
 ***********/
class GaussQuadrature {

public:
    // The number of Gauss Points
    int numGaussPoints;

    // Array containing the Gauss Point locations in the quadrature space
    const double *gaussPoints;

    // Array containing the quadrature weights
    const double *weights;

    /// Constructor, destructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numGaussPoints, number of Gauss points
     * \author Chenshen Wu
     ***********/
    GaussQuadrature(int _numGaussPoints) :
            numGaussPoints(_numGaussPoints) {
    }

    virtual ~GaussQuadrature() {
    }

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Returns the coordinates of the Gauss point
     * \param[in] _index the Gauss point to be returned
     * \author Chenshen Wu
     ***********/
    virtual const double* getGaussPoint(int _index) = 0;
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
    
    /***********************************************************************************************
     * \brief Returns the coordinates of the Gauss point
     * \param[in] _index the Gauss point to be returned
     * \author Chenshen Wu
     ***********/
    const double* getGaussPoint(int _index) {
        return &gaussPoints[_index * 2];
    }

};

/**********
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
    
    /***********************************************************************************************
     * \brief Returns the coordinates of the Gauss point
     * \param[in] _index the Gauss point to be returned
     * \author Chenshen Wu
     ***********/
    const double* getGaussPoint(int _index) {
        return &gaussPoints[_index * 2];
    }
};

/***********************************************************************************************
 * \brief This class is used to clip fluid element on the IGA patch. The fluid nodes are already
 *  				projected to the IGA in the previous step.
 *        The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice,
 *        2nd edition, Addison-Wesley" P237 - P240.
 * ***********/
class IGAPolygonClipper {

private:
    /// Use the polygon clipper in MortarMath temporarily
    MortarMath::PolygonClipper *clipper;

    /// Instersection area between the polygon and window
    double polygonWindow[12];

public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _u1, _u2, _v1, _v2 define the IGA element that clips other the fluid element.
     * \author Chenshen Wu
     ***********/
    IGAPolygonClipper(const double _u1, const double _u2, const double _v1, const double _v2);
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
bool solve2x2LinearSystem(const double* _A, double* _b, double _EPS);

/***********************************************************************************************
 * \brief Returns the scalar/vector value from the linear combination of the values on the vertices of the element with the shape functions
 * \param[in] _nNodes The number of nodes of the element
 * \param[in] _nValue Takes value 1 for a scalar, or the values 2-3 for a vector
 * \param[in] _values The values on the nodes of the element
 * \param[in] _coords The coordinates in the element
 * \param[in/out] _returnValue The resulting linear combination of the nodal values
 * \author Chenshen Wu
 ***********/
void computeLinearCombinationValueFromVerticesValues(int _nNodes, int _nValue, const double *_values,
        const double* _coords, double *_returnValue);
//computeValuesInLowOrderElement

/***********************************************************************************************
 * \brief Computes the value of a data field in the interior of an element
 * \param[in] _nNodes The number of nodes of the element
 * \param[in] _nValue Takes value 1 for a scalar, or the values 2-3 for a vector
 * \param[in] _values The values on the nodes of the element
 * \param[in] _shapeFuncs The values of the shape functions at the interior of the element
 * \param[in/out] _returnValue The resulting linear combination of the nodal values
 * \author Chenshen Wu
 ***********/
void computeLinearCombination(int _nNodes, int _nValue, const double * _values,
        const double *_shapeFuncs, double* _returnValue);

/***********************************************************************************************
 * \brief Computes the values of the low order shape functions (linear for triangle and bilinear
 *        for the quadrilateral)
 * \param[in] _nNodes The number of nodes in the element level
 * \param[in] _coords The coordinates of the point where to evaluate the shape functions
 * \param[in/out] _shapeFuncs The evaluated shape functions
 ***********/
void computeLowOrderShapeFunc(int _nNodes, const double* _coords, double* _shapeFuncs);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle in a 2D space
 * \param[in] _coordsTriangle, coordinates of the triangle. double[6].
 * \param[in] _coordsNode, coordinates of the point. double[2].
 * \param[out] _localCoords, local coordinates of the point. double[3]
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInTriangle(const double *_coordsTriangle, const double *_coordsNode,
        double* _localCoords);

/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadriliteral in a 2D space by solving a nonlinear
 *        system using the Newton-Raphson scheme
 * \param[in] _coordsQuad Coordinates of the quadriliteral. double[8].
 * \param[in] _coordsNode Coordinates of the point. double[2].
 * \param[out] _localCoords local coordinates of the point. double[2]
 * \return a boolean saying whether the point is inside the quadriliteral or not (true of false)
 * \author Chenshen Wu
 ***********/
bool computeLocalCoordsInQuad(const double *_coordsQuad, const double *_coordsNode,
        double* _localCoords);


/***********************************************************************************************
 * \brief Returns the intersection between a line and a triangle
 * \param[out] Flag whether the intersection has been found
 * \param[in] _coordsTriangle Coordinates of the triangle
 * \param[in] _coordsNode Coordinates of the point
 * \param[in] _direction Direction of the line
 * \param[in/out] _localCoords Local coordinates of the point
 * \author Chenshen Wu
 ***********/
bool computeIntersectionBetweenLineAndTriangle(const double *_coordsTriangle, const double* _coordsNode, const double* _direction, double* _localCoords);
//computeIntersectionInTri

/***********************************************************************************************
 * \brief Returns the intersection between a line and a quadrilateral
 * \param[out] Flag whether the intersection has been found
 * \param[in] _coordsQuad Coordinates of the quadrilateral
 * \param[in] _coordsNode Coordinates of the point
 * \param[in] _direction Direction of the line
 * \param[in/out] _localCoords Local coordinates of the point
 * \author Chenshen Wu
 ***********/
bool computeIntersectionBetweenLineAndQuad(const double *_coordsQuad, const double* _coordsNode, const double* _direction, double* _localCoords);
//computeIntersectionInQuad

/***********************************************************************************************
 * \brief Computes the area of a given triangle defined by two vectors
 * \param[in] _x1, _y1, _z1 Vector 1
 * \param[in] _x2, _y2, _z2 Vector 2
 * \param[out] The area of the triangle included by the two vectors
 * \author Chenshen Wu
 ***********/
double computeAreaTriangle(double x1, double y1, double z1, double x2, double y2, double z2);

/***********************************************************************************************
 * \brief Solve 3x3 Linear system, Ax = b
 * \param[in] _A Square 3x3 matirx
 * \param[in/out] _b Right hand side vector which stores also the solution
 * \param[in] _EPS tolerance up to which matrix _A is regular
 * \param[out] Flag on the solvability of the system
 * \author Chenshen Wu
 ***********/
bool solve3x3LinearSystem(const double* _A, double* _b, double _EPS);

/***********************************************************************************************
 * \brief Computes the Determinant of a given 3x3 matrix
 * \param[in] _A Matrix
 * \param[out] The Determinant of the given matrix
 * \author Chenshen Wu
 ***********/
double det3x3(const double* _A);

/***********************************************************************************************
 * \brief Computes the distance between two points
 * \param[in] _x1 Point 1
 * \param[in] _x2 Point 2
 * \param[out] The distance between the two given points
 * \author Chenshen Wu
 ***********/
double computePointDistance(double* _x1, double* _x2);

/***********************************************************************************************
 * \brief Computes the cross product of between two vectors with zero component on z direction
 * \param[in] _x1, x component of the first vector
 * \param[in] _y1, y component of the first vector
 * \param[in] _x2, x component of the second vector
 * \param[in] _y2, y component of the second vector
 * \param[out] Z component of the cross product between (x1,y1,0) and (x2,y2,0)
 * \author Chenshen Wu
 ***********/
double computeCrossProduct2D(double x1, double y1, double x2, double y2);

/***********************************************************************************************
 * \brief Clean the input polygon, remove colinear or identical point
 * \param[in/out] _polygon, vector of points coordinates
 * \author Fabien Pean
 ***********/
void cleanPolygon(std::vector<double>& _polygon);
void cleanPolygon(std::vector<std::pair<double,double> >& polygon);
void cleanPolygon(std::vector<std::pair<double,double> >& polygon,std::vector<std::pair<double,double> >& polygonSlave);
}
}

#endif /* IGAMORTARMATH_H_ */
