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
/******************************************************************************//**
 * \file MortarMath.h
 * The header file of math functions for mortar mapping.
 * \date 1/29/2013
 *********************************************************************************/
#ifndef MORTARMATH2_H_
#define MORTARMATH2_H_

#include <vector>

namespace EMPIRE {

namespace MortarMath {
//##############################################################################################
//############################ Gauss quadrature functions ######################################
//##############################################################################################
/// 3 Gauss points of a triangle
extern const double triGaussPoints3[9];
/// weights on 3 Gauss points of a triangle
extern const double triWeights3[3];
/// 6 Gauss points of a triangle
extern const double triGaussPoints6[18];
/// weights on 6 Gauss points of a triangle
extern const double triWeights6[6];
/// 7 Gauss points of a triangle
extern const double triGaussPoints7[21];
/// weights on 7 Gauss points of a triangle
extern const double triWeights7[7];
/// 12 Gauss points of a triangle
extern const double triGaussPoints12[36];
/// weights on 12 Gauss points of a triangle
extern const double triWeights12[12];
/// 1 Gauss point of a quadrilateral
extern const double quadGaussPoints1[2];
/// weights on 1 Gauss point of a quadrilateral
extern const double quadWeights1[1];
/// 4 Gauss points of a quadrilateral
extern const double quadGaussPoints4[8];
/// weights on 4 Gauss points of a quadrilateral
extern const double quadWeights4[4];
/// 9 Gauss points of a quadrilateral
extern const double quadGaussPoints9[18];
/// weights on 9 Gauss points of a quadrilateral
extern const double quadWeights9[9];

/********//**
 * \brief Class IntegrandFunction is the mother class of all integrand functions
 * \author Tianyang Wang
 ***********/
class IntegrandFunction {
public:
    IntegrandFunction() {
    }
    virtual ~IntegrandFunction() {
    }
    /***********************************************************************************************
     * \brief Compute the function value on the Gauss point
     * \param[in] gaussPoint x,y,z coordinates of the Gauss point
     * \return the function value
     ***********/
    virtual double operator()(double *gaussPoint) = 0;
};

/********//**
 * \brief Class GaussQuadratureOnTriangle performs Gauss quadrature on triangle
 * \author Tianyang Wang
 ***********/
class GaussQuadratureOnTriangle {
public:
    GaussQuadratureOnTriangle(double *_triangle, int _numGaussPoints);
    virtual ~GaussQuadratureOnTriangle();
    void setIntegrandFunc(IntegrandFunction *_integrandFunc);
    double computeIntegral();
    const int numGaussPoints;
    double *gaussPointsGlobal;
private:
    const double *triangle;
    const double *gaussPointsLocal;
    const double *weights;
    IntegrandFunction *integrandFunc;
    double area;
};

/********//**
 * \brief Class GaussQuadratureOnQuad performs Gauss quadrature on quad
 * \author Tianyang Wang
 ***********/
class GaussQuadratureOnQuad {
public:
    GaussQuadratureOnQuad(double *_quad, int _numGaussPoints);
    virtual ~GaussQuadratureOnQuad();
    void setIntegrandFunc(IntegrandFunction *_integrandFunc);
    double computeIntegral();
    const int numGaussPoints;
    double *gaussPointsGlobal;
private:
    const double *quad;
    const double *gaussPointsLocal;
    const double *weights;
    double *detJ;
    IntegrandFunction *integrandFunc;
};

/***********************************************************************************************
 * \brief Compute mass matrix of a triangle element
 * \param[in] triangle the triangle
 * \param[in] numGaussPoints number of Gauss points used in the Gauss quadrature
 * \param[in] dual whether dual or not
 * \param[out] mass matrix (3x3)
 * \author Tianyang Wang
 ***********/
void computeMassMatrixOfTrianlge(const double *triangle, int numGaussPoints, bool dual,
        double *massMatrix);
/***********************************************************************************************
 * \brief Compute mass matrix of a quad element
 * \param[in] quad the quad
 * \param[in] numGaussPoints number of Gauss points used in the Gauss quadrature
 * \param[in] dual whether dual or not
 * \param[out] mass matrix (4x4)
 * \author Tianyang Wang
 ***********/
void computeMassMatrixOfQuad(const double *quad, int numGaussPoints, bool dual, double *massMatrix);

//##############################################################################################
//############################ FE functions ####################################################
//##############################################################################################

/***********************************************************************************************
 * \brief Compute the shape function value by local coordinates in a quadrilateral
 * \param[in] xi_eta local coordinates
 * \param[out] shapeFuncValues shape function values of xi_eta
 * \author Tianyang Wang
 ***********/
void computeShapeFuncOfQuad(const double *xi_eta, double *shapeFuncValues);
/***********************************************************************************************
 * \brief Compute the determinant of Jocobian matrix by local coordinates in a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] xi_eta local coordinates
 * \return determinant of Jocobian matrix
 * \author Tianyang Wang
 ***********/
double computeDetJOfQuad(const double *quad, const double *xi_eta);
/***********************************************************************************************
 * \brief Compute global coordinates of a point in a triangle
 * \param[in] triangle the triangle
 * \param[in] localCoor local coordinates of the point
 * \param[out] globalCoor global coordinates of the point
 * \author Tianyang Wang
 ***********/
void computeGlobalCoorInTriangle(const double *triangle, const double *localCoor,
        double *globalCoor);
/***********************************************************************************************
 * \brief Compute global coordinates of a point in a quad
 * \param[in] quad the quad
 * \param[in] localCoor local coordinates of the point
 * \param[out] globalCoor global coordinates of the point
 * \author Tianyang Wang
 ***********/
void computeGlobalCoorInQuad(const double *quad, const double *localCoor, double *globalCoor);
/***********************************************************************************************
 * \brief Compute local coordinates of a point in a triangle
 * \param[in] triangle the triangle
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[in] point the point
 * \param[out] localCoor local coordinates of the point
 * \return a boolean saying whether the point is inside the triangle or not (true of false)
 * \author Tianyang Wang
 ***********/
bool computeLocalCoorInTriangle(const double *triangle, int planeToProject, const double *point,
        double *localCoor);
/***********************************************************************************************
 * \brief Compute local coordinates of a point in a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[in] point the point
 * \param[out] localCoor local coordinates of the point
 * \return a boolean saying whether the point is inside the quadrilateral or not (true of false)
 * \author Tianyang Wang
 ***********/
bool computeLocalCoorInQuad(const double *quad, int planeToProject, const double *point,
        double *localCoor);
//##############################################################################################
//############################ Geometric functions #############################################
//##############################################################################################

/***********************************************************************************************
 * \brief This class is used to cooperate with the clipping algorithm
 *        The algorithm is from the book "J.Foley, Computer Graphics, Principles and Practice,
 *        2nd edition, Addison-Wesley" P237 - P240.
 *        In practice this algorithm can have any strictly convex polygon as clipping window.
 *        This means consecutive points are not authorized to be colinear
 * ***********/
class PolygonClipper {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _polygonWindow the polygon that clips other polygons
     * \param[in] _sizePolygonWindow number of nodes/edges in this polygon
     * \param[in] _planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
     * \author Tianyang Wang
     ***********/
    PolygonClipper(const double *_polygonWindow, int _sizePolygonWindow, int _planeToProject);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~PolygonClipper();
    /***********************************************************************************************
     * \brief Clip.
     * \param[in] polygonToBeClipped polygon to be clipped by the window polygon
     * \param[in] sizePolygonToBeClipped size of the polygon
     * \param[out] the polygonResult results from the clipping
     * \return if there is a clip, return true, otherwise, return false
     * \author Tianyang Wang
     ***********/
    bool clip(const double *polygonToBeClipped, int sizePolygonToBeClipped,
            std::vector<double*> *polygonResult);
    /***********************************************************************************************
     * \brief Compute intersection between two lines. The result is the intersection with tolerance,
     *        which means it may locate outside of both lines. This error should be taken into account
     *        in the future. Make it static so that it is easily tested by some unit test.
     * \param[in] la0 the 1st point on la
     * \param[in] la1 the 2nd point on la
     * \param[in] lb0 the 1st point on lb
     * \param[in] lb1 the 2nd point on lb
     * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
     * \param[out] intersection the intersection of la and lb
     * \return boolean saying whether the intersection is on lb or not (true of false)
     * \author Tianyang Wang
     ***********/
    static bool intersect(const double *la0, const double *la1, const double *lb0, const double *lb1,
            int planeToProject, double *intersection);

private:
    /// the polygon that clips other polygons
    double *polygonWindow;
    /// number of edges of the polygon
    int sizePolygonWindow;
    /// project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
    const int planeToProject;
    /// flags used to determine inside/outside
    bool *insideFlag;
    /// slopes of all edges
    double *edgeSlopes;
    /// sign of whether to use y=f(x) or x=f(y)
    bool *reverseXY;
    // put the magic numbers here, if the number is too small, numerical error is high; if the number is
    // too large, computation error is high. Therefore, set the suitable value by experiment
    /// epsilon 1
    static const double EPS1;
    /// epsilon 2
    static const double EPS2;
    /// unit test class
    friend class TestMortarMath;
    /***********************************************************************************************
     * \brief Decides whether a point is on the "inside" side of the i-th edge of the polygon
     * \param[in] edgeID id of the edge
     * \param[in] point the point
     * \return if inside, return true, otherwise, return false
     * \author Tianyang Wang
     ***********/
    bool inside(int edgeID, double *point);
};

/***********************************************************************************************
 * \brief project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \brief The result is the plane which has smallest angle with the unitNormal
 * \param[in] planeNormal normal vector of the plane
 * \return plane to project (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \author Tianyang Wang
 ***********/
int computePlaneToProject(const double *planeNormal);
/***********************************************************************************************
 * \brief Compute the center of a polygon
 * \param[in] polygon the polygon
 * \param[in] num num of points of the polygon
 * \param[out] center the center
 * \author Tianyang Wang
 ***********/
void computePolygonCenter(const double *polygon, int num, double *center);
/***********************************************************************************************
 * \brief Compute polygon area (points should be on the same plane)
 * \param[in] polygon the polygon
 * \param[in] num num of points of the polygon
 * \return the area of the polygon
 * \author Tianyang Wang
 ***********/
double computePolygonArea(const double *polygon, int num);
/***********************************************************************************************
 * \brief Project a number of points to a plane
 * \param[in] pointOnPlane a point on the plane
 * \param[in] unitNormal unit normal of the plane
 * \param[in] points the points to be projected
 * \param[in] num number of points to be projected
 * \param[in] projections the projections
 * \author Tianyang Wang
 ***********/
void projectToPlane(const double *pointOnPlane, const double *unitNormal, const double *points,
        int num, double *projections);
/***********************************************************************************************
 * \brief Compute the normal vector of a triangle
 * \param[in] triangle the triangle
 * \param[in] unitLength whether set the length of the normal vector to 1 or not
 * \param[out] normal normal vector
 * \author Tianyang Wang
 ***********/
void computeNormalOfTriangle(const double *triangle, bool unitLength, double *normal);
/***********************************************************************************************
 * \brief Compute the normal vector of a quadrilateral
 * \param[in] quad the quadrilateral
 * \param[in] unitLength whether set the length of the normal vector to 1 or not
 * \param[out] normal normal vector
 * \author Tianyang Wang
 ***********/
void computeNormalOfQuad(const double *quad, bool unitLength, double *normal);
/***********************************************************************************************
 * \brief Compute the area of a triangle
 * \param[in] triangle the triangle
 * \return the area
 * \author Tianyang Wang
 ***********/
double computeAreaOfTriangle(const double *triangle);
/***********************************************************************************************
 * \brief Compute the length of a vector
 * \param[in] vec the vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength(const double *vec);
/***********************************************************************************************
 * \brief Compute the length of a 2D vector
 * \param[in] vec the 2D vector
 * \return the vector length
 * \author Tianyang Wang
 ***********/
double computeVectorLength2D(const double *vec);
/***********************************************************************************************
 * \brief Compute the cross product of two vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \param[out] crossProduct the crossProduct
 * \author Tianyang Wang
 ***********/
void computeVectorCrossProduct(const double *vec1, const double *vec2, double *crossProduct);
/***********************************************************************************************
 * \brief Compute the dot product of two vectors
 * \param[in] vec1 the 1st vector
 * \param[in] vec2 the 2nd vector
 * \return dot product
 * \author Tianyang Wang
 ***********/
double computeVectorDotProduct(const double *vec1, const double *vec2);
/***********************************************************************************************
 * \brief Compute the square of the distance between two points
 * \param[in] p1 the 1st point
 * \param[in] p2 the 2nd point
 * \return the square of the distance between p1 and p2
 * \author Tianyang Wang
 ***********/
double distanceSquare(const double *p1, const double *p2);
/***********************************************************************************************
 * \brief Copy a point
 * \param[in] origin to be copied
 * \param[out] copy the copy
 * \author Tianyang Wang
 ***********/
void copyPoint(const double *origin, double *copy);
/***********************************************************************************************
 * \brief Copy a point
 * \param[in] origin to be copied
 * \param[out] copy the copy
 * \author Tianyang Wang
 ***********/
void copyElem(const double *origin, int size, double *copy);
/***********************************************************************************************
 * \brief Build a triangle
 * \param[in] p0 the 1st point
 * \param[in] p1 the 2nd point
 * \param[in] p2 the 3rd point
 * \param[out] triangle the triangle built by the three pointss
 * \author Tianyang Wang
 ***********/
void buildTrianagle(const double *p0, const double *p1, const double *p2, double *triangle);
/***********************************************************************************************
 * \brief Compute the longest edge length of the element
 * \param[in] elem the element
 * \param[in] size 3--->triangle, 4--->quadrilateral
 * \return the square of the longest length of the element
 * \author Tianyang Wang
 ***********/
double longestEdgeLengthSquare(const double *elem, int size);

//##############################################################################################
//######################### Linear system solvers and linear algebra functions #################
//##############################################################################################

/***********************************************************************************************
 * \brief Solve a 2x2 linear system by close form formula
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[out] b the solution is written to b
 * \return whether the determinant is zero or not
 * \author Tianyang Wang
 ***********/
bool solve2x2LinearSystem(const double *A, double *b, double EPS);
/***********************************************************************************************
 * \brief Solve a 3x3 linear system by close form formula, the system has one row which has all entries 1.
 * \brief Therefore, it cannot solve general 3x3 linear system.
 * \param[in] A the left hand side matrix
 * \param[in] b the right hand side vector
 * \param[in] planeToProject project to plane (case: {2:x-y ; 0:y-z ;1: z-x} )
 * \param[out] b the solution is written to b
 * \author Tianyang Wang
 ***********/
void solve3x3LinearSystem(const double *A, int planeToProject, double *b);
/***********************************************************************************************
 * \brief Compute the matrix product between two matrices (for mortar, not so general)
 * \param[in] n n
 * \param[in] m m
 * \param[in] A the matrix with size n*n
 * \param[in] B the matrix with size n*m
 * \param[out] B B is overwritten by A*B (n*m)
 * \author Tianyang Wang
 ***********/
void computeMatrixProduct(int n, int m, const double *A, double *B);
/***********************************************************************************************
 * \brief My own sparse matrix (csr format, non-symmetric) vector multiplication routine (A * x = y).
 *        The interface is simplified and compatible with mkl_dcsrmv.
 *        This routine supports only one-based indexing of the input arrays.
 * \param[in] trans 'N' --- no transpose, 'T' --- transpose
 * \param[in] numRows number of rows of matrix A
 * \param[in] numCols number of columns of matrix A
 * \param[in] A sparse matrix A in CSR format
 * \param[in] JA JA of A
 * \param[in] IA IA of A
 * \param[in] x vector x
 * \param[in] y vector y
 * \author Tianyang Wang
 ***********/
void dcsrmv(char trans, int numRows, int numCols, const double *A, const int *JA, const int *IA,
        const double *x, double *y);
/***********************************************************************************************
 * \brief My own sparse matrix (csr format, symmetric) vector multiplication routine (A * x = y).
 *        The interface is simplified and compatible with mkl_dcsrsymv.
 *        This routine supports only one-based indexing of the input arrays.
 * \param[in] n size of matrix A
 * \param[in] A sparse matrix A in CSR format (symmetric)
 * \param[in] IA IA of A
 * \param[in] JA JA of A
 * \param[in] x vector x
 * \param[in] y vector y
 * \author Tianyang Wang
 ***********/
void dcsrsymv(int n, const double *A, const int *IA, const int *JA, const double *x, double *y);
//##############################################################################################
//######################### Printing functions #################################################
//##############################################################################################

/***********************************************************************************************
 * \brief Print the coordinates of an element, used in debugging
 * \param[in] elem the element
 * \param[in] size 3---triangle, 4---quadrilateral
 * \author Tianyang Wang
 ***********/
void printElem(const double *elem, int size);
/***********************************************************************************************
 * \brief Print the coordinates of a point, used in debugging
 * \param[in] p the point
 * \author Tianyang Wang
 ***********/
void printPoint(const double *p);
/***********************************************************************************************
 * \brief print a diagonal matrix
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \author Tianyang Wang
 ***********/
void printDiagonalMatrix(const double *A, int numRows);
/***********************************************************************************************
 * \brief print a general matrix
 * \param[in] A the matrix
 * \param[in] numRows number of rows of A
 * \param[in] numCols number of columns of A
 * \author Tianyang Wang
 ***********/
void printGeneralMatrix(const double *A, int numRows, int numCols);
/***********************************************************************************************
 * \brief print an unsymmetric CSR formatted matrix
 * \param[in] A A
 * \param[in] IA IA
 * \param[in] JA JA
 * \param[in] numRows number of rows of A
 * \param[in] numCols number of columns of A
 * \author Tianyang Wang
 ***********/
void printCSRMatrixUnsymmetric(const double *A, const int *IA, const int *JA, int numRows,
        int numCols);
/***********************************************************************************************
 * \brief print a symmetric CSR formatted matrix
 * \param[in] A A
 * \param[in] IA IA
 * \param[in] JA JA
 * \param[in] n nxn matrix
 * \author Tianyang Wang
 ***********/
void printCSRMatrixSymmetric(const double *A, const int *IA, const int *JA, int n);

} /* namespace MortarMath */
} /* namespace EMPIRE */
#endif /* MORTARMATH2_H_ */
