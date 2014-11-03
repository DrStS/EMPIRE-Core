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
 * \file IGAMath.h
 * The header file of math functions for mortar isogeometric mapping and for the basis functions
 * \date 1/29/2013
 *********************************************************************************/

#ifndef IGAMATH_H_
#define IGAMATH_H_

namespace EMPIRE {

/// The binomial coefficients stored in a vector format. The precomputed values are up to k=50
extern const double binomialCoefficients[2500];

/// The tolerance for the linear equation system solvers
extern const double EPS;

/***********************************************************************************************
 * \brief Computes the index function for the binomial coefficients (_i;_j)
 * \param[out] The index function for the binomial coefficients (_i;_j)
 * \param[in] _i The integer on the nominator
 * \param[in] _j The integer on the denominator
 * \author Andreas Apostolatos
 ***********/
int indexBinomialCoefficients(int, int);

/***********************************************************************************************
 * \brief Compute the square of the Euclidean distance of two points in n-D space.
 * \param[out] The square of the Euclidean distance between two points in the n-D space
 * \param[in] _length The length of the n-dimensional space
 * \param[in] _Pi The first point
 * \param[in] _Pj The second point
 * \author Andreas Apostolatos
 ***********/
double squareEuclideanDistance(int, double*, double*);

/***********************************************************************************************
 * \brief Compute the square of the 2-norm of a vector in the n-D space
 * \param[out] The square of the 2-norm of a vector in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vector A vector in the nD space
 * \author Andreas Apostolatos
 ***********/
double square2normVector(int, double*);

/***********************************************************************************************
 * \brief Compute the dot product between two vectors in the n-D space
 * \param[out] The dot product between two vectors in the n-D space
 * \param[in] _length The dimensinality of the n-D space
 * \param[in] _vecI The 1st vector
 * \param[in] _vecJ The 2nd vector
 * \author Andreas Apostolatos
 ***********/
double dotProduct(int, double*, double*);

/***********************************************************************************************
* \brief Compute the cross product between two vectors in the 3-D space
* \param[in/out] _product The product of vector1 and vector 2
* \param[in] _v1 The 1st vector
* \param[in] _v2 The 2nd vector
* \author Chenshen Wu
***********/
void crossProduct(double* _product, double* _v1, double* _v2);

/***********************************************************************************************

 * \brief Solve a 2x2 linear equation system
 * \param[out] Flag on whether the linear system is solvable up to tolerance EPS or not
 * \param[in/out] _b The right-hand side vector where the solution is also stored, _b = double[2]
 * \param[in] _A The 2x2 matrix stored in a vector format namely A[i][j] = V[2 * i + j]
 * \author Andreas Apostolatos
 ***********/
bool solve2x2linearSystem(double*, double*);

double distanceLinePlane(double* Pline,double* Uline, double* Pplane,double* Nplane);
double distancePointSegment(double* _P, double* _P1, double* _P2);
// See http://paulbourke.net/geometry/pointlineplane/
double distanceLineLine(double& _ratioA, double& _ratioB, double* _P1, double* _P2,double* P3, double* P4);


}/* namespace EMPIRE */

#endif /* IGAMATH_H_ */
