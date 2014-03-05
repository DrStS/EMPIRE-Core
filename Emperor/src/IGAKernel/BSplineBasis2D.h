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
/***********************************************************************************************//**
 * \file BSplineBasis2D.h
 * This file holds the class BSplineBasis2D.h
 * \date 11/3/2013
 **************************************************************************************************/

#ifndef BSPLINEBASIS2D_H_
#define BSPLINEBASIS2D_H_

// Inclusion of user defined libraries
#include "AbstractBSplineBasis2D.h"
#include "BSplineBasis1D.h"
#include "Message.h"

namespace EMPIRE {

/********//**
 * \brief class BSplineBasis2D is related to the BSpline basis functions in 2D case
 ***********/

class BSplineBasis2D: public AbstractBSplineBasis2D {

protected:
	/// The B-Spline basis spanning the u-direction
	BSplineBasis1D* uBSplineBasis1D;

	/// The B-Spline basis spanning the v-direction
	BSplineBasis1D* vBSplineBasis1D;

	/// The constructor, the destructor, the copy constructor and the copy assignment
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _ID The id of the basis
	 * \param[in] _pDegree The polynomial degree of the basis in u-direction
	 * \param[in] _noKnotsU The number of knots in u-direction
	 * \param[in] _KnotVectorU The knot vector in u-direction
	 * \param[in] _qDegree The polynomial degree of the basis in v-direction
	 * \param[in] _noKnotsV The number of knots in v-direction
	 * \param[in] _KnotVectorV The knot vector in v-direction
	 * \author Andreas Apostolatos
	 ***********/
	BSplineBasis2D(int, int, int, double*, int, int, double*);

	/***********************************************************************************************
	 * \brief Destructor
	 * \author Andreas Apostolatos
	 ***********/
	~BSplineBasis2D();

	/***********************************************************************************************
	 * \brief Copy constructor
	 * \param[in] _bSplineBasis2D Constant reference to an object of class BSplineBasis2D
	 * \author Andreas Apostolatos
	 ***********/
	BSplineBasis2D(const BSplineBasis2D&);

	/***********************************************************************************************
	 * \brief Copy assignment
	 * \author Andreas Apostolatos
	 ***********/
	BSplineBasis2D& operator=(const BSplineBasis2D&);

	/// On the 2D B-Spline basis functions
public:
	/***********************************************************************************************
	 * \brief Compute the index related to the uDerivIndex-th partial derivative w.r.t to u, vDerivIndex-th partial derivative w.r.t. v of the basisIndex-th basis function
	 * \param[out] The index of the functional value in the 1D pointer array containing the basis functions and their derivatives
	 * \param[in] _derivDegree The absolute derivative degree of the basis functions
	 * \param[in] _uDerivIndex The order of the partial derivative w.r.t. the u-parametric direction
	 * \param[in] _vDerivIndex The order of the partial derivative w.r.t. the v-parametric direction
	 * \param[in] _basisIndex The index of the basis function in 2D sorted as in the documentation for the computation of the basis functions in 2D
	 * \author Andreas Apostolatos
	 ***********/
	int indexDerivativeBasisFunction(int, int, int, int);

	/***********************************************************************************************
	 * \brief Compute the non-zero basis functions at the given parameter
	 * \param[in/out] _basisFcts The non-zero basis functions at the given parameter
	 * \param[in] _uPrm The parameter where the basis functions are evaluated
	 * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives in
	 * \param[in] _vPrm The parameter where the basis functions are evaluated
	 * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives in
	 * \author Andreas Apostolatos
	 ***********/
	void computeLocalBasisFunctions(double*, double, int, double, int);

	/***********************************************************************************************
	 * \brief Compute the Index of non-zero basis functions at the given knot span
	 * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives in
	 * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives in
	 * \param[in/out] _funcsIndex The index of the basis function
	 * \author Andreas Apostolatos
	 ***********/
	void getBasisFunctionsIndex(int, int, int*);

	/***********************************************************************************************
	 * \brief Compute the non-zero B-Spline basis functions and their derivatives of the  at the given parameter and stores them in an 1D pointer array of doubles
	 * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span and their _derivDegree-th derivatives
	 * \param[in] _derivDegree The maximum partial derivative degree with 0 <= _uDerivDegree + _vDerivDegree <= _derivDegree
	 * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated in u-direction
	 * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives
	 * \param[in] _vPrm The parameter where the basis functions are evaluated
	 * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives
	 * \author Andreas Apostolatos
	 ***********/
	void computeLocalBasisFunctionsAndDerivatives(double*, int, double, int,
			double, int);

	/// Get and set functions
public:
	/***********************************************************************************************
	 * \brief Get the 1D B-Spline basis in the u-direction
	 * \author Andreas Apostolatos
	 ***********/
	BSplineBasis1D* getUBSplineBasis1D() {
		return uBSplineBasis1D;
	}

	/***********************************************************************************************
	 * \brief Get the 1D B-Spline basis in the v-direction
	 * \author Andreas Apostolatos
	 ***********/
	BSplineBasis1D* getVBSplineBasis1D() {
		return vBSplineBasis1D;
	}

};

/***********************************************************************************************
 * \brief Allows for nice debug output
 * \author Chenshen Wu
 ***********/
Message &operator<<(Message &message, BSplineBasis2D &bSplineBasis2D);

}/* namespace EMPIRE */

#endif /* BSPLINEBASIS2D_H_ */
