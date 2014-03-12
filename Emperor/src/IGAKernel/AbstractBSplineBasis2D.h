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
 * \file AbstractBSplineBasis2D.h
 * This file holds the class AbstractBSplineBasis2D.h
 * \date 10/3/2013
 **************************************************************************************************/

#ifndef ABSTRACTBSPLINEBASIS2D_H_
#define ABSTRACTBSPLINEBASIS2D_H_

// Inclusion of user defined libraries
#include "AbstractIGABasis.h"

namespace EMPIRE {

/********//**
 * \brief class AbstractBSplineBasis2D is used as a base class for the BSplineBasis2D class
 ***********/

class AbstractBSplineBasis2D: public AbstractIGABasis {

    /// The constructor, the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the basis
     * \author Andreas Apostolatos
     ***********/
    AbstractBSplineBasis2D(int _ID = 0) :
            AbstractIGABasis(_ID) {
    }

    /***********************************************************************************************
     * \brief Virtual destructor
     * \author Andreas Apostolatos
     ***********/
    ~AbstractBSplineBasis2D() {
    }

    /***********************************************************************************************
     * \brief Virtual copy constructor
     * \param[in] _abstractBSplineBasis1D Constant reference to an object of class AbstractBSplineBasis1D
     * \author Andreas Apostolatos
     ***********/
    AbstractBSplineBasis2D(const AbstractBSplineBasis2D& _abstractBSplineBasis2D) :
            AbstractIGABasis(_abstractBSplineBasis2D) {
    }

    /// Functions related to the basis functions
public:
    /***********************************************************************************************
     * \brief Compute the non-zero basis functions at the given parameter
     * \param[in/out] _basisFcts The non-zero basis functions at the given parameter
     * \param[in] _uPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives in
     * \param[in] _vPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    virtual void computeLocalBasisFunctions(double*, double, int, double, int) = 0;

    /***********************************************************************************************
     * \brief Compute the non-zero basis functions and their derivatives of the  at the given parameter and stores them in a double array
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span and their _derivDegree-th derivatives
     * \param[in] _derivDegree The maximum partial derivative degree with 0 <= _uDerivDegree + _vDerivDegree <= _derivDegree
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated in u-direction
     * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives
     * \param[in] _vPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives
     * \author Andreas Apostolatos
     ***********/
    virtual void computeLocalBasisFunctionsAndDerivatives(double*, int, double, int, double,
            int) = 0;
};

}/* namespace EMPIRE */

#endif /* ABSTRACTBSPLINEBASIS1D_H_ */
