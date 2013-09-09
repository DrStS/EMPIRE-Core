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
/***********************************************************************************************//**
 * \file NurbsBasis1D.h
 * This file holds the class NurbsBasis1D.h
 * \date 12/3/2013
 **************************************************************************************************/

#ifndef NURBSBASIS1D_H_
#define NURBSBASIS1D_H_

// Inclusion of user defined libraries
#include "BSplineBasis1D.h"

namespace EMPIRE {

/********//**
 * \brief class NurbsBasis1D is related to the NURBS basis functions in 1D case
 ***********/

class NurbsBasis1D: public BSplineBasis1D {

protected:
    // The number of Control Points
    int NoControlPoints;

    // The sorted vector of the Control Point weights
    double* ControlPointWeights;

    /// The constructor and the destructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the basis
     * \param[in] _pDegree The polynomial degree of the basis
     * \param[in] _noKnots The number of knots
     * \param[in] _KnotVector The knot vector
     * \param[in] _noControlPoints The number of Control Points
     * \param[in] _controlPointWeights The Control Point weights for the 1D NURBS basis
     * \author Andreas Apostolatos
     ***********/
    NurbsBasis1D(int, int, int, double*, int, double*);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~NurbsBasis1D();

    /***********************************************************************************************
     * \brief Copy constructor (Problematic to be defined in the derived class, postponed for later if needed)
     * \author Andreas Apostolatos
     ***********/
    NurbsBasis1D(const NurbsBasis1D&);

    /// On the 1D NURBS basis functions
public:
    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions at the given parameter
     * \param[in/out] _localBasisFunctions The non-zero basis functions at the given parameter
     * \param[in] _uPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctions(double*, double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter (returns a 2D pointer array)
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivativesInefficient(double**, int, double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter and stores them in a vector of doubles
     *        however this is still not adapted by the NURBS book and handles only up to 2nd derivatives of the basis functions inefficiently
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivativesInefficient(double*, int, double, int);

    /***********************************************************************************************
     * \brief Compute the denominator function w(u) = Sum_(i=1)^n N_(i,p)*w_i and its derivatives at the given NURBS parameter u
     * \param[in/out] _denominatorFctAndDerivs The denominator function and its derivatives at u
     * \param[in] _bSplineBasisFctsAndDerivs The non zero B-Spline basis functions at the given knot span and their up to their _derivDegree-th derivatives at u
     * \param[in] _derivDegree For which degree the B-Spline derivative to be computed
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeDenominatorFunctionAndDerivatives(double*, double*, int, int);

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter (returns a 2D pointer array)
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivatives(double*, int, double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter (returns a 2D pointer array)-overloaded
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span and their up to their _derivDegree-th derivatives
     * \param[in] _bSplineBasisFctsAndDerivs The non -zero B-Spline basis functions at the given knot span and their up to their _derivDegree-th derivatives
     * \param[in] _denominatorFctAndDerivs The denominator function w(u) = Sum_(i=1)^n N_(i,p)*w_i plus up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivatives(double*, double*, double*, int, double, int);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the number of Control Points
     * \param[out] NoControlPoints
     * \author Andreas Apostolatos
     ***********/
    inline int getNoControlPoints() {
        return NoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the Control Point weights
     * \param[out] _controlPointWeights
     * \author Andreas Apostolatos
     ***********/
    double* getControlPointWeights();

    /***********************************************************************************************
     * \brief Sets the Control Point weights
     * \param[in] _noControlPoints The number of the Control Points
     * \param[in] _controlPointWeights The new Control Point weights to be assigned to the NURBS basis
     * \author Andreas Apostolatos
     ***********/
    void setControlPointWeights(int, double*);

    /// DEBUGGING functions
public:
    /***********************************************************************************************
     * \brief Outputs the Control Points weights of the net
     * \author Andreas Apostolatos
     ***********/
    void printControlPointWeights();
};

}/* namespace EMPIRE */

#endif /* NURBSBASIS1D_H_ */
