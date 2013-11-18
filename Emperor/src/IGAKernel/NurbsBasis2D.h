/***********************************************************************************************//**
 * \file NurbsBasis2D.h
 * This file holds the class BSplineBasis2D.h
 * \date 11/3/2013
 **************************************************************************************************/

#ifndef NURBSBASIS2D_H_
#define NURBSBASIS2D_H_

// Inclusion of user defined libraries
#include "BSplineBasis2D.h"
#include "IGAControlPoint.h"

namespace EMPIRE {

/********//**
 * \brief class BSplineBasis2D is related to the BSpline basis functions in 2D case
 ***********/

class NurbsBasis2D: public BSplineBasis2D {

protected:
    /// Number of Control Points in u-direction
    int uNoBasisFnc;

    /// Number of Control Points in u-direction
    int vNoBasisFnc;

    /// The Control Point net
    double* IGAControlPointWeights;

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
     * \param[in] _uNoBasisFnc The number of Control Points/Basis functions in u-direction
     * \param[in] _vNoBasisFnc The number of Control Points/Basis functions in v-direction
     * \param[in] _igaControlPointWeights The Control Point weights for the 2D NURBS basis
     * \author Andreas Apostolatos
     ***********/
    NurbsBasis2D(int, int, int, double*, int, int, double*, int, int, double*);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~NurbsBasis2D();

    /***********************************************************************************************
     * \brief Copy constructor
     * \param[in] _nurbsBasis2D Constant reference to an object of class NurbsBasis2D
     * \author Andreas Apostolatos
     ***********/
    NurbsBasis2D(const NurbsBasis2D&);

    /***********************************************************************************************
     * \brief Copy assignment
     * \author Andreas Apostolatos
     ***********/
    // NurbsBasis2D& operator=(const NurbsBasis2D&);
    /// On the 2D B-Spline basis functions
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
    void computeLocalBasisFunctions(double*, double, int, double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero basis functions and their derivatives of the  at the given parameter (Inefficient algorithm)
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span and their _derivDegree-th derivatives
     * \param[in] _maxMixDerivOrd The maximum degree of the mixed derivative
     * \param[in] _derivDegreeU Up to Which degree the B-Spline derivative to be computed in u-direction
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated in u-direction
     * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives
     * \param[in] _derivDegreeV Up to Which degree the B-Spline derivative to be computed in v-direction
     * \param[in] _vPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivativesInefficient(double**, int, int, double, int, int,
            double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero basis functions and their derivatives of the  at the given parameter and stores them in a double array but
     *        can still handle up to 3rd order derivatives and the algorithm is not adapted by the NURBS book
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span and their _derivDegree-th derivatives
     * \param[in] _maxMixDerivOrd The maximum degree of the mixed derivative
     * \param[in] _derivDegreeU Up to Which degree the B-Spline derivative to be computed in u-direction
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated in u-direction
     * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives
     * \param[in] _derivDegreeV Up to Which degree the B-Spline derivative to be computed in v-direction
     * \param[in] _vPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivativesInefficient(double*, int, int, double, int, int,
            double, int);

    /***********************************************************************************************
     * \brief Compute the denominator function w(u,v) = Sum_(i=1)^n N_i(u,v) * w_i and its derivatives at the given surface parameters u and v
     * \param[in/out] _denominatorFctAndDerivs The denominator function and its derivatives at (u,v)
     * \param[in] _bSplineBasisFctsAndDerivs The non zero B-Spline basis functions at the given knot span and their up to their _derivDegree-th derivatives at u
     * \param[in] _derivDegree For which degree the B-Spline derivative to be computed
     * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives in
     * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeDenominatorFunctionAndDerivatives(double*, double*, int, int, int);

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter and stores them in a double array
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span and their _derivDegree-th derivatives
     * \param[in] _derivDegree The maximum partial derivative degree with 0 <= _uDerivDegree + _vDerivDegree <= _derivDegree
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated in u-direction
     * \param[in] _KnotSpanIndexU The index of the knot span where _uPrm lives
     * \param[in] _vPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndexV The index of the knot span where _vPrm lives
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivatives(double*, int, double, int, double, int);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the Control Point net
     * \author Andreas Apostolatos
     ***********/
    int getUNoBasisFnc() {
        return uNoBasisFnc;
    }

    /***********************************************************************************************
     * \brief Get the Control Point net
     * \author Andreas Apostolatos
     ***********/
    int getVNoBasisFnc() {
        return vNoBasisFnc;
    }

    /***********************************************************************************************
     * \brief Get the Control Point net
     * \author Andreas Apostolatos
     ***********/
    double* getIGAControlPointWeights() {
        return IGAControlPointWeights;
    }

};

Message &operator<<(Message &message, NurbsBasis2D &nurbsBasis2D);

}/* namespace EMPIRE */

#endif /* BSPLINEBASIS2D_H_ */
