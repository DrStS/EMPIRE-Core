/***********************************************************************************************//**
 * \file AbstractBSplineBasis1D.h
 * This file holds the class AbstractBSplineBasis1D.h
 * \date 10/3/2013
 **************************************************************************************************/

#ifndef ABSTRACTBSPLINEBASIS1D_H_
#define ABSTRACTBSPLINEBASIS1D_H_

// Inclusion of user defined libraries
#include "AbstractIGABasis.h"

namespace EMPIRE {

/********//**
 * \brief class AbstractBSplineBasis1D is used as a base class for the BSplineBasis1D class
 ***********/

class AbstractBSplineBasis1D: public AbstractIGABasis {

    /// The constructor, the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the basis
     * \author Andreas Apostolatos
     ***********/
    AbstractBSplineBasis1D(int _ID) :
            AbstractIGABasis(_ID) {
    }

    /***********************************************************************************************
     * \brief destructor
     * \author Andreas Apostolatos
     ***********/
    ~AbstractBSplineBasis1D() {
    }

    /***********************************************************************************************
     * \brief Virtual copy constructor
     * \param[in] _abstractBSplineBasis1D Constant reference to an object of class AbstractBSplineBasis1D
     * \author Andreas Apostolatos
     ***********/
    AbstractBSplineBasis1D(const AbstractBSplineBasis1D& _abstractBSplineBasis1D) :
            AbstractIGABasis(_abstractBSplineBasis1D) {
    }

    /// Functions related to the basis functions
public:
    /***********************************************************************************************
     * \brief Compute the non-zero basis functions at the given parameter
     * \param[in/out] _basisFcts The non-zero basis functions at the given parameter
     * \param[in] _uPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    virtual void computeLocalBasisFunctions(double*, double, int) = 0;

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter (returns a 2D pointer array)
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    virtual void computeLocalBasisFunctionsAndDerivativesInefficient(double**, int, double,
            int) = 0;

    /***********************************************************************************************
     * \brief Compute the non-zero NURBS basis functions and their derivatives of the  at the given parameter and stores them in a vector of doubles
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _knotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    virtual void computeLocalBasisFunctionsAndDerivatives(double*, int, double, int) = 0;
};

}/* namespace EMPIRE */

#endif /* ABSTRACTBSPLINEBASIS1D_H_ */
