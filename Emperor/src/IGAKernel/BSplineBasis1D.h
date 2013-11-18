/***********************************************************************************************//**
 * \file NurbsBasis1D.h
 * This file holds the class NurbsBasis1D.h
 * \date 11/3/2013
 **************************************************************************************************/

#ifndef BSPLINEBASIS1D_H_
#define BSPLINEBASIS1D_H_

// Inclusion of user defined libraries
#include "AbstractBSplineBasis1D.h"

namespace EMPIRE {

/********//**
 * \brief class BSplineBasis1D is related to the BSpline basis functions in 1D case
 ***********/

class BSplineBasis1D: public AbstractBSplineBasis1D {

protected:
    /// The polynomial degree of the basis
    int PDegree;

    /// The number of knots
    int NoKnots;

    /// The knot vector
    double* KnotVector;

    /// The constructor, the destructor, the copy constructor and the copy assignment
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _ID The id of the basis
     * \param[in] _pDegree The polynomial degree of the basis
     * \param[in] _noKnots The number of knots
     * \param[in] _KnotVector The knot vector
     * \author Andreas Apostolatos
     ***********/
    BSplineBasis1D(int, int, int, double*);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~BSplineBasis1D() {

    }

    /***********************************************************************************************
     * \brief Copy constructor
     * \author Andreas Apostolatos
     ***********/
    BSplineBasis1D(const BSplineBasis1D&);

    /***********************************************************************************************
     * \brief Copy assignment
     * \author Andreas Apostolatos
     ***********/
    BSplineBasis1D& operator=(const BSplineBasis1D&);

    /// On the 1D B-Spline basis functions
public:
    /***********************************************************************************************
     * \brief Returns the number of basis functions
     * \author Andreas Apostolatos
     ***********/
    inline int computeNoBasisFunctions() {
        return NoKnots - PDegree - 1;
    }

    /***********************************************************************************************
     * \brief Returns the polynomial degree of the B-Spline 1D basis
     * \param[in] _uPrm The parameter on which the knot span is searced
     * \author Andreas Apostolatos
     ***********/
    int findKnotSpan(double);

    /***********************************************************************************************
     * \brief Compute the non-zero B-Spline basis functions at the given parameter
     * \param[in/out] _basisFcts The non-zero basis functions at the given parameter
     * \param[in] _uPrm The parameter where the basis functions are evaluated
     * \param[in] _KnotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctions(double*, double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero basis functions and their derivatives of the  at the given parameter (Inefficient algorithm)
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _KnotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivativesInefficient(double**, int, double, int);

    /***********************************************************************************************
     * \brief Compute the non-zero basis functions and their derivatives of the  at the given parameter and stores them in a double array
     * \param[in/out] _basisFctsAndDerivs The non-zero basis functions at the given knot span up to their _derivDegree-th derivatives
     * \param[in] _derivDegree For Which degree the B-Spline derivative to be computed
     * \param[in] _uPrm The parameter where the basis functions and their derivatives are evaluated
     * \param[in] _KnotSpanIndex The index of the knot span where _uPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeLocalBasisFunctionsAndDerivatives(double*, int, double, int);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Returns the polynomial degree of the B-Spline 1D basis
     * \author Andreas Apostolatos
     ***********/
    inline int getPolynomialDegree() {
        return PDegree;
    }

    /***********************************************************************************************
     * \brief Returns the number of knots for the 1D B-Spline basis
     * \author Andreas Apostolatos
     ***********/
    inline int getNoKnots() {
        return NoKnots;
    }

    /***********************************************************************************************
     * \brief Returns the knot vector of the B-Spline 1D basis
     * \author Andreas Apostolatos
     ***********/
    double* getKnotVector() {
        return KnotVector;
    }

    /***********************************************************************************************
     * \brief Sets the polynomial degree of the B-Spline 1D basis
     * \author Andreas Apostolatos
     ***********/
    void setPolynomialDegree(int _pDegree) {
        PDegree = _pDegree;
    }

    /***********************************************************************************************
     * \brief Sets the number of knots for the 1D B-Spline basis
     * \param[in] _noKnots The number of knots for the new knot vector
     * \param[in] _knotVector The new knot vector to be assigned
     * \author Andreas Apostolatos
     ***********/
    void setKnotVector(int, double*);

    /// DEBUGGING functions
public:
    /***********************************************************************************************
     * \brief Outputs the polynomial degree on the terminal
     * \author Andreas Apostolatos
     ***********/
    void printPolynomialDegree();

    /***********************************************************************************************
     * \brief Outputs the knot vector in the terminal
     * \author Andreas Apostolatos
     ***********/
    void printNoKnots();

    /***********************************************************************************************
     * \brief Outputs the knot vector in the terminal
     * \author Andreas Apostolatos
     ***********/
    void printKnotVector();

    /***********************************************************************************************
     * \brief Outputs the number of basis functions
     * \author Andreas Apostolatos
     ***********/
    void printNoBasisFunctions();

    /// The tolerance for accepting a very small number as interior to the knot span
    static const double EPS_ACCPETEDINTOKNOTSPAN;
};

}/* namespace EMPIRE */

#endif /* BSPLINEBASIS1D_H_ */
