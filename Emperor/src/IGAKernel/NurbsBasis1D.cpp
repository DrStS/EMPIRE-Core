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
// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Inclusion of user defined libraries
#include "NurbsBasis1D.h"
#include "IGAMath.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

NurbsBasis1D::NurbsBasis1D(int _ID = 0, int _pDegree = 0, int _noKnots = 0, double* _KnotVector =
        NULL, int _noControlPoints = 0, double* _controlPointWeights = NULL) :
        BSplineBasis1D(_ID, _pDegree, _noKnots, _KnotVector), NoControlPoints(_noControlPoints) {

    if (_noControlPoints != _noKnots - _pDegree - 1) {
        ERROR_OUT() << " in NurbsBasis1D::NurbsBasis1D" << endl;
        ERROR_OUT() << "The Number of Control Points, the number of knots and the polynomial degree do not match" << endl;
        exit(-1);
    }

    assert(_controlPointWeights!=NULL);
    ControlPointWeights = _controlPointWeights;
}

NurbsBasis1D::~NurbsBasis1D() {
    delete[] ControlPointWeights;
}

NurbsBasis1D::NurbsBasis1D(const NurbsBasis1D& _nurbsBasis1D) :
        BSplineBasis1D(_nurbsBasis1D) {

    NoControlPoints = _nurbsBasis1D.NoControlPoints;
    ControlPointWeights = new double[NoControlPoints];
    for (int i = 0; i < NoControlPoints; i++)
        ControlPointWeights[i] = _nurbsBasis1D.ControlPointWeights[i];
}

void NurbsBasis1D::computeLocalBasisFunctions(double* _localBasisFunctions, double _uPrm,
        int _KnotSpanIndex) {

    /* Initialize the output array (It must be initialized outside the function call)
     * _localBasisFunctions = double_array[Number of non-zero basis functions];
     *
     * e.g.  _localBasisFunctions = |R1 	    R2		 ... 	    Rn|
     */

    // Read input
    assert(_localBasisFunctions!=NULL);

    // Compute the non-zero B-Spline basis functions at the given knot span
    double* BSplineBasisFunctions = new double[getPolynomialDegree() + 1];
    BSplineBasis1D::computeLocalBasisFunctions(BSplineBasisFunctions, _uPrm, _KnotSpanIndex);

    // Initialize auxiliary variable
    double sum = 0;

    // Loop over all non-zero basis functions
    for (int i = 0; i <= getPolynomialDegree(); i++) {
        _localBasisFunctions[i] = BSplineBasisFunctions[i]
                * ControlPointWeights[_KnotSpanIndex - getPolynomialDegree() + i];
        sum += _localBasisFunctions[i];
    }

    // Divide through by the sum
    for (int i = 0; i <= getPolynomialDegree(); i++) {
        _localBasisFunctions[i] /= sum;
    }

    // Clear the memory on the heap
    delete[] BSplineBasisFunctions;
}

void NurbsBasis1D::computeDenominatorFunctionAndDerivatives(double* _denominatorFctAndDerivs,
        double* _bSplineBasisFctsAndDerivs, int _derivDegree, int _knotSpanIndex) {
    /*
     * Initialize the output array (This must be done outside of the function call):
     * _denominatorFctAndDerivs = new double[_derivDegree + 1]
     * The rest of the pointer input must be of the following size:
     * _bSplineBasisFctsAndDerivs = new double[(_derivDegree + 1)*(pDegree + 1)]
     *
     * The algorithm is taken from the following reference;
     * Reference: Piegl, Les and Tiller, Wayne, The NURBS Book. Springer-Verlag: Berlin 1995; p. 127.
     */

    // Read input
    assert(_denominatorFctAndDerivs!=NULL);

    // Initialize the Control Point index
    int CPIndex = 0;

    // The polynomial degree of the basis
    int pDegree = getPolynomialDegree();

    // Initialize the output to zero
    for (int i = 0; i <= _derivDegree; i++)
        _denominatorFctAndDerivs[i] = 0.0;

    // Compute the denominator function

    // Loop over all the derivatives
    for (int i = 0; i <= _derivDegree; i++) {
        // Loop over all the basis functions
        for (int j = 0; j <= pDegree; j++) {

            // Compute the Control Point index
            CPIndex = _knotSpanIndex - pDegree + j;

            // Compute the denominator function iteratively
            _denominatorFctAndDerivs[i] += _bSplineBasisFctsAndDerivs[i * (pDegree + 1) + j]
                    * ControlPointWeights[CPIndex];
        }
    }
}

void NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives(double* _basisFctsAndDerivs,
        int _derivDegree, double _uPrm, int _knotSpanIndex) {
    /*
     * Initialize the output array (This must be done outside of the function call):
     * _localBasisFctsAndDerivs = double_array[(_derivDegree + 1)*(pDegree + 1)]
     *
     * The algorithm returns an 1D pointer array which is efficient taken from the following reference.
     *                                                    ---------
     *
     * Reference: Piegl, Les and Tiller, Wayne, The NURBS Book. Springer-Verlag: Berlin 1995; p. 127.
     *
     * i.e. _localBasisFctsAndDerivs = |R1 R2 ... Rn dR1 dR2 ... dRn ddR1 ddR2 ... ddRn d..k..dR1 d..k..dR2 ... d..k..dRn |
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);

    // Initialize the Control Point index
    int CPIndex = 0;

    // The polynomial degree of the basis
    int pDegree = getPolynomialDegree();

    // Initialize the output array to zero
    for (int i = 0; i < (_derivDegree + 1) * (pDegree + 1); i++)
        _basisFctsAndDerivs[i] = 0.0;

    // Compute the non-zero B-Spline basis functions at the given knot span
    double* BSplineBasisFunctionsAndDerivs = new double[(_derivDegree + 1) * (pDegree + 1)];
    BSplineBasis1D::computeLocalBasisFunctionsAndDerivatives(BSplineBasisFunctionsAndDerivs,
            _derivDegree, _uPrm, _knotSpanIndex);

    // Initialize and compute the denominator function the w(u) = Sum_(i=1)^n N^(i,p)*w_i and its derivatives
    double* denominatorFunction = new double[_derivDegree + 1];
    computeDenominatorFunctionAndDerivatives(denominatorFunction, BSplineBasisFunctionsAndDerivs,
            _derivDegree, _knotSpanIndex);

    // Compute the derivatives of the NURBS basis functions

    // Initialize auxiliary variable
    double v = 0.0;

    // Loop over all the basis functions
    for (int i = 0; i <= pDegree; i++) {
        // Loop over all the derivatives
        for (int j = 0; j <= _derivDegree; j++) {

            // Compute the Control Point index
            CPIndex = _knotSpanIndex - pDegree + i;

            // Compute the product of the derivatives of the basis functions with the Control Point weights
            v = BSplineBasisFunctionsAndDerivs[j * (pDegree + 1) + i]
                    * ControlPointWeights[CPIndex];

            // Loop over all the involved derivatives
            for (int k = 1; k <= j; k++) {
                v -= binomialCoefficients[indexBinomialCoefficients(j, k)] * denominatorFunction[k]
                        * _basisFctsAndDerivs[(j - k) * (pDegree + 1) + i];
            }

            // Divide by the denominator function
            _basisFctsAndDerivs[j * (pDegree + 1) + i] = v / denominatorFunction[0];
        }
    }

    // Free the memory from the heap
    delete[] BSplineBasisFunctionsAndDerivs;
    delete[] denominatorFunction;
}

void NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives(double* _basisFctsAndDerivs,
        double* _bSplineBasisFctsAndDerivs, double* _denominatorFctAndDerivs, int _derivDegree,
        double _uPrm, int _knotSpanIndex) {
    /* Initialize the output array (This must be done outside of the function call). On the input/output parameters:
     *
     * The local NURBS basis functions and their derivatives: _localBasisFctsAndDerivs = double_array[(_derivDegree + 1)*(pDegree + 1)]
     * The local B-Spline basis functions and their derivatives: _bSplineBasisFctsAndDerivs = new double[(_derivDegree + 1) * (pDegree + 1)]
     * The given denominator function and its derivatives: _denominatorFctAndDerivs = new double[_derivDegree + 1]
     *
     * The algorithm returns an 1D pointer array which is efficient taken from the following reference:
     *                                                    ---------
     *
     * Reference: Piegl, Les and Tiller, Wayne, The NURBS Book. Springer-Verlag: Berlin 1995; p. 127.
     *
     * i.e. _localBasisFctsAndDerivs = |R1 R2 ... Rn dR1 dR2 ... dRn ddR1 ddR2 ... ddRn d..k..dR1 d..k..dR2 ... d..k..dRn |
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);
    assert(_bSplineBasisFctsAndDerivs!=NULL);
    assert(_denominatorFctAndDerivs!=NULL);

    // Initialize the Control Point index
    int CPIndex = 0;

    // The polynomial degree of the basis
    int pDegree = getPolynomialDegree();

    // Initialize the output array to zero
    for (int i = 0; i < (_derivDegree + 1) * (pDegree + 1); i++)
        _basisFctsAndDerivs[i] = 0.0;

    // Compute the derivatives of the NURBS basis functions

    // Initialize auxiliary variable
    double v = 0.0;

    // Loop over all the basis functions
    for (int i = 0; i <= pDegree; i++) {
        // Loop over all the derivatives
        for (int j = 0; j <= _derivDegree; j++) {

            // Compute the Control Point index
            CPIndex = _knotSpanIndex - pDegree + i;

            // Compute the product of the derivatives of the basis functions with the Control Point weights
            v = _bSplineBasisFctsAndDerivs[j * (pDegree + 1) + i] * ControlPointWeights[CPIndex];

            // Loop over all the involved derivatives
            for (int k = 1; k <= j; k++) {
                v -= binomialCoefficients[indexBinomialCoefficients(j, k)]
                        * _denominatorFctAndDerivs[k]
                        * _basisFctsAndDerivs[(j - k) * (pDegree + 1) + i];
            }

            // Divide by the denominator function
            _basisFctsAndDerivs[j * (pDegree + 1) + i] = v / _denominatorFctAndDerivs[0];
        }
    }
}

double* NurbsBasis1D::getControlPointWeights() {
    return ControlPointWeights;
}

void NurbsBasis1D::setControlPointWeights(int _noControlPoints, double* _controlPointWeights) {
    // Read input
    assert(_controlPointWeights!=NULL);

    if (_noControlPoints != this->computeNoBasisFunctions()) {

        ERROR_OUT() << "Error in NurbsBasis1D::setControlPointNet" << endl;
        ERROR_OUT() << "The assigned number of Control Points does not match with the number of basis functions!" << endl;
        exit(-1);
    }

    NoControlPoints = _noControlPoints;
    delete[] ControlPointWeights;
    assert(_controlPointWeights!=NULL);
    ControlPointWeights = _controlPointWeights;
}

Message &operator<<(Message &message, NurbsBasis1D &nurbsBasis1D) {

    message << "\t+" << "NurbsBasis1D: " << endl;
    message << "\t\t+" << "PDegree = " << nurbsBasis1D.getPolynomialDegree() << endl;
    message << "\t\t+" << "NoKnots = " << nurbsBasis1D.getNoKnots() << endl;
    message << "\t\t+" << "KnotVector = [\t";
    for (int i = 0; i < nurbsBasis1D.getNoKnots(); i++) {
        message << nurbsBasis1D.getKnotVector()[i] << "\t";
    }
    message << "\t\t+" << "ControlPointWeights = [\t";
    for (int i = 0; i < nurbsBasis1D.getNoControlPoints(); i++) {
        message << "CP " << i << ": " << nurbsBasis1D.getControlPointWeights()[i] << endl;
    }
    message << "]" << endl;
    message() << "\t+" << "---------------------------------" << endl;
    return message;
}

}/* namespace EMPIRE */
