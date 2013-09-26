// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Inclusion of user defined libraries
#include "NurbsBasis1D.h"
#include "IGAMath.h"

using namespace std;

namespace EMPIRE {

NurbsBasis1D::NurbsBasis1D(int _ID = 0, int _pDegree = 0, int _noKnots = 0, double* _KnotVector =
        NULL, int _noControlPoints = 0, double* _controlPointWeights = NULL) :
        BSplineBasis1D(_ID, _pDegree, _noKnots, _KnotVector), NoControlPoints(_noControlPoints) {

    if (_noControlPoints != _noKnots - _pDegree - 1) {
        cout << endl;
        cout << "Error in NurbsBasis1D::NurbsBasis1D" << endl;
        cout
                << "The Number of Control Points, the number of knots and the polynomial degree do not match"
                << endl;
        cout << endl;
        cout << endl;
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

void NurbsBasis1D::computeLocalBasisFunctionsAndDerivativesInefficient(double** _basisFctsAndDerivs,
        int _derivDegree, double _uPrm, int _knotSpanIndex) {

    /* Initialize the output array (This must be done outside of the function call)
     * _localBasisFctsAndDerivs = double_array[Number of derivatives to be computed][Number of non-zero basis functions].
     * The algorithm returns a 2D pointer array which is highly inefficient, use instead the following algorithm:
     * NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives   -----------
     *
     *  							   |R1		  R2		... 	   Rn|
     * 								   |dR1    	  dR2 		... 	  dRn|
     * e.g. _localBasisFctsAndDerivs = |ddR1 	  ddR2      ... 	 ddRn|
     * 								   |...		  ...	   	...    	 ... |
     * 								   |d..k..dR1 d..k..dR2 ... d..k..dRn|
     *
     * Current implementation handles NURBS derivatives up to 2nd order i.e. _derivDegree = 1,2
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);

    if (_derivDegree == 0) {
        cout << endl;
        cout << endl;
        cout << "Error in function NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives" << endl;
        cout << "Asked for NURBS derivatives of order _derivDegree = " << _derivDegree << endl;
        cout << "Alternatively use function NurbsBasis1D::computeLocalBasisFunctions" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    } else if (_derivDegree == 1 || _derivDegree == 2) {

        // Compute the non-zero B-Spline basis functions at the given knot span
        double** BSplineBasisFunctionsAndDerivs = new double*[_derivDegree + 1];
        for (int i = 0; i <= _derivDegree; i++)
            BSplineBasisFunctionsAndDerivs[i] = new double[getPolynomialDegree() + 1];
        BSplineBasis1D::computeLocalBasisFunctionsAndDerivativesInefficient(
                BSplineBasisFunctionsAndDerivs, _derivDegree, _uPrm, _knotSpanIndex);

        // Initialize auxiliary variables
        double sum = 0.0;
        double dsum = 0.0;
        double ddsum = 0.0;

        for (int i = 0; i <= getPolynomialDegree(); i++) {
            // The NURBS basis functions themselves
            _basisFctsAndDerivs[0][i] = BSplineBasisFunctionsAndDerivs[0][i]
                    * ControlPointWeights[_knotSpanIndex - getPolynomialDegree() + i];
            sum += _basisFctsAndDerivs[0][i];

            // The first derivatives of the NURBS basis functions
            _basisFctsAndDerivs[1][i] = BSplineBasisFunctionsAndDerivs[1][i]
                    * ControlPointWeights[_knotSpanIndex - getPolynomialDegree() + i];
            dsum += _basisFctsAndDerivs[1][i];

            if (_derivDegree == 2) {
                // The second derivatives of the NURBS basis functions
                _basisFctsAndDerivs[2][i] = BSplineBasisFunctionsAndDerivs[2][i]
                        * ControlPointWeights[_knotSpanIndex - getPolynomialDegree() + i];
                ddsum += _basisFctsAndDerivs[2][i];
            }
        }

        // Divide through by the sum (order of computation must be reversed due to dependencies)
        for (int i = 0; i <= getPolynomialDegree(); i++) {
            if (_derivDegree == 2) {
                // The second derivatives of the NURBS basis functions
                _basisFctsAndDerivs[2][i] = _basisFctsAndDerivs[2][i] / sum
                        - 2 * _basisFctsAndDerivs[1][i] * dsum / pow(sum, 2)
                        - _basisFctsAndDerivs[0][i] * ddsum / pow(sum, 2)
                        + 2 * _basisFctsAndDerivs[0][i] * pow(dsum, 2) / pow(sum, 3);
            }

            // The first derivatives of the NURBS basis functions
            _basisFctsAndDerivs[1][i] = _basisFctsAndDerivs[1][i] / sum
                    - _basisFctsAndDerivs[0][i] * dsum / pow(sum, 2);

            // The NURBS basis functions themselves
            _basisFctsAndDerivs[0][i] /= sum;
        }

        // Clear the memory on the heap
        for (int i = 0; i <= _derivDegree; i++)
            delete[] BSplineBasisFunctionsAndDerivs[i];
        delete[] BSplineBasisFunctionsAndDerivs;
    } else if (_derivDegree > 2) {
        cout << endl;
        cout << endl;
        cout << "Error in function NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives" << endl;
        cout << "Asked for NURBS derivatives of order _derivDegree = " << _derivDegree
                << "which have not yet implemented" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }
}

void NurbsBasis1D::computeLocalBasisFunctionsAndDerivativesInefficient(double* _basisFctsAndDerivs,
        int _derivDegree, double _uPrm, int _knotSpanIndex) {

    /* Initialize the output array (This must be done outside of the function call)
     * _localBasisFctsAndDerivs = double_array[Number of derivatives to be computed*Number of non-zero basis functions]. The algorithm returns
     * an 1D pointer array which is efficient however is not adapted from the NURBS book hence it handles only derivatives up to 2nd order yet not
     * optimally efficient
     *                              ---------
     *
     * i.e.                          _localBasisFctsAndDerivs =
     *
     *  |R1 R2 ... Rn dR1 dR2 ... dRn ddR1 ddR2 ... ddRn d..k..dR1 d..k..dR2 ... d..k..dRn |
     *
     * Current implementation handles NURBS derivatives up to 2nd order i.e. _derivDegree = 1,2
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);

    if (_derivDegree == 0) {
        cout << endl;
        cout << endl;
        cout << "Error in function NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives" << endl;
        cout << "Asked for NURBS derivatives of order _derivDegree = " << _derivDegree << endl;
        cout << "Alternatively use function NurbsBasis1D::computeLocalBasisFunctions" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    } else if (_derivDegree == 1 || _derivDegree == 2) {

        // The polynomial degree of the basis
        int pDegree = getPolynomialDegree();

        // Compute the non-zero B-Spline basis functions at the given knot span
        double* BSplineBasisFunctionsAndDerivs = new double[(_derivDegree + 1) * (pDegree + 1)];
        BSplineBasis1D::computeLocalBasisFunctionsAndDerivatives(BSplineBasisFunctionsAndDerivs,
                _derivDegree, _uPrm, _knotSpanIndex);

        // Initialize auxiliary variables
        double sum = 0.0;
        double dsum = 0.0;
        double ddsum = 0.0;

        for (int i = 0; i <= this->getPolynomialDegree(); i++) {
            // The NURBS basis functions themselves
            _basisFctsAndDerivs[i] = BSplineBasisFunctionsAndDerivs[i]
                    * ControlPointWeights[_knotSpanIndex - this->getPolynomialDegree() + i];
            sum += _basisFctsAndDerivs[i];

            // The first derivatives of the NURBS basis functions
            _basisFctsAndDerivs[pDegree + 1 + i] = BSplineBasisFunctionsAndDerivs[pDegree + 1 + i]
                    * ControlPointWeights[_knotSpanIndex - pDegree + i];
            dsum += _basisFctsAndDerivs[pDegree + 1 + i];

            if (_derivDegree == 2) {
                // The second derivatives of the NURBS basis functions
                _basisFctsAndDerivs[2 * (pDegree + 1) + i] = BSplineBasisFunctionsAndDerivs[2
                        * (pDegree + 1) + i] * ControlPointWeights[_knotSpanIndex - pDegree + i];
                ddsum += _basisFctsAndDerivs[2 * (pDegree + 1) + i];
            }
        }

        // Divide through by the sum (order of computation must be reversed due to dependencies)
        for (int i = 0; i <= pDegree; i++) {
            if (_derivDegree == 2) {
                // The second derivatives of the NURBS basis functions
                _basisFctsAndDerivs[2 * (pDegree + 1) + i] = _basisFctsAndDerivs[2 * (pDegree + 1)
                        + i] / sum - 2 * _basisFctsAndDerivs[pDegree + 1 + i] * dsum / pow(sum, 2)
                        - _basisFctsAndDerivs[i] * ddsum / pow(sum, 2)
                        + 2 * _basisFctsAndDerivs[i] * pow(dsum, 2) / pow(sum, 3);
            }

            // The first derivatives of the NURBS basis functions
            _basisFctsAndDerivs[pDegree + 1 + i] = _basisFctsAndDerivs[pDegree + 1 + i] / sum
                    - _basisFctsAndDerivs[i] * dsum / pow(sum, 2);

            // The NURBS basis functions themselves
            _basisFctsAndDerivs[i] /= sum;
        }

        // Clear the memory on the heap
        delete[] BSplineBasisFunctionsAndDerivs;

    } else if (_derivDegree > 2) {
        cout << endl;
        cout << endl;
        cout << "Error in function NurbsBasis1D::computeLocalBasisFunctionsAndDerivatives" << endl;
        cout << "Asked for NURBS derivatives of order _derivDegree = " << _derivDegree
                << "which have not yet implemented" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }
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

        cout << endl;
        cout << endl;
        cout << "Error in NurbsBasis1D::setControlPointNet" << endl;
        cout
                << "The assigned number of Control Points does not match with the number of basis functions!"
                << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }

    NoControlPoints = _noControlPoints;
    delete[] ControlPointWeights;
    assert(_controlPointWeights!=NULL);
    ControlPointWeights = _controlPointWeights;
}

void NurbsBasis1D::printControlPointWeights() {

    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class NurbsBasis1D" << endl;
    cout << "NurbsBasis1D::ControlPointWeights" << endl;
    cout << endl;
    for (int i = 0; i < NoControlPoints; i++) {
        cout << "CP " << i << ": " << ControlPointWeights[i] << endl;
    }
    cout << "_____________________________________________" << endl;
    cout << endl;
}

}/* namespace EMPIRE */
