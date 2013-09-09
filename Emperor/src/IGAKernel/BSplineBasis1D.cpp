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
// Inclusion of standard libraries
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include<assert.h>

// Inclusion of user defined libraries
#include "BSplineBasis1D.h"

using namespace std;

namespace EMPIRE {

BSplineBasis1D::BSplineBasis1D(int _ID = 0, int _pDegree = 0, int _noKnots = 0,
        double* _KnotVector = 0) :
        AbstractBSplineBasis1D(_ID), PDegree(_pDegree), NoKnots(_noKnots) {

    // Assign the given pointer to the knot vector
    assert(_KnotVector!=NULL);
    KnotVector = _KnotVector;

    for (int i = 1; i < PDegree + 1; i++) {
        if (KnotVector[i] != KnotVector[0]) {
            cout << endl;
            cout << "Error in BSplineBasis1D::BSplineBasis1D" << endl;
            cout << "p+1 first knots are not the same but knot vector is expected to be open"
                    << endl;
            cout << endl;
            cout << endl;
            exit(-1);
        }
    }
    for (int i = NoKnots - 2; i >= NoKnots - PDegree; i--) {
        if (KnotVector[i] != KnotVector[NoKnots - 1]) {
            cout << endl;
            cout << "Error in BSplineBasis1D::BSplineBasis1D" << endl;
            cout << "p+1 last knots are not the same but knot vector is expected to be open"
                    << endl;
            cout << endl;
            cout << endl;
            exit(-1);
        }
    }
}

BSplineBasis1D::BSplineBasis1D(const BSplineBasis1D& _bsplineBasis1D) :
        AbstractBSplineBasis1D(_bsplineBasis1D) {

    PDegree = _bsplineBasis1D.PDegree;
    NoKnots = _bsplineBasis1D.NoKnots;
    KnotVector = new double[NoKnots];
    for (int i = 0; i < NoKnots; i++) {
        KnotVector[i] = _bsplineBasis1D.KnotVector[i];
    }
}

BSplineBasis1D& BSplineBasis1D::operator=(const BSplineBasis1D& _bsplineBasis1D) {
    PDegree = _bsplineBasis1D.PDegree;
    if (this != &_bsplineBasis1D) {
        if (NoKnots != _bsplineBasis1D.NoKnots) {
            delete[] KnotVector;
            KnotVector = new double[NoKnots = _bsplineBasis1D.NoKnots];
        }
        for (int i = 0; i < NoKnots; i++) {
            KnotVector[i] = _bsplineBasis1D.KnotVector[i];
        }
    }
    return *this;
}

int BSplineBasis1D::findKnotSpan(double _uPrm) {
    // Check input
    if (_uPrm < KnotVector[0] || _uPrm > KnotVector[NoKnots - 1]) {
        cout << endl;
        cout << "Error in BSplineBasis1D::find_knot_span" << endl;
        cout << "Given parameter is outside of the knot span" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }

    // Compute the number of basis functions
    int n = computeNoBasisFunctions();

    // Special case, for the last knot
    if (_uPrm == KnotVector[n + 1])
        return n - 1;

    // Do binary search
    int low = PDegree;
    int high = n + 1;
    int mid = (low + high) / 2;
    while (_uPrm < KnotVector[mid] || _uPrm >= KnotVector[mid + 1]) {
        if (_uPrm < KnotVector[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return mid;
}

void BSplineBasis1D::computeLocalBasisFunctions(double* _localBasisFunctions, double _uPrm,
        int _KnotSpanIndex) {

    /* Initialize the output array (It must be initialized outside the function call)
     * _localBasisFunctions = double_array[Number of non-zero basis functions];
     *
     * e.g.  _localBasisFunctions = |N1 	    N2		 ... 	    Nn|
     */

    // Initialize the first element of the output array
    _localBasisFunctions[0] = 1.0;

    // Initialize auxiliary variables
    double saved = 0.0;
    double temp = 0.0;

    // Initialize auxiliary arrays
    double* left = new double[PDegree + 1];
    double* right = new double[PDegree + 1];

    for (int j = 1; j <= PDegree; j++) {
        left[j - 1] = _uPrm - KnotVector[_KnotSpanIndex + 1 - j];
        right[j - 1] = KnotVector[_KnotSpanIndex + j] - _uPrm;
        saved = 0.0;
        for (int r = 0; r < j; r++) {
            temp = _localBasisFunctions[r] / (right[r] + left[j - r - 1]);
            _localBasisFunctions[r] = saved + right[r] * temp;
            saved = left[j - r - 1] * temp;
        }
        _localBasisFunctions[j] = saved;
    }

    // Clean up the heap from the auxiliary pointers
    delete[] left;
    delete[] right;
}

void BSplineBasis1D::computeLocalBasisFunctionsAndDerivativesInefficient(
        double** _localBasisFctsAndDerivs, int _derivDegree, double _uPrm, int _KnotSpanIndex) {

    /* Initialize the output array (This must be done outside of the function call)
     * _localBasisFctsAndDerivs = double_array[Number of derivatives to be computed][Number of non-zero basis functions]
     * The algorithm is inefficient since an array of pointers is used for the construction of the matrix
     *                  -----------
     *
     *  							   |N1		  N2		... 	   Nn|
     * 								   |dN1    	  dN2 		... 	  dNn|
     * e.g. _localBasisFctsAndDerivs = |ddN1 	  ddN2      ... 	 ddNn|
     * 								   |...		  ...	   	...    	 ... |
     * 								   |d..k..dN1 d..k..dN2 ... d..k..dNn|
     */

    // Initialize auxiliary variables
    double saved = 0.0;
    double temp = 0.0;
    int s1 = 0;
    int s2 = 0;
    double d = 0.0;
    int rk = 0;
    int pk = 0;
    int j1 = 0;
    int j2 = 0;
    int jTemp = 0;
    int rTemp = 0;

    // Initialize auxiliary arrays
    double** ndu = new double*[PDegree + 1];
    for (int i = 0; i <= PDegree; i++)
        ndu[i] = new double[PDegree + 1];

    double* left = new double[PDegree + 1];
    double* right = new double[PDegree + 1];
    double** a = new double*[2];
    for (int i = 0; i < 2; i++)
        a[i] = new double[PDegree + 1];

    // Initialize the first element of the output array
    ndu[0][0] = 1.0;

    for (int j = 1; j <= PDegree; j++) {
        left[j - 1] = _uPrm - KnotVector[_KnotSpanIndex + 1 - j];
        right[j - 1] = KnotVector[_KnotSpanIndex + j] - _uPrm;
        saved = 0.0;
        for (int r = 0; r < j; r++) {
            // Lower triangle
            ndu[j][r] = right[r] + left[j - r - 1];
            temp = ndu[r][j - 1] / ndu[j][r];

            // Upper triangle
            ndu[r][j] = saved + right[r] * temp;
            saved = left[j - r - 1] * temp;
        }
        ndu[j][j] = saved;
    }

    // Load the basis functions
    for (int j = 0; j <= PDegree; j++) {
        _localBasisFctsAndDerivs[0][j] = ndu[j][PDegree];
    }

    // This section computes the derivatives
    for (int r = 0; r <= PDegree; r++) { // Loop over function index
        s1 = 0;
        s2 = 1;
        a[0][0] = 1.0;

        // Loop to compute the _derivDegree-th derivative
        for (int k = 1; k <= _derivDegree; k++) {
            d = 0.0;
            rk = r - k;
            pk = PDegree - k;

            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;
            if (r - 1 <= pk)
                j2 = k - 1;
            else
                j2 = PDegree - r;

            for (int j = j1; j <= j2; j++) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }

            _localBasisFctsAndDerivs[k][r] = d;

            // Switch rows
            jTemp = s1;
            s1 = s2;
            s2 = jTemp;
        }
    }

    // Multiply through by the correct factors
    rTemp = PDegree;
    for (int k = 1; k <= _derivDegree; k++) {
        for (int j = 0; j <= PDegree; j++)
            _localBasisFctsAndDerivs[k][j] *= rTemp;

        rTemp *= (PDegree - k);
    }

    // Delete auxiliary arrays
    for (int i = 0; i <= PDegree; i++)
        delete[] ndu[i];
    delete[] ndu;
    delete[] left;
    delete[] right;
    for (int i = 0; i < 2; i++)
        delete[] a[i];
    delete[] a;
}

void BSplineBasis1D::computeLocalBasisFunctionsAndDerivatives(double* _localBasisFctsAndDerivs,
        int _derivDegree, double _uPrm, int _KnotSpanIndex) {

    /* Initialize the output array (This must be done outside of the function call)
     * _localBasisFctsAndDerivs = double_array[Number of derivatives to be computed*Number of non-zero basis functions]
     * The algorithm is efficient since only a double array is used to store the complete information
     *                  ---------
     *
     * i.e.                          _localBasisFctsAndDerivs =
     *
     *  |N1 N2 ... Nn dN1 dN2 ... dNn ddN1 ddN2 ... ddNn d..k..dN1 d..k..dN2 ... d..k..dNn |
     *
     */

    // Initialize auxiliary variables
    double saved = 0.0;
    double temp = 0.0;
    int s1 = 0;
    int s2 = 0;
    double d = 0.0;
    int rk = 0;
    int pk = 0;
    int j1 = 0;
    int j2 = 0;
    int jTemp = 0;
    int rTemp = 0;

    // Initialize auxiliary arrays
    double* ndu = new double[(PDegree + 1) * (PDegree + 1)];
    double* left = new double[PDegree + 1];
    double* right = new double[PDegree + 1];
    double* a = new double[2 * (PDegree + 1)];

    // Initialize the first element of the output array
    ndu[0] = 1.0;

    for (int j = 1; j <= PDegree; j++) {
        left[j - 1] = _uPrm - KnotVector[_KnotSpanIndex + 1 - j];
        right[j - 1] = KnotVector[_KnotSpanIndex + j] - _uPrm;
        saved = 0.0;
        for (int r = 0; r < j; r++) {
            // Lower triangle
            ndu[j * (PDegree + 1) + r] = right[r] + left[j - r - 1];
            temp = ndu[r * (PDegree + 1) + j - 1] / ndu[j * (PDegree + 1) + r];

            // Upper triangle
            ndu[r * (PDegree + 1) + j] = saved + right[r] * temp;
            saved = left[j - r - 1] * temp;
        }
        ndu[j * (PDegree + 1) + j] = saved;
    }

    // Load the basis functions
    for (int j = 0; j <= PDegree; j++) {
        _localBasisFctsAndDerivs[j] = ndu[j * (PDegree + 1) + PDegree];
    }

    // This section computes the derivatives
    for (int r = 0; r <= PDegree; r++) { // Loop over function index
        s1 = 0;
        s2 = 1;
        a[0] = 1.0;

        // Loop to compute the _derivDegree-th derivative
        for (int k = 1; k <= _derivDegree; k++) {
            d = 0.0;
            rk = r - k;
            pk = PDegree - k;

            if (r >= k) {
                a[s2 * (PDegree + 1)] = a[s1 * (PDegree + 1)] / ndu[(pk + 1) * (PDegree + 1) + rk];
                d = a[s2 * (PDegree + 1)] * ndu[rk * (PDegree + 1) + pk];
            }

            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;
            if (r - 1 <= pk)
                j2 = k - 1;
            else
                j2 = PDegree - r;

            for (int j = j1; j <= j2; j++) {
                a[s2 * (PDegree + 1) + j] = (a[s1 * (PDegree + 1) + j]
                        - a[s1 * (PDegree + 1) + j - 1]) / ndu[(pk + 1) * (PDegree + 1) + rk + j];
                d += a[s2 * (PDegree + 1) + j] * ndu[(rk + j) * (PDegree + 1) + pk];
            }

            if (r <= pk) {
                a[s2 * (PDegree + 1) + k] = -a[s1 * (PDegree + 1) + k - 1]
                        / ndu[(pk + 1) * (PDegree + 1) + r];
                d += a[s2 * (PDegree + 1) + k] * ndu[r * (PDegree + 1) + pk];
            }

            _localBasisFctsAndDerivs[k * (PDegree + 1) + r] = d;

            // Switch rows
            jTemp = s1;
            s1 = s2;
            s2 = jTemp;
        }
    }

    // Multiply through by the correct factors
    rTemp = PDegree;
    for (int k = 1; k <= _derivDegree; k++) {
        for (int j = 0; j <= PDegree; j++)
            _localBasisFctsAndDerivs[k * (PDegree + 1) + j] *= rTemp;

        rTemp *= (PDegree - k);
    }

    // Delete auxiliary arrays
    delete[] ndu;
    delete[] left;
    delete[] right;
    delete[] a;
}

void BSplineBasis1D::setKnotVector(int _noKnots, double* _knotVector) {
    assert(_knotVector!=NULL);
    delete[] KnotVector;
    NoKnots = _noKnots;
    KnotVector = _knotVector;
}

void BSplineBasis1D::printPolynomialDegree() {
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis1D" << endl;
    cout << "BSplineBasis1D::PDegree = " << PDegree << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

void BSplineBasis1D::printNoKnots() {
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis1D" << endl;
    cout << "BSplineBasis1D::NoKnots = " << NoKnots << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

void BSplineBasis1D::printKnotVector() {
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis1D" << endl;
    cout << "BSplineBasis1D::KnotVector = [\t";
    for (int i = 0; i < NoKnots; i++) {
        cout << KnotVector[i] << "\t";
    }
    cout << "]" << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

void BSplineBasis1D::printNoBasisFunctions() {
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis1D" << endl;
    cout << "BSplineBasis1D::NoBasisFunctions = " << computeNoBasisFunctions() << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

}/* namespace EMPIRE */

