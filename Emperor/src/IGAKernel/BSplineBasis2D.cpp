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
#include "BSplineBasis2D.h"

using namespace std;

namespace EMPIRE {

BSplineBasis2D::BSplineBasis2D(int _ID = 0, int _pDegree = 0, int _noKnotsU = 0,
        double* _KnotVectorU = NULL, int _qDegree = 0, int _noKnotsV = 0, double* _KnotVectorV =
                NULL) :
        AbstractBSplineBasis2D(_ID) {

    // The NURBS basis functions in u-direction
    uBSplineBasis1D = new BSplineBasis1D(_ID, _pDegree, _noKnotsU, _KnotVectorU);

    // The NURBS basis functions in v-direction
    vBSplineBasis1D = new BSplineBasis1D(_ID, _qDegree, _noKnotsV, _KnotVectorV);
}

BSplineBasis2D::~BSplineBasis2D() {
    delete uBSplineBasis1D;
    delete vBSplineBasis1D;
}

BSplineBasis2D::BSplineBasis2D(const BSplineBasis2D& _bSplineBasis2D) :
        AbstractBSplineBasis2D(_bSplineBasis2D) {

    // Copy the NURBS basis functions in u-direction
    uBSplineBasis1D = new BSplineBasis1D(*(_bSplineBasis2D.uBSplineBasis1D));

    // Copy the NURBS basis functions in v-direction
    vBSplineBasis1D = new BSplineBasis1D(*(_bSplineBasis2D.vBSplineBasis1D));
}

BSplineBasis2D& BSplineBasis2D::operator=(const BSplineBasis2D& _bSplineBasis2D) {
    if (this != &_bSplineBasis2D) {
        // Copy the NURBS basis functions in u-direction
        uBSplineBasis1D = new BSplineBasis1D(*(_bSplineBasis2D.uBSplineBasis1D));

        // Copy the NURBS basis functions in v-direction
        vBSplineBasis1D = new BSplineBasis1D(*(_bSplineBasis2D.vBSplineBasis1D));
    }
    return *this;
}

int BSplineBasis2D::indexDerivativeBasisFunction(int _derivDegree, int _uDerivIndex,
        int _vDerivIndex, int _basisIndex) {
    /*
     * Returns the correct index when sorting the basis functions and the their derivatives in an 1D pointer array with the rule:
     * _basisFctsAndDerivs = new double[(_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex+1) * noBasisFcts/2 + _uDerivIndex * noBasisFcts + _basisIndex]
     */

    // Read input
    if (_uDerivIndex + _vDerivIndex > _derivDegree) {
        cout << endl;
        cout << endl;
        cout << "Error in BSplineBasis2D::indexDerivativeBasisFunction" << endl;
        cout << "It has been requested the " << _uDerivIndex
                << "-th partial derivative w.r.t. u and" << endl;
        cout << "the " << _vDerivIndex
                << "-th partial derivative w.r.t. v of the basis functions but " << endl;
        cout << "the maximum absolute derivative selected is of " << _derivDegree << "-th order"
                << endl;
        cout << endl;
        cout << endl;
    }

    // The polynomial degrees of the NURBS basis in both directions
    int pDegree = uBSplineBasis1D->getPolynomialDegree();
    int qDegree = vBSplineBasis1D->getPolynomialDegree();

    // The number of basis functions
    int noBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Compute the index of the functional value
    return (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1) * noBasisFcts / 2
            + _uDerivIndex * noBasisFcts + _basisIndex;
}

void BSplineBasis2D::computeLocalBasisFunctions(double* _basisFcts, double _uPrm,
        int _KnotSpanIndexU, double _vPrm, int _KnotSpanIndexV) {
    /*
     *  Computes the B-Spline basis functions in 2D
     *
     *  eta
     *   |
     *   |	CP(1,m) --> R((m-1)*n+1)	CP(2,m) --> R((m-1)*n+2)	...		CP(n,m) --> R(n*m)
     *   |
     *   |	...							...							...	 	...
     *   |
     *   |	CP(1,2) --> R(n+1)			CP(2,2) --> R(n+2)			...		CP(n,2) --> R(2*n)
     *   |
     *   |	CP(1,1) --> R(1)			CP(2,1) --> R(2)			...		CP(n,1) --> R(n)
     *   |______________________________________________________________________________________xi
     *
     *   The functional values are sorted in a vector as follows:
     *
     *  _basisFcts = [N(1)	N(2)	...		N(n)	N(n+1)	N(n+2)	...		N(2*n)	...		N((m-1)*n+1)	N((m-1)*n+2)	...		N(n*m)]
     *
     *  Assuming that there are n basis functions in u-direction and m basis functions in v-direction
     *
     */

    // Read input
    assert(_basisFcts!=NULL);

    // Compute the local B-Spline basis functions in u-direction
    double* uBSplinebasis1DFcts = new double[uBSplineBasis1D->getPolynomialDegree() + 1];
    uBSplineBasis1D->computeLocalBasisFunctions(uBSplinebasis1DFcts, _uPrm, _KnotSpanIndexU);

    // Compute the local B-Spline basis functions in v-direction
    double* vBSplinebasis1DFcts = new double[vBSplineBasis1D->getPolynomialDegree() + 1];
    vBSplineBasis1D->computeLocalBasisFunctions(vBSplinebasis1DFcts, _vPrm, _KnotSpanIndexV);

    // Initialize counter
    int counter = 0;

    // Loop over all the non-zero contributions of the tensor product functions
    for (int j = 0; j <= vBSplineBasis1D->getPolynomialDegree(); j++) {
        for (int i = 0; i <= uBSplineBasis1D->getPolynomialDegree(); i++) {
            // Assign the basis function value
            _basisFcts[counter] = uBSplinebasis1DFcts[i] * vBSplinebasis1DFcts[j];

            // Update the counter
            counter += 1;
        }
    }

    // Free the memory from the heap
    delete[] uBSplinebasis1DFcts;
    delete[] vBSplinebasis1DFcts;
}

void BSplineBasis2D::getBasisFunctionsIndex(int _KnotSpanIndexU, int _KnotSpanIndexV,
        int* _funcsIndex) {
    /*
     *  Returns the index of the basis function given the knot span indices.
     *
     *  Note that:
     *  The numbering scheme for basis functions shows in BSplineBasis2D::computeLocalBasisFunctions
     *  The numbering scheme for control points shows in IGAPatchSurface::computeCartesianCoordinates and
     *  are different
     */

    // Check input
    assert(_funcsIndex!=NULL);

    int counter = 0;
    int numBasisFuncsV = vBSplineBasis1D->computeNoBasisFunctions();

    for (int i = _KnotSpanIndexV - vBSplineBasis1D->getPolynomialDegree(); i <= _KnotSpanIndexV;
            i++)
        for (int j = _KnotSpanIndexU - uBSplineBasis1D->getPolynomialDegree(); j <= _KnotSpanIndexU;
                j++)
            _funcsIndex[counter++] = j * numBasisFuncsV + i;

}

void BSplineBasis2D::computeLocalBasisFunctionsAndDerivativesInefficient(
        double** _basisFctsAndDerivs, int _maxMixDerivOrd, int _derivDegreeU, double _uPrm,
        int _KnotSpanIndexU, int _derivDegreeV, double _vPrm, int _KnotSpanIndexV) {
    /*
     *  Computes the B-Spline basis functions and their derivatives of the  in 2D
     *
     *  eta
     *   |
     *   |	CP(1,m) --> (m-1)*n+1	CP(2,m) --> (m-1)*n+2	...		CP(n,m) --> n*m
     *   |
     *   |	...							...					...	 	...
     *   |
     *   |	CP(1,2) --> n+1			CP(2,2) --> n+2			...		CP(n,2) --> 2*n
     *   |
     *   |	CP(1,1) --> 1			CP(2,1) --> 2			...		CP(n,1) --> n
     *   |_______________________________________________________________________________xi
     *
     * The functional values are sorted in a matrix as follows:
     *
     *                       |N(1)            ...   N(n)            ...     N(n+1)              ...     N(2*n)              ...     N((m-1)*n+1)            ...     N(n*m)         |
     *                       |dNdu(1)         ...   dNdu(n)         ...     dNdu(n+1)           ...     dNdu(2*n)           ...     dNdu((m-1)*n+1)         ...     dNdu(n*m)      |
     *                       |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdu^k(1)       ...   dNdu^k(n)       ...     dNdu^k(n+1)         ...     dNdu^k(2*n)         ...     dNdu^k((m-1)*n+1)       ...     dNdu^k(n*m)    |
     *                       |dNdv(1)         ...   dNdv(n)         ...     dNdv(n+1)           ...     dNdv(2*n)           ...     dNdv((m-1)*n+1)         ...     dNdv(n*m)      |
     *                       |dNdudv(1)       ...   dNdudv(n)       ...     dNdudv(n+1)         ...     dNdudv(2*n)         ...     dNdudv((m-1)*n+1)       ...     dNdudv(n*m)    |
     * _basisFctsAndDerivs = |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdu^kdv(1)     ...   dNdudv(n)       ...     dNdudv(n+1)         ...     dNdudv(2*n)         ...     dNdudv((m-1)*n+1)       ...     dNdudv(n*m)    |
     *                       |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdv^l(1)       ...   dNdv^l(n)       ...     dNdv^l(n+1)         ...    dNdv^l(2*n)          ...     dNdv^l((m-1)*n+1)       ...     dNdv^l(n*m)    |
     *                       |dNdudv^l(1)     ...   dNdudv^l(n)     ...     dNdudv^l(n+1)       ...    dNdudv^l(2*n)        ...     dNdudv^l((m-1)*n+1)     ...     dNdudv^l(n*m)  |
     *                       |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdu^kdv^l(1)   ...   dNdu^kdv^l(n)   ...     dNdu^kdv^l(n+1)     ...     dNdu^kdv^l(2*n)     ...     dNdu^kdv^l((m-1)*n+1)   ...     dNdu^kdv^l(n*m)|
     *
     *  Assuming that there are n basis functions in u-direction and m basis functions in v-direction and additionally that the highest derivative in
     *  u-direction is k-th order and in v-direction l-th order.
     *
     *  The number of lines to be returned is controlled by the defined maximal order of the mixed derivative degree as well as by the indivdual partial derivatives degree.
     *
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);

    if (_maxMixDerivOrd > _derivDegreeU + _derivDegreeV) {
        cout << endl;
        cout << endl;
        cout << "Error in BSplineBasis2D::computeLocalBasisFunctionsAndDerivatives" << endl;
        cout << "Requested order of partial derivatives du: " << _derivDegreeU << " and dv: "
                << _derivDegreeV << endl;
        cout << "but requested maximal order of the mixed derivatives is dudv: " << _maxMixDerivOrd
                << endl;
        cout << "which is not possible" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }

    // Compute the derivatives of the local B-Spline basis functions in u-direction
    double** uBSplinebasis1DFctsDerivs = new double*[_derivDegreeU + 1];
    for (int i = 0; i <= _derivDegreeU; i++)
        uBSplinebasis1DFctsDerivs[i] = new double[uBSplineBasis1D->getPolynomialDegree() + 1];
    uBSplineBasis1D->computeLocalBasisFunctionsAndDerivativesInefficient(uBSplinebasis1DFctsDerivs,
            _derivDegreeU, _uPrm, _KnotSpanIndexU);

    // Compute the derivatives of the local B-Spline basis functions in v-direction
    double** vBSplinebasis1DFctsDerivs = new double*[_derivDegreeV + 1];
    for (int i = 0; i <= _derivDegreeV; i++)
        vBSplinebasis1DFctsDerivs[i] = new double[vBSplineBasis1D->getPolynomialDegree() + 1];
    vBSplineBasis1D->computeLocalBasisFunctionsAndDerivativesInefficient(vBSplinebasis1DFctsDerivs,
            _derivDegreeV, _vPrm, _KnotSpanIndexV);

    // Initialize counters
    int counterBasis = 0;
    int counterDegree = 0;

    // Loop over all the non-zero contributions of the tensor product functions
    for (int j = 0; j <= vBSplineBasis1D->getPolynomialDegree(); j++) {
        for (int i = 0; i <= uBSplineBasis1D->getPolynomialDegree(); i++) {
            for (int k = 0; k <= _derivDegreeV; k++) {
                for (int l = 0; l <= _derivDegreeU; l++) {
                    if (k + l <= _maxMixDerivOrd) {
                        // Compute the functional value as a tensor product
                        _basisFctsAndDerivs[counterDegree][counterBasis] =
                                uBSplinebasis1DFctsDerivs[l][i] * vBSplinebasis1DFctsDerivs[k][j];

                        // Update the counter related to the computed derivative
                        counterDegree += 1;
                    }
                }
            }
            // Update the counter related to the individual basis function
            counterBasis += 1;

            // Reset the counter related to the computed derivative
            counterDegree = 0;
        }
    }

    // Free the memory from the pointers on the heap
    for (int i = 0; i <= _derivDegreeU; i++)
        delete[] uBSplinebasis1DFctsDerivs[i];
    delete[] uBSplinebasis1DFctsDerivs;
    for (int i = 0; i <= _derivDegreeV; i++)
        delete[] vBSplinebasis1DFctsDerivs[i];
    delete[] vBSplinebasis1DFctsDerivs;
}

void BSplineBasis2D::computeLocalBasisFunctionsAndDerivativesInefficient(
        double* _basisFctsAndDerivs, int _maxMixDerivOrd, int _derivDegreeU, double _uPrm,
        int _KnotSpanIndexU, int _derivDegreeV, double _vPrm, int _KnotSpanIndexV) {
    /*
     *  Computes the B-Spline basis functions and their derivatives of the  in 2D and stores them in the in/out array
     *  _basisFctsAndDerivs = new double [noDerivs*noBasisFunctions]
     *
     *  eta
     *   |
     *   |  CP(1,m) --> (m-1)*n+1   CP(2,m) --> (m-1)*n+2   ...     CP(n,m) --> n*m
     *   |
     *   |  ...                         ...                 ...     ...
     *   |
     *   |  CP(1,2) --> n+1         CP(2,2) --> n+2         ...     CP(n,2) --> 2*n
     *   |
     *   |  CP(1,1) --> 1           CP(2,1) --> 2           ...     CP(n,1) --> n
     *   |_______________________________________________________________________________xi
     *
     * The functional values are sorted in a matrix as computed:
     *
     *                       |N(1)            ...   N(n)            ...     N(n+1)              ...     N(2*n)              ...     N((m-1)*n+1)            ...     N(n*m)         |
     *                       |dNdu(1)         ...   dNdu(n)         ...     dNdu(n+1)           ...     dNdu(2*n)           ...     dNdu((m-1)*n+1)         ...     dNdu(n*m)      |
     *                       |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdu^k(1)       ...   dNdu^k(n)       ...     dNdu^k(n+1)         ...     dNdu^k(2*n)         ...     dNdu^k((m-1)*n+1)       ...     dNdu^k(n*m)    |
     *                       |dNdv(1)         ...   dNdv(n)         ...     dNdv(n+1)           ...     dNdv(2*n)           ...     dNdv((m-1)*n+1)         ...     dNdv(n*m)      |
     *                       |dNdudv(1)       ...   dNdudv(n)       ...     dNdudv(n+1)         ...     dNdudv(2*n)         ...     dNdudv((m-1)*n+1)       ...     dNdudv(n*m)    |
     * _basisFctsAndDerivs = |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdu^kdv(1)     ...   dNdudv(n)       ...     dNdudv(n+1)         ...     dNdudv(2*n)         ...     dNdudv((m-1)*n+1)       ...     dNdudv(n*m)    |
     *                       |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdv^l(1)       ...   dNdv^l(n)       ...     dNdv^l(n+1)         ...    dNdv^l(2*n)          ...     dNdv^l((m-1)*n+1)       ...     dNdv^l(n*m)    |
     *                       |dNdudv^l(1)     ...   dNdudv^l(n)     ...     dNdudv^l(n+1)       ...    dNdudv^l(2*n)        ...     dNdudv^l((m-1)*n+1)     ...     dNdudv^l(n*m)  |
     *                       |...             ...   ...             ...     ...                 ...     ...                 ...     ...                     ...     ...            |
     *                       |dNdu^kdv^l(1)   ...   dNdu^kdv^l(n)   ...     dNdu^kdv^l(n+1)     ...     dNdu^kdv^l(2*n)     ...     dNdu^kdv^l((m-1)*n+1)   ...     dNdu^kdv^l(n*m)|
     *
     *  but they are stored in a vector as
     *
     *  _basisFctsAndDerivs = |_basisFctsAndDerivs[1,:] _basisFctsAndDerivs[2,:] ... _basisFctsAndDerivs[(_derivDegreeU)*(_derivDegreeV+1),:]|
     *
     *  It is also assumed that there are n basis functions in u-direction and m basis functions in v-direction and additionally that the highest derivative in
     *  u-direction is k-th order and in v-direction l-th order.
     *
     *  The number of lines to be returned is controlled by the defined maximal order of the mixed derivative degree as well as by the indivdual partial derivatives degree.
     *
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);

    if (_maxMixDerivOrd > _derivDegreeU + _derivDegreeV) {
        cout << endl;
        cout << endl;
        cout << "Error in BSplineBasis2D::computeLocalBasisFunctionsAndDerivatives" << endl;
        cout << "Requested order of partial derivatives du: " << _derivDegreeU << " and dv: "
                << _derivDegreeV << endl;
        cout << "but requested maximal order of the mixed derivatives is dudv: " << _maxMixDerivOrd
                << endl;
        cout << "which is not possible" << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }

    // Compute the number of basis functions
    int noBasisFcts = (uBSplineBasis1D->getPolynomialDegree() + 1)
            * (vBSplineBasis1D->getPolynomialDegree() + 1);

    // Compute the derivatives of the local B-Spline basis functions in u-direction
    double* uBSplinebasis1DFctsDerivs = new double[(_derivDegreeU + 1)
            * (uBSplineBasis1D->getPolynomialDegree() + 1)];
    uBSplineBasis1D->computeLocalBasisFunctionsAndDerivatives(uBSplinebasis1DFctsDerivs,
            _derivDegreeU, _uPrm, _KnotSpanIndexU);

    // Compute the derivatives of the local B-Spline basis functions in v-direction
    double* vBSplinebasis1DFctsDerivs = new double[(_derivDegreeV + 1)
            * (vBSplineBasis1D->getPolynomialDegree() + 1)];
    vBSplineBasis1D->computeLocalBasisFunctionsAndDerivatives(vBSplinebasis1DFctsDerivs,
            _derivDegreeV, _vPrm, _KnotSpanIndexV);

    // Initialize counters
    int counterBasis = 0;
    int counterDegree = 0;

    // Loop over all the non-zero contributions of the tensor product functions
    for (int j = 0; j <= vBSplineBasis1D->getPolynomialDegree(); j++) {
        for (int i = 0; i <= uBSplineBasis1D->getPolynomialDegree(); i++) {
            for (int k = 0; k <= _derivDegreeV; k++) {
                for (int l = 0; l <= _derivDegreeU; l++) {
                    if (k + l <= _maxMixDerivOrd) {
                        // Compute the functional value as a tensor product
                        _basisFctsAndDerivs[counterDegree * noBasisFcts + counterBasis] =
                                uBSplinebasis1DFctsDerivs[l
                                        * (uBSplineBasis1D->getPolynomialDegree() + 1) + i]
                                        * vBSplinebasis1DFctsDerivs[k
                                                * (vBSplineBasis1D->getPolynomialDegree() + 1) + j];

                        // Update the counter related to the computed derivative
                        counterDegree++;
                    }
                }
            }
            // Update the counter related to the individual basis function
            counterBasis++;

            // Reset the counter related to the computed derivative
            counterDegree = 0;
        }
    }

    // Free the memory from the pointers on the heap
    delete[] uBSplinebasis1DFctsDerivs;
    delete[] vBSplinebasis1DFctsDerivs;
}

void BSplineBasis2D::computeLocalBasisFunctionsAndDerivatives(double* _basisFctsAndDerivs,
        int _derivDegree, double _uPrm, int _KnotSpanIndexU, double _vPrm, int _KnotSpanIndexV) {
    /*
     *  Sorting idea for the 2D B-Spline basis functions
     *
     *  eta
     *   |
     *   |  CP(1,m) --> (m-1)*n+1   CP(2,m) --> (m-1)*n+2   ...     CP(n,m) --> n*m
     *   |
     *   |  ...                         ...                 ...     ...
     *   |
     *   |  CP(1,2) --> n+1         CP(2,2) --> n+2         ...     CP(n,2) --> 2*n
     *   |
     *   |  CP(1,1) --> 1           CP(2,1) --> 2           ...     CP(n,1) --> n
     *   |_______________________________________________________________________________xi
     *
     * On the input/output array _basisFctsAndDerivs:
     * _____________________________________________
     *
     * · The functional values are sorted in a 3D dimensional array which is in turned sorted in an 1D pointer array:
     *   _basisFctsAndDerivs = new double[(_derivDegree + 1) * (_derivDegree + 2) * noBasisFcts / 2]
     *   computing all the partial derivatives in u-direction k-th and in v-direction l-th as 0 <= k + l <= _derivDegree
     *
     * · The element _basisFctsAndDerivs[basisIndex] where index = indexDerivativeBasisFunction(_derivDegree,i,j,k) of the output array
     *   returns the i-th derivative with respect to u and the j-th derivative with respect to v (in total i+j-th partial derivative)
     *   of the k-th B-Spline basis function N_(k) sorted as in the above table
     *
     *  Reference: Piegl, Les and Tiller, Wayne, The NURBS Book. Springer-Verlag: Berlin 1995; p. 115.
     *  _________
     *
     *  It is also assumed that there are n basis functions in u-direction and m basis functions in v-direction and additionally that the highest derivative in
     *  u-direction is k-th order and in v-direction l-th order.
     */

    // Read input
    assert(_basisFctsAndDerivs!=NULL);

    // Get the polynomial degrees of the basis
    int pDegree = uBSplineBasis1D->getPolynomialDegree();
    int qDegree = vBSplineBasis1D->getPolynomialDegree();

    // Find the number of basis functions which affect the current knot span
    int noBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the index of the functional values into the 1D pointer array
    int index = 0;

    // Compute the derivatives of the local B-Spline basis functions in u-direction
    double* uBSplinebasis1DFctsDerivs = new double[(_derivDegree + 1) * (pDegree + 1)];
    uBSplineBasis1D->computeLocalBasisFunctionsAndDerivatives(uBSplinebasis1DFctsDerivs,
            _derivDegree, _uPrm, _KnotSpanIndexU);

    // Compute the derivatives of the local B-Spline basis functions in v-direction
    double* vBSplinebasis1DFctsDerivs = new double[(_derivDegree + 1) * (qDegree + 1)];
    vBSplineBasis1D->computeLocalBasisFunctionsAndDerivatives(vBSplinebasis1DFctsDerivs,
            _derivDegree, _vPrm, _KnotSpanIndexV);

    // Find the minimum between the derivative order and the polynomial degree in u-direction and store zeros for computational efficiency
    int du = min(_derivDegree, pDegree);

    // Loop over all the basis functions, over all derivatives in u-direction and then over all involved derivatives in v-direction
    for (int i = 0; i < noBasisFcts; i++)
        for (int k = pDegree + 1; k <= _derivDegree; k++)
            for (int l = 0; l <= _derivDegree - k; l++) {
                // Compute the index of the functional value
                index = indexDerivativeBasisFunction(_derivDegree, k, l, i);

                // Assign the value to zero
                _basisFctsAndDerivs[index] = 0.0;
            }

    // Find the minimum between the derivative order and the polynomial degree in v-direction and store zeros for computational efficiency
    int dv = min(_derivDegree, qDegree);

    // Loop over all the basis functions, over all derivatives in u-direction and then over all involved derivatives in v-direction
    for (int i = 0; i < noBasisFcts; i++)
        for (int l = qDegree + 1; l <= _derivDegree; l++)
            for (int k = 0; k <= _derivDegree - l; k++) {
                // Compute the index of the functional value
                index = indexDerivativeBasisFunction(_derivDegree, k, l, i);

                // Assign the value to zero
                _basisFctsAndDerivs[index] = 0.0;
            }

    // Initialize auxiliary variables
    int dd = 0;
    int counterBasis = 0;

    for (int vBasis = 0; vBasis <= qDegree; vBasis++) {
        for (int uBasis = 0; uBasis <= pDegree; uBasis++) {
            for (int k = 0; k <= du; k++) {
                dd = min(_derivDegree - k, dv);
                for (int l = 0; l <= dd; l++) {
                    // Compute the index of the functional value
                    index = indexDerivativeBasisFunction(_derivDegree, k, l, counterBasis);

                    // Compute the functional value in 2D as a tensor product of the basis functions in 1D
                    _basisFctsAndDerivs[index] = uBSplinebasis1DFctsDerivs[k * (pDegree + 1)
                            + uBasis] * vBSplinebasis1DFctsDerivs[l * (qDegree + 1) + vBasis];
                }
            }
            // Update the counter of the 2D B-Spline basis function
            counterBasis++;
        }
    }

    // Free the memory from the heap
    delete[] uBSplinebasis1DFctsDerivs;
    delete[] vBSplinebasis1DFctsDerivs;
}

void BSplineBasis2D::printPolynomialDegrees() {

    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis2D" << endl;
    cout << "BSplineBasis2D::PDegree = " << uBSplineBasis1D->getPolynomialDegree() << endl;
    cout << "BSplineBasis2D::QDegree = " << vBSplineBasis1D->getPolynomialDegree() << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

void BSplineBasis2D::printNoKnots() {

    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis2D" << endl;
    cout << "BSplineBasis2D::uNoKnots = " << uBSplineBasis1D->getNoKnots() << endl;
    cout << "BSplineBasis2D::vNoKnots = " << vBSplineBasis1D->getNoKnots() << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

void BSplineBasis2D::printKnotVectors() {
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis2D" << endl;
    cout << "BSplineBasis1D::uKnotVector = [\t";
    for (int i = 0; i < uBSplineBasis1D->getNoKnots(); i++) {
        cout << uBSplineBasis1D->getKnotVector()[i] << "\t";
    }
    cout << "]" << endl;
    cout << "BSplineBasis1D::vKnotVector = [\t";
    for (int i = 0; i < vBSplineBasis1D->getNoKnots(); i++) {
        cout << vBSplineBasis1D->getKnotVector()[i] << "\t";
    }
    cout << "]" << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

void BSplineBasis2D::printNoBasisFunctions() {
    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class BSplineBasis2D" << endl;
    cout << endl;
    cout << "BSplineBasis1D::uNoBasisFunctions = " << uBSplineBasis1D->computeNoBasisFunctions()
            << endl;
    cout << "BSplineBasis1D::vNoBasisFunctions = " << vBSplineBasis1D->computeNoBasisFunctions()
            << endl;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

}/* namespace EMPIRE */
