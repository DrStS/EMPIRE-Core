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
#include <limits>

// Inclusion of user defined libraries
#include "IGAPatchSurface.h"
#include "IGAMath.h"
#include "DataField.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

const int IGAPatchSurface::MAX_NUM_ITERATIONS = 20;

const double IGAPatchSurface::EPS_ORTHOGONALITY_CONDITION = 1e-9;

const double IGAPatchSurface::EPS_ORTHOGONALITY_CONDITION_RELAXED = 1e-7;

const double IGAPatchSurface::EPS_DISTANCE = 1e-9;

const double IGAPatchSurface::EPS_DISTANCE_RELAXED = 1e-6;


IGAPatchSurface::IGAPatchSurface(int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
        int _qDegree, int _vNoKnots, double* _vKnotVector, int _uNoControlPoints,
        int _vNoControlPoints, IGAControlPoint** _controlPointNet) :
        uNoControlPoints(_uNoControlPoints), vNoControlPoints(_vNoControlPoints) {

    // Read input
    bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
    bool vcondition = _vNoControlPoints != _vNoKnots - _qDegree - 1;

    if (ucondition || vcondition) {
        ERROR_OUT() << " in IGAPatchSurface::IGAPatchSurface" << endl;
        ERROR_OUT()
                << "Number of Control Points, number of knots and polynomial degree do not match!"
                << endl;
        exit(-1);
    }

    // Figure out whether the patch has a B-Spline or a NURBS underlying basis
    int isNurbs = 0;
    int counter = 0;
    for (int j = 0; j < vNoControlPoints; j++) {
        for (int i = 0; i < uNoControlPoints; i++) {
            if (_controlPointNet[counter]->getW() != 1.0) {
                isNurbs = 1;
                break;
            }
            // Update the counter
            counter++;
        }
    }

    // Create the NURBS or the B-Spline underlying basis
    if (!isNurbs) {
        IGABasis = new BSplineBasis2D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree,
                _vNoKnots, _vKnotVector);
    } else {
        double* controlPointWeights = new double[uNoControlPoints * vNoControlPoints];
        for (int i = 0; i < uNoControlPoints * vNoControlPoints; i++)
            controlPointWeights[i] = _controlPointNet[i]->getW();
        IGABasis = new NurbsBasis2D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector, _qDegree,
                _vNoKnots, _vKnotVector, _uNoControlPoints, _vNoControlPoints, controlPointWeights);
    }

    // On the Control Point net
    assert(_controlPointNet != NULL);
    ControlPointNet = _controlPointNet;
}

IGAPatchSurface::~IGAPatchSurface() {

    delete IGABasis;
    delete[] ControlPointNet;

}

double IGAPatchSurface::computePostprocessingScalarValue(double _u, double _v,
        double* _valuesOnCP) {
    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */
    // Read input
    assert(_valuesOnCP != NULL);

    int spanU = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
    int spanV = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
    double* localBasisFunctions = new double[noLocalBasisFunctions];

    IGABasis->computeLocalBasisFunctions(localBasisFunctions, _u, spanU, _v, spanV);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    double result = 0.0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (spanV - qDegree + j) * uNoControlPoints + (spanU - pDegree + i);

            //
            result += localBasisFunctions[counter_basis] * _valuesOnCP[CPindex];

            // Update basis function's counter
            counter_basis++;
        }
    }

    // Free the memory from the heap
    delete[] localBasisFunctions;

    return result;
}
void IGAPatchSurface::addTrimInfo(int* _knotSpanBelonging){
	int u=IGABasis->getUBSplineBasis1D()->getNoKnots();
	int v=IGABasis->getVBSplineBasis1D()->getNoKnots();
	Trimming.addTrimInfo(u,v,_knotSpanBelonging);
}
void IGAPatchSurface::addTrimLoop(int inner, int numCurves) {
    Trimming.addTrimLoop(inner, numCurves);
}

void IGAPatchSurface::addTrimCurve(int direction, int _pDegree, int _uNoKnots, double* _uKnotVector,
                  int _uNoControlPoints, double* _controlPointNet) {
                    
    int IDBasis = 0; ///???
    
    int numCPs = _uNoControlPoints;
    IGAControlPoint **cpNet;
    cpNet = new IGAControlPoint*[numCPs];
    
    for (int i = 0; i < numCPs; i++) {
            cpNet[i] = new IGAControlPoint(i, &_controlPointNet[i * 4]);
    }
    Trimming.addTrimCurve(direction, IDBasis, _pDegree, _uNoKnots, _uKnotVector,
                                               _uNoControlPoints, cpNet); 
}

void IGAPatchSurface::getUntrimmedCPindexes(std::set<int>& out) {
	const std::vector<std::vector<int> > knotSpan=Trimming.getKnotSpanInfo();

	int uNoKnots=IGABasis->getUBSplineBasis1D()->getNoKnots();
	int vNoKnots=IGABasis->getVBSplineBasis1D()->getNoKnots();

	for(int uSpan=0;uSpan<uNoKnots-1;uSpan++) {
		for(int vSpan=0;vSpan<vNoKnots-1;vSpan++) {
			int notOutside=knotSpan[uSpan][vSpan]>=0;
			if(notOutside) {
				addCPidsToSet(out,uSpan,vSpan);
			}
		}
	}
}

void IGAPatchSurface::addCPidsToSet(std::set<int>& CPids,const int uSpan, const int vSpan) {
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

	for(int p=uSpan-pDegree;p<=uSpan+1;p++ ) {
		for(int q=vSpan-qDegree;q<=vSpan+1;q++ ) {
           int CPindex = q * uNoControlPoints + p;
           int dofIndex=ControlPointNet[CPindex]->getDofIndex();
           CPids.insert(dofIndex);
		}
	}
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm,
        int _uKnotSpanIndex, double _vPrm, int _vKnotSpanIndex) {
    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Read input
    assert(_cartesianCoordinates != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);
    double* localBasisFunctions = new double[noLocalBasisFunctions];

    IGABasis->computeLocalBasisFunctions(localBasisFunctions, _uPrm, _uKnotSpanIndex, _vPrm,
            _vKnotSpanIndex);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Compute iteratively the x-coordinate of the point
            _cartesianCoordinates[0] += localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getX();
            // Compute iteratively the y-coordinate of the point
            _cartesianCoordinates[1] += localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getY();
            // Compute iteratively the z-coordinate of the point
            _cartesianCoordinates[2] += localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getZ();

            // Update basis function's counter
            counter_basis++;
        }
    }

    // Free the memory from the heap
    delete[] localBasisFunctions;
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localBasisFunctions, int _uKnotSpanIndex, int _vKnotSpanIndex) {
    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  It is also expected that the local basis functions have been precomputed outside the scope of this function and are given as arguments.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Read input
    assert(_cartesianCoordinates != NULL);
    assert(_localBasisFunctions != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Compute iteratively the x-coordinate of the point
            _cartesianCoordinates[0] += _localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getX();
            // Compute iteratively the y-coordinate of the point
            _cartesianCoordinates[1] += _localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getY();
            // Compute iteratively the z-coordinate of the point
            _cartesianCoordinates[2] += _localBasisFunctions[counter_basis]
                    * ControlPointNet[CPindex]->getZ();

            // Update basis function's counter
            counter_basis++;
        }
    }
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localCoordinates) {
    int _uKnotSpanIndex = IGABasis->getUBSplineBasis1D()->findKnotSpan(_localCoordinates[0]);
    int _vKnotSpanIndex = IGABasis->getVBSplineBasis1D()->findKnotSpan(_localCoordinates[1]);
    IGAPatchSurface::computeCartesianCoordinates(_cartesianCoordinates, _localCoordinates[0],
            _uKnotSpanIndex, _localCoordinates[1], _vKnotSpanIndex);
}

void IGAPatchSurface::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localBasisFctsAndDerivs, int _derivDegree, int _uKnotSpanIndex,
        int _vKnotSpanIndex) {
    /*
     *  Returns the Cartesian coordinates of a point on the 2D IGA patch whose surface parameters are _uPrm and _vPrm.
     *  It is also expected that the local basis functions have been precomputed outside the scope of this function and are given as arguments.
     *  The coordinates of the point are assumed on the 3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Read input
    assert(_cartesianCoordinates != NULL);
    assert(_localBasisFctsAndDerivs != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0.0;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

    // Initialize the Control Point index
    int CPindex = 0;

    // The basis function index
    int indexBasis = 0;

    // Initialize a basis functions counter
    int counter_basis = 0;

    // The derivative index so that it is exploited only the basis functions themselves
    int derivIndex = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Update the basis function index
            indexBasis = IGABasis->indexDerivativeBasisFunction(_derivDegree, derivIndex,
                    derivIndex, counter_basis);

            // Compute iteratively the x-coordinate of the point
            _cartesianCoordinates[0] += _localBasisFctsAndDerivs[indexBasis]
                    * ControlPointNet[CPindex]->getX();
            // Compute iteratively the y-coordinate of the point
            _cartesianCoordinates[1] += _localBasisFctsAndDerivs[indexBasis]
                    * ControlPointNet[CPindex]->getY();
            // Compute iteratively the z-coordinate of the point
            _cartesianCoordinates[2] += _localBasisFctsAndDerivs[indexBasis]
                    * ControlPointNet[CPindex]->getZ();

            // Update basis function's counter
            counter_basis++;
        }
    }
}

void IGAPatchSurface::computeBaseVectors(double* _baseVectors,
        double* _localBasisFunctionsAndDerivatives, int _uKnotSpanIndex, int _vKnotSpanIndex) {
    /*
     *  Returns the base vectors at a given pair of surface parameters on the NURBS 2D patch. The function expects the computation of the
     *  basis functions and their derivatives to have been carried out outside the function and given as arguments.
     *
     *  The sorting idea is:
     *
     *                 | g1x g1y g1z |
     *  _baseVectors = | g2x g2y g2z |
     *
     *  The output is sorted in an 1D pointer array as follows:
     *
     *  _baseVectors = [g1x g1y g1z g2x g2y g2z]
     *
     *  E.g. _baseVectors[i * 3 + j] returns the j-th coordinate, j=0,…,2 (x,y,z respectively) of the i-th base vector, i=0,…,1 (with respect to
     *  u-parametric line and v-parametric line respectively)
     */

    // Read input
    assert(_baseVectors != NULL);
    assert(_localBasisFunctionsAndDerivatives != NULL);

    // Number of tangent base vectors
    int noBaseVec = 2;

    // Number of coordinates
    int noCoordinates = 3;

    // Initialize the coordinates of the base vectors
    for (int i = 0; i < noBaseVec * noCoordinates; i++)
        _baseVectors[i] = 0;

    // Initialize the index of the derivative functional value
    int indexBasis = 0;

    // Order of derivatives for the basis functions needed for the computation of the base vectors
    int derivDegree = 1;

    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Control Point index
    int CPindex = 0;

    // Initialize a basis functions counter
    int counterBasis = 0;

    // Initialize the derivative counter
    int counterBaseVector = 0;

    // Loop over all the non-zero contributions
    for (int j = 0; j <= qDegree; j++) {
        for (int i = 0; i <= pDegree; i++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            CPindex = (_vKnotSpanIndex - qDegree + j) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + i);

            // Reset the base vector's counter
            counterBaseVector = 0;

            // Loop over all the first derivatives of the basis functions and sum up the contributions
            for (int l = 0; l <= derivDegree; l++) {
                for (int k = 0; k <= derivDegree - l; k++) {

                    // Skip the case when k=l=0 i.e. the case where the basis functions themselves will be issued
                    if (k == l)
                        continue;

                    // Get the index of the derivative functional value
                    indexBasis = getIGABasis()->indexDerivativeBasisFunction(derivDegree, k, l,
                            counterBasis);

                    // Compute iteratively the x-coordinate of the point
                    _baseVectors[counterBaseVector * noCoordinates + 0] +=
                            _localBasisFunctionsAndDerivatives[indexBasis]
                                    * ControlPointNet[CPindex]->getX();
                    // Compute iteratively the y-coordinate of the point
                    _baseVectors[counterBaseVector * noCoordinates + 1] +=
                            _localBasisFunctionsAndDerivatives[indexBasis]
                                    * ControlPointNet[CPindex]->getY();
                    // Compute iteratively the z-coordinate of the point
                    _baseVectors[counterBaseVector * noCoordinates + 2] +=
                            _localBasisFunctionsAndDerivatives[indexBasis]
                                    * ControlPointNet[CPindex]->getZ();

                    // Update the base vector counter
                    counterBaseVector++;
                }
            }

            // Update basis function's counter
            counterBasis++;
        }
    }
}

void IGAPatchSurface::computeBaseVectors(double* _baseVectors, double _u, int _spanU, double _v,
        int _spanV) {

    // Get the polynomial degree of the basis in each direction
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
    int noLocalBasisFunctions = (pDegree + 1) * (qDegree + 1);

    // Derivative order of the basis functions needed for the computation of the base vectors
    int derivDegree = 1;

    // Compute the local basis functions and their derivatives
    double* localBasisFunctionsAndDerivatives = new double[(derivDegree + 1) * (derivDegree + 2)
            * noLocalBasisFunctions / 2];
    IGABasis->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
            derivDegree, _u, _spanU, _v, _spanV);

    // Compute the base vectors at the given surface parameters
    computeBaseVectors(_baseVectors, localBasisFunctionsAndDerivatives, _spanU, _spanV);
}

int IGAPatchSurface::indexDerivativeBaseVector(int _derivDegree, int _uDerivIndex, int _vDerivIndex,
        int _componentIndex, int _baseVecIndex) {
    /*
     * Returns the correct index when sorting the base vectors and the their derivatives in an 1D pointer array with the rule:
     * (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1) * noCoordinates * noBaseVec / 2 +
     *  + _uDerivIndex * noCoordinates * noBaseVec + _componentIndex * noBaseVec + _baseVecIndex
     */

    // Read input
    if (_uDerivIndex + _vDerivIndex > _derivDegree) {
        ERROR_OUT() << "in IGAPatchSurface::indexDerivativeBaseVector" << endl;
        ERROR_OUT() << "It has been requested the " << _uDerivIndex
                << "-th partial derivative w.r.t. u and the " << _vDerivIndex
                << "-th partial derivative w.r.t. v of the base vectors but " << endl;
        ERROR_OUT() << "the maximum absolute derivative selected is of " << _derivDegree
                << "-th order" << endl;
        exit(-1);
    }

    // Number of tangent base vectors
    int noBaseVec = 2;

    // Number of coordinates
    int noCoordinates = 3;

    // Compute the index of the functional value
    return (_derivDegree - _vDerivIndex) * (_derivDegree - _vDerivIndex + 1) * noCoordinates
            * noBaseVec / 2 + _uDerivIndex * noCoordinates * noBaseVec + _componentIndex * noBaseVec
            + _baseVecIndex;
}

void IGAPatchSurface::computeBaseVectorsAndDerivatives(double* _baseVectorsAndDerivatives,
        double*_localBasisFunctionsAndDerivatives, int _derivDegree, int _uKnotSpanIndex,
        int _vKnotSpanIndex) {

    /*
     * Returns the base vectors and their j-th derivative w.r.t. v-parametric line and i-th derivative w.r.t. u-parametric line partial derivatives.
     * On the input arrays the following rules must be fulfiled:
     *
     * · The array _baseVectorsAndDerivatives = new double[(_derivDegree + 1) * (_derivDegree + 2) * noCoord * noBaseVec / 2] will be filled up
     *   with the components of the base vectors and their derivatives up to absolute order _derivDegree meaning that 0 <= i + j <= _derivDegree
     *
     * · The array _localBasisFunctionsAndDerivatives = new double[(_derivDegree + 2) * (_derivDegree + 3) * noBasisFcts / 2] contains the IGA basis
     *   functions and their derivatives up to absolute order _derivDegree + 1 since the base vectors contain by definition the first derivatives of
     *   IGA basis functions
     */

    // Read input
    assert(_baseVectorsAndDerivatives != NULL);
    assert(_localBasisFunctionsAndDerivatives != NULL);

    // The polynomial degrees of the patch
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Initialize the index related to the basis functions
    int indexBasis = 0;

    // Initialize the index related to the CP's
    int indexCP = 0;

    // Initialize the index related to the derivatives of the base vectors
    int indexBaseVct = 0;

    // The derivative order of the basis functions
    int derivDegreeBasis = _derivDegree + 1;

    // Initialize factor by which to multiply the derivative of the basis function
    double factor = 0.0;

    // The number of coordinates for the base vectors
    int noCoordinates = 3;

    // The number of the base vectors to be returned
    int noBaseVct = 2;

    // Initialize the output array
    int counterBaseVec = 0;
    for (int i = 0; i <= _derivDegree; i++) {
        for (int j = 0; j <= _derivDegree - i; j++) {
            for (int k = 0; k < noCoordinates; k++) {
                for (int l = 0; l < noBaseVct; l++) {
                    _baseVectorsAndDerivatives[counterBaseVec] = 0.0;
                    counterBaseVec++;
                }
            }
        }
    }

    // Initialize the counter of the basis functions
    int counterBasis = 0;

    // Loop over all the non-zero contributions in v-direction
    for (int vBasis = 0; vBasis <= qDegree; vBasis++) {
        // Loop over all the non-zero contributions in u-direction
        for (int uBasis = 0; uBasis <= pDegree; uBasis++) {

            // Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
            indexCP = (_vKnotSpanIndex - qDegree + vBasis) * uNoControlPoints
                    + (_uKnotSpanIndex - pDegree + uBasis);

            // Loop over all the derivatives in u-direction
            for (int i = 0; i <= _derivDegree; i++) {
                // Loop over all the derivatives in v-direction
                for (int j = 0; j <= _derivDegree - i; j++) {
                    // Loop over all the coordinates
                    for (int k = 0; k < noCoordinates; k++) {
                        // Loop over all the base vectors
                        for (int l = 0; l < noBaseVct; l++) {
                            // Compute the index related to the derivatives of the base vectors
                            indexBaseVct = indexDerivativeBaseVector(_derivDegree, i, j, k, l);

                            // · Case l = 0 :: --> i + 1 - l = i + 1 && j + l --> j giving the derivative of the base vector g1 = dX/du
                            // · Case l = 1 :: --> i + 1 - l = i && j + l --> j + 1 giving the derivative of the base vector g2 = dX/dv

                            // Compute the index of the basis functions and their derivatives
                            indexBasis = IGABasis->indexDerivativeBasisFunction(derivDegreeBasis,
                                    i + 1 - l, j + l, counterBasis);

                            // Factor by which to multiply the derivative of the basis function
                            if (k == 0)
                                factor = ControlPointNet[indexCP]->getX();
                            else if (k == 1)
                                factor = ControlPointNet[indexCP]->getY();
                            else
                                factor = ControlPointNet[indexCP]->getZ();

                            // Add the contribution from each basis function in the interval
                            _baseVectorsAndDerivatives[indexBaseVct] +=
                                    _localBasisFunctionsAndDerivatives[indexBasis] * factor;
                        }
                    }
                }
            }
            // Update basis function's counter
            counterBasis++;
        }
    }
}

bool IGAPatchSurface::computePointProjectionOnPatch(double& _u, double& _v, double* _P,
        bool& _flagConverge) {

    /*
     * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
     * _P = double[3]
     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
     *
     * Function layout :
     *
     * 1. Read input and initialize the data
     *
     * 2. Loop over all the Newton-Raphson iterations
     *    2i. Update the iteration counter
     *   2ii. Find the span of the given surface parameters
     *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
     *   2iv. Compute the Cartesian components of the point on the surface
     *    2v. Compute the distance vector between the vector to be projected and the estimated one
     *   2vi. Compute the 2-norm of the distance vector
     *  2vii. Compute the base vectors and their derivatives
     * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
     *   2ix. Compute the cosine of the angle with respect to u-parametric line
     *    2x. Compute the cosine of the angle with respect to v-parametric line
     *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
     *  2xii. Compute the entries of the Jacobian matrix
     * 2xiii. Compute the entries of the right-hand side vector
     *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
     *   2xv. Update the surface parameters u += du and v += dv
     *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
     *
     * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
     *
     * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
     */

    // 1. Read input and initialize the data
    // Read input
    assert(_P != NULL);

    // Initialize the flag to true
    bool flagNewtonRaphson = true;
    _flagConverge = true;

    const double epsJ = 1e-6;
    const double epsDuv = 1e-10;
    bool fixU = false;
    bool fixV = false;

    // Initialize number of spatial dimensions
    int noSpatialDimensions = 3;

    // Save the location of the point to be projected
    double point[noSpatialDimensions];
    for (int i = 0; i < noSpatialDimensions; i++)
        point[i] = _P[i];

    // Initialize the distance vector
    double distanceVector[3];

    // Initialize the indices of the base vectors and their derivatives
    int indexGu = 0;
    int indexGv = 0;
    int indexDGuDu = 0;
    int indexDGvDv = 0;
    int indexDGuDv = 0;
    int indexDGvDu = indexDGuDv;

    // Initialize the base vectors and their derivatives
    double Gu[3];
    double Gv[3];
    double DGuDu[3];
    double DGvDv[3];
    double DGuDv[3];
    double* DGvDu=DGuDv;

    // Initialize the dot products
    double distanceVector2norm = 0.0;
    double GuXdistanceVector = 0.0;
    double squareGu2norm = 0.0;
    double Gu2norm = 0.0;
    double GvXdistanceVector = 0.0;
    double squareGv2norm = 0.0;
    double Gv2norm = 0.0;

    // Initialize Jacobian matrix
    double dR[4];

    // Initialize right-hand side solution vector
    double R[2];

    // Initialize flag on the solvability of the 2x2 equation system
    bool flagLinearSystem = 1;

    // Initialize the cosines w.r.t. each parametric line
    double cosu = 0.0;
    double cosv = 0.0;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // The lengths of the knot vectors to the NURBS patch
    int lengthUKnotVct = IGABasis->getUBSplineBasis1D()->getNoKnots();
    int lengthVKnotVct = IGABasis->getVBSplineBasis1D()->getNoKnots();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the IGA basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
            * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
            * noSpatialDimensions * noBaseVcts / 2];

    // 2. Loop over all the Newton-Raphson iterations
    while (counter <= MAX_NUM_ITERATIONS) {
        // 2i. Update the iteration counter
        counter++;

        // 2ii. Find the span of the given surface parameters
        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, _u,
                uKnotSpan, _v, vKnotSpan);

        // 2iv. Compute the Cartesian components of the point on the surface
        computeCartesianCoordinates(_P, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);

        // 2v. Compute the distance vector between the vector to be projected and the estimated one
        for (int i = 0; i < noSpatialDimensions; i++)
            distanceVector[i] = _P[i] - point[i];

        // 2vi. Compute the 2-norm of the distance vector
        distanceVector2norm = square2normVector(noSpatialDimensions, distanceVector);
        distanceVector2norm = sqrt(distanceVector2norm);

        if (distanceVector2norm < EPS_DISTANCE)
            break;

        // 2vii. Compute the base vectors and their derivatives
        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
                uKnotSpan, vKnotSpan);

        // 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
        for (int i = 0; i < noSpatialDimensions; i++) {
            // On the base vector Gu = dR/du
            indexGu = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0);
            Gu[i] = baseVecAndDerivs[indexGu];

            // On the base vector Gv = dR/dv
            indexGv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1);
            Gv[i] = baseVecAndDerivs[indexGv];

            // On the derivative of the base vector Gu w.r.t. u namely dGu/du = d^2R/du^2
            indexDGuDu = indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0);
            DGuDu[i] = baseVecAndDerivs[indexDGuDu];

            // On the derivative of the base vector Gv w.r.t. u namely dGv/dv = d^2R/dv^2
            indexDGvDv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 1);
            DGvDv[i] = baseVecAndDerivs[indexDGvDv];

            // On the mixed derivative of the base vectors namely d^2Gu/dv = d^2Gv/du = d^2R/dudv
            indexDGuDv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0);
            DGuDv[i] = baseVecAndDerivs[indexDGuDv];
        }

        // 2ix. Compute the cosine of the angle with respect to u-parametric line
        GuXdistanceVector = dotProduct(noSpatialDimensions, Gu, distanceVector);
        squareGu2norm = square2normVector(noSpatialDimensions, Gu);
        Gu2norm = sqrt(squareGu2norm);
        cosu = fabs(GuXdistanceVector) / Gu2norm / distanceVector2norm;

        // 2x. Compute the cosine of the angle with respect to v-parametric line
        GvXdistanceVector = dotProduct(noSpatialDimensions, Gv, distanceVector);
        squareGv2norm = square2normVector(noSpatialDimensions, Gv);
        Gv2norm = sqrt(squareGv2norm);
        cosv = fabs(GvXdistanceVector) / Gv2norm / distanceVector2norm;

        // 2xi. Check the orthogonality condition and if it is fulfilled break the loop

        if (cosu <= EPS_ORTHOGONALITY_CONDITION && cosv <= EPS_ORTHOGONALITY_CONDITION)
            break;

        // 2xii. Compute the entries of the Jacobian matrix
        dR[0] = squareGu2norm
                + dotProduct(noSpatialDimensions, DGuDu, distanceVector);
        dR[1] = dotProduct(noSpatialDimensions, Gu, Gv)
                + dotProduct(noSpatialDimensions, DGuDv, distanceVector);
        dR[2] = dR[1];
        dR[3] = squareGv2norm
                + dotProduct(noSpatialDimensions, DGvDv, distanceVector);

        // 2xiii. Compute the entries of the right-hand side vector
        R[0] = -dotProduct(noSpatialDimensions, Gu, distanceVector);
        R[1] = -dotProduct(noSpatialDimensions, Gv, distanceVector);

        if (fabs(dR[0]) < epsJ || fixU) {
            R[0] = 0.0;
            R[1] = R[1] / dR[3];
            fixU = false;
            fixV = true;
        } else if (fabs(dR[3]) < epsJ || fixV) {
            R[0] = R[0] / dR[1];
            R[1] = 0.0;
            fixU = true;
            fixV = false;
        } else {

            // 2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved

            // Solve the equation system
            flagLinearSystem = solve2x2linearSystem(R, dR);

            // Check if the equation system has been successfully solved
            if (!flagLinearSystem) {
                ERROR_OUT() << "Error in IGAPatchSurface::computePointProjectionOnPatch" << endl;
                ERROR_OUT()
                        << "The 2x2 equation system to find the updates of the surface parameters"
                        << endl;
                ERROR_OUT()
                        << "for the orthogonal projection of a point on the NURBS patch has been"
                        << endl;
                ERROR_OUT() << "detected not solvable up to tolerance" << EPS << endl;
                exit(-1);
            }
        }
        // 2xv. Update the surface parameters u += du and v += dv
        _u += R[0];
        _v += R[1];

        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
    	IGABasis->getUBSplineBasis1D()->clampKnot(_u);
    	IGABasis->getVBSplineBasis1D()->clampKnot(_v);
    }

////     3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
    if (counter > MAX_NUM_ITERATIONS) {
        if (cosu <= EPS_ORTHOGONALITY_CONDITION_RELAXED
                && cosv <= EPS_ORTHOGONALITY_CONDITION_RELAXED)
            flagNewtonRaphson = true;
        else
            flagNewtonRaphson = false;
        if (!flagNewtonRaphson && distanceVector2norm < EPS_DISTANCE_RELAXED)
            flagNewtonRaphson = true;
        if (R[0] * R[0] + R[1] * R[1] < epsDuv)
            _flagConverge = true;
        else
            _flagConverge = false;
    } else {
    }
    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
    // Clear the memory on the heap
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return the flag
    return flagNewtonRaphson;
}

bool IGAPatchSurface::computePointProjectionOnPatch(double& _u, double& _v, double* _P) {
    /*
     *  Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
     * _P = double[3]. Its return value is a bool flag on the convergence of the Newton-Raphson iterations. Makes use of the previosuly
     * implemented function IGAPatchSurface::computePointProjectionOnPatch
     */

    // Initialize the boolean on the convergence of the Newton iterations
    bool tmp = false;

    // Compute the closest point projection using the Newton-Rapshon algorithm
    return computePointProjectionOnPatch(_u, _v, _P, tmp);
}

bool IGAPatchSurface::computePointProjectionOnPatchBoundaryOnGivenEdge_Brute(
		double& _u, double& _v, double& _ratio, double& _distance, double* _P1,
		double* _P2) {

	/*
	 * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
	 * _P = double[3]
	 * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
	 *
	 * Function layout :
	 *
	 * 1. Read input and initialize the data
	 *
	 * 2. Loop over all the Newton-Raphson iterations
	 *    2i. Update the iteration counter
	 *   2ii. Find the span of the given surface parameters
	 *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
	 *   2iv. Compute the Cartesian components of the point on the surface
	 *    2v. Compute the distance vector between the vector to be projected and the estimated one
	 *   2vi. Compute the 2-norm of the distance vector
	 *  2vii. Compute the base vectors and their derivatives
	 * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
	 *   2ix. Compute the cosine of the angle with respect to u-parametric line
	 *    2x. Compute the cosine of the angle with respect to v-parametric line
	 *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
	 *  2xii. Compute the entries of the Jacobian matrix
	 * 2xiii. Compute the entries of the right-hand side vector
	 *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
	 *   2xv. Update the surface parameters u += du and v += dv
	 *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
	 *
	 * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
	 *
	 * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
	 */

	// 1. Read input and initialize the data
	// Read input
	assert(_P1 != NULL);
	assert(_P2 != NULL);

	// Initialize the flag to true
	bool flagNewtonRaphson = true;
	// Initialize number of spatial dimensions
	int noSpatialDimensions = 3;

	// Setup initial value depending on the edge
	double u, v;

	// Initialize the knot span indices
	int uKnotSpan = 0;
	int vKnotSpan = 0;

	// The NURBS polynomial degrees
	int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
	int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

	// Local number of basis functions
	int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

	// Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
	int derivDegreeBasis = 2;

	// The number of the base vectors
	int noBaseVcts = 2;

	// Initialize the array of the NURBS basis functions and their derivatives
	double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1)
			* (derivDegreeBasis + 2) * noLocalBasisFcts / 2];

	// Number of derivatives for the base vectors
	int derivDegreeBaseVcts = derivDegreeBasis - 1;

	// Initialize the array of the base vectors and their derivatives
	double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1)
			* (derivDegreeBaseVcts + 2) * noSpatialDimensions * noBaseVcts / 2];

	double G1[3], G2[3];
	double surfaceNormal[3];
	double P_tmp[3];

	double d, d1, d2;
	double ratio1, ratio2;

	double u1, u2, v1, v2;

	double t = 0.5;
	bool isIn;

	double w, w1, w2;
	// Define extremity 1
	double UV1[2] = { _u, _v };
	double P1[3];
	//Define extremity 2
	double P2[3];

	//Define moving point
	double UV[2] = { _u, _v };
	double P[3];

	for (int i = 0; i < noSpatialDimensions; i++) {
		P1[i] = _P1[i];
		P2[i] = _P2[i];
	}

	int iteration = 0;
	DEBUG_OUT()<<"\tUV1[0]="<<UV[0]<<" , UV1[1]="<<UV[1]<<endl;
	double P_dist[3];
	double u0 = getIGABasis()->getUBSplineBasis1D()->getFirstKnot();
	double uN = getIGABasis()->getUBSplineBasis1D()->getLastKnot();
	double v0 = getIGABasis()->getVBSplineBasis1D()->getFirstKnot();
	double vN = getIGABasis()->getVBSplineBasis1D()->getLastKnot();
	bool dontStop1 = UV1[0] - u0 > EPS_DISTANCE && uN - UV1[0] > EPS_DISTANCE
			&& UV1[1] - v0 > EPS_DISTANCE && vN - UV1[1] > EPS_DISTANCE;
	bool dontStop2 = true;
	double uRel, vRel;
	do {
		iteration++;
		for (int i = 0; i < noSpatialDimensions; i++) {
			P[i] = (1 - t) * P1[i] + t * P2[i];
			P_dist[i] = P[i] - P1[i];
			P_tmp[i] = P[i];
		}

		isIn = computePointProjectionOnPatch(UV[0], UV[1], P_tmp);

		if (isIn) {
			for (int i = 0; i < noSpatialDimensions; i++) {
				P1[i] = P[i];
				P[i] = P_tmp[i];
			}
			uRel = UV1[0] - UV[0];
			vRel = UV1[1] - UV[1];
			dontStop2 = uRel > EPS_DISTANCE && vRel > EPS_DISTANCE;
			UV1[0] = UV[0];
			UV1[1] = UV[1];
			dontStop1 = UV1[0] - u0 > EPS_DISTANCE && uN - UV1[0] > EPS_DISTANCE
					&& UV1[1] - v0 > EPS_DISTANCE && vN - UV1[1] > EPS_DISTANCE;
				   DEBUG_OUT()<<"\tUV[0]="<<UV[0]<<" , UV[1]="<<UV[1]<<endl;

		} else {
			for (int i = 0; i < noSpatialDimensions; i++)
				P2[i] = P[i];
			UV[0] = UV1[0];
			UV[1] = UV1[1];
		}

	} while (sqrt(square2normVector(noSpatialDimensions, P_dist)) > EPS_DISTANCE
			&& iteration < 40);

	double P1P[noSpatialDimensions], P1P2[noSpatialDimensions];
	for (int i = 0; i < noSpatialDimensions; i++) {
		P1P[i] = P1[i] - _P1[i];
		P1P2[i] = _P2[i] - _P1[i];
	}
	if (iteration >= 40)
		return false;
	_u = UV1[0];
	_v = UV1[1];
	_distance = distancePointSegment(P, _P1, _P2);
	_ratio = (dotProduct(noSpatialDimensions, P1P2, P1P))
			/ square2normVector(noSpatialDimensions, P1P2);
	return true;

}

//bool IGAPatchSurface::computePointProjectionOnPatchBoundaryOnGivenEdge(double& _w, double& _ratio,
//        double& _distance, double* _P1, double* _P2, int _edge) {
//
//    /*
//     * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
//     * _P = double[3]
//     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
//     *
//     * Function layout :
//     *
//     * 1. Read input and initialize the data
//     *
//     * 2. Loop over all the Newton-Raphson iterations
//     *    2i. Update the iteration counter
//     *   2ii. Find the span of the given surface parameters
//     *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
//     *   2iv. Compute the Cartesian components of the point on the surface
//     *    2v. Compute the distance vector between the vector to be projected and the estimated one
//     *   2vi. Compute the 2-norm of the distance vector
//     *  2vii. Compute the base vectors and their derivatives
//     * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
//     *   2ix. Compute the cosine of the angle with respect to u-parametric line
//     *    2x. Compute the cosine of the angle with respect to v-parametric line
//     *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
//     *  2xii. Compute the entries of the Jacobian matrix
//     * 2xiii. Compute the entries of the right-hand side vector
//     *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
//     *   2xv. Update the surface parameters u += du and v += dv
//     *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
//     *
//     * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
//     *
//     * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
//     */
//
//    // 1. Read input and initialize the data
//    // Read input
//    assert(_P1 != NULL);
//    assert(_P2 != NULL);
//
//    // Initialize the flag to true
//    bool flagNewtonRaphson = true;
//    // Initialize number of spatial dimensions
//    int noSpatialDimensions = 3;
//
//	// Setup initial value depending on the edge
//	int direction;
//    double u,v;
//	switch (_edge) {
//	case 0:
//		u = _w;
//		v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
//		direction = 0;
//		break;
//	case 1:
//		u = _w;
//		v = IGABasis->getVBSplineBasis1D()->getLastKnot();
//		direction = 0;
//		break;
//	case 2:
//		u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
//		v = _w;
//		direction = 1;
//		break;
//	case 3:
//		u = IGABasis->getUBSplineBasis1D()->getLastKnot();
//		v = _w;
//		direction = 1;
//		break;
//	}
//
//	// Initialize the knot span indices
//	int uKnotSpan = 0;
//	int vKnotSpan = 0;
//
//	// The NURBS polynomial degrees
//	int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
//	int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
//
//	// Local number of basis functions
//	int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);
//
//	// Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
//	int derivDegreeBasis = 2;
//
//	// The number of the base vectors
//	int noBaseVcts = 2;
//
//	// Initialize the array of the NURBS basis functions and their derivatives
//	double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
//			* noLocalBasisFcts / 2];
//
//	// Number of derivatives for the base vectors
//	int derivDegreeBaseVcts = derivDegreeBasis - 1;
//
//	// Initialize the array of the base vectors and their derivatives
//	double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
//			* noSpatialDimensions * noBaseVcts / 2];
//
//	double G1[3],G2[3];
//	double surfaceNormal[3];
//	double P_tmp[3];
//
//    double d,d1,d2;
//    double ratio1,ratio2;
//
//    double u1,u2,v1,v2;
//
//
//    double w,w1,w2;
//    // Define extremity 1
//    w1=_w;
//    computePointMinimumDistanceToPatchBoundaryOnGivenEdge(w1,d1,_P1,_edge);
//    DEBUG_OUT()<<"\td1="<<d1<<" , w1="<<w1<<endl;
//   double UV1[2];
//   if(direction==0) {
//	   UV1[0]=w1;
//   	   UV1[1]=v;
//   } else {
//	   UV1[0]=u;
//   	   UV1[1]=w1;
//   }
//   double P1[3];
//   computeCartesianCoordinates(P1,UV1);
//   d1=distancePointSegment(P1,_P1,_P2);
//
//   DEBUG_OUT()<<"\td1="<<d1<<" , w1="<<w1<<endl;
//
//   // Define extremity 2
//   w2=_w;
//   computePointMinimumDistanceToPatchBoundaryOnGivenEdge(w2,d2,_P2,_edge);
//   DEBUG_OUT()<<"\td2="<<d2<<" , w2="<<w2<<endl;
//
//   double UV2[2];
//   if(direction==0) {
//	   UV2[0]=w2;
//   	   UV1[1]=v;
//   } else {
//	   UV2[0]=u;
//   	   UV2[1]=w2;
//   }
//   double P2[3];
//   computeCartesianCoordinates(P2,UV2);
//   d2=distancePointSegment(P2,_P1,_P2);
//
//   //Define moving point
//   double UV[2];
//   w=(w1+w2)/2;
//   if(direction==0) {
//	   UV[0]=(w1+w2)/2;
//   	   UV[1]=v;
//   } else {
//	   UV[0]=u;
//   	   UV[1]=(w1+w2)/2;
//   }
//   double P[3];
//   computeCartesianCoordinates(P,UV);
//   d=distancePointSegment(P,_P1,_P2);
//   DEBUG_OUT()<<"\td="<<d<<" , w="<<w<<endl;
//
//   int iteration=0;
//    while(fabs(d1-d2)/d1>1e-4 && iteration<40) {
//    	DEBUG_OUT()<<fabs(d1-d2)/d1<<endl;
//    	iteration++;
//    	if(d1<d2) {
//    		d2=d;
//    		w2=w;
//    		UV2[0]=UV[0];
//    		UV2[1]=UV[1];
//    	} else {
//    		d1=d;
//    		w1=w;
//    		UV1[0]=UV[0];
//    		UV1[1]=UV[1];
//    	}
//    	   w=(w1+w2)/2;
//	   if(direction==0) {
//		   UV[0]=(w1+w2)/2;
//		   UV[1]=v;
//	   } else {
//		   UV[0]=u;
//		   UV[1]=(w1+w2)/2;
//	   }
//	   computeCartesianCoordinates(P,UV);
//	   d=distancePointSegment(P,_P1,_P2);
//	   DEBUG_OUT()<<"\td="<<d<<" , w="<<w<<endl;
//    }
//
//	double P1P[noSpatialDimensions],P1P2[noSpatialDimensions];
//	for (int i = 0; i < noSpatialDimensions; i++) {
//		P1P[i] = P[i ]-_P1[i];
//		P1P2[i]=_P2[i]-_P1[i];
//	}
//	if(iteration>=40)
//		return false;
//    _w=w;
//    _distance=d;
//    _ratio=(dotProduct(noSpatialDimensions,P1P2,P1P))/sqrt(square2normVector(noSpatialDimensions,P1P2));
//    return true;
//}

//bool IGAPatchSurface::computePointProjectionOnPatchBoundaryOnGivenEdge_MinDist(double& _w, double& _ratio,
//        double& _distance, double* _P1, double* _P2, int _edge) {
//
//    /*
//     * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
//     * _P = double[3]
//     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
//     *
//     * Function layout :
//     *
//     * 1. Read input and initialize the data
//     *
//     * 2. Loop over all the Newton-Raphson iterations
//     *    2i. Update the iteration counter
//     *   2ii. Find the span of the given surface parameters
//     *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
//     *   2iv. Compute the Cartesian components of the point on the surface
//     *    2v. Compute the distance vector between the vector to be projected and the estimated one
//     *   2vi. Compute the 2-norm of the distance vector
//     *  2vii. Compute the base vectors and their derivatives
//     * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
//     *   2ix. Compute the cosine of the angle with respect to u-parametric line
//     *    2x. Compute the cosine of the angle with respect to v-parametric line
//     *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
//     *  2xii. Compute the entries of the Jacobian matrix
//     * 2xiii. Compute the entries of the right-hand side vector
//     *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
//     *   2xv. Update the surface parameters u += du and v += dv
//     *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
//     *
//     * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
//     *
//     * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
//     */
//
//    // 1. Read input and initialize the data
//    // Read input
//    assert(_P1 != NULL);
//    assert(_P2 != NULL);
//
//    // Initialize the flag to true
//    bool flagNewtonRaphson = true;
//
//
//
//    // Initialize number of spatial dimensions
//    int noSpatialDimensions = 3;
//
//    // Initialize line data
//    double P1P2[noSpatialDimensions];
//    for (int i = 0; i < noSpatialDimensions; i++)
//        P1P2[i]=_P2[i]-_P1[i];
//    double t=1.0;
//    double P[noSpatialDimensions];
//    for (int i = 0; i < noSpatialDimensions; i++)
//    	P[i]=(1-t)*_P1[i]+t*_P2[i];
//
//    // Initialize the distance vector
//    double distanceVector[3];
//    // Save the location of the point to be projected
//    double Q[noSpatialDimensions];
//    // Initialize the indices of the base vectors and their derivatives
//    int indexUBaseVec = 0;
//    int indexVBaseVec = 0;
//    int indexdUBaseVecdu = 0;
//    int indexdVBaseVecdv = 0;
//    int indexdBaseVecMix = 0;
//
//    // Initialize the base vectors and their derivatives
//    double uBaseVec[3];
//    double vBaseVec[3];
//    double uBaseVecdu[3];
//    double vBaseVecdv[3];
//    double dBaseVecMix[3];
//
//    // Initialize the dot products
//    double distanceVector2norm = 0.0;
//    double uBaseVecXdistanceVector = 0.0;
//    double squareUBaseVec2norm = 0.0;
//    double uBaseVec2norm = 0.0;
//    double vBaseVecXdistanceVector = 0.0;
//    double squareVBaseVec2norm = 0.0;
//    double vBaseVec2norm = 0.0;
//
//    // Initialize Jacobian matrix
//    double jacobianMatrix[4];
//
//    // Initialize right-hand side solution vector
//    double rightSideVct[2];
//
//    // Initialize the cosines w.r.t. each parametric line
//    double cosu = 1.0;
//    double cosv = 1.0;
//    double cosW = 1.0;
//    double cosT = 1.0;
//    double cosW_prev=1.0;
//    double cosT_prev=1.0;
//
//
//    // Initialize the knot span indices
//    int uKnotSpan = 0;
//    int vKnotSpan = 0;
//
//    // The NURBS polynomial degrees
//    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
//    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
//
//    // Local number of basis functions
//    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);
//
//    // Initialize the Newton-Raphson iteration counter
//    int counter = 0;
//
//    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
//    int derivDegreeBasis = 2;
//
//    // The number of the base vectors
//    int noBaseVcts = 2;
//
//    // Initialize the array of the IGA basis functions and their derivatives
//    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
//            * noLocalBasisFcts / 2];
//
//    // Number of derivatives for the base vectors
//    int derivDegreeBaseVcts = derivDegreeBasis - 1;
//
//    // Initialize the array of the base vectors and their derivatives
//    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
//            * noSpatialDimensions * noBaseVcts / 2];
//
//    // Setup initial value depending on the edge
//    double u;
//    double v;
//    int direction;
//    switch (_edge) {
//    case 0:
//        u = _w;
//        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
//        direction = 0;
//        break;
//    case 1:
//        u = _w;
//        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
//        direction = 0;
//        break;
//    case 2:
//        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
//        v = _w;
//        direction = 1;
//        break;
//    case 3:
//        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
//        v = _w;
//        direction = 1;
//        break;
//    }
//
//    // 2. Loop over all the Newton-Raphson iterations
//    while (counter <= 100) {
//
//        // 2i. Update the iteration counter
//        counter++;
//
//        IGABasis->getUBSplineBasis1D()->clampKnot(u);
//        IGABasis->getVBSplineBasis1D()->clampKnot(v);
//        // 2ii. Find the span of the given surface parameters
//        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
//        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);
//
//        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
//        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, u,
//                uKnotSpan, v, vKnotSpan);
//
//        // 2iv. Compute the Cartesian components of the point on the surface
//        computeCartesianCoordinates(Q, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);
//
//        // 2v. Compute the distance vector between the vector to be projected and the estimated one
//        for (int i = 0; i < noSpatialDimensions; i++) {
//        	P[i]=(1-t)*_P1[i]+t*_P2[i];
//            distanceVector[i] = P[i] - Q[i];
//        }
//
//        // 2vi. Compute the 2-norm of the distance vector
//        distanceVector2norm = square2normVector(noSpatialDimensions, distanceVector);
//        distanceVector2norm = sqrt(distanceVector2norm);
//
//        if (distanceVector2norm < EPS_DISTANCE)
//            break;
//
//        // 2vii. Compute the base vectors and their derivatives
//        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
//                uKnotSpan, vKnotSpan);
//
//        // 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
//        for (int i = 0; i < noSpatialDimensions; i++) {
//            // On the base vector Gu = dR/du
//            indexUBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0);
//            uBaseVec[i] = baseVecAndDerivs[indexUBaseVec];
//
//            // On the base vector Gv = dR/dv
//            indexVBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1);
//            vBaseVec[i] = baseVecAndDerivs[indexVBaseVec];
//
//            // On the derivative of the base vector Gu w.r.t. u namely dGu/du = d^2R/du^2
//            indexdUBaseVecdu = indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0);
//            uBaseVecdu[i] = baseVecAndDerivs[indexdUBaseVecdu];
//
//            // On the derivative of the base vector Gv w.r.t. u namely dGv/dv = d^2R/dv^2
//            indexdVBaseVecdv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 1);
//            vBaseVecdv[i] = baseVecAndDerivs[indexdVBaseVecdv];
//
//            // On the mixed derivative of the base vectors namely d^2Gu/dv = d^2Gv/du = d^2R/dudv
//            indexdBaseVecMix = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0);
//            dBaseVecMix[i] = baseVecAndDerivs[indexdBaseVecMix];
//        }
//
//        // 2ix. Compute the cosine of the angle with respect to u-parametric line
//        uBaseVecXdistanceVector = dotProduct(noSpatialDimensions, uBaseVec, distanceVector);//R1
//        squareUBaseVec2norm = square2normVector(noSpatialDimensions, uBaseVec);
//        uBaseVec2norm = sqrt(squareUBaseVec2norm);
//        cosu = fabs(uBaseVecXdistanceVector) / uBaseVec2norm / distanceVector2norm;
//
//        // 2x. Compute the cosine of the angle with respect to v-parametric line
//        vBaseVecXdistanceVector = dotProduct(noSpatialDimensions, vBaseVec, distanceVector);//R2
//        squareVBaseVec2norm = square2normVector(noSpatialDimensions, vBaseVec);
//        vBaseVec2norm = sqrt(squareVBaseVec2norm);
//        cosv = fabs(vBaseVecXdistanceVector) / vBaseVec2norm / distanceVector2norm;
//        // 2x. Compute the cosine of the angle with respect to physical line
//        cosT=fabs(dotProduct(noSpatialDimensions,distanceVector,P1P2))/
//        		sqrt(square2normVector(noSpatialDimensions,distanceVector))/
//        		sqrt(square2normVector(noSpatialDimensions,P1P2));
//        // 2xi. Check the orthogonality condition and if it is fulfilled break the loop
//        cosW_prev=cosW;
//        cosT_prev=cosT;
//        if(direction==0)
//        	cosW=cosu;
//        else
//        	cosW=cosv;
//
//        //DEBUG_OUT()<<"cosW="<<cosW<<" / cosT="<<cosT<<endl;
//
//        //if(fabs(cosW-cosW_prev)/cosW_prev<EPS_ORTHOGONALITY_CONDITION && fabs(cosT-cosT_prev)/cosT_prev<EPS_ORTHOGONALITY_CONDITION && cosW>EPS_ORTHOGONALITY_CONDITION && cosT>EPS_ORTHOGONALITY_CONDITION)
//        //	return false;
//        if (cosW <= EPS_ORTHOGONALITY_CONDITION && cosT <= EPS_ORTHOGONALITY_CONDITION)
//            break;
//
//        // 2xii. Compute the entries of the Jacobian matrix
//        if(direction==0) {
//        	jacobianMatrix[0] = square2normVector(noSpatialDimensions, P1P2);
//			jacobianMatrix[1] = dotProduct(noSpatialDimensions, uBaseVec, P1P2);
//			jacobianMatrix[2] = jacobianMatrix[1];
//			jacobianMatrix[3] = dotProduct(noSpatialDimensions, uBaseVecdu, distanceVector)
//							- square2normVector(noSpatialDimensions,uBaseVec);
//        } else {
//        	jacobianMatrix[0] = square2normVector(noSpatialDimensions, P1P2);
//			jacobianMatrix[1] = dotProduct(noSpatialDimensions, vBaseVec, P1P2);
//			jacobianMatrix[2] = jacobianMatrix[1];
//			jacobianMatrix[3] = dotProduct(noSpatialDimensions, vBaseVecdv, distanceVector)
//							- square2normVector(noSpatialDimensions,vBaseVec);
//        }
//        // 2xiii. Compute the entries of the right-hand side vector
//        rightSideVct[0] = -dotProduct(noSpatialDimensions, P1P2, distanceVector);
//        if(direction==0)
//        	rightSideVct[1] = -dotProduct(noSpatialDimensions, uBaseVec, distanceVector);
//        else
//        	rightSideVct[1] = -dotProduct(noSpatialDimensions, vBaseVec, distanceVector);
//
//
//        bool flagLinearSystem = solve2x2linearSystem(rightSideVct, jacobianMatrix);
////        DEBUG_OUT()<<"Resdiual norm="<<sqrt(square2normVector(2,rightSideVct))<<endl;
//
////        if (sqrt(square2normVector(2,rightSideVct))<=1e-2)
////        	break;
////        double dw=-R/dR;
//        // 2xv. Update the surface parameters u += du and v += dv
//        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
//		if(direction==0) {
//			u += rightSideVct[1];
//	    	IGABasis->getUBSplineBasis1D()->clampKnot(u);
//		} else {
//			v += rightSideVct[1];
//	    	IGABasis->getVBSplineBasis1D()->clampKnot(v);
//		}
//		t+=rightSideVct[0];
//		if(t<0)t=0;
//		if(t>1.0)t=1.0;
//    }
//
//////     3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
//    if (counter > 100) {
//        if (cosW <= EPS_ORTHOGONALITY_CONDITION_RELAXED
//                && cosT <= 1e-5)
//            flagNewtonRaphson = true;
//        else
//            flagNewtonRaphson = false;
//        if (!flagNewtonRaphson && distanceVector2norm < EPS_DISTANCE_RELAXED)
//            flagNewtonRaphson = true;
//    } else {
//    }
//    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
//    // Clear the memory on the heap
//    delete[] basisFctsAndDerivs;
//    delete[] baseVecAndDerivs;
//
//    // Return the flag
////    double surfaceNormal[3];
////    crossProduct(surfaceNormal,uBaseVec,vBaseVec);
////    double P_tmp[3];
////    P_tmp[0]=Q[0]+surfaceNormal[0];
////    P_tmp[1]=Q[1]+surfaceNormal[1];
////    P_tmp[2]=Q[2]+surfaceNormal[2];
////    double ratioA,ratioB;
////    distanceLineLine(ratioA,ratioB,Q,P_tmp,_P1,_P2);
////    if(ratioB>1) ratioB=1;
////    if(ratioB<0) ratioB=0;
////    _ratio=ratioB;
//    _ratio=t;
//    _distance=distanceVector2norm;
//	if(direction==0)
//		_w=u;
//	else
//		_w=v;
//    return flagNewtonRaphson;
//}

//bool IGAPatchSurface::computePointProjectionOnPatchBoundaryOnGivenEdge_Ortho(double& _w, double& _ratio,
//        double& _distance, double* _P1, double* _P2, int _edge) {
//
//    /*
//     * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
//     * _P = double[3]
//     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
//     *
//     * Function layout :
//     *
//     * 1. Read input and initialize the data
//     *
//     * 2. Loop over all the Newton-Raphson iterations
//     *    2i. Update the iteration counter
//     *   2ii. Find the span of the given surface parameters
//     *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
//     *   2iv. Compute the Cartesian components of the point on the surface
//     *    2v. Compute the distance vector between the vector to be projected and the estimated one
//     *   2vi. Compute the 2-norm of the distance vector
//     *  2vii. Compute the base vectors and their derivatives
//     * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
//     *   2ix. Compute the cosine of the angle with respect to u-parametric line
//     *    2x. Compute the cosine of the angle with respect to v-parametric line
//     *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
//     *  2xii. Compute the entries of the Jacobian matrix
//     * 2xiii. Compute the entries of the right-hand side vector
//     *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
//     *   2xv. Update the surface parameters u += du and v += dv
//     *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
//     *
//     * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
//     *
//     * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
//     */
//
//    // 1. Read input and initialize the data
//    // Read input
//    assert(_P1 != NULL);
//    assert(_P2 != NULL);
//
//    // Initialize the flag to true
//    bool flagNewtonRaphson = true;
//
//
//
//    // Initialize number of spatial dimensions
//    int noSpatialDimensions = 3;
//
//    // Save the location of the point to be projected
//    double P[noSpatialDimensions],Q[noSpatialDimensions];
//    double P1P2[noSpatialDimensions];
//    for (int i = 0; i < noSpatialDimensions; i++) {
//        P[i] = _P1[i];
//        P1P2[i]=_P2[i]-_P1[i];
//    }
//    double P1P2norm=sqrt(square2normVector(noSpatialDimensions, P1P2));
//    double P1P2unit[3];
//    for (int i = 0; i < noSpatialDimensions; i++) {
//    	P1P2unit[i]=P1P2[i]/P1P2norm;
//    }
//    // Initialize the distance vector
//    double distanceVector[3];
//    double t=1;
//
//
//    // Initialize the indices of the base vectors and their derivatives
//    int indexUBaseVec = 0;
//    int indexVBaseVec = 0;
//    int indexdUBaseVecdu = 0;
//    int indexdVBaseVecdv = 0;
//    int indexdBaseVecMix = 0;
//
//    // Initialize the base vectors and their derivatives
//    double uBaseVec[3];
//    double vBaseVec[3];
//    double uBaseVecdu[3];
//    double vBaseVecdv[3];
//    double dBaseVecMix[3];
//
//    // Initialize the dot products
//    double distanceVector2norm = 0.0;
//    double uBaseVecXdistanceVector = 0.0;
//    double squareUBaseVec2norm = 0.0;
//    double uBaseVec2norm = 0.0;
//    double vBaseVecXdistanceVector = 0.0;
//    double squareVBaseVec2norm = 0.0;
//    double vBaseVec2norm = 0.0;
//
//    // Initialize Jacobian matrix
//    double jacobianMatrix[4];
//
//    // Initialize right-hand side solution vector
//    double rightSideVct[2];
//
//    // Initialize the cosines w.r.t. each parametric line
//    double cosu = 0.0;
//    double cosv = 0.0;
//
//    // Initialize the knot span indices
//    int uKnotSpan = 0;
//    int vKnotSpan = 0;
//
//    // The NURBS polynomial degrees
//    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
//    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();
//
//    // Local number of basis functions
//    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);
//
//    // Initialize the Newton-Raphson iteration counter
//    int counter = 0;
//
//    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
//    int derivDegreeBasis = 2;
//
//    // The number of the base vectors
//    int noBaseVcts = 2;
//
//    // Initialize the array of the IGA basis functions and their derivatives
//    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
//            * noLocalBasisFcts / 2];
//
//    // Number of derivatives for the base vectors
//    int derivDegreeBaseVcts = derivDegreeBasis - 1;
//
//    // Initialize the array of the base vectors and their derivatives
//    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
//            * noSpatialDimensions * noBaseVcts / 2];
//
//
//    double u;
//    double v;
//
//    // Setup initial value depending on the edge
//    int direction;
//    switch (_edge) {
//    case 0:
//        u = _w;
//        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
//        direction = 0;
//        break;
//    case 1:
//        u = _w;
//        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
//        direction = 0;
//        break;
//    case 2:
//        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
//        v = _w;
//        direction = 1;
//        break;
//    case 3:
//        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
//        v = _w;
//        direction = 1;
//        break;
//    }
//
//    double R,dR;
//    double cosW,cosT;
//    // 2. Loop over all the Newton-Raphson iterations
//    while (counter <= 100) {
//
//        // 2i. Update the iteration counter
//        counter++;
//
//        IGABasis->getUBSplineBasis1D()->clampKnot(u);
//        IGABasis->getVBSplineBasis1D()->clampKnot(v);
//        // 2ii. Find the span of the given surface parameters
//        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
//        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);
//
//        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
//        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, u,
//                uKnotSpan, v, vKnotSpan);
//
//        // 2iv. Compute the Cartesian components of the point on the surface
//        computeCartesianCoordinates(Q, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);
//
//        // 2v. Compute the distance vector between the vector to be projected and the estimated one
//        for (int i = 0; i < noSpatialDimensions; i++) {
//        	P[i]=(1-t)*_P1[i]+t*_P2[i];
//            distanceVector[i] = P[i] - Q[i];
//        }
//
//        // 2vi. Compute the 2-norm of the distance vector
//        distanceVector2norm = square2normVector(noSpatialDimensions, distanceVector);
//        distanceVector2norm = sqrt(distanceVector2norm);
//
//        if (distanceVector2norm < EPS_DISTANCE)
//            break;
//
//        // 2vii. Compute the base vectors and their derivatives
//        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
//                uKnotSpan, vKnotSpan);
//
//        // 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
//        for (int i = 0; i < noSpatialDimensions; i++) {
//            // On the base vector Gu = dR/du
//            indexUBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0);
//            uBaseVec[i] = baseVecAndDerivs[indexUBaseVec];
//
//            // On the base vector Gv = dR/dv
//            indexVBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1);
//            vBaseVec[i] = baseVecAndDerivs[indexVBaseVec];
//
//            // On the derivative of the base vector Gu w.r.t. u namely dGu/du = d^2R/du^2
//            indexdUBaseVecdu = indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0);
//            uBaseVecdu[i] = baseVecAndDerivs[indexdUBaseVecdu];
//
//            // On the derivative of the base vector Gv w.r.t. u namely dGv/dv = d^2R/dv^2
//            indexdVBaseVecdv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 1);
//            vBaseVecdv[i] = baseVecAndDerivs[indexdVBaseVecdv];
//
//            // On the mixed derivative of the base vectors namely d^2Gu/dv = d^2Gv/du = d^2R/dudv
//            indexdBaseVecMix = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0);
//            dBaseVecMix[i] = baseVecAndDerivs[indexdBaseVecMix];
//        }
//
//        // 2ix. Compute the cosine of the angle with respect to u-parametric line
//        uBaseVecXdistanceVector = dotProduct(noSpatialDimensions, uBaseVec, distanceVector);//R1
//        squareUBaseVec2norm = square2normVector(noSpatialDimensions, uBaseVec);
//        uBaseVec2norm = sqrt(squareUBaseVec2norm);
//        cosu = fabs(uBaseVecXdistanceVector) / uBaseVec2norm / distanceVector2norm;
//
//        // 2x. Compute the cosine of the angle with respect to v-parametric line
//        vBaseVecXdistanceVector = dotProduct(noSpatialDimensions, vBaseVec, distanceVector);//R2
//        squareVBaseVec2norm = square2normVector(noSpatialDimensions, vBaseVec);
//        vBaseVec2norm = sqrt(squareVBaseVec2norm);
//        cosv = fabs(vBaseVecXdistanceVector) / vBaseVec2norm / distanceVector2norm;
//
//
//        // 2xi. Check the orthogonality condition and if it is fulfilled break the loop
//
//        if(direction==0)
//        	cosW=cosu;
//        else
//        	cosW=cosv;
//
//        cosT=dotProduct(noSpatialDimensions,distanceVector,P1P2)/
//        		sqrt(square2normVector(noSpatialDimensions,distanceVector))/
//        		sqrt(square2normVector(noSpatialDimensions,P1P2));
//
//        //DEBUG_OUT()<<"cosW="<<cosW<<" / cosT="<<cosT<<endl;
//
//        if (cosu <= EPS_ORTHOGONALITY_CONDITION && cosv <= EPS_ORTHOGONALITY_CONDITION)
//            break;
//
//        // 2xii. Compute the entries of the Jacobian matrix
//        if(direction==0) {
//			// dR11=dGu/du * QP + Gu*Gu
//        	jacobianMatrix[0] = dotProduct(noSpatialDimensions, uBaseVec, P1P2)/uBaseVec2norm;
//			// dR12=Gu * P1P2
//			jacobianMatrix[1] = dotProduct(noSpatialDimensions, uBaseVecdu, distanceVector)-squareUBaseVec2norm;
//        	// dR21=dGv/du * QP + Gv*Gu
//			jacobianMatrix[2] = dotProduct(noSpatialDimensions, vBaseVec, P1P2)/vBaseVec2norm;
//			// dR22=Gv * P1P2
//			jacobianMatrix[3] = -dotProduct(noSpatialDimensions, uBaseVec, vBaseVec)+dotProduct(noSpatialDimensions,dBaseVecMix,distanceVector);
//        } else {
//			// dR11=Gu * P1P2
//        	jacobianMatrix[0] = dotProduct(noSpatialDimensions, uBaseVec, P1P2)/uBaseVec2norm;
//			// dR12=dGu/dv * QP - Gv*Gu
//			jacobianMatrix[1] = -dotProduct(noSpatialDimensions,uBaseVec,vBaseVec)+dotProduct(noSpatialDimensions,distanceVector,dBaseVecMix);
//        	// dR21=Gv * P1P2
//			jacobianMatrix[2] = dotProduct(noSpatialDimensions, vBaseVec, P1P2)/vBaseVec2norm;
//			// dR22=dGv/dv * QP - Gv*Gv
//			jacobianMatrix[3] = dotProduct(noSpatialDimensions, vBaseVecdv, distanceVector)-squareVBaseVec2norm;
//        }
//        // 2xiii. Compute the entries of the right-hand side vector
//        rightSideVct[0] = dotProduct(noSpatialDimensions, uBaseVec, distanceVector)/uBaseVec2norm;
//        rightSideVct[1] = dotProduct(noSpatialDimensions, vBaseVec, distanceVector)/vBaseVec2norm;
//
//
//        bool flagLinearSystem = solve2x2linearSystem(rightSideVct, jacobianMatrix);
//        assert(flagLinearSystem !=0);
//        DEBUG_OUT()<<"Resdiual norm="<<sqrt(square2normVector(2,rightSideVct))<<endl;
//
////        if (sqrt(square2normVector(2,rightSideVct))<=1e-2)
////        	break;
//		//double dw=(-rightSideVct[1]+jacobianMatrix[3]/jacobianMatrix[0]*rightSideVct[0])/(jacobianMatrix[3]-jacobianMatrix[2]/jacobianMatrix[0]*jacobianMatrix[1]);
//		//double dt=(-rightSideVct[0]-jacobianMatrix[1]*dw)/jacobianMatrix[0];
//        // 2xv. Update the surface parameters u += du and v += dv
//        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
//		if(direction==0) {
//			u +=rightSideVct[1];
//	    	IGABasis->getUBSplineBasis1D()->clampKnot(u);
//		} else {
//			v += rightSideVct[1];
//	    	IGABasis->getVBSplineBasis1D()->clampKnot(v);
//		}
//		t+=rightSideVct[0];
//		DEBUG_OUT()<<"u="<<u<<"/v="<<v<<"/t="<<t<<endl;
//		if(t<0)t=0;
//		if(t>1.0)t=1.0;
//    }
//
//
//
//////     3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
//    if (counter > 100) {
//        if (cosW <= EPS_ORTHOGONALITY_CONDITION_RELAXED
//                && cosT <= EPS_ORTHOGONALITY_CONDITION_RELAXED)
//            flagNewtonRaphson = true;
//        else
//            flagNewtonRaphson = false;
//        if (!flagNewtonRaphson && distanceVector2norm < EPS_DISTANCE_RELAXED)
//            flagNewtonRaphson = true;
//    } else {
//    }
//    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
//    // Clear the memory on the heap
//    delete[] basisFctsAndDerivs;
//    delete[] baseVecAndDerivs;
//
//    // Return the flag
////    double surfaceNormal[3];
////    crossProduct(surfaceNormal,uBaseVec,vBaseVec);
////    double P_tmp[3];
////    P_tmp[0]=Q[0]+surfaceNormal[0];
////    P_tmp[1]=Q[1]+surfaceNormal[1];
////    P_tmp[2]=Q[2]+surfaceNormal[2];
////    double ratioA,ratioB;
////    distanceLineLine(ratioA,ratioB,Q,P_tmp,_P1,_P2);
////    if(ratioB>1) ratioB=1;
////    if(ratioB<0) ratioB=0;
//    _ratio=t;
//    _distance=distanceVector2norm;
//	if(direction==0)
//		_w=u;
//	else
//		_w=v;
//    return flagNewtonRaphson;
//}

bool IGAPatchSurface::computePointProjectionOnPatchBoundaryOnGivenEdge(double& _w, double& _ratio,
        double& _distance, double* _P1, double* _P2, int _edge) {

    /*
     * Returns the projection of a point _P on the NURBS patch given an initial guess for the surface parameters _u, _v via references:
     * _P = double[3]
     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
     *
     * Function layout :
     *
     * 1. Read input and initialize the data
     *
     * 2. Loop over all the Newton-Raphson iterations
     *    2i. Update the iteration counter
     *   2ii. Find the span of the given surface parameters
     *  2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
     *   2iv. Compute the Cartesian components of the point on the surface
     *    2v. Compute the distance vector between the vector to be projected and the estimated one
     *   2vi. Compute the 2-norm of the distance vector
     *  2vii. Compute the base vectors and their derivatives
     * 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
     *   2ix. Compute the cosine of the angle with respect to u-parametric line
     *    2x. Compute the cosine of the angle with respect to v-parametric line
     *   2xi. Check the orthogonality condition and if it is fulfilled break the loop
     *  2xii. Compute the entries of the Jacobian matrix
     * 2xiii. Compute the entries of the right-hand side vector
     *  2xiv. Solve the linear 2x2 equation system to get the increment of the surface parameters and check if the equation system has been successfully solved
     *   2xv. Update the surface parameters u += du and v += dv
     *  2xvi. Check and modify the surface parameters if they stay out of their knot spans
     *
     * 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
     *
     * 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
     */

    // 1. Read input and initialize the data
    // Read input
    assert(_P1 != NULL);
    assert(_P2 != NULL);

    // Initialize the flag to true
    bool flagNewtonRaphson = true;



    // Initialize number of spatial dimensions
    int noSpatialDimensions = 3;

    // Save the location of the point to be projected
    double P[noSpatialDimensions],Q[noSpatialDimensions];
    double P1P2[noSpatialDimensions];
    for (int i = 0; i < noSpatialDimensions; i++) {
        P[i] = _P1[i];
        P1P2[i]=_P2[i]-_P1[i];
    }
    double P1P2norm=sqrt(square2normVector(noSpatialDimensions, P1P2));
    double P1P2unit[3];
    for (int i = 0; i < noSpatialDimensions; i++) {
    	P1P2unit[i]=P1P2[i]/P1P2norm;
    }
    // Initialize the distance vector
    double distanceVector[3];
    double t=1;


    // Initialize the indices of the base vectors and their derivatives
    int indexUBaseVec = 0;
    int indexVBaseVec = 0;
    int indexdUBaseVecdu = 0;
    int indexdVBaseVecdv = 0;
    int indexdBaseVecMix = 0;

    // Initialize the base vectors and their derivatives
    double uBaseVec[3];
    double vBaseVec[3];
    double uBaseVecdu[3];
    double vBaseVecdv[3];
    double dBaseVecMix[3];

    // Initialize the dot products
    double distanceVector2norm = 0.0;
    double uBaseVecXdistanceVector = 0.0;
    double squareUBaseVec2norm = 0.0;
    double uBaseVec2norm = 0.0;
    double vBaseVecXdistanceVector = 0.0;
    double squareVBaseVec2norm = 0.0;
    double vBaseVec2norm = 0.0;

    // Initialize Jacobian matrix
    double jacobianMatrix[4];

    // Initialize right-hand side solution vector
    double rightSideVct[2];

    // Initialize the cosines w.r.t. each parametric line
    double cosu = 0.0;
    double cosv = 0.0;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the IGA basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
            * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
            * noSpatialDimensions * noBaseVcts / 2];


    double u;
    double v;

    // Setup initial value depending on the edge
    int direction;
    switch (_edge) {
    case 0:
        u = _w;
        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
        direction = 0;
        break;
    case 1:
        u = _w;
        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
        direction = 0;
        break;
    case 2:
        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
        v = _w;
        direction = 1;
        break;
    case 3:
        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
        v = _w;
        direction = 1;
        break;
    }

    double R,dR;
    double cosW,cosT;
    double surfaceNormal[3];
    double surfaceNormalNorm;
    double planeNormal[3];
    double planeNormalNorm;
    // 2. Loop over all the Newton-Raphson iterations
    while (counter <= 100) {

        // 2i. Update the iteration counter
        counter++;

        IGABasis->getUBSplineBasis1D()->clampKnot(u);
        IGABasis->getVBSplineBasis1D()->clampKnot(v);
        // 2ii. Find the span of the given surface parameters
        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);

        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, u,
                uKnotSpan, v, vKnotSpan);

        // 2iv. Compute the Cartesian components of the point on the surface
        computeCartesianCoordinates(Q, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);

        // 2v. Compute the distance vector between the vector to be projected and the estimated one
        for (int i = 0; i < noSpatialDimensions; i++) {
        	P[i]=(1-t)*_P1[i]+t*_P2[i];
            distanceVector[i] = _P1[i] - Q[i];
        }

        // 2vi. Compute the 2-norm of the distance vector
        distanceVector2norm = square2normVector(noSpatialDimensions, distanceVector);
        distanceVector2norm = sqrt(distanceVector2norm);

        if (distanceVector2norm < EPS_DISTANCE)
            break;

        // 2vii. Compute the base vectors and their derivatives
        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
                uKnotSpan, vKnotSpan);

        // 2viii. Filter out the base vectors and their derivatives from the continuous 1D pointer array
        for (int i = 0; i < noSpatialDimensions; i++) {
            // On the base vector Gu = dR/du
            indexUBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 0);
            uBaseVec[i] = baseVecAndDerivs[indexUBaseVec];

            // On the base vector Gv = dR/dv
            indexVBaseVec = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i, 1);
            vBaseVec[i] = baseVecAndDerivs[indexVBaseVec];

            // On the derivative of the base vector Gu w.r.t. u namely dGu/du = d^2R/du^2
            indexdUBaseVecdu = indexDerivativeBaseVector(derivDegreeBaseVcts, 1, 0, i, 0);
            uBaseVecdu[i] = baseVecAndDerivs[indexdUBaseVecdu];

            // On the derivative of the base vector Gv w.r.t. u namely dGv/dv = d^2R/dv^2
            indexdVBaseVecdv = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 1);
            vBaseVecdv[i] = baseVecAndDerivs[indexdVBaseVecdv];

            // On the mixed derivative of the base vectors namely d^2Gu/dv = d^2Gv/du = d^2R/dudv
            indexdBaseVecMix = indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 1, i, 0);
            dBaseVecMix[i] = baseVecAndDerivs[indexdBaseVecMix];
        }

        // 2ix. Compute the cosine of the angle with respect to u-parametric line
        uBaseVecXdistanceVector = dotProduct(noSpatialDimensions, uBaseVec, distanceVector);//R1
        squareUBaseVec2norm = square2normVector(noSpatialDimensions, uBaseVec);
        uBaseVec2norm = sqrt(squareUBaseVec2norm);
        cosu = fabs(uBaseVecXdistanceVector) / uBaseVec2norm / distanceVector2norm;

        // 2x. Compute the cosine of the angle with respect to v-parametric line
        vBaseVecXdistanceVector = dotProduct(noSpatialDimensions, vBaseVec, distanceVector);//R2
        squareVBaseVec2norm = square2normVector(noSpatialDimensions, vBaseVec);
        vBaseVec2norm = sqrt(squareVBaseVec2norm);
        cosv = fabs(vBaseVecXdistanceVector) / vBaseVec2norm / distanceVector2norm;


        // 2xi. Check the orthogonality condition and if it is fulfilled break the loop
        crossProduct(surfaceNormal, uBaseVec, vBaseVec);
        surfaceNormalNorm=sqrt(square2normVector(noSpatialDimensions,surfaceNormal));
        crossProduct(planeNormal, P1P2, surfaceNormal);
        // 2xiii. Compute the entries of the right-hand side vector
        rightSideVct[0] = dotProduct(noSpatialDimensions, distanceVector, planeNormal)/(P1P2norm*uBaseVec2norm*vBaseVec2norm);
        // 2xii. Compute the entries of the Jacobian matrix
        if(direction==0) {
        	jacobianMatrix[0] = dotProduct(noSpatialDimensions, uBaseVec, surfaceNormal);
        } else {
        	jacobianMatrix[0] = dotProduct(noSpatialDimensions, vBaseVec, surfaceNormal);
        }

        DEBUG_OUT()<<"Resdiual norm="<<rightSideVct[0]<<endl;

		double dw=-rightSideVct[0]/jacobianMatrix[0];
        if (dw<EPS_DISTANCE)
            break;
        // 2xv. Update the surface parameters u += du and v += dv
        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
		if(direction==0) {
			u +=dw;
	    	IGABasis->getUBSplineBasis1D()->clampKnot(u);
		} else {
			v += dw;
	    	IGABasis->getVBSplineBasis1D()->clampKnot(v);
		}
		DEBUG_OUT()<<"u="<<u<<"/v="<<v<<endl;
    }



////     3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
    if (counter > 100) {
    	flagNewtonRaphson=false;
    } else {
    	flagNewtonRaphson=true;
    }
    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
    // Clear the memory on the heap
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return the flag
    if(direction==0)
		for (int i = 0; i < noSpatialDimensions; i++) {
			vBaseVec[i]=vBaseVec[i]/vBaseVec2norm;
		}
    else
		for (int i = 0; i < noSpatialDimensions; i++) {
			uBaseVec[i]=uBaseVec[i]/uBaseVec2norm;
		}
    _ratio=distanceLinePlane(_P1,P1P2unit,Q,uBaseVec)/P1P2norm;
    _distance=distanceVector2norm;
	if(direction==0)
		_w=u;
	else
		_w=v;
    return flagNewtonRaphson;
}

bool IGAPatchSurface::computePointProjectionOnPatchBoundary(double& _u, double& _v, double& _ratio,
        double& _distance, double* _P1, double* _P2) {
    double u1 = _u;
    double v1 = _v;
    double t;
    double distance = numeric_limits<double>::max();
    double div = 0.0;
    bool isConverged = false;
    _distance = numeric_limits<double>::max();


		double u=_u;
		double v=_v;
		// Compute point projection from the line to the NURBS patch boundary
		isConverged = computePointProjectionOnPatchBoundaryOnGivenEdge_Brute(u,v, div,distance, _P1, _P2);
		DEBUG_OUT()<<"\tConverged ? "<<isConverged<<" and distance is "<<distance<<" and div is "<<div<<endl;

		if(isConverged){
			DEBUG_OUT()<<"\tu["<<u<<"], v["<<v<<"]"<<endl;

			// Fix possible numerical error
			if(fabs(div-1.0)<EPS_DISTANCE_RELAXED && div-1.0>0) div=1.0;
			if(fabs(div)	<EPS_DISTANCE_RELAXED && div<0) div=0.0;
			_ratio=div;
			_distance=distance;
			_u=u;
			_v=v;
			return true;
		}
		return false;

    // Loop over all the edges of the NURBS patch (for tensor product surfaces there are 4 edges)
    for (int edge = 0; edge < 4; edge++) {
    	double u=_u;
    	double v=_v;
    					if (edge == 0 || edge == 1)
    						t = u1;
    					else
    						t = v1;
			// Compute point projection from the line to the NURBS patch boundary
//			isConverged = computePointProjectionOnPatchBoundaryOnGivenEdge(t, div, distance, _P1, _P2, edge);
			isConverged = computePointProjectionOnPatchBoundaryOnGivenEdge_Brute(u,v, div,distance, _P1, _P2);

			// Fix possible numerical error
			if(fabs(div-1.0)<EPS_DISTANCE_RELAXED && div-1.0>0) div=1.0;
			if(fabs(div)	<EPS_DISTANCE_RELAXED && div<0) div=0.0;

			// Debug information about computePointOnPatchBoundaryOnGivenEdge result
			DEBUG_OUT()<<"\tEdge["<<edge<<"] converged ? "<<isConverged<<" and distance is "<<distance<<" and div is "<<div<<" and t is "<<t<<endl;

			if (isConverged && distance<_distance) {
				switch (edge) {
				case 0:
					IGABasis->getUBSplineBasis1D()->clampKnot(t);
					u = t;
					v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
					break;
				case 1:
					IGABasis->getUBSplineBasis1D()->clampKnot(t);
					u = t;
					v = IGABasis->getVBSplineBasis1D()->getLastKnot();
					break;
				case 2:
					IGABasis->getVBSplineBasis1D()->clampKnot(t);
					u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
					v = t;
					break;
				case 3:
					IGABasis->getVBSplineBasis1D()->clampKnot(t);
					u = IGABasis->getUBSplineBasis1D()->getLastKnot();
					v = t;
					break;
				}
				DEBUG_OUT()<<"\tu["<<u<<"], v["<<v<<"]"<<endl;
				_u=u;
				_v=v;
				_distance = distance;
				_ratio = div;
			}
    	}
    if (_distance == numeric_limits<double>::max())
        return false;
    else if (_ratio >= 0.0 && _ratio <= 1.0) {
        return true;
    } else {
        return false;
    }
}
//
//bool IGAPatchSurface::computePointProjectionOnPatchBoundary(double& _u, double& _v, double& _ratio,
//        double& _distance, double* _P1, double* _P2) {
//    double u1 = _u;
//    double v1 = _v;
//    double t;
//    double distance = numeric_limits<double>::max();
//    double div = 0.0;
//    bool isConverged = false;
//    _distance = numeric_limits<double>::max();
//
//    // Loop over all the edges of the NURBS patch (for tensor product surfaces there are 4 edges)
//    for (int edge = 0; edge < 4; edge++) {
//    	double u[3];
//    	double v[3];
//        bool hasConverged=false;
//
//    	// Do the test for every edge for each extremity of the boundary patch
//    	for(int point=0;point<3;point++) {
//			// Find the fixed and the running parameter on the patch boundary
//			if (edge == 0 || edge == 1)
//				if(point==0)
//					t=IGABasis->getUBSplineBasis1D()->getFirstKnot();
//				else
//					t=IGABasis->getUBSplineBasis1D()->getLastKnot();
//			else
//				if(point==0)
//					t=IGABasis->getVBSplineBasis1D()->getFirstKnot();
//				else
//					t=IGABasis->getVBSplineBasis1D()->getLastKnot();
//
//	        // If extremities of edge have not converged then try with initial guess
//			if(point==2 && !hasConverged) {
//				if (edge == 0 || edge == 1)
//					t = u1;
//				else
//					t = v1;
//			}
//			// Compute point projection from the line to the NURBS patch boundary
//			isConverged = computePointProjectionOnPatchBoundaryOnGivenEdge(t, div, distance, _P1, _P2, edge);
//
//			// Fix possible numerical error
//			if(fabs(div-1.0)<EPS_DISTANCE && div-1.0>0) div=1.0;
//			if(fabs(div)	<EPS_DISTANCE && div<0) div=0.0;
//
//			// Debug information about computePointOnPatchBoundaryOnGivenEdge result
//			DEBUG_OUT()<<"\tPoint["<<point<<"], Edge["<<edge<<"] converged ? "<<isConverged<<" and distance is "<<distance<<" and div is "<<div<<" and t is "<<t<<endl;
//
//			if (isConverged  || distance<_distance) {
//				switch (edge) {
//				case 0:
//					IGABasis->getUBSplineBasis1D()->clampKnot(t);
//					u[point] = t;
//					v[point] = IGABasis->getVBSplineBasis1D()->getFirstKnot();
//					break;
//				case 1:
//					IGABasis->getUBSplineBasis1D()->clampKnot(t);
//					u[point] = t;
//					v[point] = IGABasis->getVBSplineBasis1D()->getLastKnot();
//					break;
//				case 2:
//					IGABasis->getVBSplineBasis1D()->clampKnot(t);
//					u[point] = IGABasis->getUBSplineBasis1D()->getFirstKnot();
//					v[point] = t;
//					break;
//				case 3:
//					IGABasis->getVBSplineBasis1D()->clampKnot(t);
//					u[point] = IGABasis->getUBSplineBasis1D()->getLastKnot();
//					v[point] = t;
//					break;
//				}
//				// If the point is the same as the entry point
//				bool validPoint1=(u[point]!=u1 || v[point]!=v1);//Different from entry point then true else false
//				bool validPoint2=(point==1)?(u[0]==u1 && v[0]==v1 && u[1]==u1 && v[1]==v1):false;//Both are the same then true else false
//				// If it is not a point to take into account, continue
//				if(!(validPoint1 || validPoint2) && point<2) continue;
//				//Otherwise store it under following conditions
//				if (distance < _distance && div >= 0.0 && div <= 1.0) {
//					DEBUG_OUT()<<"\tu["<<u[point]<<"], v["<<v[point]<<"]"<<endl;
//					hasConverged=true;
//					_u=u[point];
//					_v=v[point];
//					_distance = distance;
//					_ratio = div;
//				}
//			}
//    	}
//    }
//    if (_distance == numeric_limits<double>::max())
//        return false;
//    else if (_ratio >= 0.0 && _ratio <= 1.0) {
//        return true;
//    } else {
//        return false;
//    }
//}

bool IGAPatchSurface::computePointMinimumDistanceToPatchBoundaryOnGivenEdge(double& _t,
        double& _distance, double* _P1, int _edge) {

    assert(_P1 != NULL);

    bool flagNewtonRaphson = true;

    // Initialize number of spatial dimensions
    int noSpatialDimensions = 3;

    // Initialize the distance vector
    double distanceVector[3];

    // Initialize the base vectors and their derivatives
    double baseVec[3];
    double dBaseVec[3];

    // Initialize Jacobian matrix
    double df;

    // Initialize right-hand side solution vector
    double f;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // The lengths of the knot vectors to the NURBS patch
    int lengthUKnotVct = IGABasis->getUBSplineBasis1D()->getNoKnots();
    int lengthVKnotVct = IGABasis->getVBSplineBasis1D()->getNoKnots();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the IGA basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
            * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
            * noSpatialDimensions * noBaseVcts / 2];

    int direction;
    double u;
    double v;
    switch (_edge) {
    case 0:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
        direction = 0;
        break;
    case 1:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
        direction = 0;
        break;
    case 2:
        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
        v = _t;
        direction = 1;
        break;
    case 3:
        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
        v = _t;
        direction = 1;
        break;
    }

    double P0[3];

	IGABasis->getUBSplineBasis1D()->clampKnot(u);
	IGABasis->getVBSplineBasis1D()->clampKnot(v);
    // 2. Loop over all the Newton-Raphson iterations
    while (counter <= MAX_NUM_ITERATIONS) {

        // 2i. Update the iteration counter
        counter++;

        // 2ii. Find the span of the given surface parameters
        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);

        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, u,
                uKnotSpan, v, vKnotSpan);

        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
                uKnotSpan, vKnotSpan);

        // 2iv. Compute the Cartesian components of the point on the surface
        computeCartesianCoordinates(P0, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);

        // 2v. Compute the distance vector between the vector to be projected and the estimated one
        for (int i = 0; i < noSpatialDimensions; i++) {
            distanceVector[i] = P0[i] - _P1[i];
            baseVec[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i,
                    direction)];
            dBaseVec[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts,
                    1 - direction, direction, i, direction)];
        }
        f = dotProduct(noSpatialDimensions, baseVec, distanceVector);
        df = dotProduct(noSpatialDimensions, dBaseVec, distanceVector)
                + dotProduct(noSpatialDimensions, baseVec, baseVec);


        if(df!=df)
        	continue;
        if(df<1e-9)
        	break;
        if (fabs(f / df) < 1e-13)
            break;

        if (direction == 0)
            u -= f / df;
        else
            v -= f / df;

        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
		IGABasis->getUBSplineBasis1D()->clampKnot(u);
		IGABasis->getVBSplineBasis1D()->clampKnot(v);
    }

    // 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
    if (counter > MAX_NUM_ITERATIONS)
        flagNewtonRaphson = false;

    if (direction == 0)
        _t = u;
    else
        _t = v;

    _distance = sqrt(square2normVector(noSpatialDimensions, distanceVector));
    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
    // Clear the memory on the heap
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return the flag
    return flagNewtonRaphson;

}

void IGAPatchSurface::computePointMinimumDistanceToPatchBoundary(double& _u, double& _v,
        double& _distance, double* _P1, int* _edge) {
	_edge[0]=-1;
	_edge[1]=-1;
    double u1 = _u;
    double v1 = _v;
    double t;
    double distance = numeric_limits<double>::max();
    bool isConverged = false;
    _distance = distance;

    // Loop over all the edges of the NURBS patch (for tensor product surfaces there are 4 edges)
    for (int edge = 0; edge < 4; edge++) {
    	// Do the test for every edge for each extremity of the boundary patch
    	for(int point=0;point<2;point++) {
			// Find the running and the fixed parameter on the NURBS boundary
			if (edge == 0 || edge == 1)
				if(point==0)
					t=IGABasis->getUBSplineBasis1D()->getFirstKnot();
				else
					t=IGABasis->getUBSplineBasis1D()->getLastKnot();
			else
				if(point==0)
					t=IGABasis->getVBSplineBasis1D()->getFirstKnot();
				else
					t=IGABasis->getVBSplineBasis1D()->getLastKnot();

			// Compute point projection from the line to the NURBS patch boundary
			isConverged = computePointMinimumDistanceToPatchBoundaryOnGivenEdge(t, distance, _P1,
					edge);
			if (distance > _distance + EPS_DISTANCE) {
				continue;
			} else if (distance < _distance) {
				_edge[1]=_edge[0];
				_edge[0]=edge;
				_distance = distance;
				switch (edge) {
				case 0:
					_u = t;
					_v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
					break;
				case 1:
					_u = t;
					_v = IGABasis->getVBSplineBasis1D()->getLastKnot();
					break;
				case 2:
					_u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
					_v = t;
					break;
				case 3:
					_u = IGABasis->getUBSplineBasis1D()->getLastKnot();
					_v = t;
					break;
				}
			} else if(fabs(distance-_distance)<EPS_DISTANCE) {
				_edge[1]=edge;
			}
    	}
    }
}


bool IGAPatchSurface::computeLineMinimumDistanceToPatchBoundaryOnGivenEdge(double& _t,
        double& _distance, double* _P1, double* _P2, int _edge) {

    assert(_P1 != NULL);
    assert(_P2 != NULL);

    bool flagNewtonRaphson = true;

    // Initialize number of spatial dimensions
    int noSpatialDimensions = 3;

    // Initialize the distance vector
    double distanceVector01[3];
    double distanceVector02[3];
    double distanceVector12[3];
    double unitVector12[3];
    for (int i = 0; i < noSpatialDimensions; i++)
    	distanceVector12[i] = _P1[i] - _P2[i];
    double distanceNorm12=sqrt(square2normVector(noSpatialDimensions,distanceVector12));
    for (int i = 0; i < noSpatialDimensions; i++)
    	unitVector12[i] = distanceVector12[i]/distanceNorm12;

    double product1[3];
    double product2[3];
    double product3[3];

    // Initialize the base vectors and their derivatives
    double baseVec[3];
    double dBaseVec[3];

    // Initialize Jacobian matrix
    double df;

    // Initialize right-hand side solution vector
    double f;

    // Initialize the knot span indices
    int uKnotSpan = 0;
    int vKnotSpan = 0;

    // The NURBS polynomial degrees
    int pDegree = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int qDegree = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // The lengths of the knot vectors to the NURBS patch
    int lengthUKnotVct = IGABasis->getUBSplineBasis1D()->getNoKnots();
    int lengthVKnotVct = IGABasis->getVBSplineBasis1D()->getNoKnots();

    // Local number of basis functions
    int noLocalBasisFcts = (pDegree + 1) * (qDegree + 1);

    // Initialize the Newton-Raphson iteration counter
    int counter = 0;

    // Number of derivatives needed for the basis functions (cause also 1st derivatives of the base vectors are needed for the Newton-Raphson iterations)
    int derivDegreeBasis = 2;

    // The number of the base vectors
    int noBaseVcts = 2;

    // Initialize the array of the IGA basis functions and their derivatives
    double* basisFctsAndDerivs = new double[(derivDegreeBasis + 1) * (derivDegreeBasis + 2)
            * noLocalBasisFcts / 2];

    // Number of derivatives for the base vectors
    int derivDegreeBaseVcts = derivDegreeBasis - 1;

    // Initialize the array of the base vectors and their derivatives
    double* baseVecAndDerivs = new double[(derivDegreeBaseVcts + 1) * (derivDegreeBaseVcts + 2)
            * noSpatialDimensions * noBaseVcts / 2];

    int direction;
    double u;
    double v;
    switch (_edge) {
    case 0:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
        direction = 0;
        break;
    case 1:
        u = _t;
        v = IGABasis->getVBSplineBasis1D()->getLastKnot();
        direction = 0;
        break;
    case 2:
        u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
        v = _t;
        direction = 1;
        break;
    case 3:
        u = IGABasis->getUBSplineBasis1D()->getLastKnot();
        v = _t;
        direction = 1;
        break;
    }
	IGABasis->getUBSplineBasis1D()->clampKnot(u);
	IGABasis->getVBSplineBasis1D()->clampKnot(v);

    double P0[3];

    // 2. Loop over all the Newton-Raphson iterations
    while (counter <= MAX_NUM_ITERATIONS) {

        // 2i. Update the iteration counter
        counter++;

        // 2ii. Find the span of the given surface parameters
        uKnotSpan = IGABasis->getUBSplineBasis1D()->findKnotSpan(u);
        vKnotSpan = IGABasis->getVBSplineBasis1D()->findKnotSpan(v);

        // 2iii. Compute the IGA basis functions and their derivatives at current (_u,_v) pair of surface parameters
        IGABasis->computeLocalBasisFunctionsAndDerivatives(basisFctsAndDerivs, derivDegreeBasis, u,
                uKnotSpan, v, vKnotSpan);

        computeBaseVectorsAndDerivatives(baseVecAndDerivs, basisFctsAndDerivs, derivDegreeBaseVcts,
                uKnotSpan, vKnotSpan);

        // 2iv. Compute the Cartesian components of the point on the surface
        computeCartesianCoordinates(P0, basisFctsAndDerivs, derivDegreeBasis, uKnotSpan, vKnotSpan);

        // 2v. Compute the distance vector between the vector to be projected and the estimated one
        for (int i = 0; i < noSpatialDimensions; i++) {
            distanceVector01[i] = P0[i] - _P1[i];
            distanceVector02[i] = P0[i] - _P2[i];
            baseVec[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts, 0, 0, i,
                    direction)];
            dBaseVec[i] = baseVecAndDerivs[indexDerivativeBaseVector(derivDegreeBaseVcts,
                    1 - direction, direction, i, direction)];
        }

        crossProduct(product1, baseVec, distanceVector12);
        crossProduct(product2, dBaseVec, distanceVector12);
        crossProduct(product3, distanceVector01, distanceVector02);
        f = dotProduct(noSpatialDimensions, product1, product3);
        df = dotProduct(noSpatialDimensions, product2, product3)
                + dotProduct(noSpatialDimensions, product1, product1);

        if(df==0)
        	continue;
        if (fabs(f / df) < 1e-13)
            break;

        if (direction == 0)
            u -= f / df;
        else
            v -= f / df;


        // 2xvi. Check and modify the surface parameters if they stay out of their knot spans
		IGABasis->getUBSplineBasis1D()->clampKnot(u);
		IGABasis->getVBSplineBasis1D()->clampKnot(v);
    }

    // 3. Check whether maximum number of iterations has been reached and if yes return 0 to the flag (non-converged iterations)
    if (counter > MAX_NUM_ITERATIONS)
        flagNewtonRaphson = false;

    if (direction == 0)
        _t = u;
    else
        _t = v;

    _distance = sqrt(square2normVector(noSpatialDimensions, product3))
            / sqrt(square2normVector(noSpatialDimensions, distanceVector12));

    // 4. Function appendix (Clear the memory from the dynamically allocated variables and return the flag on convergence)
    // Clear the memory on the heap
    delete[] basisFctsAndDerivs;
    delete[] baseVecAndDerivs;

    // Return the flag
    return flagNewtonRaphson;

}

void IGAPatchSurface::computeLineMinimumDistanceToPatchBoundary(double& _u, double& _v,
        double& _distance, double* _P1, double* _P2) {
    double u1 = _u;
    double v1 = _v;
    double t;
    double distance = numeric_limits<double>::max();
    bool isConverged = false;
    _distance = distance;

    // Loop over all the edges of the NURBS patch (for tensor product surfaces there are 4 edges)
    for (int edge = 0; edge < 4; edge++) {

        // Find the running and the fixed parameter on the NURBS boundary
        if (edge == 0 || edge == 1)
            t = u1;
        else
            t = v1;

        // Compute point projection from the line to the NURBS patch boundary
        isConverged = computeLineMinimumDistanceToPatchBoundaryOnGivenEdge(t, distance, _P1, _P2,
                edge);

        if (distance < _distance) {

            _distance = distance;
            switch (edge) {
            case 0:
                _u = t;
                _v = IGABasis->getVBSplineBasis1D()->getFirstKnot();
                break;
            case 1:
                _u = t;
                _v = IGABasis->getVBSplineBasis1D()->getLastKnot();
                break;
            case 2:
                _u = IGABasis->getUBSplineBasis1D()->getFirstKnot();
                _v = t;
                break;
            case 3:
                _u =IGABasis->getUBSplineBasis1D()->getLastKnot();
                _v = t;
                break;
            }
        }

    }
}

void IGAPatchSurface::findInitialGuess4PointProjection(double& _u, double& _v, double* _P,
        int _uDiv, int _vDiv) {

    assert(_P != NULL);

    const int noSpatialDimensions = 3;

    // Number of Basis function for U and V
    int nU = IGABasis->getUBSplineBasis1D()->computeNoBasisFunctions();
    int nV = IGABasis->getVBSplineBasis1D()->computeNoBasisFunctions();

    // The NURBS polynomial degrees for U and V
    int pU = IGABasis->getUBSplineBasis1D()->getPolynomialDegree();
    int pV = IGABasis->getVBSplineBasis1D()->getPolynomialDegree();

    // Knot vector for U and V
    double *knotU = IGABasis->getUBSplineBasis1D()->getKnotVector();
    double *knotV = IGABasis->getVBSplineBasis1D()->getKnotVector();

    double coords[3];
    double minDis = numeric_limits<double>::max();
    double Dis;

    double u0 = knotU[pU];
    double v0 = knotV[pV];
    double uEnd = knotU[nU];
    double vEnd = knotV[nV];
    double du = (uEnd - u0) / _uDiv;
    double dv = (vEnd - v0) / _vDiv;
    double uv[2];

    for (int i = 1; i < _uDiv - 1; i++)
        for (int j = 1; j < _vDiv - 1; j++) {
            uv[0] = u0 + du * i;
            uv[1] = v0 + dv * j;
            computeCartesianCoordinates(coords, uv);

            for (int k = 0; k < noSpatialDimensions; k++)
                coords[k] -= _P[k];

            Dis = square2normVector(noSpatialDimensions, coords);
            Dis = sqrt(Dis);

            if (Dis < minDis) {
                minDis = Dis;
                _u = uv[0];
                _v = uv[1];
            }

        }
}

void IGAPatchSurface::computeCartesianCoordinatesAndNormalVector(double* _coords, double* _normal,
        double _u, double _v) {

    // The knot spans
    int spanU = IGABasis->getUBSplineBasis1D()->findKnotSpan(_u);
    int spanV = IGABasis->getVBSplineBasis1D()->findKnotSpan(_v);

    // Compute the Cartesian coordinates of (_u,_v)
    computeCartesianCoordinates(_coords, _u, spanU, _v, spanV);

    // Compute the base vectors
    double baseVec[6];
    computeBaseVectors(baseVec, _u, spanU, _v, spanV);

    // Compute the cross product of the surface base vectors to get the surface normal
    _normal[0] = baseVec[1] * baseVec[5] - baseVec[2] * baseVec[4];
    _normal[1] = baseVec[2] * baseVec[3] - baseVec[0] * baseVec[5];
    _normal[2] = baseVec[0] * baseVec[4] - baseVec[1] * baseVec[3];
}

Message &operator<<(Message &message, const IGAPatchSurface &mesh) {
//  message << "\t" << "IGA Patch name: " << mesh.name << endl;

    message << "\t\tpDegree:  " << mesh.getIGABasis()->getUBSplineBasis1D()->getPolynomialDegree()
            << endl;
    message << "\t\tqDegree:  " << mesh.getIGABasis()->getVBSplineBasis1D()->getPolynomialDegree()
            << endl;

    message << "\t\tKnots Vector U: \t";
    for (int i = 0; i < mesh.getIGABasis()->getUBSplineBasis1D()->getNoKnots(); i++)
        message << mesh.getIGABasis()->getUBSplineBasis1D()->getKnotVector()[i] << "  ";
    message << endl;

    message << "\t\tKnots Vector V: \t";
    for (int i = 0; i < mesh.getIGABasis()->getVBSplineBasis1D()->getNoKnots(); i++)
        message << mesh.getIGABasis()->getVBSplineBasis1D()->getKnotVector()[i] << "  ";
    message << endl;

    message << "\t\t" << "number of control points U: " << mesh.getUNoControlPoints() << endl;
    message << "\t\t" << "number of control points V: " << mesh.getVNoControlPoints() << endl;

    message << "\t\tControl Points Net: " << endl;
    int count = 0;
    for (int j = 0; j < mesh.getVNoControlPoints(); j++) {
        message << "\t\t";
        for (int i = 0; i < mesh.getUNoControlPoints(); i++) {
            message << mesh.getControlPointNet()[count]->getX() << ", "
                    << mesh.getControlPointNet()[count]->getY() << ", "
                    << mesh.getControlPointNet()[count]->getZ() << "\t";
            count++;
        }
        message << endl;
    }
//    //Printf patch to be written in C++ unit style
//    count=0;
//    for (int j = 0; j < mesh.getVNoControlPoints(); j++) {
//        message << "\t\t";
//        for (int i = 0; i < mesh.getUNoControlPoints(); i++) {
//            message << "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+0]="<<mesh.getControlPointNet()[count]->getX() << ";" <<endl
//            		<< "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+1]="<<mesh.getControlPointNet()[count]->getY() << ";" <<endl
//            		<< "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+2]="<<mesh.getControlPointNet()[count]->getZ() << ";" <<endl
//            		<< "controlPoints[4*"<<(i+j* mesh.getUNoControlPoints())<<"+3]="<<mesh.getControlPointNet()[count]->getW() << ";" <<endl;
//            count++;
//        }
//        message << endl;
//    }
//    for(int j = 0; j < mesh.getVNoControlPoints()*mesh.getUNoControlPoints(); j++) {
//        message << "dofIndexNet["<<j<<"]="<<mesh.getControlPointNet()[j]->getDofIndex() << ";" <<endl;
//    }
    if(mesh.isTrimmed()) {
		message << mesh.getTrimming();
	}
    message << "\t" << "---------------------------------End Patch" << endl;
    return message;
}

}/* namespace EMPIRE */

