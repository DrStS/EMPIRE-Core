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
#include <assert.h>

// Inclusion of user defined libraries
#include "IGAPatch1D.h"

using namespace std;

namespace EMPIRE {

IGAPatch1D::IGAPatch1D(int _ID, int _IDBasis, int _pDegree, int _noKnots, double* _knotVector,
        int _noControlPoints, IGAControlPoint* _controlPointNet) :
        AbstractIGAPatch(_ID), NoControlPoints(_noControlPoints) {

    if (_noControlPoints != _noKnots - _pDegree - 1) {
        cout << endl;
        cout << endl;
        cout << "Error in IGAPatch1D::IGAPatch1D" << endl;
        cout << "Number of Control Points, number of knots and polynomial degree do not match!"
                << endl;
        cout << endl;
        cout << endl;
        exit(-1);
    }

    // Figure out whether the patch has a B-Spline or a NURBS underlying basis
    int isNurbs = 0;
    for (int i = 0; i < NoControlPoints; i++) {
        if (_controlPointNet[i].getW() != 1) {
            isNurbs = 1;
            break;
        }
    }

    // Create the NURBS or the B-Spline underlying basis
    if (!isNurbs) {
        IGABasis = new BSplineBasis1D(_IDBasis, _pDegree, _noKnots, _knotVector);
    } else {
        double* controlPointWeights = new double[NoControlPoints];
        for (int i = 0; i < NoControlPoints; i++)
            controlPointWeights[i] = _controlPointNet[i].getW();
        IGABasis = new NurbsBasis1D(_IDBasis, _pDegree, _noKnots, _knotVector, _noControlPoints,
                controlPointWeights);
    }

    // On the Control Point net
    assert(_controlPointNet!=NULL);
    ControlPointNet = _controlPointNet;
}

IGAPatch1D::~IGAPatch1D() {
    delete IGABasis;
    delete[] ControlPointNet;
}

IGAPatch1D::IGAPatch1D(const IGAPatch1D& _igaPatch1D) :
        AbstractIGAPatch(_igaPatch1D) {

    // Figure out whether the given patch has a B-Spline or a NURBS underlying basis
    int isNurbs = 0;
    for (int i = 0; i < NoControlPoints; i++) {
        if (_igaPatch1D.ControlPointNet[i].getW() != 1) {
            isNurbs = 1;
            break;
        }
    }

    // Copy the IGA basis member
    IGABasis = new BSplineBasis1D(*_igaPatch1D.IGABasis);

    // Copy the number of the Control Points
    NoControlPoints = _igaPatch1D.NoControlPoints;

    // Copy the Control Point net
    ControlPointNet = new IGAControlPoint[NoControlPoints];
    for (int i = 0; i < NoControlPoints; i++)
        ControlPointNet[i] = _igaPatch1D.ControlPointNet[i];
}

void IGAPatch1D::computeCartesianCoordinates(double* _cartesianCoordinates, double _uPrm,
        int _knotSpanIndex) {
    /*
     *  Returns the Cartesian coordinates of a point on the 1D IGA patch whose parameter is _uPrm. The coordinates of the point are assumed on the
     *  3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0;

    // Compute the local basis functions in the vicinity of the point
    int noLocalBasisFunctions = IGABasis->getPolynomialDegree() + 1;
    double* localBasisFunctions = new double[noLocalBasisFunctions];
    IGABasis->computeLocalBasisFunctions(localBasisFunctions, _uPrm, _knotSpanIndex);

    // Loop over all the non-zero contributions
    for (int i = 0; i <= IGABasis->getPolynomialDegree(); i++) {
        // Compute iteratively the x-coordinate of the point
        _cartesianCoordinates[0] += localBasisFunctions[i]
                * ControlPointNet[_knotSpanIndex - IGABasis->getPolynomialDegree() + i].getX();
        // Compute iteratively the y-coordinate of the point
        _cartesianCoordinates[1] += localBasisFunctions[i]
                * ControlPointNet[_knotSpanIndex - IGABasis->getPolynomialDegree() + i].getY();
        // Compute iteratively the z-coordinate of the point
        _cartesianCoordinates[2] += localBasisFunctions[i]
                * ControlPointNet[_knotSpanIndex - IGABasis->getPolynomialDegree() + i].getZ();
    }

    // Free the memory from the heap
    delete[] localBasisFunctions;
}

void IGAPatch1D::computeCartesianCoordinates(double* _cartesianCoordinates,
        double* _localBasisFunctions, double _uPrm, int _knotSpanIndex) {
    /*
     *  Returns the Cartesian coordinates of a point on the 1D IGA patch whose parameter is _uPrm and the basis functions
     *  which are assumed to have been computed outside of the function call. The coordinates of the point are assumed on the
     *  3D space that is _cartesianCoordinates = [X Y Z]
     */

    // Initialize the coordinates of the point
    for (int i = 0; i < 3; i++)
        _cartesianCoordinates[i] = 0;

    // Compute the local basis functions in the vicinity of the point
    int noLocalBasisFunctions = IGABasis->getPolynomialDegree() + 1;

    // Loop over all the non-zero contributions
    for (int i = 0; i <= IGABasis->getPolynomialDegree(); i++) {
        // Compute iteratively the x-coordinate of the point
        _cartesianCoordinates[0] += _localBasisFunctions[i]
                * ControlPointNet[_knotSpanIndex - IGABasis->getPolynomialDegree() + i].getX();
        // Compute iteratively the y-coordinate of the point
        _cartesianCoordinates[1] += _localBasisFunctions[i]
                * ControlPointNet[_knotSpanIndex - IGABasis->getPolynomialDegree() + i].getY();
        // Compute iteratively the z-coordinate of the point
        _cartesianCoordinates[2] += _localBasisFunctions[i]
                * ControlPointNet[_knotSpanIndex - IGABasis->getPolynomialDegree() + i].getZ();
    }
}

void IGAPatch1D::computeBaseVector(double* _baseVector, double _uPrm, int _knotSpanIndex) {
    /*
     * Returns the base vector g = dR/d(xi) at the given parametric coordinate _uPrm in the Cartesian space
     * _baseVector = [gx gy gz]
     */

    // Check input
    assert(_baseVector!=NULL);

    // Initialize Control Point index
    int CPindex = 0;

    // Initialize the coordinates of the base vector
    for (int i = 0; i < 3; i++)
        _baseVector[i] = 0;

    // Compute the local basis functions and their derivatives at _uPrm
    int noLocalBasisFunctions = IGABasis->getPolynomialDegree() + 1;
    int derivDegree = 1;
    double* localBasisFunctionsAndDerivatives =
            new double[noLocalBasisFunctions * (derivDegree + 1)];
    IGABasis->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
            derivDegree, _uPrm, _knotSpanIndex);

    // Loop over all the non-zero contributions
    for (int i = 0; i <= IGABasis->getPolynomialDegree(); i++) {

        // Find the correct index for the Control Point
        CPindex = _knotSpanIndex - IGABasis->getPolynomialDegree() + i;

        // Compute iteratively the x-coordinate of the base vector
        _baseVector[0] += localBasisFunctionsAndDerivatives[noLocalBasisFunctions + i]
                * ControlPointNet[CPindex].getX();

        // Compute iteratively the x-coordinate of the base vector
        _baseVector[1] += localBasisFunctionsAndDerivatives[noLocalBasisFunctions + i]
                * ControlPointNet[CPindex].getY();

        // Compute iteratively the x-coordinate of the base vector
        _baseVector[2] += localBasisFunctionsAndDerivatives[noLocalBasisFunctions + i]
                * ControlPointNet[CPindex].getZ();

    }

    // Free the memory from the heap
    delete[] localBasisFunctionsAndDerivatives;
}

void IGAPatch1D::computeBaseVector(double* _baseVector, double* _localBasisFunctionsAndDerivatives,
        double _uPrm, int _knotSpanIndex) {
    /*
     * Returns the base vector g = dR/d(xi) at the given parametric coordinate _uPrm in the Cartesian space.
     * The function expects also that the basis functions and their derivatives have been precomputed outside the scope
     * of this function. The output is sorted as follows:
     *
     * _baseVector = [gx gy gz]
     */

    // Check input
    assert(_baseVector!=NULL);

    // Initialize Control Point index
    int CPindex = 0;

    // Initialize the coordinates of the base vector
    for (int i = 0; i < 3; i++)
        _baseVector[i] = 0;

    // Number of local basis functions
    int noLocalBasisFunctions = IGABasis->getPolynomialDegree() + 1;

    // Loop over all the non-zero contributions
    for (int i = 0; i <= IGABasis->getPolynomialDegree(); i++) {

        // Find the correct index for the Control Point
        CPindex = _knotSpanIndex - IGABasis->getPolynomialDegree() + i;

        // Compute iteratively the x-coordinate of the base vector
        _baseVector[0] += _localBasisFunctionsAndDerivatives[noLocalBasisFunctions + i]
                * ControlPointNet[CPindex].getX();

        // Compute iteratively the x-coordinate of the base vector
        _baseVector[1] += _localBasisFunctionsAndDerivatives[noLocalBasisFunctions + i]
                * ControlPointNet[CPindex].getY();

        // Compute iteratively the x-coordinate of the base vector
        _baseVector[2] += _localBasisFunctionsAndDerivatives[noLocalBasisFunctions + i]
                * ControlPointNet[CPindex].getZ();

    }
}

void IGAPatch1D::computeBaseVectorAndDerivatives(double* _baseVectorAndDerivatives,
        int _derivDegree, double _uPrm, int _knotSpanIndex) {
    /*
     *  Returns the base vector and its derivatives up to order _derivDegree. The output is sorted in consecutive order that is
     *  _baseVectorAndDerivatives = [gx gy gz dgxd(xi) dgyd(xi) dgzd(xi) ... d..._derivDegree...dgxd(xi) d..._derivDegree...dgyd(xi) d..._derivDegree...dgzd(xi)]
     */

    // Check input
    assert(_baseVectorAndDerivatives!=NULL);

    // Initialize Control Point index
    int CPindex = 0;

    // Initialize the coordinates of the base vector and its derivatives
    for (int i = 0; i < 3 * (_derivDegree + 1); i++)
        _baseVectorAndDerivatives[i] = 0;

    /*
     * Compute the local basis functions and their derivatives at _uPrm
     * Important: For the computation of the n-th derivative of the base vector the (n+1)-th derivative must be used and therefore
     * we issue the 3rd derivatives of the basis functions below
     */
    int noLocalBasisFunctions = IGABasis->getPolynomialDegree() + 1;
    double* localBasisFunctionsAndDerivatives = new double[noLocalBasisFunctions
            * (_derivDegree + 2)];
    IGABasis->computeLocalBasisFunctionsAndDerivatives(localBasisFunctionsAndDerivatives,
            _derivDegree + 1, _uPrm, _knotSpanIndex);

    // Loop over all the non-zero contributions
    for (int i = 0; i <= IGABasis->getPolynomialDegree(); i++) {

        // Find the correct index for the Control Point
        CPindex = _knotSpanIndex - IGABasis->getPolynomialDegree() + i;

        // Loop over the entries of the base vector coordinates and the coordinates of its derivatives
        for (int j = 0; j <= _derivDegree; j++) {
            // Compute iteratively the x-coordinate of the base vector and its derivatives
            _baseVectorAndDerivatives[j * 3 + 0] += localBasisFunctionsAndDerivatives[(j + 1)
                    * noLocalBasisFunctions + i] * ControlPointNet[CPindex].getX();

            // Compute iteratively the x-coordinate of the base vector and its derivatives
            _baseVectorAndDerivatives[j * 3 + 1] += localBasisFunctionsAndDerivatives[(j + 1)
                    * noLocalBasisFunctions + i] * ControlPointNet[CPindex].getY();

            // Compute iteratively the x-coordinate of the base vector and its derivatives
            _baseVectorAndDerivatives[j * 3 + 2] += localBasisFunctionsAndDerivatives[(j + 1)
                    * noLocalBasisFunctions + i] * ControlPointNet[CPindex].getZ();
        }
    }

    // Free the memory from the heap
    delete[] localBasisFunctionsAndDerivatives;
}

void IGAPatch1D::computeBaseVectorAndDerivatives(double* _baseVectorAndDerivatives,
        double* _localBasisFunctionsAndDerivatives, int _derivDegree, double _uPrm,
        int _knotSpanIndex) {
    /*
     *  Returns the base vector and its derivatives up to order _derivDegree. The basis functions and their derivatives are expected to have been
     *  precomputed outside the scope of this function. The output is sorted in consecutive order that is:
     *
     *  _baseVectorAndDerivatives = [gx gy gz dgxd(xi) dgyd(xi) dgzd(xi) ... d..._derivDegree...dgxd(xi) d..._derivDegree...dgyd(xi) d..._derivDegree...dgzd(xi)]
     */

    // Check input
    assert(_baseVectorAndDerivatives!=NULL);

    // Initialize Control Point index
    int CPindex = 0;

    // Initialize the coordinates of the base vector and its derivatives
    for (int i = 0; i < 3 * (_derivDegree + 1); i++)
        _baseVectorAndDerivatives[i] = 0;

    // Number of local basis functions
    int noLocalBasisFunctions = IGABasis->getPolynomialDegree() + 1;

    // Loop over all the non-zero contributions
    for (int i = 0; i <= IGABasis->getPolynomialDegree(); i++) {

        // Find the correct index for the Control Point
        CPindex = _knotSpanIndex - IGABasis->getPolynomialDegree() + i;

        // Loop over the entries of the base vector coordinates and the coordinates of its derivatives
        for (int j = 0; j <= _derivDegree; j++) {
            // Compute iteratively the x-coordinate of the base vector and its derivatives
            _baseVectorAndDerivatives[j * 3 + 0] += _localBasisFunctionsAndDerivatives[(j + 1)
                    * noLocalBasisFunctions + i] * ControlPointNet[CPindex].getX();

            // Compute iteratively the x-coordinate of the base vector and its derivatives
            _baseVectorAndDerivatives[j * 3 + 1] += _localBasisFunctionsAndDerivatives[(j + 1)
                    * noLocalBasisFunctions + i] * ControlPointNet[CPindex].getY();

            // Compute iteratively the x-coordinate of the base vector and its derivatives
            _baseVectorAndDerivatives[j * 3 + 2] += _localBasisFunctionsAndDerivatives[(j + 1)
                    * noLocalBasisFunctions + i] * ControlPointNet[CPindex].getZ();
        }
    }
}

void IGAPatch1D::printControlPointNet() {

    cout << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Debugging information in class IGAPatch1D" << endl;
    cout << endl;
    cout << "IGAPatch2D::NoControlPoints = " << getNoControlPoints() << endl;
    cout << endl;
    for (int i = 0; i < getNoControlPoints(); i++)
        cout << "CP " << ControlPointNet[i].getId() << ": ( " << ControlPointNet[i].getX() << " , "
                << ControlPointNet[i].getY() << " , " << ControlPointNet[i].getZ()
                << " ) --Weight--> " << ControlPointNet[i].getW() << endl;
    ;
    cout << "_____________________________________________" << endl;
    cout << endl;
}

}/* namespace EMPIRE */
