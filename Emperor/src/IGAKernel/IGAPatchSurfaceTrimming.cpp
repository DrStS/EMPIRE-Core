 
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
#include "IGAPatchSurfaceTrimming.h"
#include "IGAPatchSurface.h"
#include "IGAControlPoint.h"
#include "IGAMath.h"
#include "Message.h"

using namespace std;

namespace EMPIRE {

IGAPatchSurfaceTrimming::IGAPatchSurfaceTrimming():outter(-1) {
}

IGAPatchSurfaceTrimming::~IGAPatchSurfaceTrimming() {
	for(int i=0;i<loops.size();i++)
		delete loops[i];
}
void IGAPatchSurfaceTrimming::addTrimInfo(int _uNoKnots, int _vNoKnots, int* _knotSpanBelonging) {
	knotSpanBelonging.resize(_uNoKnots);
	for(int i=0;i<knotSpanBelonging.size();i++)
		knotSpanBelonging[i].resize(_vNoKnots);
	for(int i=0;i<knotSpanBelonging.size();i++)
		for(int j=0;j<knotSpanBelonging[i].size();j++)
			knotSpanBelonging[i][j]=_knotSpanBelonging[i*_vNoKnots+j];
}
void IGAPatchSurfaceTrimming::addTrimLoop(int _inner, int _numCurves) {
    if(_inner) {
        loops.push_back(new IGAPatchSurfaceTrimmingLoop(_numCurves));
    	return;
    } else {
		if(outter>=0) {
			ERROR_OUT() << "Outter loop boundary for trimming has already been defined" << endl;
			exit(-1);
		}
		outter=loops.size();
        loops.push_back(new IGAPatchSurfaceTrimmingLoop(_numCurves));
	    return;
    }
}

void IGAPatchSurfaceTrimming::addTrimCurve(int _direction,int _IDBasis, int _pDegree, int _uNoKnots, double* _uKnotVector,
                                           int _uNoControlPoints, IGAControlPoint** _controlPointNet) {
    // Read input
    bool ucondition = _uNoControlPoints != _uNoKnots - _pDegree - 1;
    
    if (ucondition) {
        ERROR_OUT() << " in IGAPatchSurfaceTrimming::IGAPatchSurfaceTrimming" << endl;
        ERROR_OUT()
        << "Number of Control Points, number of knots and polynomial degree do not match!"
        << endl;
        exit(-1);
    }
    // Get the loop currently worked on
    IGAPatchSurfaceTrimmingLoop& loop=*loops.back();
    // Check that size is not going over allocated during instantiation
    assert(loop.IGABasis.size()<loop.IGABasis.capacity());
    assert(loop.direction.size()<loop.direction.capacity());
    assert(loop.uNoControlPoints.size()<loop.uNoControlPoints.capacity());
    assert(loop.ControlPointNet.size()<loop.ControlPointNet.capacity());
    // Add direction of the curve
    loop.direction.push_back(_direction);
    // Add number of control points of the curve
    loop.uNoControlPoints.push_back(_uNoControlPoints);
    // Figure out whether the patch has a B-Spline or a NURBS underlying basis
    int isNurbs = 0;
    int counter = 0;
    for (int i = 0; i < _uNoControlPoints; i++) {
        if (_controlPointNet[counter]->getW() != 1.0) {
            isNurbs = 1;
            break;
        }
        // Update the counter
        counter++;
    }
    // Create the NURBS or the B-Spline underlying basis
    if (!isNurbs) {
        loop.IGABasis.push_back(BSplineBasis1D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector));
    } else {
        double* controlPointWeights = new double[_uNoControlPoints];
        for (int i = 0; i < _uNoControlPoints; i++)
            controlPointWeights[i] = _controlPointNet[i]->getW();
        loop.IGABasis.push_back(NurbsBasis1D(_IDBasis, _pDegree, _uNoKnots, _uKnotVector, _uNoControlPoints, controlPointWeights));
    }
    // On the Control Point net
    assert(_controlPointNet != NULL);
    // Push back new vector of control points dedicated to this curve
    (loop.ControlPointNet).push_back(std::vector<IGAControlPoint*>());
    // Fill it
    for (int i = 0; i < _uNoControlPoints; i++) {
        (loop.ControlPointNet.back()).push_back(_controlPointNet[i]);
    }
}

void IGAPatchSurfaceTrimming::linearizeLoops() {
     for(int i=0;i<loops.size();i++) {
    	loops[i]->linearize();
    }
}

IGAPatchSurfaceTrimmingLoop::IGAPatchSurfaceTrimmingLoop(int _numCurves) {
	 IGABasis.reserve(_numCurves);
	 direction.reserve(_numCurves);
	 uNoControlPoints.reserve(_numCurves);
	 ControlPointNet.reserve(_numCurves);
	 assert(IGABasis.size()==0);
	 assert(IGABasis.capacity()!=0);
}

void IGAPatchSurfaceTrimmingLoop::linearize() {
    /*
     * Linearize approximation of the nurbs curves
     * Return value is a bool flag on the convergence of the Newton-Raphson iterations.
     *
     * Function layout :
     *
     * 1. For every NURBS curve
     * 1.1. For every control points of the curve
     * 1.1.1. Compute Greville abscissae defined Knot space
     * 1.1.2. Compute basis function at the Greville abscissae
     * 1.1.3. Compute position in parametric space
     * 1.1.4. Store point in data structure of polylines
     */
	for(int j=0;j<IGABasis.size();j++) {
		double* knotVector=IGABasis[j].getKnotVector();
		int pDegree=IGABasis[j].getPolynomialDegree();
		/// Check direction to put points in the right sequence (counter clockwise for outter loop, clockwise for inner)
		if(direction[j]) {
			for(int k=0;k<uNoControlPoints[j];k++) {
				double knotGreville=computeGrevilleAbscissae(k,pDegree,knotVector);
				double parametricCoordinates[2] = {0};
				computeCartesianCoordinates(j,parametricCoordinates,knotGreville);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
			}
		} else {
			for(int k=uNoControlPoints[j]-1;k>=0;k--) {
				double knotGreville=computeGrevilleAbscissae(k,pDegree,knotVector);
				double parametricCoordinates[2] = {0};
				computeCartesianCoordinates(j,parametricCoordinates,knotGreville);
				polylines.push_back(parametricCoordinates[0]);
				polylines.push_back(parametricCoordinates[1]);
			}
		}
	}
	cleanPolygon();
}

double IGAPatchSurfaceTrimmingLoop::computeGrevilleAbscissae(const int cp, const int pDeg, const double*const knotVector) {
	double GrevilleAbscissae=0;
	int limit=max(pDeg-1,1);
	for(int i=0;i<limit;i++) {
		GrevilleAbscissae+=knotVector[limit+cp+i];
	}
	GrevilleAbscissae/=limit;
	if(GrevilleAbscissae<knotVector[limit+cp]-EPS || GrevilleAbscissae>knotVector[limit+cp+(limit-1)]+EPS) {
		ERROR_OUT()<<"Greville abscissae "<<GrevilleAbscissae<<" is out of Knot vector bound ["
				<<knotVector[limit+cp]<<", "<<knotVector[limit+cp+(limit-1)]<<"]"<<endl;
		exit(-1);
	}
	return GrevilleAbscissae;
}

void IGAPatchSurfaceTrimmingLoop::cleanPolygon() {
	// Remove duplicated points and consecutive aligned points
	int tmp_n_pts=polylines.size()/2;
	for(int k=0;k<tmp_n_pts;k++) {
		int p1=k%tmp_n_pts;
		int p2=(k+1)%tmp_n_pts;
		int p3=(k+2)%tmp_n_pts;
		bool isSame=polylines[p1*2]==polylines[p2*2] && polylines[p1*2+1]==polylines[p2*2+1];

		double v1x=polylines[p2*2]-polylines[p1*2];
		double v1y=polylines[p2*2+1]-polylines[p1*2+1];
		double v2x=polylines[p3*2]-polylines[p1*2];
		double v2y=polylines[p3*2+1]-polylines[p1*2+1];
		double v1[3]={v1x, v1y, 0};
		double v2[3]={v2x, v2y, 0};
		double v[3];
		crossProduct(v,v1,v2);
		// Result of cross product only in Z direction
		bool isColinear=fabs(v[2])>1e-9?false:true;
		if(isSame || isColinear) {
			// Remove middle point, noted as idx2 here
			polylines.erase(polylines.begin()+p2*2+1);
			polylines.erase(polylines.begin()+p2*2);
			// Update loop to keep loop traversal consistent
			tmp_n_pts=polylines.size()/2;
			k--;
		}
	}
}

void IGAPatchSurfaceTrimmingLoop::computeCartesianCoordinates(int _curve, double* _cartesianCoordinates, double _uPrm,
        int _uKnotSpanIndex) {
    // Read input
    assert(_cartesianCoordinates != NULL);

    // Initialize the coordinates of the point
    for (int i = 0; i < 2; i++)
        _cartesianCoordinates[i] = 0;
    // Compute the local basis functions in the vicinity of the point
    int pDegree = IGABasis[_curve].getPolynomialDegree();
    int noLocalBasisFunctions = IGABasis[_curve].computeNoBasisFunctions();
    double* localBasisFunctions = new double[noLocalBasisFunctions];
    IGABasis[_curve].computeLocalBasisFunctions(localBasisFunctions,_uPrm,_uKnotSpanIndex);
    // Initialize the Control Point index
    int CPindex = 0;
    // Initialize a basis functions counter
    int counter_basis = 0;
    // Loop over all the non-zero contributions
	for (int i = 0; i <= pDegree; i++) {

		// Update the correct index for the Control Points in 2D. Pattern A[i][j] = V[j*n+i]
		CPindex =(_uKnotSpanIndex - pDegree + i);

		// Compute iteratively the x-coordinate of the point
		_cartesianCoordinates[0] += localBasisFunctions[counter_basis]
				* ControlPointNet[_curve][CPindex]->getX();
		// Compute iteratively the y-coordinate of the point
		_cartesianCoordinates[1] += localBasisFunctions[counter_basis]
				* ControlPointNet[_curve][CPindex]->getY();

		// Update basis function's counter
		counter_basis++;
	}
    // Free the memory from the heap
    delete[] localBasisFunctions;
}

void IGAPatchSurfaceTrimmingLoop::computeCartesianCoordinates(int _curve, double* _cartesianCoordinates, double _localCoordinates) {
    int _uKnotSpanIndex = IGABasis[_curve].findKnotSpan(_localCoordinates);
    IGAPatchSurfaceTrimmingLoop::computeCartesianCoordinates(_curve, _cartesianCoordinates, _localCoordinates, _uKnotSpanIndex);
}

Message &operator<<(Message &message, const IGAPatchSurfaceTrimming &trim) {
    message << "\t" << "---------------------------------Start Trimming" << endl;
    message << "\t" << "Trimming Info"<< endl;
    message << "\t\t" << "Number of loops: "<<trim.getNumOfLoops()<< endl;
    message << "\t\t" << "Outer loop index: "<<trim.getOutterLoopIndex()<< endl;
    message << "\t" << "Outter Loop["<<trim.getOutterLoopIndex()<<"]"<< endl;
    message << trim.getOutterLoop();
	for(int i=0; i<trim.getNumOfLoops();i++) {
		if(i==trim.getOutterLoopIndex()) continue;
    	message << "\t" << "InnerLoop["<<i<<"]"<< endl;
		message << trim.getLoop(i);
	}
    message << "\t" << "---------------------------------End Trimming" << endl;
	return message;
}

Message &operator<<(Message &message, const IGAPatchSurfaceTrimmingLoop &trim) {
    /// output loop
    for(int i=0;i<trim.getIGABasis().size();++i) {
        message << "\t" << "Curve["<<i<<"]"<< endl;
        message << "\t\tpDegree:  " << trim.getIGABasis(i).getPolynomialDegree()<< endl;

        message << "\t\tKnots Vector U: \t";
        for (int k = 0; k < trim.getIGABasis(i).getNoKnots(); k++)
            message << trim.getIGABasis(i).getKnotVector()[k] << "  ";
        message << endl;
        message << "\t\t" << "number of control points: " << trim.getNoControlPoints(i) << endl;

        message << "\t\tControl Points Net: " << endl;
        for (int k = 0; k < trim.getNoControlPoints(i); k++) {
            message << "\t\t";
            message << trim.getControlPointNet(i)[k]->getX() << ", "
                    << trim.getControlPointNet(i)[k]->getY() << ", "
                    << trim.getControlPointNet(i)[k]->getZ() << endl;
        }
    }
    message << "\tLinear Polygon: " << endl;
    for (int k = 0; k < trim.getPolylines().size()/2; k++) {
        message << "\t\t";
        message << trim.getPolylines()[2*k] << ", "
                << trim.getPolylines()[2*k+1]<< endl;
    }
    return message;
}

}/* namespace EMPIRE */
