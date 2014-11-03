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
/***********************************************************************************************//**
 * \file IGAPatchSurface.h
 * This file holds the class IGAPatchSurface.h
 * \date 28/5/2013
 **************************************************************************************************/

#ifndef IGAPatchSurface_H_
#define IGAPatchSurface_H_

// Inclusion of user defined libraries
#include "AbstractMesh.h"
#include "NurbsBasis2D.h"
#include "IGAControlPoint.h"
#include "IGAPatchSurfaceTrimming.h"
#include <limits>
#include <set>

namespace EMPIRE {
class DataField;
class Message;

/********//**
 * \brief class IGAPatchSurface is a specialization of the class AbstractMesh used for IGA Mesh with two parameters like shell elements
 ***********/

class IGAPatchSurface {

protected:
    /// The basis functions of the 2D NURBS patch
    BSplineBasis2D* IGABasis;

    /// Number of Control Points in u-direction
    int uNoControlPoints;

    /// Number of Control Points in v-direction
    int vNoControlPoints;

    /// The set of the Control Points of the patch
    IGAControlPoint** ControlPointNet;
    
    /// The class holding the trimming information
    IGAPatchSurfaceTrimming Trimming;

    /// The constructor and the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _IDBasis The id of the underlying basis to the IGA 2D patch
     * \param[in] _pDegree The polynomial degree of the IGA 2D patch in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
     * \param[in] _qDegree The polynomial degree of the IGA 2D patch in the v-direction
     * \param[in] _vNoKnots The number of knots for the knot vector in the v-direction
     * \param[in] _vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
     * \param[in] _vNoControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 2D NURBS patch
     * \author Andreas Apostolatos
     ***********/
    IGAPatchSurface(int, int, int, double*, int, int, double*, int, int, IGAControlPoint**);

    /***********************************************************************************************
     * \brief Destructor
     * \author Andreas Apostolatos
     ***********/
    ~IGAPatchSurface();

    
    /// Trimming related functions
public:
    void addTrimInfo(int* _knotSpanBelonging);
    /***********************************************************************************************
     * \brief Setup information about the loop soon to be received
     * \param[in] inner 0 for outter and 1 for inner
     * \param[in] numCurves Number of curves to be received for this loop 
     * \author Fabien Pean
     ***********/
    void addTrimLoop(int inner, int numCurves);
    /***********************************************************************************************
     * \brief Add a Nurbs curve for the current loop and its attached information
     * \param[in] inner Value 0 for outter and 1 for inner
     * \param[in] direction The direction of the curve if is following standard or not
     * \param[in] ID The id of the curve
     * \param[in] _pDegree The polynomial degree of the IGA 1D curve in the u-direction
     * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
     * \param[in] _uKnotVector The underlying knot vector of the IGA 1D curve in the u-direction
     * \param[in] _uNoControlPoints The number of the Control Points for the 1D NURBS patch in the u-direction
     * \param[in] _controlPointNet The set of the Control Points related to the 1D NURBS patch
     * \author Fabien Pean
     ***********/
    void addTrimCurve(int direction, int _pDegree, int _uNoKnots, double* _uKnotVector,
                      int _uNoControlPoints, double* _controlPointNet);
    
    void linearizeTrimming(){Trimming.linearizeLoops();};

    void getUntrimmedCPindexes(std::set<int>& out);
private:
    void addCPidsToSet(std::set<int>& CPids,const int uSpan, const int vSpan);

    /// Basis related functions
public:
    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters are known
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the Patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vPrm The parameter on the v-coordinate line
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double, int, double, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctions The local basis functions
     * \compute the knot span Index inside the function. Convenient but in-efficient.
     * \author Chenshen Wu
     ***********/
    void computeCartesianCoordinates(double*, double*);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctions The local basis functions
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double*, int, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of a point on a NURBS surface whose surface parameters and the local basis functions and their derivatives are given
     * \param[in/out] _cartesianCoordinates The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFctsAndDerivs The local basis functions and their derivatives
     * \param[in] _derivDegree The derivative degree up to which the basis functions have been computed contained ion array _localBasisFctsAndDerivs
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeCartesianCoordinates(double*, double*, int, int, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the base vectors at a given pair of surface parameters given the basis functions and their derivatives
     * \param[in/out] _baseVectors The Cartesian coordinates of the base vectors on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctionsAndDerivatives The local basis functions and their derivatives
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectors(double* _baseVectors, double* _localBasisFunctionsAndDerivatives,
            int _uKnotSpanIndex, int _vKnotSpanIndex);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the base vectors at a given pair of surface parameters given the basis functions and their derivatives
     * \param[in/out] _baseVectors The Cartesian coordinates of the base vectors on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vPrm The parameter on the v-coordinate line
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Chenshen Wu
     ***********/
    void computeBaseVectors(double* _baseVectors, double _uPrm, int _uKnotSpanIndex, double _vPrm,
            int _vKnotSpanIndex);

    /***********************************************************************************************
     * \brief Returns the index of the i-th partial derivative w.r.t. to u , j-th partial derivative w.r.t v of the _componentIndex-th component to the l-th base vector
     * \param[out] The index of the _uDerivIndex-th partial derivative w.r.t. to u , j-th partial derivative w.r.t v of the k-th component to the l-th base vector
     * \param[in] _derivDegree The absolute order of the partial derivatives to the base vectors
     * \param[in] _uDerivIndex The order of the partial derivative w.r.t. u-parametric coordinate, _uDerivIndex = 0 , … , _derivDegree
     * \param[in] _vDerivIndex The order of the partial derivative w.r.t. v-parametric coordinate, _uDerivIndex = 0 , … , _derivDegree - _uDerivIndex
     * \param[in] _componentIndex The component of the base vector, _componentIndex = 1, 2, 3
     * \param[in] _baseVecIndex The base vector for which to compute the derivatives _baseVecIndex = 1, 2
     * \author Andreas Apostolatos
     ***********/
    int indexDerivativeBaseVector(int, int, int, int, int);

    /***********************************************************************************************
     * \brief Returns the Cartesian Coordinates of the base vectors and their derivatives at a given pair of surface parameters given the basis functions and their derivatives
     * \param[in/out] _baseVectorsAndDerivatives The Cartesian coordinates of the point on the patch whose surface parameters are _uPrm and _vPrm
     * \param[in] _localBasisFunctionsAndDerivatives The local basis functions and their derivatives
     * \param[in] _derivDegree The derivative order with respect to both u-,v- directions
     * \param[in] _uPrm The parameter on the u-coordinate line
     * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
     * \param[in] _vPrm The parameter on the v-coordinate line
     * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
     * \author Andreas Apostolatos
     ***********/
    void computeBaseVectorsAndDerivatives(double*, double*, int, int, int);

    /// Projection related functions
public:
    /***********************************************************************************************
     * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the NURBS pacth
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _P Given the Cartesian components of the point to be projected on the NURBS patch it is returned the Cartesian components of its orthogonal projection
     * \param[in/out] _flagConverge Flag indicating whether the Newton iterations have converged true/false
     * \author Andreas Apostolatos
     ***********/
    bool computePointProjectionOnPatch(double&, double&, double*, bool&);

    /***********************************************************************************************
     * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the NURBS pacth (overloaded)
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _P Given the Cartesian components of the point to be projected on the NURBS patch it is returned the Cartesian components of its orthogonal projection
     * \author Andreas Apostolatos
     ***********/
    bool computePointProjectionOnPatch(double&, double&, double*);

    /***********************************************************************************************
     * \brief Returns the point on the given NURBS patch boundary which defines an orthogonal projection from the given line to the NURBS boundary
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _t The running parameter on the given NURBS patch boundary
     * \param[in/out] _ratio The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _edge (0,1,2,3) --> (uRunsvStart,uRunsvEnd,uStartvRuns,uEndvRuns)
     * \author Chenshen Wu
     ***********/
    bool computePointProjectionOnPatchBoundaryOnGivenEdge(double& _t, double& _ratio,
            double& _distance, double* _P1, double* _P2, int _edge);
    bool computePointProjectionOnPatchBoundaryOnGivenEdge_Brute(double& _u,double& _v, double& _ratio,
            double& _distance, double* _P1, double* _P2);
    /***********************************************************************************************
     * \brief Returns the point on the given NURBS patch boundary which defines an orthogonal projection from the given line to the NURBS boundary
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _ratio The ratio between the line segment that is projected on the NURBS patch to the complete line segment
     * \param[in/out] _distance The orthogonal distance from the NURBS surface to the line segment
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \author Chenshen Wu
     ***********/
    bool computePointProjectionOnPatchBoundary(double& _u, double& _v, double& _ratio,
            double& _distance, double* _P1, double* _P2);

    /***********************************************************************************************
     * \brief Returns the point on the NURBS patch boundary which is closest to the physical point provided
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _t The running parameter on the given NURBS patch boundary
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _distance The distance from the NURBS boundary to the physical point
     * \param[in] _P1 The physical point
     * \param[in] _edge (0,1,2,3) --> (uRunsvStart,uRunsvEnd,uStartvRuns,uEndvRuns)
     * \author Fabien Pean
     ***********/
    bool computePointMinimumDistanceToPatchBoundaryOnGivenEdge(double& _w,
            double& _distance, double* _P1, int _edge);
    /***********************************************************************************************
     * \brief Returns the point on the NURBS patch boundary which is closest to the physical point provided
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _distance The distance from the point to the physical point
     * \param[in] _P1 The physical point
     * \param[out] _edge The list of edge the point is on
     * \author Fabien Pean
     ***********/
    void computePointMinimumDistanceToPatchBoundary(double& _u, double& _v, double& _distance, double* _P1, int* _edge);
    /***********************************************************************************************
     * \brief Returns the point on the NURBS patch boundary which is closest to the line segment
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _distance The distance from the point to the line segment
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \param[in] _edge (0,1,2,3) --> (uRunsvStart,uRunsvEnd,uStartvRuns,uEndvRuns)
     * \author Chenshen Wu
     ***********/
    bool computeLineMinimumDistanceToPatchBoundaryOnGivenEdge(double& _t, double& _distance,
            double* _P1, double* _P2, int _edge);

    /***********************************************************************************************
     * \brief Returns the point on the NURBS patch boundary which is closest to the line segment
     * \param[out] The flag on whether or not the Newton-Raphson iterations have converged for the defined set of parameters
     * \param[in/out] _u Given is the initial guess for the Newton-Raphson iterations and returned value is the converged u-surface parameter
     * \param[in/out] _v Given is the initial guess for the Newton-Raphson iterations and returned value is the converged v-surface parameter
     * \param[in/out] _distance The distance from the point to the line segment
     * \param[in] _P1 The first point of the line segment
     * \param[in] _P2 The second point of the line segment
     * \author Chenshen Wu
     ***********/
    void computeLineMinimumDistanceToPatchBoundary(double& _u, double& _v, double& _distance,
            double* _P1, double* _P2);

    /***********************************************************************************************
     * \brief Find the nearest knot intersection on the patch as an initial guess for the projection
     * \param[in/out] _u Given is the u-surface parameter of the nearest knot intersection.
     * \param[in/out] _v Given is the v-surface parameter of the nearest knot intersection.
     * \param[in] _P Given the Cartesian components of the point to be projected on the NURBS patch
     * \param[in] _uDiv Given the number of division in u-direction
     * \param[in] _vDiv Given the number of division in u-direction
     * \author Chenshen Wu
     ***********/
    void findInitialGuess4PointProjection(double& _u, double& _v, double* _P, int uDiv = 5,
            int _vDiv = 5);

    /***********************************************************************************************
     * \brief Find the nearest knot intersection on the patch as an initial guess for the projection
     * \param[in/out] _coords The Cartesian coordinates of the point on the patch
     * \param[in/out] _normal The normal to the patch vector
     * \param[in] _u Given is the u-surface parameter
     * \param[in] _v Given is the v-surface parameter
     * \author Chenshen Wu
     ***********/
    void computeCartesianCoordinatesAndNormalVector(double* _coords, double* _normal, double _u,
            double _v);

    /// Postprocessing functions
public:
    /***********************************************************************************************
     * \brief Returns the approximate value of a field given its values on the Control Points at the specified parametric location
     * \param[in] _u Given is the u-surface parameter
     * \param[in] _v Given is the v-surface parameter
     * \param[in] _valuesOnCP The scalar values on the Control Points
     * \param[out] The approximate value of the scalar on (_u,_v)
     * \author Chenshen Wu
     ***********/
    double computePostprocessingScalarValue(double _u, double _v, double* _valuesOnCP);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the underlying IsoGeometric basis of the patch
     * \author Andreas Apostolatos
     ***********/
    inline BSplineBasis2D* getIGABasis() const {
        return IGABasis;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch in u-direction
     * \author Andreas Apostolatos
     ***********/
    inline int getUNoControlPoints() const {
        return uNoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch in v-direction
     * \author Andreas Apostolatos
     ***********/
    inline int getVNoControlPoints() const {
        return vNoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points of the patch
     * \author Chenshen Wu
     ***********/
    inline int getNoControlPoints() const {
        return uNoControlPoints * vNoControlPoints;
    }

    /***********************************************************************************************
     * \brief Get the Control Points of the patch
     * \author Andreas Apostolatos
     ***********/
    inline IGAControlPoint** getControlPointNet() const {
        return ControlPointNet;
    }

    /***********************************************************************************************
     * \brief Find know span on u direction
     * \author Chenshen Wu
     ***********/
    inline int findSpanU(double _u) {
        return getIGABasis()->getUBSplineBasis1D()->findKnotSpan(_u);
    }

    /***********************************************************************************************
     * \brief Find know span on v direction
     * \author Chenshen Wu
     ***********/
    inline int findSpanV(double _v) {
        return getIGABasis()->getVBSplineBasis1D()->findKnotSpan(_v);
    }

    /***********************************************************************************************
     * \brief Get Trimming class
     * \author Fabien Pean
     ***********/
    inline IGAPatchSurfaceTrimming& getTrimming() {
    	return Trimming;
    }
    inline const IGAPatchSurfaceTrimming& getTrimming() const {
        return Trimming;
    }
    /***********************************************************************************************
     * \brief Check if patch is trimmed
     * \author Fabien Pean
     ***********/
    inline bool isTrimmed() const {
    	return Trimming.isTrimmed();
    }

    /// The maximum number of Newton-Raphson iterations for the computation of the orthogonal projection of point on the NURBS patch
    static const int MAX_NUM_ITERATIONS;

    /// The tolerance for the Newton-Raphson iterations for the computation of the orthogonal projection of point on the NURBS patch
    static const double EPS_ORTHOGONALITY_CONDITION;

    /// The tolerance for the Newton-Raphson iterations for the computation of the orthogonal projection of point on the NURBS patch for points which are expected to be projected in irregular locations of the patch
    static const double EPS_ORTHOGONALITY_CONDITION_RELAXED;

    /// The tolerance for the distance of the computed point to the surface
    static const double EPS_DISTANCE;

    /// The tolerance for the distance of the computed point to the surface for points which are expected to be projected in irregular locations of the patch
    static const double EPS_DISTANCE_RELAXED;
};

/***********************************************************************************************
 * \brief Allows for nice debug output
 * \author Fabien Pean, Chenshen Wu
 ***********/
Message &operator<<(Message &message, const IGAPatchSurface &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */
