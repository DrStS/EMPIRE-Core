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
#include <limits>

namespace EMPIRE {
class DataField;
class Message;

/********//**
 * \brief class IGAPatchSurface is a specialization of the class AbstractMesh used for IGA Mesh with two parameters like shell elements
 ***********/

class IGAPatchSurface: public AbstractMesh {

protected:
	/// The basis functions of the 2D NURBS patch
	BSplineBasis2D* IGABasis;

	/// Number of Control Points in u-direction
	int uNoControlPoints;

	/// Number of Control Points in v-direction
	int vNoControlPoints;

	/// The set of the Control Points of the patch
	IGAControlPoint** ControlPointNet;

	/// The constructor and the destructor and the copy constructor
public:
	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _name name of the mesh
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
	IGAPatchSurface(std::string, int, int, int, double*, int, int, double*, int,
			int, IGAControlPoint**);

	/***********************************************************************************************
	 * \brief Constructor
	 * \param[in] _name name of the mesh
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
	IGAPatchSurface(std::string, int, int, int, double*, int, int, double*, int,
			int, double *);

	/***********************************************************************************************
	 * \brief Destructor
	 * \author Andreas Apostolatos
	 ***********/
	~IGAPatchSurface();

	/// Functions from AbstractMesh
	/***********************************************************************************************
	 * \brief Add a new data field to this mesh
	 * \param[in] dataFieldName name of the data field
	 * \param[in] location at node or at element centroid
	 * \param[in] dimension vector or scalar
	 * \param[in] typeOfQuantity field or field integral
	 * \author Tianyang Wang
	 ***********/
	void addDataField(std::string dataFieldName,
			EMPIRE_DataField_location location,
			EMPIRE_DataField_dimension dimension,
			EMPIRE_DataField_typeOfQuantity typeOfQuantity);
	/***********************************************************************************************
	 * \brief Return a pointer to the data field by its name
	 * \param[in] dataFieldName name of the data field
	 * \author Tianyang Wang
	 ***********/
	DataField *getDataFieldByName(std::string dataFieldName);
	/***********************************************************************************************
	 * \brief Compute the bounding box of the mesh
	 * \author Tianyang Wang
	 ***********/
	void computeBoundingBox();

	/// Basis related functions

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
	 * \param[in] _uPrm The parameter on the u-coordinate line
	 * \param[in] _uKnotSpanIndex The index of the knot span where the parametric coordinates _uPrm lives in
	 * \param[in] _vPrm The parameter on the v-coordinate line
	 * \param[in] _vKnotSpanIndex The index of the knot span where the parametric coordinates _vPrm lives in
	 * \author Andreas Apostolatos
	 ***********/
	void computeBaseVectors(double*, double*, int, int);

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

	/***********************************************************************************************
	 * \brief Computes the orthogonal projection of point of the 3D Euclidean space onto the NURBS pacth
	 * \param[out] The flag on whether or not the Newton-Rapson iterations have converged for the defined set of parameters
	 * \param[in/out] _u Given is the initial guess for the Newton-Rapson iterations and returned value is the converged u-surface parameter
	 * \param[in/out] _v Given is the initial guess for the Newton-Rapson iterations and returned value is the converged v-surface parameter
	 * \param[in/out] _P Given the Cartesian components of the point to be projected on the NURBS patch it is returned the Cartesian components of its orthogonal projection
	 * \author Andreas Apostolatos
	 ***********/
	bool computePointProjectionOnPatch(double&, double&, double*, bool&);
	bool computePointProjectionOnPatch(double&, double&, double*);

	/***********************************************************************************************
	 * \brief Find the nearest knot intersection on the patch as an initial guess for the projection
	 * \param[out] _u Given is the u-surface parameter of the nearest knot intersection.
	 * \param[out] _v Given is the v-surface parameter of the nearest knot intersection.
	 * \param[in] _P Given the Cartesian components of the point to be projected on the NURBS patch
	 * \author Chenshen Wu
	 ***********/
	void findNearestKnotIntersection(double&, double&, double*);

	/// Get and set functions

	/***********************************************************************************************
	 * \brief Get the underlying IsoGeometric basis of the patch
	 * \author Andreas Apostolatos
	 ***********/
	inline BSplineBasis2D* getIGABasis() {
		return IGABasis;
	}

	/***********************************************************************************************
	 * \brief Get the number of the Control Points of the patch in u-direction
	 * \author Andreas Apostolatos
	 ***********/
	inline int getUNoControlPoints() {
		return uNoControlPoints;
	}

	/***********************************************************************************************
	 * \brief Get the number of the Control Points of the patch in v-direction
	 * \author Andreas Apostolatos
	 ***********/
	inline int getVNoControlPoints() {
		return vNoControlPoints;
	}

	inline int getNoControlPoints() {
		return uNoControlPoints * vNoControlPoints;
	}

	/***********************************************************************************************
	 * \brief Get the Control Points of the patch
	 * \author Andreas Apostolatos
	 ***********/
	inline IGAControlPoint** getControlPointNet() {
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

	/// DEBUGGING functions

	/***********************************************************************************************
	 * \brief Prints the Control Point net for the 2D NURBS patch
	 * \author Andreas Apostolatos
	 ***********/
	void printControlPointNet();

	void printSelf();

	/// The maximum number of Newton-Rapson iterations for the computation of the orthogonal projection of point on the NURBS patch
	static const int MAX_NUM_ITERATIONS;

	/// Expected number of iterations for convergence for a regular problem set up with small distance between the fluid and the structural mesh
	static const int REGULAR_NUM_ITERATIONS;

	/// The tolerance for the Newton-Rapson iterations for the computation of the orthogonal projection of point on the NURBS patch
	static const double EPS_ORTHOGONALITY_CONDITION;

	/// The tolerance for the distance of the computed point to the surface
	static const double EPS_DISTANCE;

	/***********************************************************************************************
	 * \brief Allows for nice debug output
	 * \author Chenshen Wu
	 ***********/

};

Message &operator<<(Message &message, IGAPatchSurface &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */
