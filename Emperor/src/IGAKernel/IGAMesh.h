/***********************************************************************************************//**
 * \file IGAMesh.h
 * This file holds the class IGAMesh.h
 * \date 6/8/2013
 **************************************************************************************************/


#ifndef IGAMesh_H_
#define IGAMesh_H_

#include <string>
// Inclusion of user defined libraries
#include "AbstractMesh.h"

namespace EMPIRE {
class DataField;
class Message;
class IGAPatchSurface;
class IGAControlPoint;

/********//**
 * \brief class IGAMesh is a specialization of the class AbstractMesh used for IGA Mesh contains one or several IGA surface patches
 ***********/

class IGAMesh: public AbstractMesh {

protected:

	/// Array of IGA Surface Patches
	std::vector<IGAPatchSurface*> surfacePatches;

	/// Array of IGA Control Points
	std::vector<IGAControlPoint*> globalControlPoints;

	/// The array of the global ID of the Control Points
	int* controlPointID;

	/// The map from the global ID of the Control Points to the Index of its array
	std::map<int, int> cpMap;

	int numControlPoints;

public:
	/***********************************************************************************************
	 * \brief Constructor
	 * \author Chenshen Wu
	 ***********/
	IGAMesh(std::string _name, int numControlPoints, double* _globalControlPoints, int* _controlPointID);
//	IGAMesh(int numControlPoints, double* _globalControlPoints, int* _controlPointID);

	/***********************************************************************************************
	 * \brief Destructor
	 * \author Chenshen Wu
	 ***********/
	~IGAMesh();

	/***********************************************************************************************
	 * brief Add a new surface patch to this mesh
	 * \param[in] _pDegree The polynomial degree of the IGA 2D patch in the u-direction
	 * \param[in] _uNoKnots The number of knots for the knot vector in the u-direction
	 * \param[in] _uKnotVector The underlying knot vector of the IGA 2D patch in the u-direction
	 * \param[in] _qDegree The polynomial degree of the IGA 2D patch in the v-direction
	 * \param[in] _vNoKnots The number of knots for the knot vector in the v-direction
	 * \param[in] _vKnotVector The underlying knot vector of the IGA 2D patch in the v-direction
	 * \param[in] _uNoControlPoints The number of the Control Points for the 2D NURBS patch in the u-direction
	 * \param[in] _vNoControlPoints The number of the Control Points for the 2D NURBS patch in the v-direction
	 * \param[in] _controlPointNet The set of the Control Points related to the 2D NURBS patch
	 * \author Chenshen Wu
	 ***********/
	void addPatch(int _pDegree,	int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
			double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints, int* _controlPointNetID);

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


	inline std::vector<IGAPatchSurface*> getSurfacePatches(){ return surfacePatches;}
	inline int getNumControlPoints(){return numControlPoints;}
	inline std::map<int, int> getCpMap(){return cpMap;}
	/***********************************************************************************************
	 * \brief Return a pointer to the data field by its name
	 * \param[in] dataFieldName name of the data field
	 * \author Tianyang Wang
	 ***********/
	DataField *getDataFieldByName(std::string dataFieldName);
	/***********************************************************************************************
	 * \brief Compute the bounding box of the mesh
	 * \author Chenshen Wu
	 ***********/
	void computeBoundingBox();

};
Message &operator<<(Message &message, IGAMesh &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */
