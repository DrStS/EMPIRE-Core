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
 * \brief class IGAMesh is a specialization of the class AbstractMesh used for IGA Mesh containing number of IGA surface patches
 ***********/

class IGAMesh: public AbstractMesh {

protected:
    /// Array of IGA Surface Patches
    std::vector<IGAPatchSurface*> surfacePatches;

    /// Array of IGA Control Points
    std::vector<IGAControlPoint*> globalControlPoints;

    /// The global IDs of the Control Points in the IGAMesh
    int* controlPointID;

    /// The map from the the global IDs of the Control Points to the elements in the array controlPointID
    std::map<int, int> mapControlPointIDToIndex;

    /// The number of the Control Points in the IGAMesh
    int numControlPoints;

    /// The constructor, the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the IGA mesh
     * \param[in] _numControlPoints The number of the Control Points
     * \param[in] _globalControlPoints The coordinates and the weights of the Control Points sorted in an array
     * \param[in] _controlPointID The Control Point IDs sorted in an array
     * \author Chenshen Wu
     ***********/
    IGAMesh(std::string _name, int _numControlPoints, double* _globalControlPoints,
            int* _controlPointID);

    /***********************************************************************************************
     * \brief Destructor
     * \author Chenshen Wu
     ***********/
    ~IGAMesh();

    /***********************************************************************************************
     * brief Add a new surface patch to the IGA mesh
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
    void addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
            double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
            int* _controlPointNetID);

    /// Specializing abstract functions from AbstractMesh class
public:
    /***********************************************************************************************
     * \brief Add a new data field to this mesh
     * \param[in] _dataFieldName name of the data field
     * \param[in] _location at node or at element centroid
     * \param[in] _dimension vector or scalar
     * \param[in] _typeOfQuantity field or field integral
     * \author Chenshen Wu
     ***********/
    void addDataField(std::string _dataFieldName, EMPIRE_DataField_location _location,
            EMPIRE_DataField_dimension _dimension, EMPIRE_DataField_typeOfQuantity _typeOfQuantity);

    /***********************************************************************************************
     * \brief Compute the bounding box of the mesh
     * \author Chenshen Wu
     ***********/
    void computeBoundingBox();

    /// Postprocessing
public:
    /***********************************************************************************************
     * \brief Returns the displacement component at the specified parametric location
     * \param[in] _dataFieldName name of the data field
     * \param[in] _patchid The ID of the patch
     * \param[in] _u The parametric coordinate u
     * \param[in] _v The parametric coordinate v
     * \param[in] _component The component of the displacement field (0,1,2)=(x,y,z)
     * \author Chenshen Wu
     ***********/
    double computeDisplacementComponent(std::string _dataFieldName, int _patchid, double _u, double _v,
            int _component);

    /// Get and set functions
public:
    /***********************************************************************************************
     * \brief Get the surface patches
     * \param[out] A container vector of type std::vector<IGAPatchSurface*>
     * \author Chenshen Wu
     ***********/
    inline std::vector<IGAPatchSurface*> getSurfacePatches() {
        return surfacePatches;
    }

    /***********************************************************************************************
     * \brief Get the map of the global ID of the Control Points to the Index of its array controlPointID
     * \param[out] The map of the global ID of the Control Points to the Index of its array
     * \author Chenshen Wu
     ***********/
    inline std::map<int, int> getMapControlPointIDToIndex() {
        return mapControlPointIDToIndex;
    }

    /***********************************************************************************************
     * \brief Get the number of the Control Points
     * \param[out] The number of the Control Points
     * \author Chenshen Wu
     ***********/
    inline int getNumControlPoints() {
        return numControlPoints;
    }

    /// DeBug
public:
    /***********************************************************************************************
     * \brief Prints all the patches of the IGAMesh
     ***********/
    void print();
};

/***********************************************************************************************
 * \brief Allows for nice debug output later
 * \author Chenshen Wu
 ***********/
Message &operator<<(Message &message, IGAMesh &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */
