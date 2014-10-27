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

    /// The number of the Control Points in the IGAMesh
    int numNodes;
    /// The number of NOT TRIMMED OUT Control Points in the IGAMesh
    int untrimmedNumNodes;

    /// The constructor, the destructor and the copy constructor
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name The name of the IGA mesh
     * \param[in] _numControlPoints The number of the Control Points
     * \author Chenshen Wu
     ***********/
    IGAMesh(std::string _name, int _numNodes);

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
     * \param[in] _dofIndexNet The index of the dof of the each Control Points related to
     * \return The pointer to the patch just created
     * \author Chenshen Wu
     ***********/
    IGAPatchSurface* addPatch(int _pDegree, int _uNoKnots, double* _uKnotVector, int _qDegree, int _vNoKnots,
                  double* _vKnotVector, int _uNoControlPoints, int _vNoControlPoints,
                  double* controlPointNet, int* _dofIndexNet);

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
    inline const std::vector<IGAPatchSurface*>& getSurfacePatches() const {
        return surfacePatches;
    }

    /***********************************************************************************************
     * \brief Get the number of the Nodes
     * \param[out] The number of the Nodes
     * \author Chenshen Wu
     ***********/
    inline int getNumNodes() const {
        return numNodes;
    }

    int getUntrimmedNumNodes();

};

/***********************************************************************************************
 * \brief Allows for nice debug output
 * \author Fabien Pean, Chenshen Wu
 ***********/
Message &operator<<(Message &message, const IGAMesh &mesh);

}/* namespace EMPIRE */

#endif /* IGAPatchSurface_H_ */
