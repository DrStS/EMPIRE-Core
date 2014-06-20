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
/***********************************************************************************************//**
 * \file FEMesh.h
 * This file holds the class FEMesh
 * \date 3/5/2012
s **************************************************************************************************/

#ifndef FEMESH_H_
#define FEMESH_H_

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"
#include "AbstractMesh.h"

namespace EMPIRE {
class DataField;
class Message;
/********//**
 * \brief Class FEMesh has all data w.r.t. a finite element mesh
 ***********/
class FEMesh : public AbstractMesh {
public:
    /***********************************************************************************************
     * \brief Constructor, allocating the storage of the mesh
     * \param[in] _name name of the mesh
     * \param[in] _numNodes number of nodes
     * \param[in] _numElems number of elements
     * \param[in] _triangulateAll triangulate all elements
     * \author Tianyang Wang
     ***********/
    FEMesh(std::string _name, int _numNodes, int _numElems, bool _triangulateAll = false);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~FEMesh();
    /***********************************************************************************************
     * \brief The array elems can only be initialized after numNodesPerElem is defined.
     * Determine mesh type.
     * \author Tianyang Wang
     ***********/
    void initElems();
    /***********************************************************************************************
     * \brief Add a new data field to this mesh
     * \param[in] dataFieldName name of the data field
     * \param[in] location at node or at element centroid
     * \param[in] dimension vector or scalar
     * \param[in] typeOfQuantity field or field integral
     * \author Tianyang Wang
     ***********/
    void addDataField(std::string dataFieldName, EMPIRE_DataField_location location,
            EMPIRE_DataField_dimension dimension, EMPIRE_DataField_typeOfQuantity typeOfQuantity);
    /***********************************************************************************************
     * \brief Triangulate the mesh and return the triangulated mesh
     * \return the triangualted mesh, NULL if it is not a to-be-triangulated mesh.
     * \author Tianyang Wang
     ***********/
    FEMesh *triangulate();
    /***********************************************************************************************
     * \brief Compute the bounding box of the mesh
     * \author Tianyang Wang
     ***********/
    void computeBoundingBox();

    /// triangulate all elments
    bool triangulateAll;
    /// number of nodes
    const int numNodes;
    /// coordinates of all nodes
    double *nodes;
    /// IDs of all nodes
    int *nodeIDs;
    /// number of elements
    const int numElems;
    /// number of nodes of each element
    int *numNodesPerElem;
    /// nodes connectivity inside all elements
    int *elems;
    /// size of array elems
    int elemsArraySize;
    /// IDs of all elements (now it is not received from clients, therefore it is fixed as 1,2,3...)
    int *elemIDs;
private:
    /// to be triangualted or not
    bool tobeTriangulated;
    /// the triangulated mesh
    FEMesh *triangulatedMesh;
    /// unit test class
    friend class TestFEMesh;
};

/***********************************************************************************************
 * \brief Revert the surface normal of the mesh
 * \param[in/out] mesh the mesh
 * \author Tianyang Wang
 ***********/
void revertSurfaceNormalOfFEMesh(FEMesh *mesh);
/***********************************************************************************************
 * \brief Allows for nice debug output later
 * \author Stefan Sicklinger
 ***********/
Message &operator<<(Message &message, FEMesh &mesh);

} /* namespace EMPIRE */
#endif /* FEMESH_H_ */
