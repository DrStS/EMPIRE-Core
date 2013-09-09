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
#ifndef ABSTRACTMESH_H_
#define ABSTRACTMESH_H_

#include <string>
#include <map>
#include <vector>
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class DataField;
class Message;
/********//**
 * \brief Class AbstractMesh is the superclass of FEMesh and IGAMesh
 ***********/
class AbstractMesh {
public:
    /***********************************************************************************************
     * \brief Constructor, allocating the storage of the mesh
     * \param[in] _name name of the mesh
     * \param[in] _type type of the mesh
     * \author Tianyang Wang
     ***********/
    AbstractMesh(std::string _name, EMPIRE_Mesh_type _type);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~AbstractMesh();
    /***********************************************************************************************
     * \brief Add a new data field to this mesh
     * \param[in] dataFieldName name of the data field
     * \param[in] location at node or at element centroid
     * \param[in] dimension vector or scalar
     * \param[in] typeOfQuantity field or field integral
     * \author Tianyang Wang
     ***********/
    virtual void addDataField(std::string dataFieldName, EMPIRE_DataField_location location,
            EMPIRE_DataField_dimension dimension,
            EMPIRE_DataField_typeOfQuantity typeOfQuantity) = 0;
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
    virtual void computeBoundingBox() = 0;

    /// name of the mesh
    const std::string name;
    /// the map of data fields
    std::map<std::string, DataField*> nameToDataFieldMap;
    /// type of the mesh
    const EMPIRE_Mesh_type type;

    /// the boundingBox of the mesh
    struct structBoundingBox {
        double xmax;
        double xmin;
        double ymax;
        double ymin;
        double zmax;
        double zmin;
        bool isComputed;
    } boundingBox;
};

/***********************************************************************************************
 * \brief Output bounding box
 * \author Tianyang Wang
 ***********/
Message &operator<<(Message &message, AbstractMesh::structBoundingBox &boundingBox);

} /* namespace EMPIRE */

#endif /* ABSTRACTMESH_H_ */
