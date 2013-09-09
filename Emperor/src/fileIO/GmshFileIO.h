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
 * \file GmshFileIO.h
 * This file holds the class GmshFileIO
 * \date 01/31/2013
 **************************************************************************************************/
#ifndef GMSHFILEIO_H_
#define GMSHFILEIO_H_

#include <string>

namespace EMPIRE {

/********//**
 * \brief Class GmshFileIO holds a number of static functions for Gmsh file formats
 ***********/
class GmshFileIO {
public:
    /***********************************************************************************************
     * \brief Read a .msh file, initialize mesh data
     * \param[in] fileName name of the mesh file
     * \param[out] numberOfMeshNodes number of nodes belonging to elements
     * \param[out] numberOfElements number of elements in the mesh
     * \param[out] meshNodeCoordinates coordinates of mesh nodes (x,y,z)
     * \param[out] MeshNodeIds index/id of each mesh node
     * \param[out] numberOfNodesPerElement number of nodes in each element
     * \param[out] elementNodeTables "element tables" or "connectivity tables"
     * \param[out] elementIds index/id of each element
     * \author Michael Andre
     ***********/
    static void readDotMsh(std::string fileName, int &numberOfMeshNodes, int &numberOfElements, double *&meshNodeCoordinates,
			   int *&meshNodeIds, int *&numberOfNodesPerElement, int *&elementNodeTables, int *&elementIds);
    /***********************************************************************************************
     * \brief Write a .msh file
     * \param[in] fileName name of the mesh file
     * \param[out] numberOfMeshNodes number of nodes belonging to elements
     * \param[out] numberOfElements number of elements in the mesh
     * \param[out] meshNodeCoordinates coordinates of mesh nodes (x,y,z)
     * \param[out] MeshNodeIds index/id of each mesh node
     * \param[out] numberOfNodesPerElement number of nodes in each element
     * \param[out] elementNodeTables "element tables" or "connectivity tables"
     * \param[out] elementIds index/id of each element
     * \author Michael Andre
     ***********/
     static void writeDotMsh(std::string fileName, int numberOfMeshNodes, int numberOfElements, double *&meshNodeCoordinates,
			     int *meshNodeIds, int *numberOfNodesPerElement, int *elementNodeTables, int *elementIds);
};
} /* namespace EMPIRE */

#endif
