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
 * \file GiDFileIO.h
 * This file holds functions for GiD file formats
 * \date 3/15/2013
 **************************************************************************************************/
#ifndef GIDFILEIO_H_
#define GIDFILEIO_H_

#include <string>

namespace GiDFileIO {
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
 * \author Michael Andre, Tianyang Wang
 ***********/
void readDotMsh(std::string fileName, int &numberOfMeshNodes, int &numberOfElements,
        double *&meshNodeCoordinates, int *&meshNodeIds, int *&numberOfNodesPerElement,
        int *&elementNodeTables, int *&elementIds);
/***********************************************************************************************
 * \brief Write a .msh file
 * \param[in] fileName name of the mesh file
 * \param[in] numberOfMeshNodes number of nodes belonging to elements
 * \param[in] numberOfElements number of elements in the mesh
 * \param[in] meshNodeCoordinates coordinates of mesh nodes (x,y,z)
 * \param[in] MeshNodeIds index/id of each mesh node
 * \param[in] numberOfNodesPerElement number of nodes in each element
 * \param[in] elementNodeTables "element tables" or "connectivity tables"
 * \param[in] elementIds index/id of each element
 * \author Michael Andre, Tianyang Wang
 ***********/
void writeDotMsh(std::string fileName, int numberOfMeshNodes, int numberOfElements,
        const double *meshNodeCoordinates, const int *meshNodeIds,
        const int *numberOfNodesPerElement, const int *elementNodeTables, const int *elementIds);
/***********************************************************************************************
 * \brief Initialize a .res file
 * \param[in] fileName name of the result file
 * \param[in] numNodesPerElem number of nodes per element, used to output Gauss points information
 * \author Tianyang Wang
 ***********/
void initDotRes(std::string fileName);
/***********************************************************************************************
 * \brief Append data of a certain time step to a .res file
 * \param[in] fileName name of the result file
 * \param[in] resultName name of the data
 * \param[in] analysisName name of the analysis
 * \param[in] stepNum step number
 * \param[in] type "Scalar" or "Vector"
 * \param[in] numberOfNodes number of nodes
 * \param[in] nodeIds Ids of nodes
 * \param[in] data data of this step
 * \author Tianyang Wang
 ***********/
void appendNodalDataToDotRes(std::string fileName, std::string resultName, std::string analysisName,
        int stepNum, std::string type, int numberOfNodes, const int *nodeIds, const double *data);
/***********************************************************************************************
 * \brief Append data of a certain time step to a .res file
 * \param[in] fileName name of the result file
 * \param[in] resultName name of the data
 * \param[in] analysisName name of the analysis
 * \param[in] stepNum step number
 * \param[in] type "Scalar" or "Vector"
 * \param[in] numberOfElements number of elements
 * \param[in] elemIds Ids of elements
 * \param[in] numberOfNodesPerElement number of nodes in each element
 * \param[in] data data of this step
 * \author Tianyang Wang
 ***********/
void appendElementalDataToDotRes(std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type, int numberOfElements,
        const int *elemIds, const int *numberOfNodesPerElement, const double *data);
/***********************************************************************************************
 * \brief Read nodal data in a .res file given the parameters of the data
 * \param[in] fileName name of the result file
 * \param[in] resultName name of the data
 * \param[in] analysisName name of the analysis
 * \param[in] stepNum step number
 * \param[in] type "Scalar" or "Vector"
 * \param[in] numberOfNodes number of nodes
 * \param[in] nodeIds Ids of nodes
 * \param[out] data data of this step
 * \author Tianyang Wang
 ***********/
void readNodalDataFromDotRes(std::string fileName, std::string resultName, std::string analysisName,
        int stepNum, std::string type, int numberOfNodes, const int *nodeIds, double *data);
/***********************************************************************************************
 * \brief Read nodal data in a .res file in a fast way by using the existing file stream
 * \param[in] dotResFile stream of the result file
 * \param[in] fileName name of the result file
 * \param[in] resultName name of the data
 * \param[in] analysisName name of the analysis
 * \param[in] stepNum step number
 * \param[in] type "Scalar" or "Vector"
 * \param[in] numberOfNodes number of nodes
 * \param[in] nodeIds Ids of nodes
 * \param[out] data data of this step
 * \author Tianyang Wang
 ***********/
void readNodalDataFromDotResFast(std::ifstream &dotResFile, std::string fileName, std::string resultName, std::string analysisName,
        int stepNum, std::string type, int numberOfNodes, const int *nodeIds, double *data);
/***********************************************************************************************
 * \brief Read elemental data in a .res file given the parameters of the data
 * \param[in] fileName name of the result file
 * \param[in] resultName name of the data
 * \param[in] analysisName name of the analysis
 * \param[in] stepNum step number
 * \param[in] type "Scalar" or "Vector"
 * \param[in] numberOfElements number of elements
 * \param[in] elemIds Ids of elements
 * \param[in] numberOfNodesPerElement number of nodes in each element
 * \param[out] data data of this step
 * \author Tianyang Wang
 ***********/
void readElementalDataFromDotRes(std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type, int numberOfElements,
        const int *elemIds, const int *numberOfNodesPerElement, double *data);
/***********************************************************************************************
 * \brief Read elemental data in a .res in a fast way by using the existing file stream
 * \param[in] dotResFile stream of the result file
 * \param[in] fileName name of the result file
 * \param[in] resultName name of the data
 * \param[in] analysisName name of the analysis
 * \param[in] stepNum step number
 * \param[in] type "Scalar" or "Vector"
 * \param[in] numberOfElements number of elements
 * \param[in] elemIds Ids of elements
 * \param[in] numberOfNodesPerElement number of nodes in each element
 * \param[out] data data of this step
 * \author Tianyang Wang
 ***********/
void readElementalDataFromDotResFast(std::ifstream &dotResFile, std::string fileName, std::string resultName,
        std::string analysisName, int stepNum, std::string type, int numberOfElements,
        const int *elemIds, const int *numberOfNodesPerElement, double *data);

} /* namespace GiDFileIO */

#endif /* GIDFILEIO_H_ */
