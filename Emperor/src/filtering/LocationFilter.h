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
 * \file LocationFilter.h
 * This file holds the class LocationFilter
 * \date 9/3/2012
 **************************************************************************************************/
#ifndef LOCATIONFILTER_H_
#define LOCATIONFILTER_H_

#include <string>
#include <vector>
#include "AbstractFilter.h"

namespace EMPIRE {

class AbstractMesh;
class FEMesh;
/********//**
 * \brief Class LocationFilter filters the data between element centroids and nodes
 ***********/
class LocationFilter: public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Tianyang Wang
     ***********/
    LocationFilter();
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~LocationFilter();
    /***********************************************************************************************
     * \brief Filtering
     * \author Tianyang Wang
     ***********/
    void filtering();
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Tianyang Wang
     ***********/
    void init();

private:
    /// the mesh
    AbstractMesh *mesh;
    /// cast mesh to feMesh
    FEMesh *feMesh;
    /// table that links a node position (instead of node ID) to all elements containing it
    std::vector<int> **nodePosToElemTable;
    /// case number --- 1. field from element centroids to nodes; 2. to be implemented
    int caseNum;
    /***********************************************************************************************
     * \brief Filter the data field in case 1
     * \author Tianyang Wang
     ***********/
    void filterDataFieldCase1();
    /***********************************************************************************************
     * \brief Compute the member nodePosToElemTable
     * \author Tianyang Wang
     ***********/
    void computeNodePosToElemTable();
};

} /* namespace EMPIRE */
#endif /* LOCATIONFILTER_H_ */
