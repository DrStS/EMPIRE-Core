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
 * \file BarycentricInterpolationMapper.h
 * This file holds the class BarycentricInterpolationMapper
 * \date 5/2/2013
 **************************************************************************************************/

#ifndef BARYCENTRICINTERPOLATIONMAPPER_H_
#define BARYCENTRICINTERPOLATIONMAPPER_H_

#include <string>
#include "AbstractMapper.h"

namespace EMPIRE {
/********//**
 * \brief Class BarycentricInterpolationMapper performs barycentric interpolation mapping
 ***********/
class BarycentricInterpolationMapper: public AbstractMapper {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numNodesA number of nodes of A
     * \param[in] _nodesA nodes of A
     * \param[in] _numNodesB number of nodes of B
     * \param[in] _nodesB nodes of B
     * \author Tianyang Wang
     ***********/
    BarycentricInterpolationMapper(int _numNodesA, const double *_nodesA, int _numNodesB, const double *_nodesB);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~BarycentricInterpolationMapper();
    /***********************************************************************************************
     * \brief Do consistent mapping on fields (e.g. displacements or tractions)
     * \param[in] fieldA the field of mesh A (e.g. x-displacements on all structure nodes)
     * \param[out] fieldB the field of mesh B (e.g. x-displacements on all fluid nodes)
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const double *fieldA, double *fieldB);
    /***********************************************************************************************
     * \brief Do conservative mapping on integrated fields (e.g. forces)
     * \param[in] fieldB the field of mesh B (e.g. x-forces on all fluid nodes)
     * \param[out] fieldA the field of mesh A (e.g. x-forces on all structure nodes)
     * \author Tianyang Wang
     * ***********/
    void conservativeMapping(const double *fieldB, double *fieldA);
private:
    /// number of nodes of A
    int numNodesA;
    /// nodes of A
    const double *nodesA;
    /// number of nodes of B
    int numNodesB;
    /// nodes of B
    const double *nodesB;
    /// neighbors of nodes in B
    int *neighborsTable;
    /// weights of the neighbors
    double *weightsTable;
    /// number of neighbors to search (more than 3 is needed, because sometimes 3 nodes are on the same line)
    static const int MAX_NUM_NEIGHBORS_TO_SEARCH;
    /***********************************************************************************************
     * \brief Compute the neighbors
     * \author Tianyang Wang
     ***********/
    void computeNeighbors();
    /***********************************************************************************************
     * \brief Compute the weights
     * \author Tianyang Wang
     ***********/
    void computeWeights();
    /***********************************************************************************************
     * \brief whether the three nodes are on the same line or not
     * \return true if not on the same line, otherwise, false
     * \author Tianyang Wang
     ***********/
    static bool areNotOnTheSameLine(const double *p1, const double *p2, const double *p3);
};

} /* namespace EMPIRE */
#endif /* BARYCENTRICINTERPOLATIONMAPPER_H_ */
