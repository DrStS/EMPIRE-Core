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
 * \file NearestElementMapper.h
 * This file holds the class NearestElementMapper
 * \date 8/14/2013
 **************************************************************************************************/
#ifndef NEARESTELEMENTMAPPER_H_
#define NEARESTELEMENTMAPPER_H_

#include <vector>
#include "AbstractMapper.h"

namespace EMPIRE {
/********//**
 * \brief Class NearestElementMapper performs nearest element mapping
 ***********/
class NearestElementMapper: public AbstractMapper {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _numNodesA number of nodes (A)
     * \param[in] _numElemsA number of elements (A)
     * \param[in] _numNodesPerElemA number of nodes for each element (only 3 or 4, enable hybrid mesh) (A)
     * \param[in] _nodesA x,y,z coordinates of all slave nodes (A)
     * \param[in] _nodeIDsA the id of each slave node (A)
     * \param[in] _elemTableA how the elements are constructed by the nodes (A)
     *
     * \param[in] _numNodesB number of nodes (B)
     * \param[in] _numElemsB number of elements (B)
     * \param[in] _numNodesPerElemB number of nodes for each element (only 3 or 4, enable hybrid mesh) (B)
     * \param[in] _nodesB x,y,z coordinates of all slave nodes (B)
     * \param[in] _nodeIDsB the id of each slave node (B)
     * \param[in] _elemTableB how the elements are constructed by the nodes (B)
     *
     * \author Tianyang Wang
     ***********/
    NearestElementMapper(int _numNodesA, int _numElemsA, const int *_numNodesPerElemA,
            const double *_nodesA, const int *_nodeIDsA, const int *_elemTableA, int _numNodesB,
            int _numElemsB, const int *_numNodesPerElemB, const double *_nodesB,
            const int *_nodeIDsB, const int *_elemTableB);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~NearestElementMapper();
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
    /// defines number of threads used for mapper routines
    static int mapperSetNumThreads;

private:
    /// number of nodes (A)
    int numNodesA;
    ///  number of elements (A)
    int numElemsA;
    ///  number of nodes for each element (only 3 or 4, enable hybrid mesh) (A)
    const int *numNodesPerElemA;
    ///  x,y,z coordinates of all slave nodes (A)
    const double *nodesA;
    /// the id of each slave node (A)
    const int *nodeIDsA;
    /// how the elements are constructed by the nodes (A)
    const int *elemTableA;
    /// number of nodes (B)
    int numNodesB;
    /// number of elements (B)
    int numElemsB;
    /// number of nodes for each element (only 3 or 4, enable hybrid mesh) (B)
    const int *numNodesPerElemB;
    /// x,y,z coordinates of all slave nodes (B)
    const double *nodesB;
    /// the id of each slave node (B)
    const int *nodeIDsB;
    /// how the elements are constructed by the nodes (B)
    const int *elemTableB;

    /// number of nodes per neighbor element
    int *numNodesPerNeighborElem;
    /// neighbors of nodes in B
    std::vector<int*> *neighborsTable;
    /// weights of the neighbors
    std::vector<double*> *weightsTable;
    /// directElemTable means the entries is not the node number, but the position in nodeCoors
    std::vector<int> **directElemTableA;
    /// number of neighbors to search
    static const int MAX_NUM_NEIGHBORS_TO_SEARCH;
    /***********************************************************************************************
     * \brief Compute the neighbors and the weights
     * \author Tianyang Wang
     ***********/
    void computeNeighborsAndWeights();
    /***********************************************************************************************
     * \brief Given the element index/id, return the element
     * \param[in] elemIndex the element index/id
     * \param[out] elem the element corresponds to the element index
     * \author Tianyang Wang
     ***********/
    void getElemCoorInA(int elemIndex, double *elem);
    /***********************************************************************************************
     * \brief Determine whether a node is inside the element or not
     * \param[in] numNodesThisElem number of nodes of this element (3 or 4)
     * \param[out] localCoor local coordinates of this node
     * \author Tianyang Wang
     ***********/
    bool insideElement(int numNodesThisElem, double *localCoor);
};

} /* namespace EMPIRE */
#endif /* NEARESTELEMENTMAPPER_H_ */
