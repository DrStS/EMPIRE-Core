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
 * \file MapperAdapter.h
 * This file holds the class MapperAdapter
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef MAPPERADAPTER_H_
#define MAPPERADAPTER_H_

#include <string>
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class MortarMapper;
class AbstractMesh;
class DataField;
class AbstractMapper;

/********//**
 * \brief Class MapperAdapter is the adaptor of the mapper.
 ***********/
class MapperAdapter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name name of the mapper
     * \param[in] _meshA meshA
     * \param[in] _meshB meshB
     * \author Tianyang Wang
     ***********/
    MapperAdapter(std::string _name, AbstractMesh *_meshA, AbstractMesh *_meshB);
    /***********************************************************************************************
     * \brief Initialize MortarMapper
     * \param[in] oppositeSurfaceNormal whether the surface normal of A and B are opposite
     * \param[in] dual use dual mortar or not
     * \param[in] enforceConsistency enforce consistensy in consistent mapping
     * \author Tianyang Wang
     ***********/
    void initMortarMapper(bool oppositeSurfaceNormal, bool dual, bool enforceConsistency);
    /***********************************************************************************************
     * \brief Initialize IGA Mortar Mapper
     * \param[in] _tolProjectionDistance tolerance of distance for projection
     * \param[in] _numGPsTriangle number of Gauss points when performs integration on triangle
     * \param[in] _numGPsQuad number of Gauss points when performs integration on quadrilateral
     * \author Chenshen Wu
     ***********/
    void initIGAMortarMapper(double tolProjectionDistance, int numGPsTriangle, int numGPsQuad);
    /***********************************************************************************************
     * \brief Initialize NearestNeighborMapper
     * \author Tianyang Wang
     ***********/
    void initNearestNeighborMapper();
    /***********************************************************************************************
     * \brief Initialize BarycentricInterpolationMapper
     * \author Tianyang Wang
     ***********/
    void initBarycentricInterpolationMapper();
    /***********************************************************************************************
     * \brief Initialize NearestElementMapper
     * \author Tianyang Wang
     ***********/
    void initNearestElementMapper();
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~MapperAdapter();
    /***********************************************************************************************
     * \brief Do consistent mapping from A to B (map displacements)
     * \param[in] fieldA is the input data
     * \param[out] fieldB is the output data
     * \author Tianyang Wang
     ***********/
    void consistentMapping(const DataField *fieldA, DataField *fieldB);
    /***********************************************************************************************
     * \brief Do conservative mapping from B to A (map forces)
     * \param[in] fieldB is the input data
     * \param[out] fieldA is the output data
     * \author Tianyang Wang
     ***********/
    void conservativeMapping(const DataField *fieldB, DataField *fieldA);
    /***********************************************************************************************
     * \brief is it meshA or not
     * \param[in] mesh mesh
     * \return ture if it is meshA
     * \author Tianyang Wang
     ***********/
    bool isMeshA(AbstractMesh *mesh) {
        return meshA == mesh;
    }
    /***********************************************************************************************
     * \brief is it meshB or not
     * \param[in] mesh mesh
     * \return ture if it is meshB
     * \author Tianyang Wang
     ***********/
    bool isMeshB(AbstractMesh *mesh) {
        return meshB == mesh;
    }
private:
    /// the adapted mapper
    AbstractMapper *mapperImpl;
    /// name of the mapper
    std::string name;
    /// mesh A
    AbstractMesh *meshA;
    /// mesh B
    AbstractMesh *meshB;
};

} /* namespace EMPIRE */
#endif /* MAPPERADAPTER_H_ */
