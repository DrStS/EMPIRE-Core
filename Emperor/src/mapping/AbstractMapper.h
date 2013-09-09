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
 * \file AbstractMapper.h
 * This file holds the class AbstractMapper
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef ABSTRACTMAPPER_H_
#define ABSTRACTMAPPER_H_

#include <string>

namespace EMPIRE {
class DataField;
class AbstractMesh;

/********//**
 * \brief Class AbstractMapper is the superclass of all mappers
 ***********/
class AbstractMapper {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name name of the mapper
     * \author Tianyang Wang
     ***********/
    AbstractMapper() {
    }
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~AbstractMapper() {
    }
    /***********************************************************************************************
     * \brief Do consistent mapping on fields (e.g. displacements or tractions)
     * \param[in] fieldA the field of mesh A (e.g. x-displacements on all structure nodes)
     * \param[out] fieldB the field of mesh B (e.g. x-displacements on all fluid nodes)
     * \author Tianyang Wang
     ***********/
    virtual void consistentMapping(const double *fieldA, double *fieldB) = 0;
    /***********************************************************************************************
     * \brief Do conservative mapping on integrated fields (e.g. forces)
     * \param[in] fieldB the field of mesh B (e.g. x-forces on all fluid nodes)
     * \param[out] fieldA the field of mesh A (e.g. x-forces on all structure nodes)
     * \author Tianyang Wang
     * ***********/
    virtual void conservativeMapping(const double *fieldB, double *fieldA) = 0;
};

} /* namespace EMPIRE */
#endif /* ABSTRACTMAPPER_H_ */
