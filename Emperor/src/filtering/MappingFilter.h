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
 * \file MappingFilter.h
 * This file holds the class MappingFilter
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef MAPPINGFILTER_H_
#define MAPPINGFILTER_H_

#include "AbstractFilter.h"

namespace EMPIRE {

class MapperAdapter;
class DataFieldIntegrationFilter;
class DataField;

/********//**
 * \brief Class MappingFilter filters the data by calling a mapper
 * \author Tianyang Wang
 ***********/
class MappingFilter: public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _mapper pointer to the mapper
     * \param[in] _function use the conservativeMapping or the consistentMapping
     * \author Tianyang Wang
     ***********/
    MappingFilter(MapperAdapter *_mapper);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~MappingFilter();
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
    /// The mapper
    MapperAdapter *mapper;
    /// true if do consistent mapping, false if do conservative mapping
    bool consistentMapping;
    /// modify the input datafield
    DataFieldIntegrationFilter *inputModifier;
    /// modify the output datafield
    DataFieldIntegrationFilter *outputModifier;
    /// temporary input for mapper
    DataField *tmpInput;
    /// temporary output for mapper
    DataField *tmpOutput;
};

} /* namespace EMPIRE */
#endif /* MAPPINGFILTER_H_ */
