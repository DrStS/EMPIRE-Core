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
 * \file CopyFilter.h
 * This file holds the class CopyFilter
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef COPYFILTER_H_
#define COPYFILTER_H_

#include "AbstractFilter.h"

namespace EMPIRE {

/********//**
 * \brief Class CopyFilter copy the input to the output, if not enough input data, fill by 0
 ***********/
class CopyFilter: public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _factor for scaling
     * \author Tianyang Wang, Stefan Sicklinger
     ***********/
    CopyFilter(int _signalOffset);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~CopyFilter();
    /***********************************************************************************************
     * \brief Copy the input to the output, if not enough input data, fill by 0
     * \author Tianyang Wang
     ***********/
    void filtering();
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Tianyang Wang
     ***********/
    void init();
private:
/// Signal offset value 0--> no offset
int signalOffset;
};


} /* namespace EMPIRE */
#endif /* COPYFILTER_H_ */
