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
 * \file AbstractFilter.h
 * This file holds the class AbstractFilter
 * \date 3/5/2012
 **************************************************************************************************/

#ifndef ABSTRACTFILTER_H_
#define ABSTRACTFILTER_H_

#include <assert.h>
#include <stdlib.h>
#include <vector>
#include "EMPEROR_Enum.h"

namespace EMPIRE {
class ConnectionIO;
/********//**
 * \brief Class AbstractFilter is the superclass of all data field filters
 ***********/
class AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Tianyang Wang
     ***********/
    AbstractFilter();
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~AbstractFilter();
    /***********************************************************************************************
     * \brief Filtering
     * \author Tianyang Wang
     ***********/
    virtual void filtering() = 0;
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Tianyang Wang
     ***********/
    virtual void init() = 0;

    void addInput(ConnectionIO *input);

    void addOutput(ConnectionIO *output);

protected:
    /// inputs
    std::vector<ConnectionIO*> inputVec;
    /// outputs
    std::vector<ConnectionIO*> outputVec;

private:
    AbstractFilter(const AbstractFilter&);
};

} /* namespace EMPIRE */
#endif /* ABSTRACTFILTER_H_ */
