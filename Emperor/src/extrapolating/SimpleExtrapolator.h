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
 * \file SimpleExtrapolator.h
 * This file holds the class SimpleExtrapolator
 * \date 3/5/2012
 **************************************************************************************************/
#ifndef SIMPLEEXTRAPOLATOR_H_
#define SIMPLEEXTRAPOLATOR_H_

#include "AbstractExtrapolator.h"

namespace EMPIRE {
/********//**
 * \brief Class SimpleExtrapolator does nothing
 ***********/
class SimpleExtrapolator: public AbstractExtrapolator {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name
     * \author Tianyang Wang
     ***********/
    SimpleExtrapolator(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~SimpleExtrapolator();
    /***********************************************************************************************
     * \brief Do extrapolation
     * \author Tianyang Wang
     ***********/
    void extrapolate();
};

} /* namespace EMPIRE */
#endif /* SIMPLEEXTRAPOLATOR_H_ */
