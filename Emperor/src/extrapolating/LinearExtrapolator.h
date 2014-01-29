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
 * \file LinearExtrapolator.h
 * This file holds the class LinearExtrapolator
 * \date 1/17/2014
 **************************************************************************************************/
#ifndef LINEAREXTRAPOLATOR_H_
#define LINEAREXTRAPOLATOR_H_

#include "AbstractExtrapolator.h"
#include <vector>

namespace EMPIRE {
/********//**
 * \brief Class LinearExtrapolator does linear extrapolation at the beginning of a time step, as a guess
 ***********/
class LinearExtrapolator : public AbstractExtrapolator {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name
     * \author Tianyang Wang
     ***********/
    LinearExtrapolator(std::string _name);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~LinearExtrapolator();
    /***********************************************************************************************
     * \brief initialize the extrapolator after all the connectionIOs are added
     * \author Tianyang Wang
     ***********/
    void init();
    /***********************************************************************************************
     * \brief Do extrapolation
     * \author Tianyang Wang
     ***********/
    void extrapolate();
private:
    // data of the previous' previous time step
    std::vector<double*> data00;
    // data of the previous time step
    std::vector<double*> data0;
};

} /* namespace EMPIRE */
#endif /* LINEAREXTRAPOLATOR_H_ */
