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
 * \file AdditionFilter.h
 * This file holds the class AdditionFilter
 * \date 3/11/2014
 **************************************************************************************************/
#ifndef ADDITIONFILTER_H_
#define ADDITIONFILTER_H_

#include "AbstractFilter.h"

namespace EMPIRE {
/********//**
 * \brief Class AdditionFilter z = a*x + b*y
 ***********/
class AdditionFilter : public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _a a
     * \param[in] _b b
     * \author Tianyang Wang
     ***********/
    AdditionFilter(double _a, double _b);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~AdditionFilter();
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
    /// a
    const double a;
    /// b
    const double b;
};

} /* namespace EMPIRE */
#endif /* ADDITIONFILTER_H_ */
