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
 * \file CouplingAlgorithmFilter.h
 * This file holds the class CouplingAlgorithmFilter
 * \date 3/5/2012
 **************************************************************************************************/
#ifndef EXTRAPOLATINGFILTER_H_
#define EXTRAPOLATINGFILTER_H_

#include "AbstractFilter.h"

namespace EMPIRE {

class AbstractExtrapolator;
/********//**
 * \brief Class ExtrapolatingFilter filters the data by calling an extrapolator
 ***********/
class ExtrapolatingFilter : public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _extrapolator the extrapolator which the filter calls
     * \author Tianyang Wang
     ***********/
    ExtrapolatingFilter(AbstractExtrapolator *_extrapolator);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~ExtrapolatingFilter();
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
    /// the extrapolator
    AbstractExtrapolator *extrapolator;
};

} /* namespace EMPIRE */
#endif /* EXTRAPOLATINGFILTER_H_ */
