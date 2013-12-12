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
#ifndef SETFILTER_H_
#define SETFILTER_H_

#include "AbstractFilter.h"
#include <assert.h>
#include <stdlib.h>

namespace EMPIRE {
/********//**
 * \brief Class SetFilter sets signals and fields to a constant value
 * \author Stefan Sicklinger
 ***********/
class SetFilter : public AbstractFilter {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _value to set the signal/datafield to
     * \author Stefan Sicklinger
     ***********/
    SetFilter(std::vector<double> _value);
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    virtual ~SetFilter();
    /***********************************************************************************************
     * \brief Filtering
     * \author Stefan Sicklinger
     ***********/
    void filtering();
    /***********************************************************************************************
     * \brief Initialize data according to the inputs and outputs
     * \author Stefan Sicklinger
     ***********/
    void init();
private:
    /// values for all entries
    std::vector<double> value;
};

} /* namespace EMPIRE */
#endif /* SETFILTER_H_ */
