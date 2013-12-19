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
 * \file ConstantRelaxation.h
 * This file holds the class ConstantRelaxation
 * \date 12/19/2013
 **************************************************************************************************/
#ifndef CONSTANTRELAXATION_H_
#define CONSTANTRELAXATION_H_

#include "AbstractCouplingAlgorithm.h"

namespace EMPIRE {
/********//**
 * \brief Class ConstantRelaxation does a constant relaxation
 ***********/
class ConstantRelaxation: public AbstractCouplingAlgorithm {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \param[in] _name the name of the coupling algorithm
     * \param[in] _relaxationFactor the constant relaxation factor
     * \author Tianyang Wang
     ***********/
    ConstantRelaxation(std::string _name, double _relaxationFactor);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~ConstantRelaxation();
    /***********************************************************************************************
     * \brief Calculate the new value of the output
     * \author Tianyang Wang
     ***********/
    void calcNewValue();
private:
    /// relaxation factor
    const double RELAXATION_FACTOR;
    /// friend class in unit test
    friend class TestConstantRelaxation;
    /// whether output numbers or not
    bool debugMe;
};

} /* namespace EMPIRE */
#endif /* CONSTANTRELAXATION_H_ */
