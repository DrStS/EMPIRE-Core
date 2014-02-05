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
 * \file TimeStepLoop.h
 * This file holds the class TimeStepLoop
 * \date 5/15/2012
 **************************************************************************************************/
#ifndef TIMESTEPLOOP_H_
#define TIMESTEPLOOP_H_

#include <vector>

#include "AbstractCouplingLogic.h"
#include "EMPEROR_Enum.h"

namespace EMPIRE {

class AbstractExtrapolator;

/********//**
 * \brief Class TimeStepLoop performs time step loop on the sequence of coupling logics
 * \author Tianyang Wang
 ***********/
class TimeStepLoop: public AbstractCouplingLogic {
public:
    /***********************************************************************************************
     * \brief Constructor
     * \author Tianyang Wang
     ***********/
    TimeStepLoop(int _numTimeSteps);
    /***********************************************************************************************
     * \brief Destructor
     * \author Tianyang Wang
     ***********/
    virtual ~TimeStepLoop();
    /***********************************************************************************************
     * \brief Do time step coupling
     * \author Tianyang Wang
     ***********/
    void doCoupling();
    /***********************************************************************************************
     * \brief Set the extrapolator which will do extrapolation at the beginning of the time step
     * \param[in] _extrapolator the extrapolator
     * \author Tianyang Wang
     ***********/
    void setExtrapolator(AbstractExtrapolator *_extrapolator);

private:
    /// number of time steps
    int numTimeSteps;
    /// Extrapolator
    AbstractExtrapolator *extrapolator;
    /// output counter
    int outputCounter;
    /// the unit test classes
    friend class TestLoops;
    friend class TestEmperor;
};

} /* namespace EMPIRE */
#endif /* TIMESTEPLOOP_H_ */
