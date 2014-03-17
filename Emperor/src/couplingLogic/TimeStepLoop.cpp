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
#include <assert.h>
#include <iostream>
#include <sstream>

#include "TimeStepLoop.h"
#include "Connection.h"
#include "AbstractExtrapolator.h"
#include "Message.h"
#include "DataOutput.h"

namespace EMPIRE {

using namespace std;

TimeStepLoop::TimeStepLoop(int _numTimeSteps) :
        AbstractCouplingLogic(), numTimeSteps(_numTimeSteps), extrapolator(NULL) {
    assert(numTimeSteps > 0);
    outputCounter = 0;
}

TimeStepLoop::~TimeStepLoop() {
}

void TimeStepLoop::doCoupling() {
    // initialize output files
    outputCounter++;
    stringstream rearPart;
    rearPart << "_" << outputCounter;
    for (int i = 0; i < dataOutputVec.size(); i++)
        dataOutputVec[i]->init(rearPart.str());

    for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {
        // output to shell
        stringstream ss;
        ss << "time step: " << timeStep;
        HEADING_OUT(3, "TimeStepLoop", ss.str(), infoOut);

        // set extrapolation
        if (extrapolator != NULL)
            extrapolator->extrapolate();

        // do the coupling
        for (int i = 0; i < couplingLogicSequence.size(); i++)
            couplingLogicSequence[i]->doCoupling();

        // write data field at current time step
        for (int i = 0; i < dataOutputVec.size(); i++)
            dataOutputVec[i]->writeCurrentStep(timeStep);

    }
}

void TimeStepLoop::setExtrapolator(AbstractExtrapolator *_extrapolator) {
    extrapolator = _extrapolator;
}

} /* namespace EMPIRE */
