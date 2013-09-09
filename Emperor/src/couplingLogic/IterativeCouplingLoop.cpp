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

#include "IterativeCouplingLoop.h"
#include "ConvergenceChecker.h"
#include "DataField.h"
#include "ClientCode.h"
#include "AbstractCouplingAlgorithm.h"
#include "Message.h"
#include "DataOutput.h"

using namespace std;

namespace EMPIRE {

IterativeCouplingLoop::IterativeCouplingLoop() :
        AbstractCouplingLogic(), convergenceObserverVec(), couplingAlgorithmVec(), convergenceChecker(
                NULL) {
    outputCounter = 0;
}

IterativeCouplingLoop::~IterativeCouplingLoop() {
    delete convergenceChecker;
}

void IterativeCouplingLoop::initConvergenceChecker(ConvergenceChecker *_convergenceChecker) {
    assert(_convergenceChecker!=NULL);
    assert(convergenceChecker==NULL);
    convergenceChecker = _convergenceChecker;
}

void IterativeCouplingLoop::doCoupling() {
    assert(convergenceChecker!=NULL);
    int count = 0;
    // notify the coupling algorithms the start of the iteration
    for (int i = 0; i < couplingAlgorithmVec.size(); i++) {
        couplingAlgorithmVec[i]->setNewLoop();
    }

    // initialize output files, in each iterative coupling the previous files are overwritten
    outputCounter++;
    stringstream rearPart;
    rearPart << outputCounter;
    for (int i = 0; i < dataOutputVec.size(); i++)
        dataOutputVec[i]->init(rearPart.str());

    while (true) {
        count++;
        stringstream ss;
        ss << "iteration step: " << count;
        HEADING_OUT(4, "IterativeCouplingLoop", ss.str(), infoOut);

        // do coupling
        for (int i = 0; i < couplingLogicSequence.size(); i++)
            couplingLogicSequence[i]->doCoupling();

        // write data field at this iteration
        for (int i = 0; i < dataOutputVec.size(); i++)
            dataOutputVec[i]->writeCurrentStep(count);

        // broadcast convergence
        if (convergenceChecker->isConvergent()) {
            broadcastConvergenceToClients(true);
            break;
        } else {
            broadcastConvergenceToClients(false);
        }
        assert(count == convergenceChecker->getcurrentNumOfIterations());
    }
    //std::cout << "number of iterative coupling loops: " << count << std::endl;
}

void IterativeCouplingLoop::addConvergenceObserver(ClientCode *clientCode) {
    convergenceObserverVec.push_back(clientCode);
}

void IterativeCouplingLoop::broadcastConvergenceToClients(bool convergent) {
    for (int i = 0; i < convergenceObserverVec.size(); i++) {
        convergenceObserverVec[i]->sendConvergenceSignal(convergent);
    }
}

void IterativeCouplingLoop::addCouplingAlgorithm(AbstractCouplingAlgorithm *couplingAlgorithm) {
    couplingAlgorithmVec.push_back(couplingAlgorithm);
}

void IterativeCouplingLoop::addDataOutput(DataOutput *dataOutput) {
    dataOutputVec.push_back(dataOutput);
}

} /* namespace EMPIRE */
