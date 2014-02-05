#include "OptimizationLoop.h"
#include "ClientCode.h"
#include "Message.h"
#include "DataOutput.h"

#include <sstream>
#include <assert.h>

using namespace std;

namespace EMPIRE {

OptimizationLoop::OptimizationLoop(int _maxNumOfIterations) :
        AbstractCouplingLogic(), maxNumOfIterations(_maxNumOfIterations), convergenceSignalSender(
                NULL) {
    currentNumOfIterations = 0;
}

OptimizationLoop::~OptimizationLoop() {
}

void OptimizationLoop::doCoupling() {
    // initialize output files
    for (int i = 0; i < dataOutputVec.size(); i++)
        dataOutputVec[i]->init("");

    while (true) {
        currentNumOfIterations++;
        if (currentNumOfIterations == maxNumOfIterations) {
            ERROR_OUT("OptimizationLoop::doCoupling(): reach maximum number of iterations!");
            assert(false);
        }
        stringstream ss;
        ss << "optimization step: " << currentNumOfIterations;
        HEADING_OUT(2, "OptimizationLoop", ss.str(), infoOut);

        // do coupling
        for (int i = 0; i < couplingLogicSequence.size(); i++)
            couplingLogicSequence[i]->doCoupling();

        // write data field at this iteration
        for (int i = 0; i < dataOutputVec.size(); i++)
            dataOutputVec[i]->writeCurrentStep(currentNumOfIterations);

        // receive/send convergence signal
        assert(convergenceSignalSender != NULL);
        bool convergent = convergenceSignalSender->recvConvergenceSignal();
        for (int i = 0; i < convergenceSignalReceivers.size(); i++) {
            convergenceSignalReceivers[i]->sendConvergenceSignal(convergent);
        }

        if (convergent)
            break;
    }
}

void OptimizationLoop::setConvergenceSignalSender(ClientCode *_convergenceSignalSender) {
    assert(convergenceSignalSender == NULL);
    convergenceSignalSender = _convergenceSignalSender;
}

void OptimizationLoop::addConvergenceSignalReceiver(ClientCode *convergenceSignalReceiver) {
    convergenceSignalReceivers.push_back(convergenceSignalReceiver);
}

} /* namespace EMPIRE */
