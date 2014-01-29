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
#include <math.h>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <sstream>
#include <fstream>

#include "ConvergenceChecker.h"
#include "DataField.h"
#include "Message.h"
#include "AbstractCouplingAlgorithm.h"
#include "Signal.h"

using namespace std;

namespace EMPIRE {

ConvergenceChecker::CheckResidual::CheckResidual(double _absoluteTolerance,
        double _relativeTolerance, AbstractCouplingAlgorithm *_couplingAlgorithm,
        int _residualIndex) :
        ABS_TOL(_absoluteTolerance), REL_TOL(_relativeTolerance), couplingAlgorithm(
                _couplingAlgorithm), residualIndex(_residualIndex) {
}

ConvergenceChecker::CheckResidual::~CheckResidual() {

}

void ConvergenceChecker::CheckResidual::updateInitialResidual() {
    initialResidual = couplingAlgorithm->getResidualL2Norm(residualIndex);
}

double ConvergenceChecker::CheckResidual::getAbsoluteResidual() {
    return couplingAlgorithm->getResidualL2Norm(residualIndex);
}

double ConvergenceChecker::CheckResidual::getRelativeResidual() {
    if (couplingAlgorithm->getResidualL2Norm(residualIndex) >= initialResidual * 1E10)
        return 0.0;
    return couplingAlgorithm->getResidualL2Norm(residualIndex) / initialResidual;
}

bool ConvergenceChecker::CheckResidual::isConvergent() {
    double absoluteResidual = getAbsoluteResidual();
    double relativeResidual = getRelativeResidual();
    if (absoluteResidual < ABS_TOL || relativeResidual < REL_TOL) {
        return true;
    }
    return false;
}

void ConvergenceChecker::CheckResidual::writeResidualToShell() {
    stringstream ss;
    ss << scientific;
    ss << "ConvergenceChecker::CheckResidual(" << couplingAlgorithm->getName() << ", " << residualIndex
            << "): " << "(" << getRelativeResidual() << ", " << getAbsoluteResidual() << ")"
            << endl;
    INDENT_OUT(1, ss.str(), infoOut);
}

ConvergenceChecker::ConvergenceChecker(double maxNumOfIters) :
        MAX_NUM_ITERATIONS(maxNumOfIters) {
    currentNumOfIterations = 0;
    debugResidual = true;
    timeStepNumber = 1;

    residualFileName = "convergenceChecker";
    residualFileName.append(".log");
    fstream residualFile;
    residualFile.open(residualFileName.c_str(), ios_base::out);
    residualFile << scientific;
    residualFile << "timeStep" << '\t' << "iteration" << '\t' << "relativeResidual" << '\t'
            << "absoluteResidual" << endl;
    assert(!residualFile.fail());
}

ConvergenceChecker::~ConvergenceChecker() {
    for (int i = 0; i < checkResiduals.size(); i++) {
        delete checkResiduals[i];
    }
}

bool ConvergenceChecker::isConvergent() {
    currentNumOfIterations++;
    fstream residualFile;
    residualFile.open(residualFileName.c_str(), ios_base::out | ios_base::app);

    // update the initial residual
    if (currentNumOfIterations == 1) {
        for (int i = 0; i < checkResiduals.size(); i++) {
            checkResiduals[i]->updateInitialResidual();
        }
    }

    // output residuals to shell
    if (debugResidual) {
        for (int i = 0; i < checkResiduals.size(); i++) {
            checkResiduals[i]->writeResidualToShell();
        }
    }

    // write residuals to file
    residualFile << timeStepNumber << '\t' << currentNumOfIterations;
    residualFile << scientific;
    for (int i = 0; i < checkResiduals.size(); i++) {
        residualFile << '\t' << "(" << checkResiduals[i]->getRelativeResidual() << '\t'
                << checkResiduals[i]->getAbsoluteResidual() << ")";
    }
    residualFile << endl;

    // 3. set convergence when limit or maxNumOfIters is satisfied
    bool reachMaxNumOfIters = (currentNumOfIterations == MAX_NUM_ITERATIONS);
    bool isConvergent = true;
    for (int i = 0; i < checkResiduals.size(); i++) {
        isConvergent = isConvergent && checkResiduals[i]->isConvergent();
    }

    if (isConvergent || reachMaxNumOfIters) {
        if (reachMaxNumOfIters) {
            WARNING_BLOCK_OUT("ConvergenceChecker", "isConvergent()",
                    "reach maximum number of iterations!");
        }
        currentNumOfIterations = 0;
        timeStepNumber++;
        return true;
    }

    return false;
}

void ConvergenceChecker::addCheckResidual(double _absoluteTolerance, double _relativeTolerance,
        AbstractCouplingAlgorithm *_couplingAlgorithm, int _residualIndex) {
    checkResiduals.push_back(
            new CheckResidual(_absoluteTolerance, _relativeTolerance, _couplingAlgorithm,
                    _residualIndex));
}

int ConvergenceChecker::getcurrentNumOfIterations() {
    return currentNumOfIterations;
}

} /* namespace EMPIRE */
