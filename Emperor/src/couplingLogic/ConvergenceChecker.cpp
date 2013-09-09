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

const double ConvergenceChecker::DEFAULT_ABS_TOL = -1.0;
const double ConvergenceChecker::DEFAULT_REL_TOL = -1.0;
const double ConvergenceChecker::DEFAULT_MAX_NUM_ITERATIONS = 1e100;

ConvergenceChecker::ConvergenceChecker(double absTol, double relTol, double maxNumOfIters) :
        ABS_TOL(absTol), REL_TOL(relTol), MAX_NUM_ITERATIONS(maxNumOfIters) {
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

    dataField = NULL;
    signal = NULL;
    couplingAlgorithm = NULL;

    dataLastInnerLoopStep = NULL;
}

ConvergenceChecker::~ConvergenceChecker() {
    if (dataLastInnerLoopStep != NULL)
        delete[] dataLastInnerLoopStep;
}

void ConvergenceChecker::setDataField(DataField *_dataField) {
    assert(_dataField != NULL);
    assert(dataField == NULL);
    assert(signal == NULL);
    assert(couplingAlgorithm == NULL);
    dataField = _dataField;
    int size = dataField->numLocations * dataField->dimension;
    dataLastInnerLoopStep = new double[size];
    whichRef = EMPIRE_ConvergenceChecker_dataFieldRef;
}

void ConvergenceChecker::setSignal(Signal *_signal) {
    assert(_signal != NULL);
    assert(dataField == NULL);
    assert(signal == NULL);
    assert(couplingAlgorithm == NULL);
    signal = _signal;
    int size = signal->size;
    dataLastInnerLoopStep = new double[size];
    whichRef = EMPIRE_ConvergenceChecker_signalRef;
}

void ConvergenceChecker::setCouplingAlgorithm(AbstractCouplingAlgorithm *_couplingAlgorithm) {
    assert(_couplingAlgorithm != NULL);
    assert(dataField == NULL);
    assert(signal == NULL);
    assert(couplingAlgorithm == NULL);
    couplingAlgorithm = _couplingAlgorithm;
    whichRef = EMPIRE_ConvergenceChecker_couplingAlgorithmRef;
}

bool ConvergenceChecker::isConvergent() {
    currentNumOfIterations++;
    fstream residualFile;
    residualFile.open(residualFileName.c_str(), ios_base::out | ios_base::app);
    const double *dataCurrent = NULL;
    int size = -1;
    { // 1. compute residuals
        if (whichRef == EMPIRE_ConvergenceChecker_dataFieldRef) {
            dataCurrent = dataField->data;
            size = dataField->numLocations * dataField->dimension;
            if (currentNumOfIterations == 1) {
                for (int i = 0; i < size; i++)
                    dataLastInnerLoopStep[i] = dataCurrent[i];
                return false;
            } else if (currentNumOfIterations == 2) {
                initialResidual = calcDifferenceL2Norm(dataLastInnerLoopStep, dataCurrent, size);
                for (int i = 0; i < size; i++)
                    dataLastInnerLoopStep[i] = dataCurrent[i];
                currentResidual = initialResidual;
            } else {
                currentResidual = calcDifferenceL2Norm(dataLastInnerLoopStep, dataCurrent, size);
            }
        } else if (whichRef == EMPIRE_ConvergenceChecker_signalRef) {
            dataCurrent = signal->array;
            size = signal->size;
            if (currentNumOfIterations == 1) {
                for (int i = 0; i < size; i++)
                    dataLastInnerLoopStep[i] = dataCurrent[i];
                return false;
            } else if (currentNumOfIterations == 2) {
                initialResidual = calcDifferenceL2Norm(dataLastInnerLoopStep, dataCurrent, size);
                for (int i = 0; i < size; i++)
                    dataLastInnerLoopStep[i] = dataCurrent[i];
                currentResidual = initialResidual;
            } else {
                currentResidual = calcDifferenceL2Norm(dataLastInnerLoopStep, dataCurrent, size);
            }
        } else if (whichRef == EMPIRE_ConvergenceChecker_couplingAlgorithmRef) {
            if (currentNumOfIterations == 1) {
                return false;
            } else if (currentNumOfIterations == 2) {
                initialResidual = couplingAlgorithm->getInitialResidual();
                currentResidual = initialResidual;
            } else {
                currentResidual = couplingAlgorithm->getCurrentResidual();
            }
        } else {
            assert(false);
        }
    }
    // 2. output residuals
    if (debugResidual) {
        {
            stringstream ss;
            ss << "ConvergenceChecker absolute residual: " << currentResidual << endl;
            INDENT_OUT(1, ss.str(), infoOut);
        }
        {
            stringstream ss;
            ss << "ConvergenceChecker relative residual: " << currentResidual / initialResidual;
            INDENT_OUT(1, ss.str(), infoOut);
        }
    }

    residualFile << timeStepNumber << '\t' << currentNumOfIterations << '\t'
            << currentResidual / initialResidual << '\t' << currentResidual << endl;

    // 3. set convergence when limit or maxNumOfIters is satisfied
    bool reachMaxNumOfIters = (currentNumOfIterations == MAX_NUM_ITERATIONS);
    if ((currentResidual < ABS_TOL) || (currentResidual / initialResidual < REL_TOL)
            || reachMaxNumOfIters) {
        if (reachMaxNumOfIters) {
            WARNING_BLOCK_OUT("ConvergenceChecker", "isConvergent()",
                    "reach maximum number of iterations!");
        }
        currentNumOfIterations = 0;
        timeStepNumber++;
        return true;
    }

    // 4. update dataLastInnerLoopStep
    if (dataField != NULL || signal != NULL) {
        for (int i = 0; i < size; i++)
            dataLastInnerLoopStep[i] = dataCurrent[i];
    }

    return false;
}

int ConvergenceChecker::getcurrentNumOfIterations() {
    return currentNumOfIterations;
}

double ConvergenceChecker::calcDifferenceL2Norm(const double *array1, const double *array2,
        int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double tmp = array1[i] - array2[i];
        sum += tmp * tmp;
    }
    sum /= size;
    sum = sqrt(sum);
    return sum;
}
} /* namespace EMPIRE */
