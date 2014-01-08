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
#include <string>
#include <assert.h>
#include <math.h>

#include "AbstractCouplingAlgorithm.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include "Residual.h"
#include "EMPEROR_Enum.h"
#include <iostream>

using namespace std;

namespace EMPIRE {

AbstractCouplingAlgorithm::CouplingAlgorithmOutput::CouplingAlgorithmOutput(
        ConnectionIO *_reference) :
        reference(_reference) {
    if (reference->type == EMPIRE_ConnectionIO_DataField) {
        size = reference->dataField->dimension * reference->dataField->numLocations;
        outputCopyAtIterationBeginning = new double[size];
    } else if (reference->type == EMPIRE_ConnectionIO_Signal) {
        size = reference->signal->size;
        outputCopyAtIterationBeginning = new double[size];
    } else {
        assert(false);
    }
}

AbstractCouplingAlgorithm::CouplingAlgorithmOutput::~CouplingAlgorithmOutput() {
    delete reference;
    delete[] outputCopyAtIterationBeginning;
}

void AbstractCouplingAlgorithm::CouplingAlgorithmOutput::updateAtIterationBeginning() {
    if (reference->type == EMPIRE_ConnectionIO_DataField) {
        for (int i = 0; i < size; i++)
            outputCopyAtIterationBeginning[i] = reference->dataField->data[i];
    } else if (reference->type == EMPIRE_ConnectionIO_Signal) {
        for (int i = 0; i < size; i++)
            outputCopyAtIterationBeginning[i] = reference->signal->array[i];
    } else {
        assert(false);
    }
}

void AbstractCouplingAlgorithm::CouplingAlgorithmOutput::overwrite(double *newData) {
    if (reference->type == EMPIRE_ConnectionIO_DataField) {
        for (int i = 0; i < size; i++)
            reference->dataField->data[i] = newData[i];
    } else if (reference->type == EMPIRE_ConnectionIO_Signal) {
        for (int i = 0; i < size; i++)
            reference->signal->array[i] = newData[i];
    } else {
        assert(false);
    }
}

AbstractCouplingAlgorithm::AbstractCouplingAlgorithm(std::string _name) :
        name(_name), newTimeStep(false) {
}

AbstractCouplingAlgorithm::~AbstractCouplingAlgorithm() {
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++)
        delete it->second;
    for (map<int, CouplingAlgorithmOutput*>::iterator it = outputs.begin(); it != outputs.end();
            it++)
        delete it->second;
}

void AbstractCouplingAlgorithm::addResidual(Residual *residual, int index) {
    residuals.insert(residuals.begin(), pair<int, Residual*>(index, residual));
}

void AbstractCouplingAlgorithm::addOutput(ConnectionIO *output, int index) {
    outputs.insert(outputs.begin(),
            pair<int, CouplingAlgorithmOutput*>(index, new CouplingAlgorithmOutput(output)));
}

void AbstractCouplingAlgorithm::updateAtIterationBeginning() {
    // update residuals
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++)
        it->second->updateAtIterationBeginning();
    // update output copys
    for (map<int, CouplingAlgorithmOutput*>::iterator it = outputs.begin(); it != outputs.end();
            it++) {
        it->second->updateAtIterationBeginning();
    }
}

void AbstractCouplingAlgorithm::updateAtIterationEnd() {
    // update residuals
    for (map<int, Residual*>::iterator it = residuals.begin(); it != residuals.end(); it++)
        it->second->updateAtIterationEnd();
}

void AbstractCouplingAlgorithm::setNewTimeStep() {
	newTimeStep = true;
}

double AbstractCouplingAlgorithm::getResidualL2Norm(int index) {
    return residuals[index]->residualVectorL2Norm;
}

std::string AbstractCouplingAlgorithm::getName() {
    return name;
}

double AbstractCouplingAlgorithm::vecL2Norm(const double *vec, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += vec[i] * vec[i];
    }
    sum /= size;
    sum = sqrt(sum);
    return sum;
}

} /* namespace EMPIRE */
