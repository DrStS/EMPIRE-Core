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
#include <math.h>

#include "Residual.h"
#include "ConnectionIO.h"
#include "DataField.h"
#include "Signal.h"

namespace EMPIRE {

Residual::Component::Component(double _coefficient, std::string _timeToUpdate,
        ConnectionIO *_reference) :
        coefficient(_coefficient), timeToUpdate(_timeToUpdate), reference(_reference) {
    if (reference->type == EMPIRE_ConnectionIO_DataField) {
        size = reference->dataField->dimension * reference->dataField->numLocations;
        dataCopy = new double[size];
    } else if (reference->type == EMPIRE_ConnectionIO_Signal) {
        size = reference->signal->size;
        dataCopy = new double[size];
    } else {
        assert(false);
    }
}

Residual::Component::~Component() {
    delete[] dataCopy;
    delete reference;
}

void Residual::Component::updateAtIterationBeginning() {
    if (timeToUpdate == "iterationBeginning") {
        if (reference->type == EMPIRE_ConnectionIO_DataField) {
            for (int i = 0; i < size; i++)
                dataCopy[i] = reference->dataField->data[i];
        } else if (reference->type == EMPIRE_ConnectionIO_Signal) {
            for (int i = 0; i < size; i++)
                dataCopy[i] = reference->signal->array[i];
        } else {
            assert(false);
        }
    } else {
        assert(timeToUpdate == "iterationEnd");
    }
}

void Residual::Component::updateAtIterationEnd() {
    if (timeToUpdate == "iterationEnd") {
        if (reference->type == EMPIRE_ConnectionIO_DataField) {
            for (int i = 0; i < size; i++)
                dataCopy[i] = reference->dataField->data[i];
        } else if (reference->type == EMPIRE_ConnectionIO_Signal) {
            for (int i = 0; i < size; i++)
                dataCopy[i] = reference->signal->array[i];
        } else {
            assert(false);
        }
    } else {
        assert(timeToUpdate == "iterationBeginning");
    }
}

Residual::Residual(int _index) :
        index(_index) {
    size = 0;
    residualVector = NULL;
    residualVectorL2Norm = 0.0;
}

Residual::~Residual() {
    for (int i = 0; i < components.size(); i++)
        delete components[i];
    delete[] residualVector;
}

void Residual::init() {
    size = components[0]->size;
    for (int i = 1; i < components.size(); i++) {
        assert(components[i]->size == size);
    }
    residualVector = new double[size];
}

void Residual::updateAtIterationBeginning() {
    for (int i = 0; i < components.size(); i++) {
        components[i]->updateAtIterationBeginning();
    }

}
void Residual::updateAtIterationEnd() {
    for (int i = 0; i < components.size(); i++) {
        components[i]->updateAtIterationEnd();
    }
}

void Residual::addComponent(double _coefficient, std::string _timeToUpdate,
        ConnectionIO *_reference) {
    components.push_back(new Component(_coefficient, _timeToUpdate, _reference));
}

void Residual::computeCurrentResidual() {
    // compute the residual vector
    for (int i=0; i<size; i++)
        residualVector[i] = 0.0;
    for (int i = 0; i < components.size(); i++) {
        for (int j=0; j<size; j++) {
            residualVector[j] += components[i]->dataCopy[j] * components[i]->coefficient;
        }
    }
    // compute the L2 norm of the residual vector
    residualVectorL2Norm = 0.0;
    for (int i = 0; i < size; i++) {
        residualVectorL2Norm += residualVector[i] * residualVector[i];
    }
    residualVectorL2Norm /= size;
    residualVectorL2Norm = sqrt(residualVectorL2Norm);
}

} /* namespace EMPIRE */
