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
#include "ScalingFilter.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"

namespace EMPIRE {

ScalingFilter::ScalingFilter(double _factor) :
        AbstractFilter(), factor(_factor) {
}

ScalingFilter::~ScalingFilter() {
}

void ScalingFilter::filtering() {
    EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        DataField *inDataField = inputVec[0]->dataField;
        for (int i = 0; i < inDataField->numLocations * inDataField->dimension; i++) {
            if (factor == 0.0)
                inDataField->data[i] = 0.0;
            else
                inDataField->data[i] *= factor;
        }
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        Signal *inSignal = inputVec[0]->signal;
        for (int i = 0; i < inSignal->size; i++) {
            inSignal->array[i] *= factor;
        }
    } else {
        assert(false);
    }
}

void ScalingFilter::init() {
    assert(inputVec.size() == 1);
    assert(outputVec.size() == 1);
    assert(inputVec[0]->type == outputVec[0]->type);
    EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        assert(inputVec[0]->dataField == outputVec[0]->dataField);
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        assert(inputVec[0]->signal == outputVec[0]->signal);
    } else {
        assert(false);
    }
}

} /* namespace EMPIRE */
