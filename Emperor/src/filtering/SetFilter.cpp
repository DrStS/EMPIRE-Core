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
#include "SetFilter.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"

namespace EMPIRE {

SetFilter::SetFilter(std::vector<double> _value) :
        AbstractFilter(), value(_value) {
}

SetFilter::~SetFilter() {
}

void SetFilter::filtering() {
    EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
    	assert(false); // Not yet implemented for fields
        DataField *inDataField = inputVec[0]->dataField;
        /*for (int i = 0; i < inDataField->numLocations * inDataField->dimension; i++) {
            inDataField->data[i] = value;
        }*/
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        Signal *inSignal = inputVec[0]->signal;
        assert(inSignal->size==value.size());
        for (int i = 0; i < inSignal->size; i++) {
            inSignal->array[i] = value[i];
        }
    } else {
        assert(false);
    }
}

void SetFilter::init() {
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
