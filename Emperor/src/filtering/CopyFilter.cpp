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
#include "CopyFilter.h"
#include <assert.h>
#include <stdlib.h>
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"

namespace EMPIRE {



CopyFilter::CopyFilter(int _signalOffset) :
        AbstractFilter(), signalOffset(_signalOffset) {

	std::cout << "SERVUS: " << signalOffset <<std::endl;
}
CopyFilter::~CopyFilter() {
}
void CopyFilter::filtering() {
    EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        const DataField *inDataField = inputVec[0]->dataField;
        DataField *outDataField = outputVec[0]->dataField;
        assert(inDataField!=NULL);
        assert(outDataField!=NULL);

        int sizeIn = inDataField->numLocations * inDataField->dimension;
        int sizeOut = outDataField->numLocations * outDataField->dimension;

        if (sizeIn >= sizeOut) {
            for (int i = 0; i < sizeOut; i++)
                outDataField->data[i] = inDataField->data[i];
        } else {
            for (int i = 0; i < sizeIn; i++)
                outDataField->data[i] = inDataField->data[i];
            for (int i = sizeIn; i < sizeOut; i++)
                outDataField->data[i] = 0.0;
        }
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        const Signal *inSignal = inputVec[0]->signal;
        Signal *outSignal = outputVec[0]->signal;
        assert(inSignal!=NULL);
        assert(outSignal!=NULL);

        int sizeIn = inSignal->size;
        int sizeOut = outSignal->size;

        if (sizeIn >= sizeOut) {
            for (int i = 0; i < sizeOut; i++)
                outSignal->array[i] = inSignal->array[i+signalOffset];
        } else {
            for (int i = 0; i < sizeIn; i++)
                outSignal->array[i] = inSignal->array[i];
            for (int i = sizeIn; i < sizeOut; i++)
                outSignal->array[i] = 0.0;
        }
    } else {
        assert(false);
    }
}

void CopyFilter::init() {
    assert(inputVec.size() == 1);
    assert(outputVec.size() == 1);
    assert(inputVec[0]->type == outputVec[0]->type);
}

} /* namespace EMPIRE */

