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
#include "ExtrapolatingFilter.h"
#include "ConnectionIO.h"
#include "AbstractExtrapolator.h"
#include <assert.h>

namespace EMPIRE {

ExtrapolatingFilter::ExtrapolatingFilter(AbstractExtrapolator *_extrapolator) :
        extrapolator(_extrapolator) {
    assert(_extrapolator != NULL);
}

ExtrapolatingFilter::~ExtrapolatingFilter() {
}

void ExtrapolatingFilter::filtering() {
    extrapolator->extrapolate();
}

void ExtrapolatingFilter::init() {
  if (inputVec.size() == 1 && outputVec.size() == 1) {
    assert(inputVec[0]->type == outputVec[0]->type);
    extrapolator->setInputAndOutput(inputVec[0], outputVec[0]);
  }
  else {
    for (int i=0; i<inputVec.size(); i++) {
      extrapolator->addInput(inputVec[i]);
    }
    for (int i=0; i<outputVec.size(); i++) {
      extrapolator->addOutput(outputVec[i]);
    }
  }
    // do the following if different input and output are forbidden
    /*EMPIRE_ConnectionIO_Type IOType = inputVec[0]->type;
    if (IOType == EMPIRE_ConnectionIO_DataField) {
        assert(inputVec[0]->dataField == outputVec[0]->dataField);
    } else if (IOType == EMPIRE_ConnectionIO_Signal) {
        assert(inputVec[0]->signal == outputVec[0]->signal);
    } else {
        assert(false);
    }*/
}

} /* namespace EMPIRE */
