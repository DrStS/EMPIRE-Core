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

#include "CouplingAlgorithmFilter.h"
#include "AbstractCouplingAlgorithm.h"
#include "DataField.h"
#include "Signal.h"
#include "ConnectionIO.h"

namespace EMPIRE {

CouplingAlgorithmFilter::CouplingAlgorithmFilter(AbstractCouplingAlgorithm *_couplingAlgorithm) :
        AbstractFilter(), couplingAlgorithm(_couplingAlgorithm) {
    assert(couplingAlgorithm!=NULL);
}

void CouplingAlgorithmFilter::init() {
    assert(inputVec.size() == 1);
    assert(outputVec.size() == 1);
    assert(inputVec[0]->type == outputVec[0]->type);
    couplingAlgorithm->setInputAndOutput(inputVec[0], outputVec[0]);
}

CouplingAlgorithmFilter::~CouplingAlgorithmFilter() {
}

void CouplingAlgorithmFilter::filtering() {
    couplingAlgorithm->calcNewValue();
}

} /* namespace EMPIRE */
