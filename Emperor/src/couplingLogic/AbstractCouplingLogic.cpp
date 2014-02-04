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
#include "AbstractCouplingLogic.h"
#include "Connection.h"
#include "DataOutput.h"

namespace EMPIRE {

using namespace std;

AbstractCouplingLogic::AbstractCouplingLogic() :
        couplingLogicSequence() {
}

AbstractCouplingLogic::~AbstractCouplingLogic() {
}

void AbstractCouplingLogic::addCouplingLogic(AbstractCouplingLogic *couplingLogic) {
    couplingLogicSequence.push_back(couplingLogic);
}

int AbstractCouplingLogic::size() {
    return couplingLogicSequence.size();
}

void AbstractCouplingLogic::addDataOutput(DataOutput *dataOutput) {
    dataOutputVec.push_back(dataOutput);
}


} /* namespace EMPIRE */
