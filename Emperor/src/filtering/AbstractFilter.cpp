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
#include "AbstractFilter.h"
#include <assert.h>
#include "ConnectionIO.h"

namespace EMPIRE {

AbstractFilter::AbstractFilter() :
        inputVec(), outputVec() {
}

AbstractFilter::~AbstractFilter() {
    for (int i = 0; i < inputVec.size(); i++)
        delete inputVec[i];
    for (int i = 0; i < outputVec.size(); i++)
        delete outputVec[i];
}

void AbstractFilter::addInput(ConnectionIO *input) {
    assert(input!=NULL);
    inputVec.push_back(input);
}

void AbstractFilter::addOutput(ConnectionIO *output) {
    assert(output!=NULL);
    outputVec.push_back(output);
}

} /* namespace EMPIRE */
