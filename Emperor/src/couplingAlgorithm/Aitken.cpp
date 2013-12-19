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
#include <sstream>
#include <math.h>

#include "Aitken.h"
#include "DataField.h"
#include "Message.h"
#include "Signal.h"
#include "ConnectionIO.h"
#include "EMPEROR_Enum.h"

using namespace std;

namespace EMPIRE {

const double Aitken::LIMIT = 1E6;

Aitken::Aitken(std::string _name, double _initialRelaxationFactor) :
        AbstractCouplingAlgorithm(_name), INIT_AITKEN_FACTOR(_initialRelaxationFactor) {
    assert(INIT_AITKEN_FACTOR >= 0.0);
}

Aitken::~Aitken() {
}

void Aitken::calcNewValue() {
    assert(false);
}

void Aitken::setZeroState() {

}

} /* namespace EMPIRE */
