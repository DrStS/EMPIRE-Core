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
#include <iostream>
#include <assert.h>

#include "Connection.h"
#include "ClientCode.h"
#include "DataField.h"
#include "AbstractFilter.h"
#include "ConnectionIO.h"

using namespace std;

namespace EMPIRE {
Connection::Connection(std::string _name) :
        AbstractCouplingLogic(), name(_name), inputVec(), outputVec(), filterVec() {
}

Connection::~Connection() {
    for (int i = 0; i < inputVec.size(); i++)
        delete inputVec[i];
    for (int i = 0; i < outputVec.size(); i++)
        delete outputVec[i];
    for (int i = 0; i < filterVec.size(); i++)
        delete filterVec[i];
}

void Connection::doCoupling() {
    transferData();
}

void Connection::transferData() {
    for (unsigned i = 0; i < inputVec.size(); i++)
        inputVec[i]->receive();
    for (unsigned i = 0; i < filterVec.size(); i++)
        filterVec[i]->filtering();
    for (unsigned i = 0; i < outputVec.size(); i++)
        outputVec[i]->send();
}

void Connection::addInput(ConnectionIO *input) {
    assert(input!=NULL);
    inputVec.push_back(input);
}

void Connection::addOutput(ConnectionIO *output) {
    assert(output!=NULL);
    outputVec.push_back(output);
}

void Connection::addFilter(AbstractFilter *filter) {
    assert(filter!=NULL);
    filterVec.push_back(filter);
}

} /* namespace EMPIRE */
