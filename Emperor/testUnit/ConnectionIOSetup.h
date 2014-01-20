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
/***********************************************************************************************//**
 * \file ConnectionIOSetup.h
 * This file holds the class ConnectionIOSetup
 * \date 11/21/2012
 **************************************************************************************************/
#ifndef CONNECTIONIOSETUP_H_
#define CONNECTIONIOSETUP_H_

#include "AbstractFilter.h"
#include "AbstractCouplingAlgorithm.h"
#include "AbstractExtrapolator.h"
#include "ConnectionIO.h"
//#include <iostream>
//using namespace std;

namespace EMPIRE {

class ConnectionIOSetup {
public:
    static void setupIOForFilter(AbstractFilter *filter, AbstractMesh *meshIn, DataField *in,
            AbstractMesh *meshOut, DataField *out) {
        ConnectionIO *connectionInput = new ConnectionIO();
        connectionInput->type = EMPIRE_ConnectionIO_DataField;
        connectionInput->mesh = meshIn;
        connectionInput->dataField = in;
        ConnectionIO *connectionOutput = new ConnectionIO();
        connectionOutput->type = EMPIRE_ConnectionIO_DataField;
        connectionOutput->mesh = meshOut;
        connectionOutput->dataField = out;
        filter->addInput(connectionInput);
        filter->addOutput(connectionOutput);
        filter->init();
    }
    static void setupIOForFilter(AbstractFilter *filter, Signal *in, Signal *out) {
        ConnectionIO *connectionInput = new ConnectionIO();
        connectionInput->type = EMPIRE_ConnectionIO_Signal;
        connectionInput->signal = in;
        ConnectionIO *connectionOutput = new ConnectionIO();
        connectionOutput->type = EMPIRE_ConnectionIO_Signal;
        connectionOutput->signal = out;
        filter->addInput(connectionInput);
        filter->addOutput(connectionOutput);
        filter->init();
    }
    static ConnectionIO *constructDummyConnectionIO(DataField *dataField) {
        ConnectionIO *dummy = new ConnectionIO();
        dummy->type = EMPIRE_ConnectionIO_DataField;
        dummy->dataField = dataField;
        return dummy;
    }
    static ConnectionIO *constructDummyConnectionIO(Signal *signal) {
        ConnectionIO *dummy = new ConnectionIO();
        dummy->type = EMPIRE_ConnectionIO_Signal;
        dummy->signal = signal;
        return dummy;
    }
};

} /* namespace EMPIRE */
#endif /* CONNECTIONIOSETUP_H_ */
