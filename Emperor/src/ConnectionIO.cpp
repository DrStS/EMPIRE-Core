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
#include "ConnectionIO.h"
#include "ClientCode.h"
#include "Signal.h"
#include "DataField.h"
#include "MetaDataStructures.h"
#include "AbstractMesh.h"
#include <assert.h>

namespace EMPIRE {

ConnectionIO::ConnectionIO(std::map<std::string, ClientCode*> &nameToClientCodeMap
        , std::string clientCodeName, std::string meshName, std::string dataFieldName) {
    type = EMPIRE_ConnectionIO_DataField;
    clientCode = nameToClientCodeMap[clientCodeName];
    mesh = clientCode->getMeshByName(meshName);
    dataField = mesh->getDataFieldByName(dataFieldName);
    signal = NULL;
}

ConnectionIO::ConnectionIO(std::map<std::string, ClientCode*> &nameToClientCodeMap
        , std::string clientCodeName, std::string signalName) {
    type = EMPIRE_ConnectionIO_Signal;
    clientCode = nameToClientCodeMap[clientCodeName];
    signal = clientCode->getSignalByName(signalName);
    mesh = NULL;
    dataField = NULL;
}

ConnectionIO::ConnectionIO(ClientCode *_clientCode, AbstractMesh *_mesh, DataField *_dataField) {
    type = EMPIRE_ConnectionIO_DataField;
    clientCode = _clientCode;
    mesh = _mesh;
    dataField = _dataField;
    signal = NULL;
}

ConnectionIO::ConnectionIO(ClientCode *_clientCode, Signal *_signal) {
    type = EMPIRE_ConnectionIO_Signal;
    clientCode = _clientCode;
    signal = _signal;
    mesh = NULL;
    dataField = NULL;
}

ConnectionIO::ConnectionIO() {
    clientCode = NULL;
    signal = NULL;
    mesh = NULL;
    dataField = NULL;
}

ConnectionIO::~ConnectionIO() {
}

void ConnectionIO::send() {
    if (type == EMPIRE_ConnectionIO_DataField)
        clientCode->sendDataField(mesh->name, dataField->name);
    else if (type == EMPIRE_ConnectionIO_Signal)
        clientCode->sendSignal(signal->name);
    else
        assert(false);
}

void ConnectionIO::receive() {
    if (type == EMPIRE_ConnectionIO_DataField)
        clientCode->recvDataField(mesh->name, dataField->name);
    else if (type == EMPIRE_ConnectionIO_Signal)
        clientCode->recvSignal(signal->name);
    else
        assert(false);
}

} /* namespace EMPIRE */
